suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

input_file <- "work/orthofinder_plots/new_master_PAV.csv"
out_dir <- "outputs/pav_concordance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c("minigraph", "pggb", "orthofinder")
method_labels <- c(
  minigraph = "Minigraph",
  pggb = "PGGB",
  orthofinder = "OrthoFinder"
)
accessions <- c("a55015", "Bristol", "DSV_1", "DSV_2", "ERA3543")

master <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_cols <- c("geneID", as.vector(outer(methods, accessions, paste, sep = "_")))
missing_cols <- setdiff(expected_cols, colnames(master))
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

pav_cols <- setdiff(expected_cols, "geneID")

to_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

raw_numeric <- master
raw_numeric[pav_cols] <- lapply(raw_numeric[pav_cols], to_numeric)

value_summary <- lapply(methods, function(method) {
  method_cols <- paste(method, accessions, sep = "_")
  values <- unlist(raw_numeric[method_cols], use.names = FALSE)
  data.frame(
    method = method_labels[[method]],
    n_cells = length(values),
    n_0 = sum(values == 0, na.rm = TRUE),
    n_1 = sum(values == 1, na.rm = TRUE),
    n_gt1 = sum(values > 1, na.rm = TRUE),
    n_missing_or_non_numeric = sum(is.na(values)),
    n_other_numeric = sum(!is.na(values) & !(values %in% c(0, 1)) & values <= 1),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

write.table(
  value_summary,
  file = file.path(out_dir, "new_master_PAV_value_cleanup_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

binary <- raw_numeric
binary[pav_cols] <- lapply(binary[pav_cols], function(x) {
  ifelse(x %in% c(0, 1), as.integer(x), NA_integer_)
})

write.csv(
  binary,
  file = file.path(out_dir, "new_master_PAV.binary_0_1_NA.csv"),
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)

method_matrix <- function(df, method) {
  as.matrix(df[paste(method, accessions, sep = "_")])
}

is_valid_pattern <- function(mat) {
  rowSums(is.na(mat)) == 0
}

same_pattern <- function(mat_a, mat_b) {
  rowSums(mat_a == mat_b) == ncol(mat_a)
}

same_pattern_allowing_matching_na <- function(mat_a, mat_b) {
  same_non_na <- mat_a == mat_b
  same_non_na[is.na(same_non_na)] <- FALSE
  same_na <- is.na(mat_a) & is.na(mat_b)
  rowSums(same_non_na | same_na) == ncol(mat_a)
}

pairwise <- list()
pair_index <- 1

for (method_a in methods) {
  mat_a <- method_matrix(binary, method_a)
  valid_a <- is_valid_pattern(mat_a)

  for (method_b in methods) {
    mat_b <- method_matrix(binary, method_b)
    valid_b <- is_valid_pattern(mat_b)
    comparable <- valid_a & valid_b
    concordant <- comparable & same_pattern(mat_a, mat_b)

    pairwise[[pair_index]] <- data.frame(
      method_a = method_labels[[method_a]],
      method_b = method_labels[[method_b]],
      n_total_genes = nrow(binary),
      n_valid_method_a = sum(valid_a),
      n_valid_method_b = sum(valid_b),
      n_comparable_genes = sum(comparable),
      n_concordant_genes = sum(concordant),
      n_discordant_genes = sum(comparable) - sum(concordant),
      concordance_rate = ifelse(sum(comparable) > 0, sum(concordant) / sum(comparable), NA_real_),
      stringsAsFactors = FALSE
    )
    pair_index <- pair_index + 1
  }
}

pairwise_summary <- bind_rows(pairwise) %>%
  mutate(
    method_a = factor(method_a, levels = method_labels[methods]),
    method_b = factor(method_b, levels = method_labels[methods]),
    concordance_percent = 100 * concordance_rate
  )

write.table(
  pairwise_summary,
  file = file.path(out_dir, "method_pairwise_exact_pattern_concordance_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rate_matrix <- pairwise_summary %>%
  select(method_a, method_b, concordance_rate) %>%
  pivot_wider(names_from = method_b, values_from = concordance_rate)

write.table(
  rate_matrix,
  file = file.path(out_dir, "method_pairwise_exact_pattern_concordance_rate_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

count_matrix <- pairwise_summary %>%
  select(method_a, method_b, n_comparable_genes, n_concordant_genes, n_discordant_genes) %>%
  pivot_wider(
    names_from = method_b,
    values_from = c(n_comparable_genes, n_concordant_genes, n_discordant_genes)
  )

write.table(
  count_matrix,
  file = file.path(out_dir, "method_pairwise_exact_pattern_concordance_count_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

all_mats <- lapply(methods, function(method) method_matrix(binary, method))
names(all_mats) <- methods
all_valid <- Reduce(`&`, lapply(all_mats, is_valid_pattern))
all_same <- all_valid &
  same_pattern(all_mats[["minigraph"]], all_mats[["pggb"]]) &
  same_pattern(all_mats[["minigraph"]], all_mats[["orthofinder"]])

all_method_summary <- data.frame(
  n_total_genes = nrow(binary),
  n_all_three_methods_comparable = sum(all_valid),
  n_all_three_methods_concordant = sum(all_same),
  n_all_three_methods_discordant = sum(all_valid) - sum(all_same),
  all_three_concordance_rate = ifelse(sum(all_valid) > 0, sum(all_same) / sum(all_valid), NA_real_)
)

write.table(
  all_method_summary,
  file = file.path(out_dir, "all_three_methods_exact_pattern_concordance_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

heatmap_df <- pairwise_summary %>%
  mutate(
    method_a_index = as.integer(method_a),
    method_b_index = as.integer(method_b),
    label = ifelse(
      is.na(concordance_rate),
      "NA",
      sprintf("%.1f%%\n(n=%s)", concordance_percent, format(n_comparable_genes, big.mark = ","))
    ),
    label_no_repeat_bottom_right = ifelse(method_b_index > method_a_index, "", label)
  )

p_heat <- ggplot(heatmap_df, aes(x = method_b, y = method_a, fill = concordance_percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 4.2, lineheight = 0.9) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08519c",
    limits = c(0, 100),
    na.value = "grey90",
    name = "Concordance (%)"
  ) +
  coord_equal() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Exact PAV-pattern concordance between methods") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "method_pairwise_exact_pattern_concordance_heatmap.pdf"),
  p_heat,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

p_heat_no_repeat <- ggplot(heatmap_df, aes(x = method_b, y = method_a, fill = concordance_percent)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label_no_repeat_bottom_right), size = 4.2, lineheight = 0.9) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08519c",
    limits = c(0, 100),
    na.value = "grey90",
    name = "Concordance (%)"
  ) +
  coord_equal() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Exact PAV-pattern concordance between methods") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "method_pairwise_exact_pattern_concordance_heatmap_no_repeat_bottom_right.pdf"),
  p_heat_no_repeat,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

fixed_denominator_pairs <- list()
fixed_pair_index <- 1

for (method_a in methods) {
  mat_a <- method_matrix(binary, method_a)
  valid_a <- is_valid_pattern(mat_a)

  for (method_b in methods) {
    mat_b <- method_matrix(binary, method_b)
    valid_b <- is_valid_pattern(mat_b)
    comparable <- valid_a & valid_b
    concordant_fixed <- same_pattern_allowing_matching_na(mat_a, mat_b)

    fixed_denominator_pairs[[fixed_pair_index]] <- data.frame(
      method_a = method_labels[[method_a]],
      method_b = method_labels[[method_b]],
      n_denominator_genes = nrow(binary),
      n_valid_method_a = sum(valid_a),
      n_valid_method_b = sum(valid_b),
      n_comparable_genes = sum(comparable),
      n_concordant_genes = sum(concordant_fixed),
      n_nonconcordant_genes = nrow(binary) - sum(concordant_fixed),
      fixed_denominator_concordance_rate = sum(concordant_fixed) / nrow(binary),
      fixed_denominator_concordance_percent = 100 * sum(concordant_fixed) / nrow(binary),
      stringsAsFactors = FALSE
    )
    fixed_pair_index <- fixed_pair_index + 1
  }
}

fixed_denominator_summary <- bind_rows(fixed_denominator_pairs) %>%
  mutate(
    method_a = factor(method_a, levels = method_labels[methods]),
    method_b = factor(method_b, levels = method_labels[methods])
  )

write.table(
  fixed_denominator_summary,
  file = file.path(out_dir, "method_pairwise_exact_pattern_concordance_fixed_denominator_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

fixed_rate_matrix <- fixed_denominator_summary %>%
  select(method_a, method_b, fixed_denominator_concordance_rate) %>%
  pivot_wider(names_from = method_b, values_from = fixed_denominator_concordance_rate)

write.table(
  fixed_rate_matrix,
  file = file.path(out_dir, "method_pairwise_exact_pattern_concordance_fixed_denominator_rate_matrix.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

fixed_heatmap_df <- fixed_denominator_summary %>%
  mutate(
    method_a_index = as.integer(method_a),
    method_b_index = as.integer(method_b),
    label = sprintf("%.1f%%", fixed_denominator_concordance_percent),
    label_no_repeat_bottom_right = ifelse(method_b_index > method_a_index, "", label)
  )

p_heat_fixed <- ggplot(
  fixed_heatmap_df,
  aes(x = method_b, y = method_a, fill = fixed_denominator_concordance_percent)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 4.8) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08519c",
    limits = c(0, 100),
    na.value = "grey90",
    name = "Concordance (%)"
  ) +
  coord_equal() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Exact PAV-pattern concordance between methods") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "method_pairwise_exact_pattern_concordance_heatmap_fixed_denominator.pdf"),
  p_heat_fixed,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

p_heat_fixed_no_repeat <- ggplot(
  fixed_heatmap_df,
  aes(x = method_b, y = method_a, fill = fixed_denominator_concordance_percent)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label_no_repeat_bottom_right), size = 4.8) +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#08519c",
    limits = c(0, 100),
    na.value = "grey90",
    name = "Concordance (%)"
  ) +
  coord_equal() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Exact PAV-pattern concordance between methods") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "method_pairwise_exact_pattern_concordance_heatmap_fixed_denominator_no_repeat_bottom_right.pdf"),
  p_heat_fixed_no_repeat,
  width = 6.5,
  height = 5.5,
  dpi = 300
)

message("Wrote concordance outputs to: ", out_dir)
