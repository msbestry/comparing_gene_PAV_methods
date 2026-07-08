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

to_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

for (method in methods) {
  for (accession in accessions) {
    col <- paste(method, accession, sep = "_")
    master[[col]] <- to_numeric(master[[col]])
  }
}

make_core_calls <- function(df, mode) {
  out <- data.frame(geneID = df$geneID, stringsAsFactors = FALSE)
  for (method in methods) {
    cols <- paste(method, accessions, sep = "_")
    mat <- as.matrix(df[, cols])
    if (mode == "strict_binary") {
      out[[method]] <- rowSums(mat == 1, na.rm = FALSE) == length(accessions)
    } else if (mode == "presence_any_count") {
      out[[method]] <- rowSums(mat > 0, na.rm = FALSE) == length(accessions)
    } else {
      stop("Unknown mode: ", mode)
    }
    out[[method]][is.na(out[[method]])] <- FALSE
  }
  out
}

plot_mode <- function(mode, mode_label) {
  calls <- make_core_calls(master, mode)
  calls$n_methods_calling_core <- rowSums(calls[, methods])

  by_n <- calls %>%
    count(n_methods_calling_core, name = "n_genes") %>%
    complete(n_methods_calling_core = 0:length(methods), fill = list(n_genes = 0))

  write.table(
    by_n,
    file = file.path(out_dir, paste0("core_call_overlap_by_number_of_methods_", mode, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  combination <- calls %>%
    mutate(
      Minigraph = minigraph,
      PGGB = pggb,
      OrthoFinder = orthofinder,
      combination = case_when(
        Minigraph & PGGB & OrthoFinder ~ "All three",
        Minigraph & PGGB & !OrthoFinder ~ "Minigraph + PGGB",
        Minigraph & !PGGB & OrthoFinder ~ "Minigraph + OrthoFinder",
        !Minigraph & PGGB & OrthoFinder ~ "PGGB + OrthoFinder",
        Minigraph & !PGGB & !OrthoFinder ~ "Minigraph only",
        !Minigraph & PGGB & !OrthoFinder ~ "PGGB only",
        !Minigraph & !PGGB & OrthoFinder ~ "OrthoFinder only",
        TRUE ~ "Not core by any method"
      )
    ) %>%
    count(Minigraph, PGGB, OrthoFinder, n_methods_calling_core, combination, name = "n_genes")

  combination_order <- c(
    "All three",
    "Minigraph + PGGB",
    "Minigraph + OrthoFinder",
    "PGGB + OrthoFinder",
    "Minigraph only",
    "PGGB only",
    "OrthoFinder only",
    "Not core by any method"
  )
  combination$combination <- factor(combination$combination, levels = rev(combination_order))

  write.table(
    combination %>% arrange(desc(n_methods_calling_core), desc(n_genes)),
    file = file.path(out_dir, paste0("core_call_overlap_combinations_", mode, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  p_by_n <- by_n %>%
    filter(n_methods_calling_core > 0) %>%
    mutate(n_methods_calling_core = factor(n_methods_calling_core, levels = 1:3)) %>%
    ggplot(aes(x = n_methods_calling_core, y = n_genes)) +
    geom_col(fill = "#2b8cbe", width = 0.68) +
    geom_text(aes(label = format(n_genes, big.mark = ",")), vjust = -0.35, size = 3.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    xlab("Number of methods calling gene as core") +
    ylab("Number of genes") +
    ggtitle(mode_label) +
    cowplot::theme_cowplot()

  p_combo <- combination %>%
    filter(n_methods_calling_core > 0) %>%
    ggplot(aes(x = combination, y = n_genes)) +
    geom_col(fill = "#41ab5d", width = 0.68) +
    geom_text(aes(label = format(n_genes, big.mark = ",")), hjust = -0.08, size = 3.4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip(clip = "off") +
    xlab(NULL) +
    ylab("Number of genes") +
    cowplot::theme_cowplot() +
    theme(plot.margin = margin(5.5, 28, 5.5, 5.5))

  combined <- cowplot::plot_grid(p_by_n, p_combo, labels = c("A", "B"), ncol = 1, rel_heights = c(1, 1.35))

  ggsave(
    file.path(out_dir, paste0("core_call_overlap_", mode, ".pdf")),
    combined,
    width = 8.8,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    file.path(out_dir, paste0("core_call_overlap_", mode, ".png")),
    combined,
    width = 8.8,
    height = 8,
    dpi = 220,
    bg = "white"
  )

  pairwise <- expand.grid(method_a = methods, method_b = methods, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      method_a_label = method_labels[[method_a]],
      method_b_label = method_labels[[method_b]],
      n_genes_both_core = sum(calls[[method_a]] & calls[[method_b]]),
      percent_total_genes = 100 * n_genes_both_core / nrow(calls)
    ) %>%
    ungroup() %>%
    mutate(
      method_a_label = factor(method_a_label, levels = method_labels[methods]),
      method_b_label = factor(method_b_label, levels = method_labels[methods]),
      label = sprintf("%s\n%.1f%%", format(n_genes_both_core, big.mark = ","), percent_total_genes)
    )

  write.table(
    pairwise,
    file = file.path(out_dir, paste0("core_call_pairwise_overlap_", mode, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  p_heat <- ggplot(pairwise, aes(x = method_b_label, y = method_a_label, fill = percent_total_genes)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 4.1, lineheight = 0.9) +
    scale_fill_gradient(
      low = "#f7fbff",
      high = "#08519c",
      limits = c(0, 100),
      name = "% of all genes"
    ) +
    coord_equal() +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(paste0("Pairwise overlap of core gene calls: ", mode_label)) +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.position = "right"
    )

  ggsave(
    file.path(out_dir, paste0("core_call_pairwise_overlap_heatmap_", mode, ".pdf")),
    p_heat,
    width = 6.5,
    height = 5.5,
    dpi = 300,
    bg = "white"
  )
}

plot_mode("strict_binary", "Strict binary core calls")
plot_mode("presence_any_count", "Presence-based core calls")

message("Wrote core-call overlap outputs to: ", out_dir)
