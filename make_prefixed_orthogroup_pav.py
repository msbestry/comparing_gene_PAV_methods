#!/usr/bin/env python3
"""Prefix OrthoFinder genes with individual names and build OG-level PAV."""

from __future__ import print_function

import argparse
import csv
import re
from collections import OrderedDict


MRNA_SUFFIX = re.compile(r"-mRNA-\d+$")


def split_genes(cell):
    if not cell or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def clean_gene_name(gene):
    return MRNA_SUFFIX.sub("", gene.strip())


def make_prefixed_gene(individual, gene):
    return "%s|%s" % (individual, clean_gene_name(gene))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite OrthoFinder gene names as individual|gene_without_mRNA_suffix "
            "and generate an orthogroup-level gene PAV matrix."
        )
    )
    parser.add_argument("orthogroups_tsv", help="Input Orthogroups.tsv")
    parser.add_argument(
        "--renamed-orthogroups",
        default="Orthogroups.individual_gene.tsv",
        help="Output Orthogroups table with renamed gene IDs",
    )
    parser.add_argument(
        "--pav",
        default="gene_presence_absence.by_orthogroup.individual_gene.tsv",
        help="Output OG-level PAV matrix with renamed gene IDs",
    )
    parser.add_argument(
        "--with-orthogroup",
        default="gene_presence_absence.by_orthogroup.individual_gene.with_orthogroup.tsv",
        help="Output one row per renamed gene and orthogroup",
    )
    parser.add_argument(
        "--mapping",
        default="individual_gene_to_orthogroup.tsv",
        help="Output renamed gene-to-orthogroup mapping",
    )
    parser.add_argument(
        "--summary",
        default="gene_presence_absence.by_orthogroup.individual_gene.summary.txt",
        help="Output summary text file",
    )
    args = parser.parse_args()

    rows_for_renamed_orthogroups = []
    gene_masks = OrderedDict()
    gene_ogs = OrderedDict()
    gene_og_rows = []
    duplicate_same_og = 0
    repeated_gene_new_og = 0
    gene_occurrences = 0
    orthogroups = 0

    with open(args.orthogroups_tsv, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise SystemExit("Input needs one Orthogroup column plus at least one individual column")

        individuals = header[1:]
        n_individuals = len(individuals)

        for row in reader:
            orthogroups += 1
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            orthogroup = row[0]
            cells = row[1 : 1 + n_individuals]
            og_mask = 0
            renamed_cells = []
            renamed_genes_in_row = []
            seen_in_row = set()

            for index, cell in enumerate(cells):
                individual = individuals[index]
                genes = split_genes(cell)
                if genes:
                    og_mask |= 1 << index

                renamed_genes = []
                for gene in genes:
                    renamed_gene = make_prefixed_gene(individual, gene)
                    renamed_genes.append(renamed_gene)
                    gene_occurrences += 1

                    if renamed_gene in seen_in_row:
                        duplicate_same_og += 1
                    else:
                        seen_in_row.add(renamed_gene)
                        renamed_genes_in_row.append(renamed_gene)

                renamed_cells.append(", ".join(renamed_genes))

            rows_for_renamed_orthogroups.append([orthogroup] + renamed_cells)

            for renamed_gene in renamed_genes_in_row:
                if renamed_gene in gene_masks and orthogroup not in gene_ogs[renamed_gene]:
                    repeated_gene_new_og += 1

                gene_masks[renamed_gene] = gene_masks.get(renamed_gene, 0) | og_mask
                gene_ogs.setdefault(renamed_gene, [])
                if orthogroup not in gene_ogs[renamed_gene]:
                    gene_ogs[renamed_gene].append(orthogroup)
                gene_og_rows.append((renamed_gene, orthogroup, og_mask))

    with open(args.renamed_orthogroups, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows_for_renamed_orthogroups)

    with open(args.pav, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene"] + individuals)
        for renamed_gene, mask in gene_masks.items():
            writer.writerow([renamed_gene] + [1 if mask & (1 << index) else 0 for index in range(n_individuals)])

    with open(args.with_orthogroup, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "orthogroup"] + individuals)
        for renamed_gene, orthogroup, mask in gene_og_rows:
            writer.writerow([renamed_gene, orthogroup] + [1 if mask & (1 << index) else 0 for index in range(n_individuals)])

    with open(args.mapping, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "orthogroup"])
        for renamed_gene, orthogroups_for_gene in gene_ogs.items():
            writer.writerow([renamed_gene, ",".join(orthogroups_for_gene)])

    multi_og_genes = sum(1 for orthogroups_for_gene in gene_ogs.values() if len(orthogroups_for_gene) > 1)

    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("input\t%s\n" % args.orthogroups_tsv)
        handle.write("renamed_orthogroups\t%s\n" % args.renamed_orthogroups)
        handle.write("pav\t%s\n" % args.pav)
        handle.write("with_orthogroup\t%s\n" % args.with_orthogroup)
        handle.write("mapping\t%s\n" % args.mapping)
        handle.write("presence_definition\torthogroup-level; every renamed gene in an OG receives the OG's individual presence pattern\n")
        handle.write("gene_name_definition\tindividual|gene_name, with trailing -mRNA-N removed\n")
        handle.write("orthogroups\t%d\n" % orthogroups)
        handle.write("individual_columns\t%d\n" % n_individuals)
        handle.write("individual_names\t%s\n" % "\t".join(individuals))
        handle.write("gene_occurrences\t%d\n" % gene_occurrences)
        handle.write("unique_renamed_gene_names\t%d\n" % len(gene_masks))
        handle.write("renamed_gene_names_seen_in_multiple_orthogroups\t%d\n" % multi_og_genes)
        handle.write("duplicate_renamed_gene_occurrences_within_same_orthogroup\t%d\n" % duplicate_same_og)
        handle.write("repeated_renamed_gene_new_orthogroup_events\t%d\n" % repeated_gene_new_og)

    print("Wrote %s" % args.renamed_orthogroups)
    print("Wrote %s" % args.pav)
    print("Wrote %s" % args.with_orthogroup)
    print("Wrote %s" % args.mapping)
    print("Wrote %s" % args.summary)
    print("Unique renamed gene names: %d" % len(gene_masks))
    print("Renamed gene names seen in multiple orthogroups: %d" % multi_og_genes)


if __name__ == "__main__":
    main()
