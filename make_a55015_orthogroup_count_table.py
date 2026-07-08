#!/usr/bin/env python3
"""Create an a55015-anchored orthogroup count table from Orthogroups.tsv."""

from __future__ import print_function

import argparse
import csv
import re


MRNA_SUFFIX = re.compile(r"-mRNA-\d+$")


def split_genes(cell):
    if not cell or not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def clean_gene_name(gene):
    return MRNA_SUFFIX.sub("", gene.strip())


def prefixed_gene(individual, gene):
    return "%s|%s" % (individual, clean_gene_name(gene))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build an anchor-accession orthogroup count table. One row is one "
            "orthogroup containing at least one anchor gene; the second column "
            "lists all anchor genes in that orthogroup; subsequent columns are "
            "per-individual gene counts in that orthogroup."
        )
    )
    parser.add_argument("orthogroups_tsv", help="Input Orthogroups.tsv")
    parser.add_argument("--anchor", default="a55015", help="Anchor individual column")
    parser.add_argument(
        "-o",
        "--output",
        default="a55015_orthogroup_gene_counts.tsv",
        help="Output anchor orthogroup count table",
    )
    parser.add_argument(
        "--summary",
        default="a55015_orthogroup_gene_counts.summary.txt",
        help="Output summary text file",
    )
    args = parser.parse_args()

    orthogroups_total = 0
    anchor_orthogroups = 0
    anchor_gene_total = 0
    anchor_ogs_with_multiple_anchor_genes = 0
    max_anchor_genes_in_one_og = 0
    max_anchor_genes_og = ""
    max_count = 0
    max_count_og = ""
    max_count_individual = ""

    with open(args.orthogroups_tsv, "r", newline="", encoding="utf-8", errors="replace") as in_handle, \
            open(args.output, "w", newline="", encoding="utf-8") as out_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")

        header = next(reader)
        if len(header) < 2:
            raise SystemExit("Input needs one Orthogroup column plus individual columns")

        individuals = header[1:]
        if args.anchor not in individuals:
            raise SystemExit(
                "Anchor %s was not found in header: %s" % (args.anchor, ", ".join(individuals))
            )

        anchor_index = individuals.index(args.anchor)
        writer.writerow(["Orthogroup", "%s_genes" % args.anchor] + individuals)

        for row in reader:
            orthogroups_total += 1
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            orthogroup = row[0]
            cells = row[1 : 1 + len(individuals)]
            genes_by_individual = [split_genes(cell) for cell in cells]
            anchor_genes_raw = genes_by_individual[anchor_index]

            if not anchor_genes_raw:
                continue

            anchor_orthogroups += 1
            anchor_gene_total += len(anchor_genes_raw)
            if len(anchor_genes_raw) > 1:
                anchor_ogs_with_multiple_anchor_genes += 1
            if len(anchor_genes_raw) > max_anchor_genes_in_one_og:
                max_anchor_genes_in_one_og = len(anchor_genes_raw)
                max_anchor_genes_og = orthogroup

            counts = [len(genes) for genes in genes_by_individual]
            for individual, count in zip(individuals, counts):
                if count > max_count:
                    max_count = count
                    max_count_og = orthogroup
                    max_count_individual = individual

            anchor_genes = [prefixed_gene(args.anchor, gene) for gene in anchor_genes_raw]
            writer.writerow([orthogroup, ",".join(anchor_genes)] + counts)

    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("input\t%s\n" % args.orthogroups_tsv)
        handle.write("output\t%s\n" % args.output)
        handle.write("anchor\t%s\n" % args.anchor)
        handle.write("row_definition\torthogroups containing at least one anchor gene\n")
        handle.write("gene_column_definition\tall anchor genes in that orthogroup, formatted as individual|gene_name without trailing -mRNA-N\n")
        handle.write("count_definition\tnumber of genes from each individual in the orthogroup\n")
        handle.write("orthogroups_total\t%d\n" % orthogroups_total)
        handle.write("anchor_orthogroups\t%d\n" % anchor_orthogroups)
        handle.write("anchor_genes_total\t%d\n" % anchor_gene_total)
        handle.write("anchor_orthogroups_with_multiple_anchor_genes\t%d\n" % anchor_ogs_with_multiple_anchor_genes)
        handle.write("max_anchor_genes_in_one_orthogroup\t%d\n" % max_anchor_genes_in_one_og)
        handle.write("max_anchor_genes_orthogroup\t%s\n" % max_anchor_genes_og)
        handle.write("max_single_individual_count_in_anchor_orthogroups\t%d\n" % max_count)
        handle.write("max_single_individual_count_orthogroup\t%s\n" % max_count_og)
        handle.write("max_single_individual_count_individual\t%s\n" % max_count_individual)

    print("Wrote %s" % args.output)
    print("Wrote %s" % args.summary)
    print("Anchor orthogroups: %d" % anchor_orthogroups)
    print("Anchor genes total: %d" % anchor_gene_total)


if __name__ == "__main__":
    main()
