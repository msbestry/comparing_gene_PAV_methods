#!/usr/bin/env python3
"""Build OrthoFinder orthogroup gene-count matrices."""

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


def renamed_gene(individual, gene):
    return "%s|%s" % (individual, clean_gene_name(gene))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create count matrices from Orthogroups.tsv. Values are the number "
            "of genes from each accession in the orthogroup, not binary 0/1."
        )
    )
    parser.add_argument("orthogroups_tsv", help="Input Orthogroups.tsv")
    parser.add_argument(
        "--anchor",
        default="a55015",
        help="Anchor accession for the anchor-gene count table",
    )
    parser.add_argument(
        "--orthogroup-counts",
        default="orthogroup_gene_counts.tsv",
        help="Output OG x accession gene-count matrix",
    )
    parser.add_argument(
        "--anchor-counts",
        default="a55015_anchor_gene_counts.by_orthogroup.tsv",
        help="Output anchor-gene x accession gene-count matrix",
    )
    parser.add_argument(
        "--all-gene-counts",
        default="gene_copy_number.by_orthogroup.individual_gene.tsv",
        help="Output all-renamed-gene x accession gene-count matrix",
    )
    parser.add_argument(
        "--mapping",
        default="individual_gene_to_orthogroup.count_matrix.tsv",
        help="Output renamed gene-to-orthogroup mapping",
    )
    parser.add_argument(
        "--summary",
        default="orthogroup_gene_counts.summary.txt",
        help="Output summary text file",
    )
    args = parser.parse_args()

    og_count_rows = []
    anchor_rows = []
    all_gene_rows = []
    mapping_rows = []
    unique_renamed_genes = set()

    orthogroups = 0
    gene_occurrences = 0
    ogs_with_any_count_gt1 = 0
    ogs_with_any_count_gt5 = 0
    ogs_with_any_count_gt10 = 0
    max_count = 0
    max_count_og = ""
    max_count_individual = ""
    per_individual_total_genes = OrderedDict()

    with open(args.orthogroups_tsv, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise SystemExit("Input needs one Orthogroup column plus at least one individual column")

        individuals = header[1:]
        if args.anchor not in individuals:
            raise SystemExit("Anchor %s was not found in header: %s" % (args.anchor, ", ".join(individuals)))

        anchor_index = individuals.index(args.anchor)
        for individual in individuals:
            per_individual_total_genes[individual] = 0

        for row in reader:
            orthogroups += 1
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            orthogroup = row[0]
            cells = row[1 : 1 + len(individuals)]
            cell_genes = []
            count_row = []

            for index, cell in enumerate(cells):
                individual = individuals[index]
                genes = split_genes(cell)
                cell_genes.append(genes)
                count = len(genes)
                count_row.append(count)
                gene_occurrences += count
                per_individual_total_genes[individual] += count

                if count > max_count:
                    max_count = count
                    max_count_og = orthogroup
                    max_count_individual = individual

            if any(count > 1 for count in count_row):
                ogs_with_any_count_gt1 += 1
            if any(count > 5 for count in count_row):
                ogs_with_any_count_gt5 += 1
            if any(count > 10 for count in count_row):
                ogs_with_any_count_gt10 += 1

            og_count_rows.append([orthogroup] + count_row)

            for index, genes in enumerate(cell_genes):
                individual = individuals[index]
                for gene in genes:
                    rgene = renamed_gene(individual, gene)
                    unique_renamed_genes.add(rgene)
                    mapping_rows.append([rgene, orthogroup])
                    all_gene_rows.append([rgene] + count_row)

            for gene in cell_genes[anchor_index]:
                rgene = renamed_gene(args.anchor, gene)
                anchor_rows.append([rgene, orthogroup] + count_row)

    with open(args.orthogroup_counts, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup"] + individuals)
        writer.writerows(og_count_rows)

    with open(args.anchor_counts, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "orthogroup"] + individuals)
        writer.writerows(anchor_rows)

    with open(args.all_gene_counts, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene"] + individuals)
        writer.writerows(all_gene_rows)

    with open(args.mapping, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "orthogroup"])
        writer.writerows(mapping_rows)

    with open(args.summary, "w", encoding="utf-8") as handle:
        handle.write("input\t%s\n" % args.orthogroups_tsv)
        handle.write("orthogroup_counts\t%s\n" % args.orthogroup_counts)
        handle.write("anchor_counts\t%s\n" % args.anchor_counts)
        handle.write("all_gene_counts\t%s\n" % args.all_gene_counts)
        handle.write("mapping\t%s\n" % args.mapping)
        handle.write("value_definition\tnumber of genes from each accession in the orthogroup\n")
        handle.write("gene_name_definition\tindividual|gene_name, with trailing -mRNA-N removed\n")
        handle.write("anchor\t%s\n" % args.anchor)
        handle.write("orthogroups\t%d\n" % orthogroups)
        handle.write("individual_columns\t%d\n" % len(individuals))
        handle.write("individual_names\t%s\n" % "\t".join(individuals))
        handle.write("gene_occurrences\t%d\n" % gene_occurrences)
        handle.write("unique_renamed_gene_names\t%d\n" % len(unique_renamed_genes))
        handle.write("anchor_gene_rows\t%d\n" % len(anchor_rows))
        handle.write("all_gene_rows\t%d\n" % len(all_gene_rows))
        handle.write("orthogroups_with_any_accession_count_gt1\t%d\n" % ogs_with_any_count_gt1)
        handle.write("orthogroups_with_any_accession_count_gt5\t%d\n" % ogs_with_any_count_gt5)
        handle.write("orthogroups_with_any_accession_count_gt10\t%d\n" % ogs_with_any_count_gt10)
        handle.write("max_single_accession_count\t%d\n" % max_count)
        handle.write("max_single_accession_count_orthogroup\t%s\n" % max_count_og)
        handle.write("max_single_accession_count_individual\t%s\n" % max_count_individual)
        for individual, total in per_individual_total_genes.items():
            handle.write("total_genes_%s\t%d\n" % (individual, total))

    print("Wrote %s" % args.orthogroup_counts)
    print("Wrote %s" % args.anchor_counts)
    print("Wrote %s" % args.all_gene_counts)
    print("Wrote %s" % args.mapping)
    print("Wrote %s" % args.summary)
    print("Orthogroups with any accession count >1: %d" % ogs_with_any_count_gt1)


if __name__ == "__main__":
    main()
