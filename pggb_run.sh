# Combine genome assemblies from each accession into a single input assembly for PGGB:

cat *_Genome.updated.fasta > brassica_pggb_input.fasta

# Rename contig names to match PanSN convention:

sed 's/_/#1#/g' brassica_pggb_input.fasta > brassica_pggb_input.renamed.fasta

# Generate an assembly index:

samtools faidx brassica_pggb_input.renamed.fasta

# Run PGGB v0.7.5:

module load singularity/4.1.0-nompi
singularity run pggb_latest.sif pggb -i brassica_pggb_input.renamed.fasta  -s 10000 -p 95 -n 16 -t 96 -o brassica_pggb_out

Generate input BED files for ODGI pav using provided annotation GFF files from each accession:

for gff in *.gff; do
 acc=$(basename "$gff" .gff)

 awk -v acc="$acc" '
 BEGIN{OFS="\t"}

 $3=="gene" {

     len=$5-$4
     if(len < 100 || len > 200000) next

     if (match($9, /ID=([^;]+)/, a))
         id=a[1]
     else
         id="gene_" NR

     path=acc"#1#"$1

     print path, $4-1, $5, acc"|"id
 }' "$gff" > "${acc}.bed"

done

# Generate PAV tables from ODGI (.og) graph file and BED files using ODGI pav:

BEDS=(*.bed)
BED=${BEDS[$SLURM_ARRAY_TASK_ID]}

name=$(basename $BED .bed)

odgi pav -i brassica*.og -b $BED -B 0.8 -M -S --threads=108 > ${name}_genepav.tsv

Concatenate to generate gene PAV across all accessions:

cat *_genepav.tsv > all_genepav.tsv

Summarise statistics from all_genepav.tsv:

# Create summarise_presence_classes.py from nano:

import sys
import pandas as pd

## load table
pav = pd.read_csv(sys.argv[1], sep="\t")

## keep only accession columns
pav = pav.iloc[:, 4:]

## force numeric
pav = pav.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

# number of accessions present per gene
presence = pav.sum(axis=1)

total = len(presence)

summary = []

for n in sorted(presence.unique(), reverse=True):

   count = (presence == n).sum()

   pct = round(count / total * 100, 2)

   if n == pav.shape[1]:
       label = f"{n} accessions (core)"
   elif n == 1:
       label = f"{n} accession (unique)"
   else:
       label = f"{n} accessions"

   summary.append([label, count, pct])

summary_df = pd.DataFrame(
   summary,
   columns=["Category", "Genes", "Percent"]
)

print(summary_df.to_string(index=False))

# Run in login node:

python summarise_presence_classes.py all_genepav.tsv
