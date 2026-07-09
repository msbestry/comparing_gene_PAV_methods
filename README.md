# comparing_gene_PAV_methods
Comparing orthogroup-based analysis and graph pangenome methods for determining gene PAV in crop species. Commands were run on setonix HPC using work partition unless stated otherwise.

## Obtain statistics of assemblies using assembly-stats v1.0.1:
assembly-stats *fasta > assembly_statistics.txt

## Alternatively, organise the statistics into a table:
assembly-stats -t *fasta > assembly_statistics.txt

## Overview of files in repository

build_gene_pav_master_table.py:	Combining gene PAV tables for master table
make_a55015_orthogroup_count_table.py:	Generating Orthofinder gene counts for genes present in a55015 accession
make_orthogroup_count_matrices.py:	Build Orthofinder gene count matrices
make_prefixed_orthogroup_pav.py:	Generate gene PAV table from Orthofinder gene counts
minigraph_pipeline.qmd:	Assemble pangenome graph with Minigraph and generate gene PAV table
orthofinder_run.sh:	Run Orthofinder for orthogroup-based analysis
pggb_run.sh:	Generate de novo pangenome graph with PGGB and perform ODGI post-processing to generate gene PAV table
plot_core_call_overlap.R:	Generate heatmap of core gene concordance between methods
plot_pav_method_concordance.R:	Generate heatmap of gene PAV concordance between methods
<img width="1394" height="262" alt="image" src="https://github.com/user-attachments/assets/f2811408-78a0-476d-bca4-5cc5352123e7" />
