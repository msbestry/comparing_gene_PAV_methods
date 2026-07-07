# Combine genome assemblies from each accession into a single input assembly for PGGB:

cat *_Genome.updated.fasta > brassica_pggb_input.fasta

# Rename contig names to match PGGB convention:

sed 's/_/#1#/g' brassica_pggb_input.fasta > brassica_pggb_input.renamed.fasta

# Generate an assembly index:

samtools faidx brassica_pggb_input.renamed.fasta

# Run PGGB v0.7.5:

singularity run pggb_latest.sif pggb -i brassica_pggb_input.renamed.fasta  -s 10000 -p 95 -n 16 -t 96 -o brassica_pggb_out
