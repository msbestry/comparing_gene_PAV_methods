#!/bin/bash --login

#SBATCH --job-name=orthofinder
#SBATCH --time=72:00:00
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=8G
#SBATCH --output=out.log
#SBATCH --error=err.log
#SBATCH --export=NONE
#SBATCH --mail-user=mitchell.bestry@uwa.edu.au  
#SBATCH --mail-type=BEGIN,END

# This slurm job was run on the UWA Kaya HPC, as it was unsuitable to run on setonix
# All five protein assembly fasta files were inside the ./input_fast_protein directory

orthofinder -f ./input_fasta_protein -t 48 -a 12
