#!/bin/bash -l
#SBATCH --job-name=Bnapus_minigraph
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --account=pawsey0149
#SBATCH --time=1:00:00
#SBATCH --array=0-4
#SBATCH --output=/scratch/pawsey0149/falbornoz/brassica2/logs/Bnapus_minigraph%A_%a.out
#SBATCH --error=/scratch/pawsey0149/falbornoz/brassica2/logs/Bnapus_minigraph%A_%a.err
#SBATCH --partition=work

########CHANGE THESE TO FIT YOUR FILES AND FOLDERS
#your working directory
WORKDIR="$MYSCRATCH/brassica2"
#where are your assemblies in .fasta format
GENOMES_DIR="$WORKDIR/genomes"
#number of threads to use for commands
THREADS=20
#Species name. this is for naming of files
species="Bnapus"
#genome/assembly that will be used as reference in minigraph



####these ones below do not need to be edited, unless you want a different folder naming system. but they should be copy/pasted to certain steps. but just to be safe copy/paste them all at each step. it won't hurt.
#folder for some logs
LOGDIR="$WORKDIR/logs"
#folder for some QC stats and summaries
QCDIR="$WORKDIR/QC_stats"
#folder for .fasta assemblies after renaming contigs
GENOMES_RENAMED_DIR="$GENOMES_DIR/renamed"
#folder for .fasta assemblies after renaming contigs back to their original
GENOMES_RENAMED_AGAIN_DIR="$GENOMES_RENAMED_DIR/again"
#folder for the outputs of minigraph rounds as .gfa and after transforming to .bed file
MINIDIR="$WORKDIR/minigraph_output"
#folder for the outputs of .gfa graph after sorting it
MINIDIR_SORTED="$MINIDIR/sorted"
#folder for .gaf file output of mapping of each genome agains the .sorted.bed graph file
MAPDIR="$WORKDIR/graph_map"
#folder for the transformed .gaf files into .bed files
MAPDIR_BED="$MAPDIR/gaf_to_bed"
#folder for after sorting the .bed files for each .gaf mapping
MAPDIR_BED_SORTED="$MAPDIR_BED/sorted"
#folder for the PAV matrix
PAV_DIR="$WORKDIR/PAV_OUT"
#folder for the .gff annotations of each genome.
GFF_DIR="$WORKDIR/genome_annotation"
REF="$GENOMES_RENAMED_DIR/a55015_Genome_renamed.fasta"
#set working directories
mkdir -p "$MAPDIR"

ASMS=("$GENOMES_RENAMED_AGAIN_DIR"/*_renamed_again.fasta)
graphGFA=$(<"$MINIDIR/${species}_graphGFA_file_name.txt")

#get the final round, simply for file name calling
round=$(basename "$graphGFA" .gfa)
round=${round#${species}_}

#extract number of genomes for the array
NUM_ASMS=${#ASMS[@]}
echo "Found $NUM_ASMS genome(s) in $GENOMES_RENAMED_AGAIN_DIR"

# check array index is valid
if [ -z "${SLURM_ARRAY_TASK_ID+x}" ]; then
  echo "SLURM_ARRAY_TASK_ID is not set" >&2
  exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 0 )) || (( SLURM_ARRAY_TASK_ID >= NUM_ASMS )); then
  echo "Array index ${SLURM_ARRAY_TASK_ID} out of range (0..$((NUM_ASMS-1)))" >&2
  exit 2
fi

#create base names for each genome
INPATH="${ASMS[$SLURM_ARRAY_TASK_ID]}"
BASENAME="$(basename "$INPATH")"
NAME="${BASENAME%.*}"        

echo "Task $SLURM_ARRAY_TASK_ID -> input: $INPATH"
echo "Output will be: $MAPDIR/${NAME}_vs_${species}_${round}.gaf"

# grabs our filename from a directory listing and run minigraph of each genome against the graph
minigraph -cxasm -t "${SLURM_CPUS_PER_TASK}" -f.1 "$graphGFA" "$INPATH" > "${MAPDIR}/${NAME}_vs_${species}_${round}.gaf"
