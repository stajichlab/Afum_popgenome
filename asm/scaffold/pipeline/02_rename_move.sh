#!/usr/bin/bash
#SBATCH -p batch -N 1 -n 1
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
ASM=genomes
IN=scaffolds
OUT=genomes_scaffolded
mkdir -p $OUT
CTGS=$(ls $ASM/*.sorted.fasta | sed -n ${N}p)
NAME=$(basename $CTGS .sorted.fasta)
SCAFFOLDS=$IN/$NAME/ragoo_output/ragoo.fasta
module load AAFTF
AAFTF sort -i $SCAFFOLDS -o $OUT/$NAME.scaffolds.fasta

