#!/bin/bash
#SBATCH --mem 2G --ntasks 1 --nodes 1 -p batch -J fetch --out logs/fetch.%A_%a.log

module load aspera
module load sratoolkit

ASCP=$(which ascp)
OUTDIR=fastq
mkdir -p /scratch/$USER/cache
if [ ! -e ~/ncbi ]; then
	ln -s /scratch/$USER/cache ~/ncbi
fi

mkdir -p $OUTDIR
N=1
if [ ! -z ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ ! -z $1 ]; then
 N=$1
fi

SRA=sra_samples.tsv
tail -n +2 $SRA | sed -n ${N}p | cut -f1 | while read SRARUN
do
# echo "$SRARUN fetching"
 if [ ! -f $OUTDIR/${SRARUN}_1.fastq.gz ]; then
# prefetch -a "$ASCP|$ASPERAKEY" --ascp-options "-k1 -Tr -l800m" $SRARUN
	echo "($N) $SRARUN"
	fastq-dump $SRARUN --gzip --split-files -O $OUTDIR

 fi
done
