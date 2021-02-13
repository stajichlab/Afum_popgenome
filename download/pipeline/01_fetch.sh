#!/bin/bash
#SBATCH --mem 16G --ntasks 8 --nodes 1 -p batch -J fetch --out logs/fetch.%A_%a.log
CPU=8
module load aspera
module unload perl
module load miniconda3
source activate mosdepth
module load sratoolkit
TMP=/scratch
ASCP=$(which ascp)
OUTDIR=fastq
#mkdir -p /scratch/$USER/cache
#if [ ! -e ~/ncbi ]; then
#	ln -s /scratch/$USER/cache ~/ncbi
#fi

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
#	prefetch -a "$ASCP|$ASPERAKEY" --ascp-options "-k1 -Tr -l800m" $SRARUN
	echo $SRARUN > /tmp/run.$$
	prefetch -t fasp --ascp-path "$ASCP|$ASPERAKEY" --option-file=/tmp/run.$$
	unlink /tmp/run.$$
	echo "($N) $SRARUN"
#	fastq-dump $SRARUN --gzip --split-files -O $OUTDIR
	parallel-fastq-dump --tmpdir $TMP --gzip  --sra-id $SRARUN --threads $CPU -O $OUTDIR/$SRARUN --split-files
 fi
done

