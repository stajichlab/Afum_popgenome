#!/bin/bash
#SBATCH -p short

OUTDIR=fastq
SRA=sra_samples.tsv

N=1
tail -n +2 $SRA | cut -f1 | while read SRARUN
do
 if [ ! -f $OUTDIR/${SRARUN}_1.fastq.gz ]; then
	echo "($N) $SRARUN"
 fi
 N=$(expr $N + 1)
done
