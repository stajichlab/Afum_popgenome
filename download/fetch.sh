#!/bin/bash
#SBATCH --mem 2G --ntasks 1 --nodes 1 -p batch -J prefetch  --out fetch.%A_%a.log

module load aspera
module load sratoolkit

ASCP=$(which ascp)
OUTDIR=fastq

mkdir -p $OUTDIR
N=1
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ "$1" != "" ]; then
 N=$1
fi

FILE=$(ls runs/*.tab | sed -n ${N}p)

RunCol=$(head -n 1 $FILE | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
         if($i ~ /Run/) { print i }
      }
}')

tail -n +2 $FILE | cut -f${RunCol} | while read SRARUN 
do
 echo "$SRARUN fetching"
 prefetch -a "$ASCP|$ASPERAKEY" --ascp-options "-k1 -Tr -l800m" $SRARUN
 fastq-dump  --gzip  --split-files --qual-filter-1 -O $outdir $SRARUN
done
