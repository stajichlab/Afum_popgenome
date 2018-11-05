#!/bin/bash
#SBATCH --ntasks 4 --nodes 1 --mem 24G --time 4:00:00  -J trimmomatic --out logs/trimmomatic.%A_%a.log
CPU=2
module load trimmomatic

N=1
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ "$1" != "" ]; then
 N=$1
fi
INDIR=fastq
FILE=$(ls $INDIR/*_1.fastq.gz | sed -n ${N}p)
BASE=$(basename $FILE _1.fastq.gz)
OUT=trim
MINLEN=36
if [ ! -f $INDIR/${BASE}_2.fastq.gz ]; then 
 # SE
 if [ ! -f $OUT/${BASE}_1U.fq.gz ]; then
  java -jar $TRIMMOMATIC SE \
   -threads $CPU -phred33 -trimlog logs/$BASE.log $FILE  \
   $OUT/${BASE}_SE.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MAXINFO:50:0.8 MINLEN:$MINLEN
 fi
else 
 if [ ! -f $OUT/${BASE}_1P.fq.gz ]; then
 java -jar $TRIMMOMATIC PE \
  -threads $CPU -phred33 -trimlog logs/$BASE.log -validatePairs $FILE $INDIR/${BASE}_2.fastq.gz \
  -baseout $OUT/${BASE}.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN
 fi
fi
