#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 32G --time 1-0:0:0 -J realign --out logs/realign_novogene.%a.log

module unload java
module load java/8
module load gatk/3.7
module load picard
module load samtools

MEM=64g
BAMDIR=aln
FINALBAMDIR=bam
TEMP=/scratch
mkdir -p $FINALBAMDIR
b=$(basename $GENOMEIDX .fasta)
dir=$(dirname $GENOMEIDX)
if [ ! -f $dir/$b.dict ]; then
 java -jar $PICARD CreateSequenceDictionary R=$GENOMEIDX O=$dir/$b.dict SPECIES="Aspergillus_fumigatus" TRUNCATE_NAMES_AT_WHITESPACE=true
fi

N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=Novogene_samples.csv
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN PREFIX LEFT RIGHT
do
 echo "SAMPLE=$SAMPLE"
done
