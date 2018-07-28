#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 32G --time 1-0:0:0 -J realign --out logs/realign.%a.log

module unload java
module load java/8
module load gatk/3.7
module load picard
module load samtools

MEM=64g
GENOMEIDX=genome/Af293.fasta
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

 SAMPLEORIG=$STRAIN
 if [ -e $BAMDIR/${STRAIN}.PE.bam ]; then
  SAMPLE=${STRAIN}.PE
 elif [ -e $BAMDIR/${STRAIN}.SE.bam ]; then
  SAMPLE=${STRAIN}.SE
 else 
   echo "Cannot find SE or PE for $STRAIN"
   exit
 fi

 echo "SAMPLE=$SAMPLE"
 if [ ! -f $FINALBAMDIR/$SAMPLEORIG.bam ]; then
  if [ ! -f $BAMDIR/$SAMPLE.DD.bam ]; then
    time java -jar $PICARD MarkDuplicates I=$BAMDIR/$SAMPLE.bam O=$BAMDIR/$SAMPLE.DD.bam METRICS_FILE=logs/$SAMPLE.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  fi
  if [ ! -f $BAMDIR/$SAMPLE.DD.bai ]; then
    time java -jar $PICARD BuildBamIndex I=$BAMDIR/$SAMPLE.DD.bam TMP_DIR=/scratch
  fi
  if [ ! -f $BAMDIR/$SAMPLE.intervals ]; then 
   time java -Xmx$MEM -jar $GATK \
     -T RealignerTargetCreator \
     -R $GENOMEIDX \
     -I $BAMDIR/$SAMPLE.DD.bam \
     -o $BAMDIR/$SAMPLE.intervals
  fi
  time java -Xmx$MEM -jar $GATK \
      -T IndelRealigner \
      -R $GENOMEIDX \
      -I $BAMDIR/$SAMPLE.DD.bam \
      -targetIntervals $BAMDIR/$SAMPLE.intervals \
      -o $BAMDIR/$SAMPLE.realign.bam

  mv $BAMDIR/$SAMPLE.realign.bam $FINALBAMDIR/$SAMPLEORIG.bam
 fi
 if [ ! -f $FINALBAMDIR/$SAMPLEORIG.bai ]; then
    time java -jar $PICARD BuildBamIndex I=$FINALBAMDIR/$SAMPLEORIG.bam TMP_DIR=$TEMP
 fi
done
