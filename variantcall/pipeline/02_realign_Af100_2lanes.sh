#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 32G --time 1-0:0:0 -J realign --out logs/realign.%a.log

module unload java
module load java/8
module load gatk/3.7
module load picard
module load samtools

MEM=32g
GENOMEIDX=genome/Af293.fasta
TOPBAMDIR=aln
TEMP=/scratch

b=`basename $GENOMEIDX .fasta`
dir=`dirname $GENOMEIDX`
if [ ! -f $dir/$b.dict ]; then
 java -jar $PICARD CreateSequenceDictionary R=$GENOMEIDX O=$dir/$b.dict SPECIES="Aspergillus_fumigatus" TRUNCATE_NAMES_AT_WHITESPACE=true
fi

N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=Af100_samples.csv
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
sed -n ${N}p $SAMPFILE | while read SAMPLE FWD REV;
do
 SAMPLEORIG=$SAMPLE
 if [ $REV ]; then
  SAMPLE=$SAMPLE.PE
 else
  SAMPLE=$SAMPLE.SE
 fi

 hostname
 for LANE in lane1 lane2; 
 do
  BAMDIR=$TOPBAMDIR/$LANE
  echo "SAMPLE=$SAMPLE"
  if [ ! -f $BAMDIR/$SAMPLE.realign.bam ]; then
   if [ ! -f $BAMDIR/$SAMPLE.DD.bam ]; then
     time java -jar $PICARD MarkDuplicates I=$BAMDIR/$SAMPLE.bam O=$BAMDIR/$SAMPLE.DD.bam METRICS_FILE=logs/$SAMPLE.$LANE.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
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
   if [ ! -f $BAMDIR/$SAMPLE.realign.bam ]; then
    time java -Xmx$MEM -jar $GATK \
      -T IndelRealigner \
      -R $GENOMEIDX \
      -I $BAMDIR/$SAMPLE.DD.bam \
      -targetIntervals $BAMDIR/$SAMPLE.intervals \
      -o $BAMDIR/$SAMPLE.realign.bam
   fi
  fi
 done
 # This can be the merging step
 if [ ! -f $TOPBAMDIR/$SAMPLEORIG.bam ]; then
  samtools merge -O bam $TOPBAMDIR/$SAMPLEORIG.bam aln/lane*/$SAMPLE*realign.bam
 fi
 if [ ! -f $TOPBAMDIR/$SAMPLEORIG.bai ]; then
    time java -jar $PICARD BuildBamIndex I=$TOPBAMDIR/$SAMPLEORIG.bam TMP_DIR=$TEMP
 fi
done
