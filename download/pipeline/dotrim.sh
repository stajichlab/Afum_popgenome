#!/bin/bash
#SBATCH --ntasks 1 --nodes 1 --mem 4G --time 4:00:00 -J sickle --out sickle.%A_%a.log
module load sickle

INPUT=fastq
OUTPUT=trim
LEN=30
QUAL=20

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "need a num from --array or cmdline"
 exit
fi
FILE1=$(ls $INPUT/*_1.fastq.gz | sed -n ${N}p)
BASE=`basename $FILE1 _1.fastq.gz`

FILE2=$INPUT/$BASE"_2.fastq.gz"
OUTFILE1=$OUTPUT/$BASE.1.trim.fq.gz
OUTFILE2=$OUTPUT/$BASE.2.trim.fq.gz
OUTFILES=$OUTPUT/$BASE.s.trim.fq.gz

echo "$FILE1 $FILE2 -- $OUTFILE1 $OUTFILE2 $OUTFILES"
if [ ! -f $OUTFILE1 ]; then
# only run if we haven't run this before
 if [ -f $FILE2 ]; then
  # paired end data
echo  sickle pe -f $FILE1 -r $FILE2 -t sanger -l $LEN -q $QUAL \
  -o $OUTFILE1 -p $OUTFILE2 -s $OUTFILES
 else
 echo  sickle se -f $FILE1 -t sanger -l $LEN -q $QUAL -o $OUTFILES
 fi

fi
