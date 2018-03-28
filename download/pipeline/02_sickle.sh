#!/bin/bash
#SBATCH --ntasks 1 --nodes 1 --mem 4G --time 2:00:00 -J sickle --out logs/sickle.%a.log -p short
module load sickle

INPUT=fastq
OUTPUT=trim
LEN=30
QUAL=20
if [ ! $TEMPDIR ]; then
 TEMPDIR=/tmp
fi
 
mydir=$(mktemp -dt XXXXXXXXXXXX)
if [ ! -d $mydir ]; then 
 echo "cannot create tempdir"
 exit
fi

function finish {
  rm -rf "$mydir"
}
trap finish EXIT

ALLACC=$mydir/all_SRA_acc.tab
rm -f $ALLACC
touch $ALLACC

for f in runs/*.tab
do
 RunCol=$(head -n 1 $f | awk -F'\t' ' {
      for(i=1;i < NF;i++) {
         if($i ~ /Run/) { print i }
      }
}')
 tail -n +2 $f | cut -f${RunCol} | while read SRARUN
 do
  echo $SRARUN >> $ALLACC
 done
done
acccount=$(wc -l $ALLACC | awk '{print $1}')


N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "need a num from --array or cmdline"
 exit
fi
BASE=$(sed -n ${N}p $ALLACC)

FILE1=$INPUT/$BASE"_1.fastq.gz"
FILE2=$INPUT/$BASE"_2.fastq.gz"
OUTFILE1=$OUTPUT/$BASE.1.trim.fq
OUTFILE2=$OUTPUT/$BASE.2.trim.fq
OUTFILES=$OUTPUT/$BASE.s.trim.fq

echo "$FILE1 $FILE2 -- $OUTFILE1 $OUTFILE2 $OUTFILES"
if [ ! -f $OUTFILE1.gz ]; then
# only run if we haven't run this before
 if [ -f $FILE2 ]; then
  # paired end data
  sickle pe -f $FILE1 -r $FILE2 -t sanger -l $LEN -q $QUAL \
  -o $OUTFILE1 -p $OUTFILE2 -s $OUTFILES
  pigz $OUTFILE1 $OUTFILE2 $OUTFILES
 elif [ ! -f $OUTFILES ]; then
  sickle se -f $FILE1 -t sanger -l $LEN -q $QUAL -o $OUTFILES
  pigz $OUTFILES
 fi
fi

