#!/usr/bin/bash
#SBATCH --mem 8G --ntasks 8 --nodes 1 -J bwa.Afum --out logs/Afum.bwa.%A_%a.log --time 8:00:00

module load bwa/0.7.15
module unload java
module load java/8
module load picard


GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input
TOPOUTDIR=aln
mkdir -p $TOPOUTDIR

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=SRA_samples.csv
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array or cmdline"
 exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
echo "$N $MAX for $SAMPFILE"
if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUN STRAIN SAMPLE CENTER EXP PROJ
do
  OUTDIR=$TOPOUTDIR
  PAIR1=$INDIR/${RUN}.1.trim.fq.gz
  PAIR2=$INDIR/${RUN}.2.trim.fq.gz
  SINGLE=$INDIR/${RUN}.s.trim.fq.gz
  
  SAMFILE=NULL
   
  if [ -e $PAIR1 ]; then      
   SAMFILE=$OUTDIR/${RUN}.PE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $SAMFILE ]; then
	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$SAMPLE\tLB:$RUN\tPL:illumina\tCN:$CENTER" $GENOME $PAIR1 $PAIR2 > $SAMFILE
   fi 
   if [ ! -f $OUTDIR/${RUN}.PE.bam ]; then
	samtools fixmate -O bam $SAMFILE $TEMP/${RUN}.fixmate.bam
	samtools sort -O bam -o  $OUTDIR/${RUN}.PE.bam -T $TEMP $TEMP/${RUN}.fixmate.bam
	/usr/bin/rm $TEMP/${RUN}.fixmate.bam
   fi
  elif [ -e $SINGLE ]; then
   SAMFILE=$OUTDIR/${RUN}.SE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $SAMFILE ]; then
    	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$SAMPLE\tLB:$RUN\tPL:illumina\tCN:$CENTER" $GENOME $SINGLE > $SAMFILE
   fi
   if [ ! -f $OUTDIR/${RUN}.SE.bam ]; then
	samtools view -b $SAMFILE > $TEMP/${RUN}.unsrt.bam	
	samtools sort -O bam -o $OUTDIR/${RUN}.SE.bam -T $TEMP $TEMP/${RUN}.unsrt.bam
	/usr/bin/rm $TEMP/${RUN}.unsrt.bam
   fi
  else
      echo "NO $PAIR1 and no $SINGLE?"
  fi
done
