#!/usr/bin/bash
#SBATCH --mem 8G --ntasks 8 --nodes 1 -J bwa.Afum
#SBATCH --out logs/Afum_sra.bwa.%A_%a.log --time 8:00:00

module load bwa/0.7.17

GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input/SRA
TOPOUTDIR=aln
mkdir -p $TOPOUTDIR

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=SRA_samples.csv
if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then 
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
  PAIR1=$INDIR/${RUN}_1.fastq.gz
  PAIR2=$INDIR/${RUN}_2.fastq.gz
  
  SAMFILE=NULL
   
  if [ -e $PAIR1 ]; then      
   SAMFILE=$OUTDIR/${RUN}.PE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $SAMFILE ]; then
	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$SAMPLE\tLB:$RUN\tPL:illumina\tCN:$CENTER" $GENOME $PAIR1 $PAIR2 > $SAMFILE
   fi 
   if [ ! -f $OUTDIR/${RUN}.PE.bam ]; then
	samtools fixmate -@ $CPU -O bam $SAMFILE $TEMP/${RUN}.fixmate.bam
	samtools sort -@ $CPU -O bam -o  $OUTDIR/${RUN}.PE.bam -T $TEMP $TEMP/${RUN}.fixmate.bam
	unlink $TEMP/${RUN}.fixmate.bam $SAMFILE
   fi
  elif [ -e $SINGLE ]; then
   SAMFILE=$OUTDIR/${RUN}.SE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $OUTDIR/${RUN}.SE.bam ]; then
    	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$SAMPLE\tLB:$RUN\tPL:illumina\tCN:$CENTER" $GENOME $SINGLE | samtools sort -@ $CPU -O bam -T $TEMP -o $OUTDIR/${RUN}.SE.bam -
   fi
  else
      echo "NO $PAIR1 and no $SINGLE?"
  fi
done
