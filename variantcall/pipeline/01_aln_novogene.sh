#!/usr/bin/bash
#SBATCH --mem 8G --ntasks 8 --nodes 1 -J bwa.Afum --out logs/Afum.bwa_novogene.%A_%a.log --time 8:00:00

module load bwa/0.7.15
module unload java
module load java/8
module load picard
CENTER=Novogene
GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input/Novogene
TOPOUTDIR=aln
mkdir -p $TOPOUTDIR

TEMP=/scratch
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
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN PREFIX LEFT RIGHT
do
  OUTDIR=$TOPOUTDIR
  PAIR1=$INDIR/$LEFT
  PAIR2=$INDIR/$RIGHT
  #SINGLE=$INDIR/${PREFIX}.s.trim.fq.gz
  
  SAMFILE=NULL
   
  if [ -e $PAIR1 ]; then      
   SAMFILE=$OUTDIR/$STRAIN.PE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $SAMFILE ]; then
	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$PREFIX\tPL:illumina\tCN:$CENTER" $GENOME $PAIR1 $PAIR2 > $SAMFILE
   fi 
   if [ ! -f $OUTDIR/${STRAIN}.PE.bam ]; then
	samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${STRAIN}.fixmate.bam
	samtools sort --threads $CPU -O bam -o  $OUTDIR/${STRAIN}.PE.bam -T $TEMP $TEMP/${STRAIN}.fixmate.bam
	/usr/bin/rm $TEMP/${STRAIN}.fixmate.bam
        /usr/bin/rm $SAMFILE
   fi
 fi
done
