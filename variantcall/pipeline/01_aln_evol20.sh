#!/usr/bin/bash
#SBATCH --mem 8G --ntasks 8 --nodes 1 -J evol20.Afum --out logs/Afum_Evol20.bwa.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard

CENTER=UCR

GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input/UCR_FC548
TOPOUTDIR=aln
mkdir -p $TOPOUTDIR

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=UCR_FC548_samples.csv
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
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX STRAIN LEFT RIGHT BARCODE DESC
do
  OUTDIR=$TOPOUTDIR
  PAIR1=$INDIR/$LEFT
  PAIR2=$INDIR/$RIGHT
  
  SAMFILE=NULL
   
  if [ -e $PAIR1 ]; then      
   SAMFILE=$OUTDIR/${PREFIX}.PE.unsrt.sam
   echo "SAMFILE is $SAMFILE"
   if [ ! -f $SAMFILE ]; then
	bwa mem -t $CPU -R "@RG\tID:$PREFIX\tSM:$PREFIX\tLB:$PREFIX\tPL:illumina\tCN:$CENTER" $GENOME $PAIR1 $PAIR2 > $SAMFILE
   fi 
   if [ ! -f $OUTDIR/${PREFIX}.PE.bam ]; then
	samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${PREFIX}.fixmate.bam
	samtools sort --threads $CPU -O bam -o  $OUTDIR/${PREFIX}.PE.bam -T $TEMP $TEMP/${PREFIX}.fixmate.bam
	/usr/bin/rm $TEMP/${PREFIX}.fixmate.bam
        /usr/bin/rm $SAMFILE
   fi
 fi
done
