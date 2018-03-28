#!/usr/bin/bash
#SBATCH -J GATK.HTC --out logs/GATK_HTC.%a.log --ntasks 8 --nodes 1 --mem 16G

module unload java
module load java/8
module load gatk/3.8
module load picard

MEM=32g
GENOMEIDX=genome/Af293.fasta
BAMDIR=aln
OUTDIR=gvcf

mkdir -p $OUTDIR

b=$(basename $GENOMEIDX .fasta)
dir=$(dirname $GENOMEIDX)

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

MAX=$(wc -l $SAMPFILE | awk '{print $1}')

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
sed -n ${N}p $SAMPFILE | while read SAMPLE FWD REV;
do
 hostname
 echo "SAMPLE=$SAMPLE"
 BAMFILE=$BAMDIR/$SAMPLE.bam

 if [ ! -f $OUTDIR/$SAMPLE.g.vcf ]; then
  java -Xmx${MEM} -jar $GATK \
  -T HaplotypeCaller \
  -ERC GVCF \
  -ploidy 1 \
  -I $BAMFILE -R $GENOMEIDX \
  -forceActive -disableOptimizations \
  -o $OUTDIR/$SAMPLE.g.vcf -nct $CPU
 fi

done
