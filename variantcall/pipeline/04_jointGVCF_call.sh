#!/usr/bin/bash
#SBATCH --mem 32G --nodes 1 --ntasks 24 -J GATK.GVCFGeno --out logs/GVCFGeno.log --time 1-0:00:00

MEM=32g
module unload java
module load java/8
module load picard
module load gatk/3.8

GENOME=genome/Af293.fasta
OUTDIR=Variants
mkdir -p $OUTDIR
INDIR=gvcf
OUT=$OUTDIR/A_fumigiatus_Af293.Popgen2.vcf

CPU=$SLURM_CPUS_ON_NODE

if [ ! $CPU ]; then
 CPU=2
fi

N=$(ls $INDIR/*.g.vcf | sort | perl -p -e 's/\n/ /; s/(\S+)/-V $1/')
if [ ! -f $OUT ]; then
java -Xmx$MEM -jar $GATK \
    -T GenotypeGVCFs \
    -R $GENOME \
    $N \
    -o $OUT \
    -nt $CPU
fi

#OUT=$OUTDIR/A_fumigiatus_Af293.EVOL_Obar.vcf
#N=$(ls $INDIR/Obar-*.g.vcf | sort | perl -p -e 's/\n/ /; s/(\S+)/-V $1/')

#if [ ! -f $OUT ]; then
#java -Xmx$MEM -jar $GATK \
#    -T GenotypeGVCFs \
#    -R $GENOME \
#    $N \
#    -o $OUT \
#    -nt $CPU
#fi


