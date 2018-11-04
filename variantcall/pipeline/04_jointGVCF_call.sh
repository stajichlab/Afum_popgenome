#!/usr/bin/bash
#SBATCH --mem 32G --nodes 1 --ntasks 24 -J GATK.GVCFGeno --out logs/GVCFGeno.log --time 1-0:00:00

MEM=32g
module unload java
module load java/8
module load picard
module load gatk/3.8
module load tabix
module load parallel

GENOME=genome/Af293.fasta
FINALVCF=vcf
VARIANTFOLDER=gvcf
mkdir -p $FINALVCF
if [ -f config.txt ]; then
	source config.txt
fi
OUT=$FINALVCF/$PREFIX.all.vcf

CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
if [[ $(ls $VARIANTFOLDER/*.g.vcf | wc -l | awk '{print $1}') -gt "0" ]]; then
	parallel -j $CPU bgzip {} ::: $VARIANTFOLDER/*.g.vcf
	parallel -j $CPU tabix -f {} ::: $VARIANTFOLDER/*.g.vcf.gz
fi

N=$(ls $VARIANTFOLDER/*.g.vcf.gz | sort | perl -p -e 's/\n/ /; s/(\S+)/-V $1/')

if [ ! -f $OUT.gz ]; then
 if [ ! -f $OUT ]; then
	java -Xmx$MEM -jar $GATK \
    -T GenotypeGVCFs \
    -R $GENOME \
    $N \
    -o $OUT \
    -nt $CPU
  fi
  bgzip $OUT
  tabix $OUT
fi
