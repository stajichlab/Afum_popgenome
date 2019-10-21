#!/usr/bin/bash
#SBATCH --nodes 1 -p batch
#SBATCH --ntasks 8
#SBATCH --mem 32G
#SBATCH --job-name=GATK.gatk4.select_filter
#SBATCH --output=logs/GATK.gatk4_bcftools.filter.%A.log

module load gatk/4
module unload java
module load java/8
module load bcftools
hostname


CONFIG=config.txt
if [ -f $CONFIG ]; then
    source $CONFIG
fi
GENOME=$GENOMEFOLDER/$GENOMEFASTA
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
    CPU=1
fi

if [ -z $FINALVCF ]; then
	echo "Need FINALVCF for output"
	exit
fi
TMP=/scratch/$USER.$$.snpfilter
mkdir -p $TMP
echo "$PREFIX $INFILE"
INFILE=$FINALVCF/$PREFIX.gatk4.all.vcf.gz
PREFIX=$PREFIX.gatk4_bcftools_2

SNPONLY=$FINALVCF/$PREFIX.selected.SNP.vcf.gz
INDELONLY=$FINALVCF/$PREFIX.selected.INDEL.vcf.gz

SNPFILONLY=$FINALVCF/$PREFIX.gatk_filtered.SNP.vcf
INDELFILONLY=$FINALVCF/$PREFIX.gatk_filtered.INDEL.vcf

SNPSELONLY=$FINALVCF/$PREFIX.gatk_selected.SNP.vcf
INDELSELONLY=$FINALVCF/$PREFIX.gatk_selected.INDEL.vcf

if [ ! -f $SNPONLY ]; then
    echo " running bcftools filter to SNPs"
   time bcftools view --threads 8 --max-alleles 2 --max-af 0.99 -v snps -e 'QD<2.0 || MQ<40.0 || MQRankSum < -12.5 || SOR > 3.0 || FS > 60.0 || ReadPosRankSum < -8.0' -Oz -o $SNPONLY $INFILE
   tabix $SNPONLY
fi

if [ ! -f $SNPFILONLY.gz ]; then
    rsync -av $SNPONLY $TMP
    tabix $TMP/$(basename $SNPONLY)
   gatk VariantFiltration --output $SNPFILONLY \
   --variant $TMP/$(basename $SNPONLY) -R $REFGENOME \
   --cluster-window-size 10  \
   --filter-expression "QD < 2.0" --filter-name QualByDepth \
   --filter-expression "MQ < 40.0" --filter-name MapQual \
   --filter-expression "QUAL < 100" --filter-name QScore \
   --filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
   --filter-expression "SOR > 3.0" --filter-name StrandOddsRatio \
   --filter-expression "FS > 60.0" --filter-name FisherStrandBias \
   --filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \
   --missing-values-evaluate-as-failing
   bgzip $SNPFILONLY
   tabix $SNPFILONLY.gz
fi

if [ ! -f $SNPSELONLY.gz ]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $SNPFILONLY \
	--output $SNPSELONLY \
	--exclude-filtered
    bgzip $SNPSELONLY
    tabix $SNPSELONLY.gz
fi

rm -rf $TMP
exit
if [ ! -f $INDELONLY ]; then
    time bcftools view --threads 8 --max-af 0.99 --exclude-types snps -e 'QD<2.0 || MQ<40.0 || MQRankSum < -12.5 || SOR > 3.0 || FS > 60.0 || ReadPosRankSum < -8.0' -Oz -o $INDELONLY $INFILE
    tabix $INDELONLY
fi

if [ ! -f $INDELFILONLY.gz ]; then
    gatk VariantFiltration --output $INDELFILONLY \
	--variant $INDELONLY -R $REFGENOME \
	--cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
	-filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \
	-filter "SOR > 4.0" --filter-name StrandOddsRatio \
	-filter "FS > 200.0" --filter-name FisherStrandBias \
	-filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \
	-filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank 
    
    bgzip $INDELFILONLY
    tabix $INDELFILONLY.gz
fi

if [ ! -f $INDELSELONLY.gz ]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $INDELFILONLY \
	--output $INDELSELONLY \
	--exclude-filtered
    bgzip $INDELSELONLY
    tabix $INDELSELONLY.gz
fi


