#!/usr/bin/bash
#SBATCH --nodes 1 -p batch
#SBATCH --ntasks 4
#SBATCH --mem 16G
#SBATCH --job-name=GATK.select_filter
#SBATCH --output=logs/GATK.select_filter.log

module load gatk/4
module unload java
module load java/8
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
echo "$PREFIX $INFILE"
INFILE=$FINALVCF/$PREFIX.all.vcf.gz
INSNP=$FINALVCF/$PREFIX.SNP.vcf
ININDEL=$FINALVCF/$PREFIX.INDEL.vcf
FILTEREDSNP=$FINALVCF/$PREFIX.filtered.SNP.vcf
FILTEREDINDEL=$FINALVCF/$PREFIX.filtered.INDEL.vcf
SNPONLY=$FINALVCF/$PREFIX.selected.SNP.vcf
INDELONLY=$FINALVCF/$PREFIX.selected.INDEL.vcf

FILTEREDFIXEDSNP=$FINALVCF/$PREFIX.filteredfixed.SNP.vcf
FILTEREDFIXEDINDEL=$FINALVCF/$PREFIX.filteredfixed.INDEL.vcf
SNPNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.SNP.vcf
INDELNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.INDEL.vcf

if [[ ! -f $INSNP.gz || $INFILE -nt $INSNP.gz ]]; then
    gatk SelectVariants \
	-R $REFGENOME \
	--variant $INFILE \
	-O $INSNP \
	--restrict-alleles-to BIALLELIC \
	--select-type-to-include SNP --create-output-variant-index false
    bgzip $INSNP
    tabix $INSNP.gz
    unlink $INSNP.idx
fi

if [[ ! -f $ININDEL.gz || $INFILE -nt $ININDEL.gz ]]; then
    gatk SelectVariants \
	-R $REFGENOME \
	--variant $INFILE \
	--output $ININDEL \
	--select-type-to-include INDEL --select-type-to-include MIXED \
	--select-type-to-include MNP --create-output-variant-index false
    bgzip $ININDEL
    tabix $ININDEL.gz
fi
    
if [[ ! -f $FILTEREDSNP.gz || $INSNP.gz -nt $FILTEREDSNP.gz ]]; then
    gatk VariantFiltration --output $FILTEREDSNP \
	--variant $INSNP.gz -R $REFGENOME \
	--cluster-window-size 10  \
	--filter-expression "QD < 2.0" --filter-name QualByDepth \
	--filter-expression "MQ < 40.0" --filter-name MapQual \
	--filter-expression "QUAL < 100" --filter-name QScore \
	--filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
	--filter-expression "SOR > 4.0" --filter-name StrandOddsRatio \
	--filter-expression "FS > 60.0" --filter-name FisherStrandBias \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \
	--missing-values-evaluate-as-failing --create-output-variant-index false
    bgzip $FILTEREDSNP
    tabix $FILTEREDSNP.gz
fi

if [[ ! -f $FILTEREDFIXEDSNP.gz || $INSNP.gz -nt $FILTEREDFIXEDSNP.gz ]]; then
    gatk VariantFiltration --output $FILTEREDFIXEDSNP \
	--variant $INSNP.gz -R $REFGENOME \
	--cluster-window-size 10  \
	--filter-expression "QD < 2.0" --filter-name QualByDepth \
	--filter-expression "MQ < 40.0" --filter-name MapQual \
	--filter-expression "QUAL < 100" --filter-name QScore \
	--filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
	--filter-expression "SOR > 4.0" --filter-name StrandOddsRatio \
	--filter-expression "FS > 60.0" --filter-name FisherStrandBias \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \
	--filter-expression "AF > 0.99" --filter-name FixedAllele \
	--missing-values-evaluate-as-failing --create-output-variant-index false
    bgzip $FILTEREDFIXEDSNP
    tabix $FILTEREDFIXEDSNP.gz
fi

if [[ ! -f $FILTEREDINDEL.gz || $FILTEREDINDEL.gz -nt $ININDEL.gz ]]; then
    gatk VariantFiltration --output $FILTEREDINDEL \
	--variant $ININDEL.gz -R $REFGENOME \
	--cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
	-filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \
	-filter "SOR > 10.0" --filter-name StrandOddsRatio \
	-filter "FS > 200.0" --filter-name FisherStrandBias \
	-filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \
	-filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank \
	--create-output-variant-index false
    bgzip $FILTEREDINDEL
    tabix $FILTEREDINDEL.gz
fi

if [[ ! -f $FILTEREDFIXEDINDEL.gz || $INDEL.gz -nt $FILTEREDFIXEDINDEL.gz ]]; then
    gatk VariantFiltration --output $FILTEREDFIXEDINDEL \
	--variant $ININDEL.gz -R $REFGENOME \
	--cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
	-filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \
	-filter "SOR > 10.0" --filter-name StrandOddsRatio \
	-filter "FS > 200.0" --filter-name FisherStrandBias \
	-filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \
	-filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank \
	--filter-expression "AF > 0.99" --filter-name FixedAllele \
	--create-output-variant-index false
    bgzip $FILTEREDFIXEDINDEL
    tabix $FILTEREDFIXEDINDEL.gz
fi

if [[ ! -f $SNPONLY.gz || $FILTEREDSNP -nt $SNPONLY.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTEREDSNP \
	--output $SNPONLY \
	--exclude-filtered --create-output-variant-index false
    bgzip $SNPONLY
    tabix $SNPONLY.gz
fi

if [[ ! -f $INDELONLY.gz || $FILTEREDINDEL -nt $INDELONLY.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTEREDINDEL \
	--output $INDELONLY \
	--exclude-filtered --create-output-variant-index false
    bgzip $INDELONLY
    tabix $INDELONLY.gz
fi

if [[ ! -f $SNPNOFIXED.gz || $FILTEREDFIXEDSNP -nt $SNPNOFIXED.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTEREDFIXEDSNP \
	--output $SNPNOFIXED \
	--exclude-filtered --create-output-variant-index false
    bgzip $SNPNOFIXED
    tabix $SNPNOFIXED.gz
fi

if [[ ! -f $INDELNOFIXED.gz || $FILTEREDFIXEDINDEL -nt $INDELNOFIXED.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTEREDFIXEDINDEL \
	--output $INDELNOFIXED \
	--exclude-filtered --create-output-variant-index false
    bgzip $INDELNOFIXED
    tabix $INDELNOFIXED.gz
fi
