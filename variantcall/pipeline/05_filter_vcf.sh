#!/usr/bin/bash
#SBATCH --nodes 1 -p short
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=GATK.select_filter
#SBATCH --output=logs/GATK.select_filter.%A.log

module load gatk/3.8
module unload java
module load java/8
INDIR=Variants
OUTDIR=Variants_filter
GENOME=genome/Af293.fasta
mkdir -p $OUTDIR
for PREFIX in A_fumigiatus_Af293.Popgen1 A_fumigiatus_Af293.Popgen2
do
    echo $PREFIX
    INFILE=$INDIR/$PREFIX.vcf
    INSNP=$OUTDIR/$PREFIX.SNP.vcf
    ININDEL=$OUTDIR/$PREFIX.INDEL.vcf
    FILTEREDSNP=$OUTDIR/$PREFIX.filtered.SNP.vcf
    FILTEREDINDEL=$OUTDIR/$PREFIX.filtered.INDEL.vcf
    SNPONLY=$OUTDIR/$PREFIX.selected.SNP.vcf
    INDELONLY=$OUTDIR/$PREFIX.selected.INDEL.vcf

 if [ ! -f $INSNP ]; then
  java -Xmx3g -jar $GATK \
  -T SelectVariants \
  -R $GENOME \
  --variant $INFILE \
  -o $INSNP \
  -env \
  -ef \
  -restrictAllelesTo BIALLELIC \
  -selectType SNP
 fi

 if [ ! -f $ININDEL ]; then
  java -Xmx3g -jar $GATK \
  -T SelectVariants \
  -R $GENOME \
  --variant $INFILE \
  -o $ININDEL \
  -env \
  -ef \
  -selectType INDEL -selectType MIXED -selectType MNP
 fi

 if [ ! -f $FILTEREDSNP ]; then
   java -Xmx3g -jar $GATK \
   -T VariantFiltration -o $FILTEREDSNP \
   --variant $INSNP -R $GENOME \
   --clusterWindowSize 10  -filter "QD<2.0" -filterName QualByDepth \
   -filter "MQ<40.0" -filterName MapQual \
   -filter "QUAL<100" -filterName QScore \
   -filter "MQRankSum < -12.5" -filterName MapQualityRankSum \
   -filter "SOR > 4.0" -filterName StrandOddsRatio \
   -filter "FS>60.0" -filterName FisherStrandBias \
   -filter "ReadPosRankSum<-8.0" -filterName ReadPosRank \
   --missingValuesInExpressionsShouldEvaluateAsFailing 

#-filter "HaplotypeScore > 13.0" -filterName HaplotypeScore
#-filter "MQ0>=10 && ((MQ0 / (1.0 * DP)) > 0.1)" -filterName MapQualRatio \
 fi

 if [ ! -f $FILTEREDINDEL ]; then
  java -Xmx3g -jar $GATK \
  -T VariantFiltration -o $FILTEREDINDEL \
  --variant $ININDEL -R $GENOME \
  --clusterWindowSize 10  -filter "QD<2.0" -filterName QualByDepth \
  -filter "MQRankSum < -12.5" -filterName MapQualityRankSum \
  -filter "SOR > 4.0" -filterName StrandOddsRatio \
  -filter "FS>200.0" -filterName FisherStrandBias \
  -filter "InbreedingCoeff<-0.8" -filterName InbreedCoef \
  -filter "ReadPosRankSum<-20.0" -filterName ReadPosRank 
 fi

 if [ ! -f $SNPONLY ]; then
  java -Xmx16g -jar $GATK \
   -R $GENOME \
   -T SelectVariants \
   --variant $FILTEREDSNP \
   -o $SNPONLY \
   -env \
   -ef \
   --excludeFiltered
 fi

 if [ ! -f $INDELONLY ]; then
  java -Xmx16g -jar $GATK \
   -R $GENOME \
   -T SelectVariants \
   --variant $FILTEREDINDEL \
   -o $INDELONLY \
   --excludeFiltered 
 fi
done
