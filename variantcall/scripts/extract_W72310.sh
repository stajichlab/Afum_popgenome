#!/usr/bin/bash
#SBATCH --mem 8G --nodes 1 -p short --out logs/subset_W72310.%A.log --ntasks 2
module load bcftools/1.9
IN=vcf/A_fumigiatus_Af293.Popgen5.selected_nofixed
#SNP.vcf.gz
STRAINS="W72310,Af293"
OUT=W72310/vcf/A_fumigiatus_Af293.Popgen5.selected_nofixed_fromBig
#SNP.vcf.gz

#bcftools view -Ou -s $STRAINS $IN.SNP.vcf.gz | bcftools +fill-tags -- -t AF | bcftools view -Oz --exclude 'AF==1 || AF==0' > $OUT.SNP.vcf.gz
#bcftools view -Ou -s $STRAINS $IN.INDEL.vcf.gz | bcftools +fill-tags -- -t AF | bcftools view -Oz --exclude 'AF==1 || AF==0' > $OUT.INDEL.vcf.gz

bcftools view -Oz --private -s $STRAINS $IN.SNP.vcf.gz > $OUT.v2.SNP.vcf.gz
bcftools view -Oz --private -s $STRAINS $IN.INDEL.vcf.gz > $OUT.v2.INDEL.vcf.gz

