#!/usr/bin/bash
#SBATCH --mem 2gb 

module load plink
plink --bcf ../vcf/A_fumigiatus_Af293.Popgen3.selected.SNP.bcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out AfumAf293.Run3


