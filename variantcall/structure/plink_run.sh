#!/usr/bin/bash
#SBATCH --mem 2gb 

module load plink
plink --bcf ../Variants_filter/A_fumigiatus_Af293.Popgen2.selected.bcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out AfumAf293.Run2


