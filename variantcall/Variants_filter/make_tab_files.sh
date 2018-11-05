module load bcftools/1.9

bcftools query -l A_fumigiatus_Af293.Popgen2.selected.SNP.vcf.gz | grep -v COH > COH_exclude.txt
bcftools view -S COH_exclude.txt A_fumigiatus_Af293.Popgen2.selected.SNP.vcf.gz > A_fumigiatus_Af293.Popgen2.selected.SNP.noCOH.vcf
bcftools +fill-tags A_fumigiatus_Af293.Popgen2.selected.SNP.noCOH.vcf > A_fumigiatus_Af293.Popgen2.selected.SNP.noCOH.retag.vcf

for file in *.vcf; do 
	b=$(basename $file .vcf); 
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' $file > $b.tab 
done

