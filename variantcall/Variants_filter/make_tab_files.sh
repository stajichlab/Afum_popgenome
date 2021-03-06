module load bcftools/1.9

bcftools query -l A_fumigiatus_Af293.Popgen2.selected.SNP.vcf.gz | grep -v COH > COH_exclude.txt
bcftools concat -Ou -a A_fumigiatus_Af293.Popgen2.selected.SNP.vcf.gz A_fumigiatus_Af293.Popgen2.selected.INDEL.vcf.gz  | bcftools sort -Ob -o A_fumigiatus_Af293.Popgen2.selected.bcf

bcftools view -Ou -S COH_exclude.txt A_fumigiatus_Af293.Popgen2.selected.bcf | bcftools +fill-tags -Ob -o A_fumigiatus_Af293.Popgen2.selected.noCOH.bcf 

for file in *.bcf; do 
	b=$(basename $file .bcf); 
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' $file > $b.tab 
	pigz $b.tab
done

