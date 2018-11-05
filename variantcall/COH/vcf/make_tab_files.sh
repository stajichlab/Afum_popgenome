module load bcftools/1.9


bcftools concat -Ou -a COH.Run1.selected.SNP.vcf.gz COH.Run1.selected.INDEL.vcf.gz | bcftools sort -Ob -o COH.Run1.selected.bcf
bcftools concat -Ou -a COH.Run1.selected_nofixed.SNP.vcf.gz COH.Run1.selected_nofixed.INDEL.vcf.gz | bcftools sort -Ob -o COH.Run1.selected_nofixed.bcf

for file in *.bcf; do 
	b=$(basename $file .bcf); 
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' $file > $b.tab 
	pigz $b.tab
done

