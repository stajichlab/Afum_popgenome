module load bcftools/1.9

for file in *.selected*.vcf.gz; do 
	b=$(basename $file .vcf.gz); 
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' $file > $b.tab 
done

