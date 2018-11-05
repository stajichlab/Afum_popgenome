module load snpEff

java  -jar $SNPEFFJAR eff -dataDir `pwd`/data -v FungiDB-39_AfumigatusAf293  ../0000.vcf > COH_specific_vars.vcf

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' COH_specific_vars.vcf >  COH_specific_vars.tab
grep -v downstream COH_specific_vars.tab | grep -v upstream  | grep -v intergenic > genic_only.tab
