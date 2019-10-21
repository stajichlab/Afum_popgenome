#!/usr/bin/bash

#SBATCH --mem=24gb --ntasks 16 --nodes 1
#SBATCH --time=48:00:00
#SBATCH -J makeTree --out logs/make_tree_outgroup.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [[ -f config.txt ]]; then
	source config.txt
else
	echo "Need a config.txt"
	exit
fi

if [[ -z $REFNAME ]]; then
	REFNAME=REF
fi
module load bcftools/1.9
module load samtools/1.9
module load IQ-TREE
module load fasttree
module load tabix

query_vcf_fasta() {
    printf ">%s\n%s\n" $1 $(bcftools query -R -s $1 -f '[%TGT]' $2)
}
export -f query_vcf_fasta

mkdir -p $TREEDIR
for TYPE in SNP 
do
    root=$FINALVCF/${PREFIX}_mergeasm.selected_nofixed.$TYPE
    FAS=$TREEDIR/${PREFIX}_mergeasm.nofixed.$TYPE.fasaln
    if [ -f $root.vcf ]; then
	module load tabix
	bgzip $root.vcf
	tabix $root.vcf.gz
    fi
    vcf=$root.vcf.gz
    highfreqvcf=$root.highfreq.vcf.gz
    tab=$root.bcftools.tab
    if [ ! -f $tab ]; then
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $tab
	bcftools query -H -e 'INFO/AF < 0.1' -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $root.highfreq.bcftools.tab
    fi
    if [ ! -f $highfreqvcf ]; then
	bcftools filter -Oz -e 'INFO/AF < 0.1' -o $highfreqvcf ${vcf}
	if [ -f $highfreqvcf ]; then
	    tabix $highfreqvcf
	    bcftools index $highfreqvcf
	fi
    fi

    if [ ! -f $FAS ]; then
	printf ">%s\n%s\n" $REFNAME $(bcftools query -f '%REF' ${highfreqvcf}) > $FAS
	for str in $(bcftools query -l ${vcf});
	do
	    printf ">%s\n%s\n" $str $(bcftools query -s $str -f '[%TGT]' ${highfreqvcf}) >> $FAS
	done
	perl -i -p -e 'if (/^>/) { s/[\(\)#]/_/g; s/_+/_/g; }' $FAS
    fi

    if [ ! -f $TREEDIR/${PREFIX}_mergeasm.fasttree.tre ]; then
	FastTreeMP -gtr -gamma -nt < $FAS > $TREEDIR/${PREFIX}_mergeasm.fasttree.tre
    fi
done
#iqtree-omp -nt $CPU -s $FAS -m GTR+ASC -b 100
