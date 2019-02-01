#!/usr/bin/bash

#SBATCH --mem=24gb --ntasks 16 --nodes 1
#SBATCH --time=48:00:00
#SBATCH -J makeTree --out logs/make_tree.log

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

mkdir -p $TREEDIR
for TYPE in SNP INDEL
do
    root=$FINALVCF/$PREFIX.selected_nofixed.$TYPE
    FAS=$TREEDIR/$PREFIX.nofixed.$TYPE.fasaln
    if [ -f $root.vcf ]; then
	module load tabix
	bgzip $root.vcf
	tabix $root.vcf.gz
    fi
    vcf=$root.vcf.gz
    tab=$root.bcftools.tab
    if [ ! -f $tab ]; then
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $tab
	bcftools query -H -e 'INFO/AF < 0.1' -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $root.highfreq.bcftools.tab
    fi
    if [ ! -f $FAS ]; then
	printf ">%s\n%s\n" $REFNAME $(bcftools query -e 'INFO/AF < 0.1' -f '%REF' ${vcf}) > $FAS
	for str in $(bcftools query -l ${vcf});
	do
	    printf ">%s\n%s\n" $str $(bcftools query -e 'INFO/AF < 0.1' -s $str -f '[%TGT]' ${vcf}) >> $FAS
	done
#	parallel $(bcftools query -e 'INFO/AF < 0.1' -s {} -f '[%TGT]' ${vcf}) ::: $(bcftools query -l ${vcf})
	perl -i -p -e 'if (/^>/) { s/[\(\)#]/_/g; s/_+/_/g; }' $FAS
    fi

    if [ ! -f $TREEDIR/$PREFIX.fasttree.tre ]; then
	FastTreeMP -gtr -gamma -nt < $FAS > $TREEDIR/$PREFIX.fasttree.tre
    fi
done
#iqtree-omp -nt $CPU -s $FAS -m GTR+ASC -b 100
