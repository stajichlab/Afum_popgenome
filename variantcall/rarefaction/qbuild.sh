#!/usr/bin/bash
#SBATCH -p batch --mem 16gb -n 16

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
module load parallel
module load bcftools/1.9
module unload perl

VCF=A_fumigiatus_Af293.Popgen6.selected_nofixed.SNP.vcf.gz
MAX=$(bcftools query -l ${VCF} | wc -l | awk '{print $1}')
SAMPLIST=names.txt
bcftools query -l $VCF > $SAMPLIST
for n in $(seq 1 2)
do
	ODIR=subsets/iter$n
	mkdir -p $ODIR
	parallel -j 8 shuf -n {} -o $ODIR/samp{}.txt $SAMPLIST ::: $(seq 1 10 $MAX)
	echo -e "SAMPLES\tVARIANTS" > $ODIR.counts.tsv
	for x in $(seq 1 10 $MAX)
	do
		COUNT=$(bcftools view --samples-file $ODIR/samp${x}.txt --min-ac=1 $VCF -H -G --threads $CPU | wc -l | awk '{print $1}')
		echo -e "$x\t$COUNT" >> $ODIR.counts.tsv
	done
done

