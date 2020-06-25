#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --time 2:00:00 -p short --mem 64G --out logs/mosdepth_gene.parallel.log
#SBATCH -J modepth
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
if [ -f config.txt ]; then
	source config.txt
fi
module unload perl
module load mosdepth
OUTDIR=coverage/mosdepth_gene
mkdir -p $OUTDIR

GENEBED=genome/Af293.genes.bed
echo $REFGENOME
if [ -z $REFGENOME ]; then
	echo "need a REFGENOME file for CRAM"
	exit
fi

parallel --jobs $CPU mosdepth -f $REFGENOME -T 1,10,50,100,200 -n --by $GENEBED -t 2 \
 "{= s:${HTCFOLDER}:${OUTDIR}:; s:\.cram:: =}" {} ::: ${HTCFOLDER}/*.cram

bash scripts/mosdepth_prep_ggplot_genebed.sh
bash scripts/mosdepth_prep_ggplot_genebed_nonorm.sh
#mkdir -p plots
#Rscript Rscripts/plot_mosdepth_CNV_genebed.R
