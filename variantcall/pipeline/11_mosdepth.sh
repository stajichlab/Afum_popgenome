#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --time 2:00:00 -p short --mem 64G --out mosdepth.parallel.log
#SBATCH -J modepth
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
if [ -f config.txt ]; then
	source config.txt
fi
module load miniconda3
mkdir -p coverage/mosdepth

echo $REFGENOME
if [ -z $REFGENOME ]; then
	echo "need a REFGENOME file for CRAM"
	exit
fi
WINDOW=10000
parallel --jobs $CPU mosdepth -f $REFGENOME -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:${ALNFOLDER}\/:coverage/mosdepth/:; s:\.cram:.${WINDOW}bp: =}" {} ::: ${ALNFOLDER}/*.cram

bash scripts/mosdepth_prep_ggplot.sh
mkdir -p plots
Rscript Rscripts/plot_mosdepth_CNV.R
