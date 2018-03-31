#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 3 --mem 8G --time 12:00:00  --out logs/strain_tree.log

module load IQ-TREE
module load vcftools
TREEDIR=strain_tree
VCF=Variants_filter/A_fumigiatus_Af293.Popgen1.selected.SNP.vcf
VCFTAB=Variants_filter/A_fumigiatus_Af293.Popgen1.selected.SNP.tab
OUTFAS=$(basename $VCFTAB .tab)".fasaln"
BOOTSTRAPS=100

if [ ! -f $VCFTAB ]; then
 vcf-to-tab < $VCF > $VCFTAB
fi

mkdir -p $TREEDIR

if [ ! -f $TREEDIR/$OUTFAS ]; then
 perl scripts/vcftab_to_fasta.pl -o $TREEDIR/$OUTFAS $VCFTAB
fi

pushd $TREEDIR

iqtree-omp -nt 3 -s $OUTFAS -b $BOOTSTRAPS -m GTR+ASC

