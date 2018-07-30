#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 8G --time 12:00:00  --out logs/strain_tree.%A.log

module load IQ-TREE
module load vcftools
SUBSET=2000
TREEDIR=strain_tree
BOOTSTRAPS=100
REPS=4 # make sure ntasks == 2*REPS
mkdir -p $TREEDIR
for EXT in Popgen2
do
    VCF=Variants_filter/A_fumigiatus_Af293.$EXT.selected.SNP.vcf
    VCFTAB=Variants_filter/A_fumigiatus_Af293.$EXT.selected.SNP.tab
    if [ ! -f $VCFTAB ]; then
        vcf-to-tab < $VCF > $VCFTAB
    fi
    for N in $(seq $REPS)
    do
        TABSUBSET=Variants_filter/A_fumigiatus_Af293.$EXT.rand${SUBSET}_${N}.SNP.tab
        OUTFAS=$(basename $TABSUBSET .tab)".fasaln"
        echo "running $EXT $VCF SET $N"
        # select a random subset of SNPs
        if [ ! -f $TREEDIR/$OUTFAS ]; then
            (head -n 1 $VCFTAB && tail -n +2 $VCFTAB | sort -R | head -n $SUBSET) > $TABSUBSET
            perl scripts/vcftab_to_fasta.pl -o $TREEDIR/$OUTFAS $VCFTAB
        fi
    done
    parallel -j $REPS iqtree -nt 4 -s {} -bb 1000 -m GTR+ASC ::: $TREEDIR/*.rand*.SNP.fasaln.varsites.phy
done
