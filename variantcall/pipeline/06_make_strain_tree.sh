#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 3 --mem 8G --time 12:00:00  --out logs/strain_tree.%A.log

module load IQ-TREE
module load vcftools
SUBSET=10000
TREEDIR=strain_tree
BOOTSTRAPS=50
mkdir -p $TREEDIR
for EXT in Popgen2
do
    VCF=Variants_filter/A_fumigiatus_Af293.$EXT.selected.SNP.vcf
    VCFTAB=Variants_filter/A_fumigiatus_Af293.$EXT.selected.SNP.tab
    if [ ! -f $VCFTAB ]; then
        vcf-to-tab < $VCF > $VCFTAB
    fi
    for N in $(seq 5)
    do
        TABSUBSET=Variants_filter/A_fumigiatus_Af293.$EXT.rand$SUBSET_$N.SNP.tab
        OUTFAS=$(basename $TABSUBSET .tab)".fasaln"
        echo "running $EXT $VCF SET $N"
        # select a random subset of SNPs
        (head -n 1 $VCFTAB && tail -n +2 $VCFTAB | sort -R | head -n $SUBSET) > $TABSUBSET
        if [ ! -f $TREEDIR/$OUTFAS ]; then
            perl scripts/vcftab_to_fasta.pl -o $TREEDIR/$OUTFAS $VCFTAB
        fi

        if [ ! -f $TREEDIR/$OUTFAS.contree ]; then
	    pushd $TREEDIR
#	    iqtree -nt 3 -s $OUTFAS -b $BOOTSTRAPS -m GTR+ASC
	    popd
        fi
    done
done
