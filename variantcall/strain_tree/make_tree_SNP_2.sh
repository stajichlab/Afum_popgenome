#!/bin/bash
#SBATCH --ntasks 16 --nodes 1 --out fasttree.%A.outlog --mem 8G --time 72:00:00 -p intel -J AfumSNPtree
module load IQ-TREE
module load fasttree
IN=A_fumigiatus_Af293.Popgen1.selected.SNP.fasaln
IN=A_fumigiatus_Af293.Popgen2.rand2000_1.SNP.fasaln
OUT=$(basename $IN .fasaln)".FT.tre"
#iqtree -nt AUTO -s $IN -bb 1000 -m GTR+ASC

FastTreeMP -gtr -gamma < $IN > $OUT 
