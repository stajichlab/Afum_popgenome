#!/bin/bash
#SBATCH --ntasks 5 --nodes 1 --out iqtree.%A.outlog --mem 2G --time 72:00:00 -p intel
module load IQ-TREE
IN=A_fumigiatus_Af293.Popgen1.selected.INDEL.noREF.fasaln
iqtree-omp -nt 4 -s $IN -st MORPH -b 100 -m MK 
