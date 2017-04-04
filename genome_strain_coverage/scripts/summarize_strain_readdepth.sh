#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 1

module load bedtools
module load samtools
BAMLIST=all.bams.list
SRA=sra_strains.tab

mkdir -p depth

samtools depth -f $BAMLIST | perl scripts/depth_alllib_sum.pl -b $BAMLIST -s $SRA > depth/strain.depths.tab
