#!/bin/bash
#SBATCH --out logs/mosdepth_nonorm.log

for file in $(ls coverage/mosdepth_gene/*.regions.bed.gz)
do
 b=$(basename $file .regions.bed.gz)
 pigz -dc $file | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,'\"$b\"'}'
done  > coverage/mosdepth_gene.gg_nonorm.tab
pigz -f coverage/mosdepth_gene.gg_nonorm.tab
