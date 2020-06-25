#!/bin/bash

for file in $(ls coverage/mosdepth_gene/*.regions.bed.gz)
do
 b=$(basename $file .regions.bed.gz)
 mean=$(pigz -dc $file | awk '{total += $5} END { print total/NR}')
 pigz -dc $file | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5/'$mean','\"$b\"'}'
done  > coverage/mosdepth_gene.gg.tab
pigz -f coverage/mosdepth_gene.gg.tab
