#!/bin/bash
#SBATCH --out logs/mosdepth_nonorm.log
# 5000 20000
for WINDOW in 10000 ;
do
 for file in $(ls coverage/mosdepth/*.${WINDOW}bp.regions.bed.gz)
 do
     # b=$(basename $file .${WINDOW}bp.regions.bed.gz | perl -p -e 's/(ctl([12]))\.AL1B/A$1.L1B-$2/; s/(\S+)\.(\S+)/$2/')
     b=$(basename $file .${WINDOW}bp.regions.bed.gz)
     GROUP="Clinical"
     #$(echo $b | perl -p -e 'my $in = $_; chomp($in); open(my $fh => "samples.csv") || die $!; while(my $n = <$fh>) { chomp($n); my @row = split(",",$n); $lookup{$row[0]} = [@row];} 
#my $new = $lookup{$in}->[-1]; 
#$_ = $new."\n"')
 # echo "$GROUP $b"
 mean=$(zcat $file | awk '{total += $4} END { print total/NR}') 
 zcat $file | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,'\"$GROUP\"','\"$b\"'}' 
 done  > coverage/mosdepth.${WINDOW}bp.gg_nonorm.tab
 pigz -f coverage/mosdepth.${WINDOW}bp.gg_nonorm.tab
done


