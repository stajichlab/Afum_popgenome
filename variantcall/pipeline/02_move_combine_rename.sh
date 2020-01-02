#!/usr/bin/bash
#SBATCH -p batch --ntasks 48 --mem 48g --time 12:00:00

module load samtools/1.9
module unload perl
module load parallel
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
JOBS=$(expr $CPU / 2)

module load python/3
HTCFOLDER=cram
if [ -f config.txt ]; then
    source config.txt
fi
./scripts/cram_to_combine.py  | parallel -j $JOBS
for a in $(find $HTCFOLDER -name "*.cram"); do if [ $a -nt $a.crai ]; then  echo $a; fi; done | parallel -j $CPU samtools index {}
