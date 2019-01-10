#!/usr/bin/bash
#SBATCH -p batch --ntasks 48 --mem 48g --time 12:00:00

module load samtools/1.9

CPU=$SLURM_CPUS_ON_NODE
JOBS=$(expr $CPU / 2)

module load python/3
HTCFOLDER=cram
if [ -f config.txt ]; then
    source config.txt
fi
./scripts/cram_to_combine.py  | parallel -j $JOBS
parallel -j $CPU samtools index {} ::: $HTCFOLDER/*.cram
