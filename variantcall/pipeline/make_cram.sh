#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 8gb -p short 

module load samtools/1.8
module load parallel

parallel -j 12 samtools view --threads 2 -O cram -o {.}.cram -T genome/FungiDB-39_AfumigatusAf293_Genome.fasta {} ::: bam/*.bam
parallel -j 24 samtools index {} ::: bam/*.cram
