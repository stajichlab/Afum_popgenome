#!/usr/bin/bash

#SBATCH -p stajichlab --mem 8gb -p short -n 32 --out logs/mosdepth.log

module load parallel
module unload perl
module load miniconda3

parallel -j 32 mosdepth --fasta genome/FungiDB-39_AfumigatusAf293_Genome.fasta -x -n --by lib/chroms.bed -t 2 reports/{/.} {} ::: cram/*.cram
