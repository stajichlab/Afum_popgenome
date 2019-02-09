#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 2 -p short --out logs/shred.%a.log --mem 16G

mkdir -p logs
OUTDIR=shred_fastq
TOTALREADS=6000000
module load BBMap
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [ -f config.txt ]; then
    source config.txt
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "need to provide a number by --array or cmdline"
     exit
 fi
fi

mkdir -p $OUTDIR
refgenome=$(ls *.fasta | sed -n ${N}p)
base=$(basename $refgenome .fasta)
bbmap.sh ref=$refgenome build=$N
randomreads.sh paired=t gaussian=t insrate=0 delrate=0 maxsnps=0 adderrors=false \
	out1=$OUTDIR/${base}_1.fq.gz out2=$OUTDIR/${base}_2.fq.gz reads=$TOTALREADS minlength=99 \
	maxlength=101 seed=13 build=$N 
