#!/usr/bin/bash

#SBATCH --mem=2gb --ntasks 2 --nodes 1
#SBATCH --time=2:00:00 -p short
#SBATCH -J speedyslice --out logs/speedyslice.%a.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [[ -f config.txt ]]; then
	source config.txt
else
	echo "Need a config.txt"
	exit
fi

if [[ -z $REFNAME ]]; then
	REFNAME=REF
fi
module load bcftools/1.9
module load samtools/1.9
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then 
 echo "need to provide a number by --array or cmdline"
 exit
fi

TMPOUT=tmp_out
mkdir -p $TMPOUT
for TYPE in SNP INDEL
do
    root=$FINALVCF/$PREFIX.selected_nofixed.$TYPE
    if [ -f $root.vcf ]; then
	bgzip $root.vcf
	tabix $root.vcf.gz
    fi
    vcf=$root.vcf.gz
    MAX=$(bcftools query -l ${vcf} | wc -l | awk '{print $1}')

    echo "$N $MAX"
    if [ $N -gt $MAX ]; then 
	echo "$N is too big, only $MAX lines in $SAMPFILE"
	exit
    fi
    outstrain=$(bcftools query -l ${vcf} | sed -n ${N}p)
    echo $outstrain
    printf ">%s\n%s\n" $outstrain $(bcftools query -e 'INFO/AF < 0.1' -s "$outstrain" -f '[%TGT]' ${vcf}) > $TMPOUT/$PREFIX.nofixed.$TYPE.$outstrain.fas_seq
    perl -i -p -e 'if (/^>/) { s/[\(\)#]/_/g; s/_+/_/g; } else { s/[\*\.]/-/g; }' $TMPOUT/$PREFIX.nofixed.$TYPE.$outstrain.fas_seq

done
