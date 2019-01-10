#!/usr/bin/bash
#SBATCH -J GATK.HTC --out logs/GATK_HTC.%a.log --ntasks 8 --nodes 1 --mem 16G

module unload java
module load java/8
module load gatk/3.8
module load bcftools/1.9
module load samtools/1.9
module load picard
module load tabix

MEM=32g

ALNFOLDER=aln
VARIANTFOLDER=gvcf
HTCFORMAT=cram #default but may switch back to bam
HTCFOLDER=cram # default
if [ -f config.txt ]; then
    source config.txt
fi
GENOMEIDX=$GENOMEFOLDER/$GENOMEFASTA

mkdir -p $VARIANTFOLDER

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array slurm or on the cmdline"
 exit
fi

MAX=$(ls $HTCFOLDER/*.$HTCEXT | wc -l | awk '{print $1}')

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPLESINFO"
 exit
fi
hostname
ALNFILE=$(ls $HTCFOLDER/*.$HTCEXT | sed -n ${N}p)
echo "ALNFILE=$ALNFILE"
if [[ $ALNFILE == "" ]]; then
    echo "cannot find samples in the folder $HTCFOLDER/*.$HTCEXT, exiting ($N)"
    exit
fi
SAMPLE=$(basename $ALNFILE .$HTCEXT)

if [ ! -e $ALNFILE ]; then
    echo "Cannot find $ALNFILE"
    exit
fi
if [ ! -f $VARIANTFOLDER/$SAMPLE.g.vcf.gz ]; then
    if [ ! -f $VARIANTFOLDER/$SAMPLE.g.vcf ]; then
	java -Xmx${MEM} -jar $GATK \
	    -T HaplotypeCaller \
	    -ERC GVCF \
	    -ploidy 1 \
	    -I $ALNFILE -R $GENOMEIDX \
	    -o $VARIANTFOLDER/$SAMPLE.g.vcf -nct $CPU
    fi
    if [ -f $VARIANTFOLDER/$SAMPLE.g.vcf ]; then
	bgzip $VARIANTFOLDER/$SAMPLE.g.vcf
	tabix $VARIANTFOLDER/$SAMPLE.g.vcf.gz
    fi
fi
