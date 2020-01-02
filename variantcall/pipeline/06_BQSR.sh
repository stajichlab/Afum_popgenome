#!/usr/bin/bash
#SBATCH -p intel -J BQSR --out logs/BQSR.%a.log
#SBATCH --time 48:00:00 --mem 16gb -N 1 -n 4

MEM=128g
module load samtools/1.9

HTCFOLDER=cram # default
HTCEXT=cram
if [ -f config.txt ]; then
    source config.txt
fi
DICT=$(echo $REFGENOME | sed 's/fasta$/dict/')

if [ ! -f $DICT ]; then
    module load samtools
    module load picard
    picard CreateSequenceDictionary R=$GENOMEIDX O=$DICT
fi

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then 
	echo "need to provide a number by --array slurm or on the cmdline"
	exit
    fi
fi

MAX=`ls $HTCFOLDER/*.$HTCEXT | wc -l`

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

SNPONLY=$FINALVCF/$PREFIX.selected.SNP.vcf.gz
INDELONLY=$FINALVCF/$PREFIX.selected.INDEL.vcf.gz
SNPNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.SNP.vcf.gz
INDELNOFIXED=$FINALVCF/$PREFIX.selected_nofixed.INDEL.vcf.gz

if [ ! -f $REFGENOME ]; then
    samtools faidx $REFGENOME
fi
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi

mkdir -p $BQSRALN
RECAL_BAM=$BQSRALN/${SAMPLE}_recal.bam 
RECAL_CRAM=$BQSRALN/${SAMPLE}_recal.cram

if [ ! -f $GVCFFOLDER/$SAMPLE.bqsr.g.vcf.gz ]; then
    if [ ! -f $RECAL_CRAM ]; then
	module load gatk/4
	#Base Quality Score Recalibration (BQSR) #1
	if [ ! -f $BQSRALN/${SAMPLE}_recal_data.table ]; then
	    gatk BaseRecalibrator --reference ${REFGENOME} -I ${ALNFILE} \
		--knownSites $SNPNOFIXED --knownSites $INDELNOFIXES \
		-O $BQSRALN/${SAMPLE}_recal_data.table
	fi
	
	#Base Quality Score Recalibration (BQSR) #2
	if [ ! -f $BQSRALN/${SAMPLE}_post_recal_data.table ]; then
	    gatk BaseRecalibrator --reference ${REFGENOME} -I ${ALNFILE} \
		--knownSites $SNPNOFIXED --knownSites $INDELNOFIXES \
		-BQSR $BQSRALN/${SAMPLE}_recal_data.table \
		-O $BQSRALN/${SAMPLE}_post_recal_data.table
	fi
	
	gatk ApplyBQSR --reference ${REFERENCE} -I ${ALNFILE} -O ${RECAL_BAM} \
	    --bqsr-recal-file $TMPOUT/${SAMPLE}_recal_data.table
	
	samtools index $RECAL_BAM
	
	if [ -f $RECAL_CRAM ]; then
	    time samtools view -O cram --threads $CPU \
		--reference $REFGENOME -o $RECAL_CRAM $RECAL_BAM
	    samtools index $RECAL_CRAM
	fi
	
	module unload gatk/4
    elif [ ! -f $RECAL_BAM ]; then
	time samtools view -O bam --threads $CPU \
	    --reference $REFGENOME -o $RECAL_BAM $RECAL_CRAM
	samtools index $RECAL_BAM
    fi
    
    module load gatk/3.8

    time java -Xmx${MEM} -jar $GATK \
	-T HaplotypeCaller \
	-ERC GVCF \
	-ploidy 1 \
	-I ${RECAL_BAM} -R ${REFGENOME} \
	-o $GVCFFOLDER/$SAMPLE.bqsr.g.vcf \
	-nct $CPU

    module load bcftools/1.9
    bgzip $GVCFFOLDER/$SAMPLE.bqsr.g.vcf
    tabix $GVCFFOLDER/$SAMPLE.bqsr.g.vcf.gz
    
    unlink ${RECAL_BAM}
    unlink ${RECAL_BAM}.bai
fi


