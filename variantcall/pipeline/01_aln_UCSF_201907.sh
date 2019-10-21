#!/usr/bin/bash
#SBATCH --mem 96G --ntasks 8 --nodes 1 -J alnUCSFAfum
#SBATCH --out logs/UCSF_201907.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard
module load gatk/3.7

MEM=96g
CENTER=UCR_UCSF
INDIR=input/UCSF_20190725
TOPOUTDIR=tmp
ALNFOLDER=aln
HTCEXT=cram
HTCFORMAT=cram
GENOMESTRAIN=Af293

if [ -f config.txt ]; then
    source config.txt
fi
if [ -z $REFGENOME ]; then
    echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
    exit
fi

if [ ! -f $REFGENOME.dict ]; then
    echo "NEED a $REFGENOME.dict - make sure 00_index.sh is run"
fi
mkdir -p $TOPOUTDIR
SAMPFILE=UCSF_20190725_samples.csv

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi


if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "need to provide a number by --array or cmdline"
     exit
 fi
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN PREFIX LEFT RIGHT
do
    PAIR1=$INDIR/$LEFT
    PAIR2=$INDIR/$RIGHT
    echo "$PAIR1 $PAIR2"

    SAMFILE=$TOPOUTDIR/$STRAIN.unsrt.sam
    SRTED=$TOPOUTDIR/${STRAIN}.srt.bam
    DDFILE=$TOPOUTDIR/$STRAIN.DD.bam
    REALIGN=$TOPOUTDIR/$STRAIN.realign.bam
    INTERVALS=$TOPOUTDIR/$STRAIN.intervals
    FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
    READGROUP="@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$PREFIX\tPL:illumina\tCN:$CENTER"

    if [ ! -f $FINALFILE ]; then
	if [ ! -f $DDFILE ]; then
	    if [ ! -f $SRTED ]; then
		if [ -e $PAIR1 ]; then
		    echo "SAMFILE is $SAMFILE"
		    if [ ! -f $SAMFILE ]; then
			bwa mem -t $CPU -R $READGROUP $REFGENOME $PAIR1 $PAIR2 > $SAMFILE
		    fi
		else
		    echo "Cannot find $PAIR1, skipping $STRAIN"
		    exit
		fi
		samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${STRAIN}.fixmate.bam
		samtools sort --threads $CPU -O bam -o  $SRTED -T $TEMP $TEMP/${STRAIN}.fixmate.bam
		if [ -f $SRTED ]; then
		    rm -f $TEMP/${STRAIN}.fixmate.bam $SAMFILE
		fi

	    fi # SRTED file exists or was created by this block

	    time java -jar $PICARD MarkDuplicates I=$SRTED O=$DDFILE \
		METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
	    if [ -f $DDFILE ]; then
		rm -f $SRTED
	    fi
	fi # DDFILE is created after this or already exists

	#	if [ ! -f $TOPOUTDIRDIR/$STRAIN.DD.bai ]; then
	#	    time java -jar $PICARD BuildBamIndex I=$TOPOUTDIR/$STRAIN.DD.bam TMP_DIR=/scratch
	#	fi
	if [ ! -f $INTERVALS ]; then
	    time java -Xmx$MEM -jar $GATK \
		-T RealignerTargetCreator \
		-R $REFGENOME \
		-I $DDFILE \
		-o $INTERVALS
	fi

	if [ ! -f $REALIGN ]; then
	    time java -Xmx$MEM -jar $GATK \
		-T IndelRealigner \
		-R $REFGENOME \
		-I $DDFILE \
		-targetIntervals $INTERVALS \
		-o $REALIGN
	fi # REALIGN created or already existed

	samtools view -O $HTCFORMAT --threads $CPU \
	    --reference $REFGENOME -o $FINALFILE $REALIGN
	samtools index $FINALFILE

	if [ -f $FINALFILE ]; then
	    rm -f $DDFILE $REALIGN
	    rm -f $(echo $REALIGN | sed 's/bam$/bai/')
	    rm -f $(echo $DDFILE | sed 's/bam$/bai/')
	    rm -f $INTERVALS
	fi
    fi #FINALFILE created or already exists
done
