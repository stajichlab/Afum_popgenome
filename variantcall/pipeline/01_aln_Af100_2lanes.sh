#!/usr/bin/bash
#SBATCH --mem 24G --ntasks 8 --nodes 1 -J Af100
#SBATCH --out logs/AF100.bwa.%a.log --time 2:00:00 -p short

module load bwa/0.7.17
module unload java
module load java/8
module load picard
module load samtools/1.9
module load gatk/3.7

MEM=24g
GENOMESTRAIN=Af293
INDIR=input
TOPOUTDIR=tmp
ALNFOLDER=aln
HTCEXT=cram
HTCFORMAT=cram


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
SAMPFILE=Af100_samples.csv

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
ct=0
IFS=,
sed -n ${N}p $SAMPFILE | while read STRAIN FWD REV
do
    LIBRARY=$(basename $FWD _R1_001.fastq.gz)
    FINALMERGE=$ALNFOLDER/$STRAIN.$HTCEXT
    FINALLIST=()

    for LANE in DA002_lane1 DA002_lane2
    do
	PAIR1=$INDIR/$LANE/$FWD
	PAIR2=$INDIR/$LANE/$REV
	
	SAMFILE=$TOPOUTDIR/$STRAIN.$LANE.unsrt.sam
	SRTED=$TOPOUTDIR/${STRAIN}.$LANE.srt.bam
	DDFILE=$TOPOUTDIR/$STRAIN.$LANE.DD.bam
	REALIGN=$TOPOUTDIR/$STRAIN.$LANE.realign.bam
	INTERVALS=$TOPOUTDIR/$STRAIN.$LANE.intervals
	FINALFILE=$TOPOUTDIR/$STRAIN.$LANE.$HTCEXT    
	FINALLIST[$ct]=$FINALFILE
	ct=$(expr $ct + 1)
	READGROUP="@RG\tID:$STRAIN.$LANE\tSM:$STRAIN\tLB:$LIBRARY.$LANE\tPL:illumina\tCN:Seqmatic"
	if [ ! -f $FINALFILE ]; then
	    if [ ! -f $DDFILE ]; then
		if [ ! -f $SRTED ]; then
		    if [ -e $PAIR2 ]; then      	
			echo "SAMFILE is $SAMFILE"
			if [ ! -f $SAMFILE ]; then
			    echo "bwa mem -t $CPU -R $READGROUP -o $SAMFILE $REFGENOME $PAIR1 $PAIR2"
			    bwa mem -t $CPU -R $READGROUP -o $SAMFILE $REFGENOME $PAIR1 $PAIR2
			fi 
		    else
			echo "Cannot find $PAIR2, skipping $STRAIN"
			exit
		    fi
		    if [ ! -f $SAMFILE ]; then
			echo "no $SAMFILE exiting"
			exit
		    fi
		    samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${STRAIN}.fixmate.bam
		    samtools sort --threads $CPU -O bam -o  $SRTED -T $TEMP $TEMP/${STRAIN}.fixmate.bam
		    if [ $SRTED ]; then
			rm -f $TEMP/${STRAIN}.fixmate.bam $SAMFILE
		    fi
		fi # SRTED file exists or was created by this block

		time java -jar $PICARD MarkDuplicates I=$SRTED O=$DDFILE \
		    METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
		if [ -f $DDFILE ]; then
		    rm -f $SRTED
		fi
	    fi # DDFILE is created after this or already exists
	    
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
	    
	    samtools view -O $HTCFORMAT --threads $CPU --reference $REFGENOME -o $FINALFILE $REALIGN
	    samtools index $FINALFILE
	    if [ -f $FINALFILE ]; then
		rm -f $DDFILE $REALIGN		
		rm -f $(echo $REALIGN | sed 's/bam$/bai/')
		rm -f $(echo $DDFILE | sed 's/bam$/bai/')
		rm -f $INTERVALS
	    fi
	fi
    done
    # there should be a merging now?
    echo "$FINALMERGE $FINALLIST"
    LIST=$(printf ",%s" "${FINALLIST[@]}")
    if [ ! -f $FINALMERGE ]; then
	CMD="samtools merge --reference $REFGENOME --threads $CPU -O $HTCFORMAT $FINALMERGE $LIST"
	eval $CMD
	samtools index $FINALMERGE
	for file in ${FINALLIST[@]}; do
	    rm -f $file
	    rm -f $file.crai
	done
    fi
done
