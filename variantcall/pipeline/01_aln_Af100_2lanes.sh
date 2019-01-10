#!/usr/bin/bash
#SBATCH --mem 24G --ntasks 8 --nodes 1 -J bwa.Afum --out logs/Afum.Af100.bwa.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard
module load samtools/1.9

MEM=24g
GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input
TOPOUTDIR=aln
ALNFOLDER=bam
HTCEXT=cram
HTCFORMAT=cram

if [ -f config.txt ]; then
    source config.txt
fi

mkdir -p $TOPOUTDIR

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=Af100_samples.csv
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
sed -n ${N}p $SAMPFILE | while read STRAIN FWD REV
do
    LIBRARY=$(basename $FWD _R1_001.fastq.gz)
    FINALLIST=""

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
	
	READGROUP="@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$LIBRARY\tPL:illumina\tCN:Seqmatic"
	if [ ! -f $FINALFILE ]; then
	    if [ ! -f $DDFILE ]; then
		if [ ! -f $SRTED ]; then
		    if [ -e $PAIR2 ]; then      	
			echo "SAMFILE is $SAMFILE"
			if [ ! -f $SAMFILE ]; then
			    echo "bwa mem -t $CPU -R $READGROUP -o $SAMFILE $GENOME $PAIR1 $PAIR2"
			    bwa mem -t $CPU -R $READGROUP -o $SAMFILE $GENOME $PAIR1 $PAIR2
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
		    /usr/bin/rm $TEMP/${STRAIN}.fixmate.bam $SAMFILE
		fi # SRTED file exists or was created by this block
		
		time java -jar $PICARD MarkDuplicates I=$SRTED O=$DDFILE \
		    METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
		
	    fi # DDFILE is created after this or already exists
	    
	    #	if [ ! -f $TOPOUTDIRDIR/$STRAIN.DD.bai ]; then
	    #	    time java -jar $PICARD BuildBamIndex I=$TOPOUTDIR/$STRAIN.DD.bam TMP_DIR=/scratch
	    #	fi
	    if [ ! -f $INTERVALS ]; then 
		time java -Xmx$MEM -jar $GATK \
		    -T RealignerTargetCreator \
		    -R $GENOMEIDX \
		    -I $DDFILE \
		    -o $INTERVALS
	    fi
	    
	    if [ ! -f $REALIGN ]; then
		time java -Xmx$MEM -jar $GATK \
		    -T IndelRealigner \
		    -R $GENOMEIDX \
		    -I $TOPOUTDIR/$STRAIN.DD.bam \
		    -targetIntervals $INTERVALS \
		    -o $REALIGN
	    fi # REALIGN created or already existed
	    
	    samtools view -O $HTCFORMAT --threads $CPU --reference $REFGENOME -o $FINALFILE $REALIGN
	    samtools index $FINALFILE
	    FINALLIST="$FINALCRAM $FINALFILE"
	    if [ -f $FINALFILE ]; then
		rm -f $DDFILE $REALIGN
		rm $(basename $REALIGN .bam)".bai"
	    fi
	fi
    done
    # there should be a merging now?
    if [ ! -f $ALNFOLDER/$STRAIN.$HTCEXT ]; then
	samtools merge --reference $REFGENOME --threads $CPU -O $HTCFORMAT $FINALLIST
	rm -f $FINALLIST
    fi
done
