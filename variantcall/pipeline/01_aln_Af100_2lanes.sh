#!/usr/bin/bash
#SBATCH --mem 4G --ntasks 8 --nodes 1 -J bwa.Afum --out logs/Afum.Af100.bwa.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard

GENOME=genome/Af293
GENOMESTRAIN=Af293
INDIR=input
TOPOUTDIR=aln
mkdir -p $TOPOUTDIR

TEMP=/scratch
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then 
 CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=Af100_samples.csv
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then 
 echo "need to provide a number by --array or cmdline"
 exit
fi

MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
sed -n ${N}p $SAMPFILE | while read STRAIN FWD REV;
do
    LIBRARY=$(basename $FWD .fastq.gz)
    LIBRARY1=$(basename $FWD .fastq.gz)
    LIBRARY2=$(basename $REV .fastq.gz)

    for LANE in DA002_lane1 DA002_lane2
    do
	OUTDIR=$TOPOUTDIR/$LANE
	mkdir -p $OUTDIR
	PAIR1=${INDIR}/$LANE/${LIBRARY1}_val_1.fq.gz
	PAIR2=${INDIR}/$LANE/${LIBRARY2}_val_2.fq.gz
	
	echo "... files are $PAIR1 $PAIR2 $LIBRARY"
	SAMFILE=NULL
	
	if [ -f $PAIR2 ]; then
	    SAMFILE=$OUTDIR/${STRAIN}.PE.unsrt.sam
	    echo "SAMFILE is $SAMFILE"
	    if [ ! -f $SAMFILE ]; then
		bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$LIBRARY\tPL:illumina\tCN:Seqmatic" $GENOME $PAIR1 $PAIR2 > $SAMFILE
	    fi 
	    if [ ! -f $OUTDIR/${STRAIN}.PE.bam ]; then
		samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${STRAIN}.fixmate.bam
		samtools sort --threads $CPU -O bam -o $OUTDIR/${STRAIN}.PE.bam -T $TEMP $TEMP/${STRAIN}.fixmate.bam
		/usr/bin/rm $TEMP/${STRAIN}.fixmate.bam $SAMFILE
	    fi
	else
	    echo "unsure how to deal with unpaired multi-lane results"
	    echo "skipping..."
	    exit
#	    SAMFILE=$OUTDIR/${ID}.SE.unsrt.sam
	    # echo "SAMFILE is $SAMFILE"
	    # if [ ! -f $SAMFILE ]; then
    	    # 	bwa mem -t $CPU -R "@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$LIBRARY\tPL:illumina\tCN:Seqmatic" $GENOME $PAIR1 > $SAMFILE
	    # fi
	    # if [ ! -f $OUTDIR/${STRAIN}.SE.bam ]; then
	    # 	samtools view --threads $CPU -b $SAMFILE > $TEMP/${STRAIN}.unsrt.bam	
	    # 	samtools sort --threads $CPU -O bam -o $OUTDIR/${STRAIN}.SE.bam -T $TEMP $TEMP/${STRAIN}.unsrt.bam
	    # 	/usr/bin/rm $TEMP/${STRAIN}.unsrt.bam $SAMFILE
	    # fi
	fi
    done
    # there should be a merging now?
    samtools merge ${TOPOUTDIR}/$STRAIN.PE.bam ${TOPOUTDIR}/DA002_lane?/${STRAIN}.PE.bam
    echo "unlink ${TOPOUTDIR}/DA002_lane?/${STRAIN}.PE.bam"
done
