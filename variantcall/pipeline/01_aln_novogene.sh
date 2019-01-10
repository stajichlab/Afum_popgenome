#!/usr/bin/bash
#SBATCH --mem 24G --ntasks 8 --nodes 1 -J bwa.Afum
#SBATCH --out logs/Afum_novogene.bwa.%a.log --time 8:00:00

module load bwa/0.7.17
module unload java
module load java/8
module load picard
module load gatk/3.7

MEM=24g
CENTER=Novogene
GENOME=genome/Af293
GENOMESTRAIN=Af293
GENOMEIDX=$GENOME.fasta
REFGENOME=genome/FungiDB-39_AfumigatusAf293_Genome.fasta
INDIR=input/Novogene
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

SAMPFILE=Novogene_samples.csv
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
			bwa mem -t $CPU -R $READGROUP $GENOME $PAIR1 $PAIR2 > $SAMFILE
		    fi 
		else
		    echo "Cannot find $PAIR1, skipping $STRAIN"
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

	if [ -f $FINALFILE ]; then
	    rm -f $DDFILE $REALIGN
	    rm $(basename $REALIGN .bam)".bai"
	fi
    fi #FINALFILE created or already exists  
done
