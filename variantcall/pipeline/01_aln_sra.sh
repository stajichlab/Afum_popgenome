#!/usr/bin/bash
#SBATCH --mem 8G --ntasks 8 --nodes 1 -J sraAfum
#SBATCH --out logs/Afum_sra.bwa.%a.log --time 8:00:00
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
INDIR=input/SRA
TOPOUTDIR=aln
ALNFOLDER=bam
HTCEXT=cram
HTCFORMAT=cram
SAMPFILE=SRA_samples.csv

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


if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then 
 echo "need to provide a number by --array or cmdline"
 exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
echo "$N $MAX for $SAMPFILE"
if [ $N -gt $MAX ]; then 
 echo "$N is too big, only $MAX lines in $SAMPFILE"
 exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUN STRAIN SAMPLE CENTER EXP PROJ
do
  PAIR1=$INDIR/${RUN}_1.fastq.gz
  PAIR2=$INDIR/${RUN}_2.fastq.gz

  SAMFILE=$TOPOUTDIR/$RUN.unsrt.sam
  SRTED=$TOPOUTDIR/${RUN}.srt.bam
  DDFILE=$TOPOUTDIR/${RUN}.DD.bam
  REALIGN=$TOPOUTDIR/${RUN}.realign.bam
  INTERVALS=$TOPOUTDIR/${RUN}.intervals
  FINALFILE=$ALNFOLDER/${RUN}.$HTCEXT    
  READGROUP="@RG\tID:$RUN\tSM:$STRAIN\tLB:$PREFIX\tPL:illumina\tCN:$CENTER"
    
  if [ ! -f $FINALFILE ]; then
      if [ ! -f $DDFILE ]; then
	  if [ ! -f $SRTED ]; then
	      if [ -e $PAIR2 ]; then
		  echo "RUNNING paired-end bwa"
		  bwa mem -t $CPU -R $READGROUP $GENOME $PAIR1 $PAIR2 > $SAMFILE
		  samtools fixmate --threads $CPU -O bam $SAMFILE $TEMP/${RUN}.fixmate.bam
		  samtools sort --threads $CPU -O bam -o  $SRTED -T $TEMP $TEMP/${RUN}.fixmate.bam
		  /usr/bin/rm $TEMP/${RUN}.fixmate.bam $SAMFILE
	      elif [ -e $PAIR1 ]; then
		  echo "RUNNING unpaired bwa"
    		  bwa mem -t $CPU -R "@RG\tID:$RUN\tSM:$STRAIN\tLB:$RUN\tPL:illumina\tCN:$CENTER" $GENOME $PAIR1 | samtools sort -@ $CPU -O bam -T $TEMP -o $SRTED
	      else
		  echo "NO $PAIR1 and no $PAIR2?"
		  exit
	      fi
	  fi # SRTED file exists or was created by this block
	  time picard MarkDuplicates I=$SRTED O=$DDFILE \
	      METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
	   if [ ! -f $DDFILE.bai ]; then
		picard BuildBamIndex I=$DDFILE TMP_DIR=/scratch
	   fi	    
      fi # DDFILE is created after this or already exists
      
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
	      -I $DDFILE \
	      -targetIntervals $INTERVALS \
	      -o $REALIGN
      fi
      
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
