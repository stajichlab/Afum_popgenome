#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 48gb -J Af.AAFTF --out logs/AAFTF_full.%A_%a.log -p batch --time 48:00:00

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

module load AAFTF

OUTDIR=input
SAMPLEFILE=samples.dat
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
if [[ "$BASE" =~ ^\# ]]; then
 echo "skipping $BASE"
 exit
fi
PHYLUM=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
ASM=asm

mkdir -p $ASM

if [ -z $CPU ]; then
    CPU=1
fi

ASMFILE=$ASM/${BASE}.spades.fasta
WORKDIR=working_AAFTF
VECCLEAN=$ASM/${BASE}.vecscreen.fasta
PURGE=$ASM/${BASE}.sourpurge.fasta
CLEANDUP=$ASM/${BASE}.rmdup.fasta
PILON=$ASM/${BASE}.pilon.fasta
SORTED=$ASM/${BASE}.sorted.fasta
STATS=$ASM/${BASE}.sorted.stats.txt

mkdir -p $WORKDIR
LEFT=
echo "$BASE"
if [ ! -f $WORKDIR/${BASE}_cleaned_1.fastq.gz ]; then
    echo "$OUTDIR/${BASE}_R1.fq.gz $OUTDIR/${BASE}_R2.fq.gz"

    if [ ! -f $WORKDIR/${BASE}_1P.fastq ]; then
	AAFTF trim \
	    --left $OUTDIR/${BASE}_R1.fq.gz --right $OUTDIR/${BASE}_R2.fq.gz \
	    -c $CPU -o $WORKDIR/${BASE} --trimmomatic $TRIMMOMATIC 
    fi
    if [ ! -f $WORKDIR/${BASE}_cleaned_1.fastq.gz ]; then
	AAFTF filter -w $WORKDIR -c $CPU --left $WORKDIR/${BASE}_1P.fastq.gz --right $WORKDIR/${BASE}_2P.fastq.gz --aligner bowtie2 -o $WORKDIR/${BASE}
    fi
fi
if [ ! -f $ASMFILE  ]; then
	AAFTF assemble -w $WORKDIR/spades_${BASE}  -c $CPU --debug \
	-o $ASMFILE --left $WORKDIR/${BASE}_cleaned_1.fastq.gz --right $WORKDIR/${BASE}_cleaned_2.fastq.gz --spades_tmpdir /scratch/$USER/$SLURM_JOBID
fi
if [ -s $ASMFILE ]; then
	#rm -rf $WORKDIR/spades_${BASE}/K?? $WORKDIR/spades_${BASE}/tmp
	#rm -rf $WORKDIR/${BASE}_[12][UP].fastq
	echo "rm -rf $WORKDIR/${BASE}_[12][UP].fastq"
fi

if [ ! -f $ASMFILE ]; then
 echo "SPADES must have failed, exiting"
 exit
fi
if [ ! -f $VECCLEAN ]; then
    AAFTF vecscreen -i $ASMFILE -c $CPU -o $VECCLEAN
fi

if [ ! -f $PURGE ]; then
    AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $WORKDIR/${BASE}_cleaned_1.fastq.gz --right ${WORKDIR}/${BASE}_cleaned_1.fastq.gz
fi

COUNT=$(grep -c ">" $PURGE)
if [ "$COUNT" -gt 3000 ]; then
	echo "too many contigs to run rmdup ($COUNT) skipping that step and jumping to Pilon"
	if [ ! -f $PILON ]; then
		AAFTF pilon -i $PURGE -o $PILON -c $CPU -w $WORKDIR --left $WORKDIR/${BASE}_cleaned_1.fastq.gz --right ${WORKDIR}/${BASE}_cleaned_1.fastq.gz
	fi
else
	if [ ! -f $CLEANDUP ]; then
    		AAFTF rmdup -i $PURGE -o $CLEANDUP -p $BASE -c $CPU
	fi
	if [ ! -f $PILON ]; then
		echo "pilon -i $CLEANDUP -o $PILON -p $BASE -w $WORKDIR"
    		AAFTF pilon -i $CLEANDUP -o $PILON -c $CPU -w $WORKDIR --left $WORKDIR/${BASE}_cleaned_1.fastq.gz --right ${WORKDIR}/${BASE}_cleaned_1.fastq.gz
	fi
fi

if [ ! -f $SORTED ]; then
    AAFTF sort -i $PILON -o $SORTED
fi

if [ ! -f $STATS ]; then
	AAFTF assess -i $SORTED -r $STATS
fi
