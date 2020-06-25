#!/usr/bin/bash
#SBATCH -p short

INDIR=input/SRA
TOPOUTDIR=aln

SAMPFILE=SRA_samples.csv

N=0
IFS=,
while read RUN STRAIN SAMPLE CENTER EXP PROJ
do
	if [ ! -e $TOPOUTDIR/${RUN}.cram ]; then
		echo "($N) $RUN $STRAIN"
		RUNME="$RUNME,$N"
	fi
	N=$(expr $N + 1)
done < $SAMPFILE

echo "--array=$RUNME pipeline/01_aln_sra.sh"
