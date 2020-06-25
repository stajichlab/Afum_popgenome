#!/usr/bin/bash

HTCFOLDER=cram
VARIANTFOLDER=Variants
HTCEXT=cram
if [ -f config.txt ]; then
	source config.txt
fi
N=1
RUN=""
for file in $HTCFOLDER/*.$HTCEXT
do
	b=$(basename $file .$HTCEXT)
	if [ ! -f $VARIANTFOLDER/$b.g.vcf.gz ]; then
		echo "($N) $file"
		if [ -z $RUN ]; then
			RUN=$N
		else
			RUN="$RUN,$N"
		fi
	fi
	N=$(expr $N + 1)
done
echo "sbatch --array=$RUN pipeline/03_GATK_HTC_gatk4.sh"
