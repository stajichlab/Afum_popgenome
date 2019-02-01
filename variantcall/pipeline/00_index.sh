#!/usr/bin/bash
module load samtools/1.9
module load bwa/0.7.17
module load picard
if [ -f config.txt ]; then
	source config.txt
fi
pushd genome
RELEASE=39
echo "expected $REFGENOME"
FASTAFILE=FungiDB-${RELEASE}_AfumigatusAf293_Genome.fasta
echo "working off $FASTAFILE - check if these don't match may need to update config/init script"

if [ ! -f $FASTAFILE ] ; then
	curl -O http://fungidb.org/common/downloads/release-${RELEASE}/AfumigatusAf293/fasta/data/FungiDB-${RELEASE}_AfumigatusAf293_Genome.fasta
fi

if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.bwt || $FASTAFILE -nt $FASTAFILE.bwt ]]; then
	bwa index $FASTAFILE
fi

DICT=$(basename $FASTAFILE .fasta)".dict"

if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	rm -f $DICT
	picard CreateSequenceDictionary R=$FASTAFILE O=$DICT
	ln -s $DICT $FASTAFILE.dict 
fi

popd
