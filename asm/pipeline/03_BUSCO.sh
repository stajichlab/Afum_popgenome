#!/bin/bash
#SBATCH --nodes 1 --ntasks 4 --mem 16G --time 36:00:00 --out logs/busco.%a.log -J busco

module load busco

# for augustus training
export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config

CPU=${SLURM_CPUS_ON_NODE}
N=${SLURM_ARRAY_TASK_ID}
if [ ! $CPU ]; then
     CPU=2
fi

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

if [ -z ${SLURM_ARRAY_JOB_ID} ]; then
	SLURM_ARRAY_JOB_ID=$$
fi
GENOMEFOLDER=genomes
EXT=sorted.fasta
LINEAGE=/srv/projects/db/BUSCO/v9/eurotiomycetes_odb9
OUTFOLDER=BUSCO
TEMP=/scratch/${SLURM_ARRAY_JOB_ID}_${N}
mkdir -p $TEMP
SAMPLEFILE=samples.dat
SEED_SPECIES=aspergillus_fumigatus
sed -n ${N}p $SAMPLEFILE | while read BASE NAME PHYLUM
do
	GENOMEFILE=$(realpath $GENOMEFOLDER/${NAME}.${EXT})
	if [ -d "$OUTFOLDER/run_${NAME}" ];  then
    		echo "Already have run $NAME in folder busco - do you need to delete it to rerun?"
   		exit
	else
    		pushd $OUTFOLDER
    		busco.py -i $GENOMEFILE -l $LINEAGE -o $NAME -m geno --cpu $CPU --tmp $TEMP -sp $SEED_SPECIES --tarzip
    		popd
	fi
	rm -rf $TEMP
done
