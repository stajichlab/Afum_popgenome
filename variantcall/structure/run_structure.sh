#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 16gb --out structure.%a.log
module load faststructure
K=${SLURM_ARRAY_TASK_ID}
structure.py -K $K --input=AfumAf293.Run2 --output=AfumAf293.Run2
