n=0
cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -f BUSCO/run_${strain}/single_copy_busco_sequences.tar.gz ]; then   echo "$n $strain"; fi; done 
n=0
m=$(cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -s BUSCO/run_${strain}/single_copy_busco_sequences.tar.gz ]; then   echo "$n"; fi; done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')
echo "sbatch --array=$m pipeline/03_BUSCO.sh"
