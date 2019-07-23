n=0
cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -f genomes/$strain.sorted.fasta ]; then   echo "$n $strain"; fi; done 
n=0
m=$(cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -f genomes/$strain.sorted.fasta ]; then   echo "$n"; fi; done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')
echo "sbatch --array=$m pipeline/AAFTF.sh"
