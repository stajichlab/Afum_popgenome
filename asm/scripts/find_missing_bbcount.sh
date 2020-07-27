n=0
cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -f mapping_report/$strain.bbmap_summary.txt ]; then   echo "$n $strain"; fi; done 
n=0
m=$(cut -f2 samples.dat | while read strain ; do n=$(expr $n + 1); if [ ! -f mapping_report/$strain.bbmap_summary.txt ]; then   echo "$n"; fi; done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')
echo "sbatch --array=$m -p batch --time 48:00:00 -n 16 -N 1 pipeline/02_read_count.sh"
