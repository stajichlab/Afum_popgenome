#PBS -l nodes=1:ppn=16,mem=500gb,walltime=168:00:00 -N GATK.GVCFGeno 

module load java
module load gatk/3.4-0
#module load picard

GENOME=genomes/A_fumigatus_A1163.fa
CPU=1
if [ $PBS_NP ]; then
 CPU=$PBS_NP
fi

b=`basename $GENOME .fa`
dir=`dirname $GENOME`

N=`ls project_varfiles/*.g.vcf.gz | perl -p -e 's/\n/ /; s/(\S+)/-V $1/'`
#echo $N
#N="-V project_varfiles/SRP006193.g.vcf.gz"
java -Xmx495g -jar $GATK \
    -T GenotypeGVCFs \
    --max_alternate_alleles 3 \
    -R $GENOME \
    $N \
    -o Afum.A1163.vcf \
    -nt $CPU >& GATK.genotypeall.log
