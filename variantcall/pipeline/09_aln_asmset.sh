#!/usr/bin/bash
#SBATCH --ntasks 8 -p short --nodes 1 --mem 32G --out logs/aln_asmset.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
MEM=32g
TEMP=/scratch
module load last
module load bcftools/1.9
module load parallel
module unload java
module load java/8
module load picard
module load tabix
module load gatk/3.8
ASMALNFOLDER=aln_asm
ASMFOLDER=input_asm
if [ -f config.txt ]; then
    source config.txt
else
    echo "need a config.txt with REFGENOME defined"
    exit
fi

if [ -z $REFGENOME ]; then
    echo "need a config.txt with REFGENOME defined"
    exit
fi
if [[ -z $ASMFOLDER || ! -d $ASMFOLDER ]]; then
    echo "Need assemblies ending .fasta in ASMFOLDER ($ASMFOLDER)"
    exit
fi

mkdir -p $ASMALNFOLDER/aln $ASMALNFOLDER/gvcf $ASMALNFOLDER/vcf
DB=$(echo $REFGENOME | sed 's/\.fasta/.last/')
VCF=$FINALVCF/$PREFIX.selected_nofixed.SNP.vcf.gz
if [ ! -f $VCF ]; then
    echo "need final vcf in place - run steps 1-5 first"
fi

samtools faidx $REFGENOME

if [ ! -f $DB.prj ]; then
	lastdb -P$CPU -R01 $DB $REFGENOME
fi

for asm in $ASMFOLDER/*.$ASMEXT
do
    BASE=$(basename $asm .fasta)
    MAF=$ASMALNFOLDER/aln/$BASE.maf
    
    if [ ! -s $MAF ]; then
	lastal -P$CPU $DB $asm > $MAF	
    fi
    BAM=$ASMALNFOLDER/aln/$BASE.bam
    TEMPBAM=$TEMP/$BASE.ng.bam
    if [ ! -s $BAM ]; then 
	maf-convert sam $MAF | samtools view -bt $REFGENOME.fai - |
	samtools sort --threads $CPU -o $TEMPBAM
	
	picard AddOrReplaceReadGroups I=$TEMPBAM RGID=$BASE RGSM=$BASE RGLB=$BASE RGCN=NCBI RGPL=Assembly RGPU=$BASE O=$BAM
	samtools index $BAM
    fi
    GVCF=$ASMALNFOLDER/gvcf/$BASE.g.vcf
    if [ ! -f $GVCF.gz ]; then
	if [ ! -f $GVCF ]; then 
	    java -Xmx${MEM} -jar $GATK \
		-T HaplotypeCaller \
		-ERC GVCF \
		-ploidy 1 \
		-I $BAM -R $REFGENOME \
		-L $VCF -rf ReassignMappingQuality -DMQ 60 --defaultBaseQualities 20 \
		-o $GVCF -nct $CPU
	fi    
	bgzip $GVCF
	tabix $GVCF.gz
#	bcftools view -Oz --threads $CPU -o $GVCF.gz $GVCF
#	bcftools index --threads $CPU $GVCF.gz
    fi
done

N=$(ls $ASMALNFOLDER/gvcf/*.g.vcf.gz | sort | perl -p -e 's/\n/ /; s/(\S+)/-V $1/')

OUT=$ASMALNFOLDER/vcf/assemblies.vcf

if [ ! -f $OUT.gz ]; then
 if [ ! -f $OUT ]; then
     java -Xmx$MEM -jar $GATK \
	 -T GenotypeGVCFs \
	 -R $REFGENOME \
	 -L $VCF \
	 $N \
	 -o $OUT 
 fi
 bgzip $OUT
 tabix $OUT.gz
# bcftools view -Oz --threads $CPU -o $OUT.gz $OUT
# bcftools index --threads $CPU $OUT.gz
fi
OUTMERGE=$FINALVCF/${PREFIX}_mergeasm.selected_nofixed.SNP.vcf

if [ ! -f $OUTMERGE ]; then
	bcftools merge -o $OUTMERGE $VCF $OUT.gz
	bgzip $OUTMERGE
fi
