#!/usr/bin/bash 
#SBATCH --mem=16G --nodes 1 --ntasks 2 --out logs/snpEff_alt_tab.log

SNPEFFOUT=snpEff
SNPEFFGENOME=AfumigatusAf293_FungiDB_39
snpEffConfig=snpEff.config
GFFGENOME=FungiDB-39_AfumigatusAf293.gff
MEM=16g

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
	source config.txt
fi
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA

pushd $SNPEFFOUT
INVCF=$PREFIX.comb_selected.SNP.vcf
OUTVCF=$PREFIX.snpEff.vcf
OUTTAB=$PREFIX.snpEff.tab

source activate gen220

../scripts/map_snpEff2domains.py --vcf $OUTVCF --domains ../genome/FungiDB-39_AfumigatusAf293_InterproDomains.txt --output A_fumigiatus_Af293.Popgen8.snpEf.domain_variants.tsv

../scripts/snpEff_2_tab.py A_fumigiatus_Af293.Popgen8.snpEff.vcf > A_fumigiatus_Af293.Popgen8.snpEff.matrix.tab
