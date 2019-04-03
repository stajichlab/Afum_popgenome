#!/usr/bin/env python
import sys, warnings
import vcf,re
from Bio import SeqIO
# this script will convert a VCF from snpEff into a useful table

if len(sys.argv) < 2:
        warnings.warn("Usage snpEff_to_table.py snpeff.vcf genome")
        sys.exit()

newfilename = sys.argv[1].strip('vcf')
genome = "/bigdata/stajichlab/shared/projects/Afum_popgenome/variantcall/genome/FungiDB-39_AfumigatusAf293_Genome.fasta"

if len(sys.argv) > 2:
        genome = sys.argv[2]
count=0
vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))


chrs={}
for seq in SeqIO.parse(genome, "fasta"):
	chrs[seq.id] = seq.seq

for record in vcf_reader:
        if count == 0:
                sampname = []
                title = ["CHROM","POS","FLANKING","TYPE","IMPACT","GENE",
                         "CHANGEDNA","CHANGEPEP","REF","ALT"]
                for sample in record.samples:
                        title.append(str(sample.sample))
                title.append("ANN")
                print("\t".join(title))
                count = 1
        if 'ANN' not in record.INFO:
           sys.stderr.write("Cannot find ANN in %s\n"%(record))
           continue
        anns = record.INFO['ANN']
        arrayout = [record.CHROM,record.POS]
        flanking_seq = chrs[record.CHROM][record.POS-10:record.POS+10]
        arrayout.append(flanking_seq)
        annarr = anns[0].split('|')
        dnachg = re.sub("c.","",annarr[9])
        if ( annarr[1] == 'upstream_gene_variant' or
             annarr[1] == 'downstream_gene_variant' or
             annarr[1] == 'intergenic_region'):
                arrayout.extend(('intergenic',annarr[2],annarr[3],
                         dnachg,""))
        else:
                pepchg = re.sub('p.','',annarr[10])
                arrayout.extend((annarr[1],annarr[2],annarr[3],
                                 dnachg,pepchg))

        arrayout.extend((record.REF,record.ALT))
	#print arrayout
        for sample in record.samples:
                if sample.gt_bases:
                        arrayout.append(sample.gt_bases+"/")
                else:
                        arrayout.append('./')
        for ann in anns:
                annarr=ann.split('|')
                readableann= "(%s,%s,%s)" % (annarr[3],annarr[0],annarr[1])
                arrayout.append(readableann)

        print('\t'.join(map(str,arrayout)))
