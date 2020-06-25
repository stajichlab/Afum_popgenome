library(ggplot2)
library(RColorBrewer)
library(vcfR)

vcfile = "vcf/AF100.Run2.2019-01-16.selected_nofixed.SNP.vcf.gz"
args = commandArgs(trailingOnly=TRUE)
if (length(args)> 0) {
    vcfile = args[1]
    print(args[1])
}
gff_file = "genome/FungiDB-39_AfumigatusAf293.gff"
dna_file = "genome/FungiDB-39_AfumigatusAf293_Genome.simpleid.fa"
vcf <- read.vcfR(vcfile, verbose = FALSE )
gff <- read.table(gff_file, sep="\t", quote="")
colnames(gff) <- c("CHROM","SOURCE","TYPE","START",
                   "END","SCORE","STRAND","PHASE","DESCRIPTION")
gff <- subset(gff,gff$TYPE == "gene")
dna <- ape::read.dna(dna_file, format = "fasta")

chrnames = names(dna)
nchr = length(chrnames)
pdf("SNP_qc.pdf")
for (i in 1:nchr ) {
    name = chrnames[i]
    print(name)
    chrom <- create.chromR(name=name, vcf=vcf[name], seq=dna[name],
                           ann=subset(gff,gff$CHROM == name))
    chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700,
                    min_MQ = 50,  max_MQ = 70)
    plot(chrom)
    chrom <- proc.chromR(chrom, verbose=TRUE)
    plot(chrom)

    chromoqc(chrom, dp.alpha=20)
}
