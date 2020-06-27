library(gridExtra)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

MIN_COVERAGE = 0.1
chrlist = c(1:8)

secondaryclusters = "lib/Af293.clusters_names.tab"
clusterNames <- read.table(secondaryclusters,sep=",",header=F)
colnames(clusterNames) = c("GeneLeft","GeneRight")
for (row in clusterNames) {
}
#colnames(clusters) = c("Chr", "Start", "End", "Cluster_locus", "Overlap", "Strand")
#clusters$CHR <- sub("Chr([0-9])_A_fumigatus_Af293", "\\1", clusters$Chr, perl = TRUE)

#nchr = length(unique(clusters$CHR))
#cl = clusters[clusters$CHR %in% chrlist, ]
#cl <- cl[order(cl$CHR, cl$Start), ]
#cl$index = rep.int(seq_along(unique(cl$CHR)), times = tapply(cl$Start, cl$CHR, length))



indepth = "coverage/mosdepth_gene.gg.tab.gz"
bedwindows = read.table(indepth, header = F)
colnames(bedwindows) = c("Chr", "Start", "End", "Gene", "Depth", "Strain")

bedwindows$CHR <- sub("Chr([0-9])_A_fumigatus_Af293", "\\1", bedwindows$Chr, perl = TRUE)

nchr = length(unique(chrlist))
d = bedwindows[bedwindows$CHR %in% chrlist, ]
d <- d[order(d$CHR, d$Start), ]
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start, d$CHR, length))


StrainCount = length(unique(d$Strain))

locations <- d %>% 
  mutate(Chrom = index) %>%
  select(Gene, Chrom, Start, End) %>%
  distinct()

gene_count_locations <- d %>% 
  filter(Depth > MIN_COVERAGE ) %>%
  group_by(Gene) %>%
  summarize( NumMissing = n() ) %>% 
  mutate(GeneRetainFreq = 100 * (NumMissing / StrainCount)) %>%
  inner_join(locations)

clustersL = clusterNames %>% 
  mutate(Gene = GeneLeft) %>%
  left_join(locations,by="Gene") %>%
  mutate(label = "[")

clustersR = clusterNames %>% 
    mutate(Gene = GeneRight) %>%
  left_join(locations,by="Gene") %>%
  mutate(label ="]")

dat_text <- bind_rows(clustersL,clustersR)

# + scale_color_brewer(palette = "RdYlBu", type = "seq") +
#scale_x_continuous(name = "Chromosome",expand = c(0, 0), breaks = ticks, labels = (unique(d$CHR))) + scale_y_continuous(name = "Normalized Read Depth",
#  geom_vline(mapping = NULL, xintercept = minorB, alpha = 0.5, size = 0.1, colour = "grey15") + expand = c(0, 0), limits = c(0, 3))
p <- ggplot(gene_count_locations, aes(x = Start, y = GeneRetainFreq,color=factor(Chrom))) + 
  geom_point(alpha = 0.8, size = 1) + scale_color_brewer(palette = "Set2", type = "seq")  +
  labs(title = sprintf("Gene Retention across %d strains",StrainCount), xlab = "Chrom Position", y = "Gene Retention Percent") + 
  geom_text (
    data    = dat_text,
    size    = 4,
    color   = "black",
    mapping = aes(x = Start, y = 50, label = label),
  ) +
  facet_grid(rows=vars(Chrom)) +
  theme_classic() + guides(fill = guide_legend(keywidth = 3,keyheight = 1))

ggsave(sprintf("plots/GeneRetention.pdf"),p, width = 20, height = 10)

