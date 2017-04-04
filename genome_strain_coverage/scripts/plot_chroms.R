library(gplots)
library(gplots)
library(fastcluster)
library(RColorBrewer)
library(colorRamps)
library(pheatmap)
library(ape)
palette <- colorRampPalette(c('blue','white','red'))(50)
for ( i in c(1:12,14,15,20,32,45) ) {
    filename=sprintf("plot/scf_%06d_A_fumigatus_A1163.gene_cov_norm.tab",i)
    chr <- read.table(filename,header=T,sep="\t",row.names=1)
    chr <- as.matrix(chr)
    ch <- 4
    cw <- 4
    
    fs_row = 5
    fs_col = 5

    pdffile=sprintf("plot/scf_%d.pdf",i)
    pdf(pdffile,height=100,width=10)
    res_t <- pheatmap(chr, fontsize_row = fs_row,
                      fontsize_col = fs_col,
                      cluster_cols = TRUE, cluster_rows = FALSE,
                      col = palette, scale="none",
                      cellheight = ch,
                      cellwidth  = cw,
                      legend = T,main=sprintf("SCF %d plot",i),
                      );

        res_t <- pheatmap(chr, fontsize_row = fs_row,
                      fontsize_col = fs_col,
                      cluster_cols = TRUE, cluster_rows = FALSE,
                      col = palette, scale="column",
                      cellheight = ch,
                      cellwidth  = cw,
                      legend = T,main=sprintf("SCF %d plot",i),
                      );

        res_t <- pheatmap(chr, fontsize_row = fs_row,
                      fontsize_col = fs_col,
                      cluster_cols = TRUE, cluster_rows = FALSE,
                      col = palette, scale="row",
                      cellheight = ch,
                      cellwidth  = cw,
                      legend = T,main=sprintf("SCF %d plot",i),
                      );

}
#heatmap.2(chr4)
