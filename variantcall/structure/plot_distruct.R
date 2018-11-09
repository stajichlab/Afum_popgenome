library(pophelper)
library(grid)
library(gridExtra)
clist <- list(
"shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
"strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
"oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
"keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
"vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
"muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
"teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
"merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
"funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
"retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
"cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
"cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
"morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
"wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
"krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))

group <- read.csv("Afum.popset_inferred.tab",sep="\t", header=T,stringsAsFactors = F)
onelabsetrep <- group[,2,drop=FALSE]
colnames(group) <- c("Strain","GroupAssign")
inds <- read.table("Afum.Run2.popset",header=FALSE,stringsAsFactors=F)
#inds$V1
ffiles <- list.files(path=".",pattern="*.meanQ",full.names=T)
flist <- readQ(files=ffiles)

rownames(flist[[1]]) <- (inds$V1)
if(length(unique(sapply(flist,nrow)))==1) flist <- lapply(flist,"rownames<-",inds$V1)
# show row names of all runs and all samples
#lapply(flist, rownames)

tr1 <- tabulateQ(qlist=flist)
summariseQ(tr1, writetable=TRUE)
plotQ(flist,
            imgoutput="join",returnplot=T,exportplot=T,linesize=0.8,pointsize=4,
            quiet=T,basesize=8,showlegend=T,sortind="Cluster1",sharedindlab=F,
            clustercol=clist$shiny,splab=paste0("K=",sapply(flist,ncol)),
            outputfilename="Afum.Run2.joinedplot",imgtype="pdf",
            useindlab=T,showindlab=T,
            width=100)


plotQ(flist[c(5)],returnplot=T,exportplot=T,quiet=T,basesize=8,ordergrp=T,
            grplab=onelabsetrep,grplabsize=4,linesize=0.8,pointsize=3,outputfilename="Afum.Run2.K6", imgtype="pdf",
            sharedindlab=F,showindlab=F,width=100)

plotQ(flist[c(4:8)],imgoutput="join",returnplot=T,exportplot=T,quiet=T,basesize=8,ordergrp=T,
      grplab=onelabsetrep,grplabsize=3,linesize=0.8,pointsize=3,outputfilename="Afum.Run2.combined.groups", imgtype="pdf",
      sharedindlab=F,showindlab=F,width=100)


p <- plotQMultiline(flist[2], returnplot=T,spl=100,useindlab=T,grplab=onelabsetrep,grplabsize=3,ordergrp=T,
               imgtype="pdf",exportplot=T,sortind="Cluster1",showlegend=T,
               outputfilename="Afum.Run2.joined_multiline.K_3")
#grid.arrange(p$plot[[1]][[1]])

p <- plotQMultiline(flist[3], returnplot=T,spl=100,useindlab=T,showlegend=T,
               imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
               outputfilename="Afum.Run2.joined_multiline.K_4")
p <- plotQMultiline(flist[4], returnplot=T,spl=100,useindlab=T,showlegend=T,
                    imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
                    outputfilename="Afum.Run2.joined_multiline.K_5")
p <- plotQMultiline(flist[5], returnplot=T,spl=100,useindlab=T,showlegend=T,
                    imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
                    outputfilename="Afum.Run2.joined_multiline.K_6")
p <- plotQMultiline(flist[6], returnplot=T,spl=100,useindlab=T,showlegend=T,
                    imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
                    outputfilename="Afum.Run2.joined_multiline.K_7")
p <- plotQMultiline(flist[7], returnplot=T,spl=100,useindlab=T,showlegend=T,
                    imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
                    outputfilename="Afum.Run2.joined_multiline.K_8")
p <- plotQMultiline(flist[8], returnplot=T,spl=100,useindlab=T,showlegend=T,
                    imgtype="pdf",exportplot=T,sortind="Cluster1",grplab=onelabsetrep,grplabsize=3,ordergrp=T,
                    outputfilename="Afum.Run2.joined_multiline.K_9")

