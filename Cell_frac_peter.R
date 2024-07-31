rm(list=ls())
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(DESeq2)
library(ashr)
#FLEXI count input
dat<-read.table(gzfile("Cellular_frac/raw_counts.csv.gz"),header=T,sep=",",row.names=1)

coldata_DE<- data.frame(Fraction=rep(c("Total", "Nuc", "Cyto"), 8), 
                        Cell=rep(c("Hela","K562","MDA","MCF7"), each=6), 
                        Rep=rep(c(rep("1", 3), rep("2", 3)),4))

coldata_DE$Fraction<- factor(coldata_DE$Fraction)
coldata_DE$Cell<- factor(coldata_DE$Cell)
rownames(coldata_DE)<- colnames(dat)[5:28]

rownames(dat)<- dat$derowname
dat_DE<- dat[dat$Type!="ERCC"&dat$Type!="SP"&dat$Type!="rRNA",c(5:28)]
dat_DE<- dat_DE[,1:24] + 1

#All samples
dds1<- DESeqDataSetFromMatrix(countData= dat_DE,
                              colData= coldata_DE,
                              design= ~ Fraction + Cell )


dds1<- DESeq(dds1, betaPrior = F)

#rld transform
rld1<- rlog(dds1, blind=F)
#PCA
pcaData <- plotPCA(rld1, intgroup=c("Fraction", "Cell"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("Figures/FigS20A.pdf",height=3,width=4.5)
ggplot(pcaData, aes(PC1, PC2, color=Cell, shape=Fraction)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-80, 120)) + 
  scale_y_continuous(limits=c(-80,80)) +
  theme_classic()
dev.off()

##HeLa
coldata_DE_HeLa<- coldata_DE %>% filter(Cell=="Hela")
dat_DE_HeLa<- dat_DE[,1:6]

dds1_HeLa<- DESeqDataSetFromMatrix(countData= dat_DE_HeLa,
                                   colData= coldata_DE_HeLa,
                                   design= ~ Fraction )


dds1_HeLa<- DESeq(dds1_HeLa, betaPrior = F)

#rld transform
rld1_HeLa<- rlog(dds1_HeLa, blind=F)

pdf("Figures/FigS20B_1.pdf",height=2,width=3)
pcaData <- plotPCA(rld1_HeLa, intgroup=c("Fraction"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction)) +
  geom_point(size=2, alpha=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-70, 70)) + 
  scale_y_continuous(limits=c(-70, 70)) +
  ggtitle("Hela")+
  theme_classic()
dev.off()


##K562
coldata_DE_K562<- coldata_DE %>% filter(Cell=="K562")
dat_DE_K562<- dat_DE[,7:12]

dds1_K562<- DESeqDataSetFromMatrix(countData= dat_DE_K562,
                                   colData= coldata_DE_K562,
                                   design= ~ Fraction )


dds1_K562<- DESeq(dds1_K562, betaPrior = F)

#rld transform
rld1_K562<- rlog(dds1_K562, blind=F)

pdf("Figures/FigS20B_2.pdf",height=2,width=3)
pcaData <- plotPCA(rld1_K562, intgroup=c("Fraction"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction)) +
  geom_point(size=2, alpha=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-70, 70)) + 
  scale_y_continuous(limits=c(-70, 70)) +
  ggtitle("K562")+
  theme_classic()
dev.off()

##MDA
coldata_DE_MDA<- coldata_DE %>% filter(Cell=="MDA")
dat_DE_MDA<- dat_DE[,13:18]

dds1_MDA<- DESeqDataSetFromMatrix(countData= dat_DE_MDA,
                                  colData= coldata_DE_MDA,
                                  design= ~ Fraction )


dds1_MDA<- DESeq(dds1_MDA, betaPrior = F)

#rld transform
rld1_MDA<- rlog(dds1_MDA, blind=F)

pdf("Figures/FigS20B_3.pdf",height=2,width=3)
pcaData <- plotPCA(rld1_MDA, intgroup=c("Fraction"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction)) +
  geom_point(size=2, alpha=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-70, 70)) + 
  scale_y_continuous(limits=c(-70, 70)) +
  ggtitle("MDA")+
  theme_classic()
dev.off()

##MCF7
coldata_DE_MCF7<- coldata_DE %>% filter(Cell=="MCF7")
dat_DE_MCF7<- dat_DE[,19:24]

dds1_MCF7<- DESeqDataSetFromMatrix(countData= dat_DE_MCF7,
                                   colData= coldata_DE_MCF7,
                                   design= ~ Fraction )


dds1_MCF7<- DESeq(dds1_MCF7, betaPrior = F)

#rld transform
rld1_MCF7<- rlog(dds1_MCF7, blind=F)

pdf("Figures/FigS20B_4.pdf",height=2,width=3)
pcaData <- plotPCA(rld1_MCF7, intgroup=c("Fraction"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Fraction)) +
  geom_point(size=2, alpha=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_x_continuous(limits=c(-70, 70)) + 
  scale_y_continuous(limits=c(-70, 70)) +
  ggtitle("MCF7")+
  theme_classic()
dev.off()

#Z_prep
gldata0<- assay(rld1)
gldata0<-gldata0[,c(2,5,1,4,3,6,8,11,7,10,9,12,14,17,13,16,15,18,20,23,19,22,21,24)]
Z0<- t(scale(t(gldata0)))

#Samples distance
sampleDists <- dist(t(gldata0))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("Figures/FigS20C.pdf",height=3.5,width=5)
pheatmap(sampleDistMatrix,
         cluster_rows=F,
         cluster_cols=F,
         #clustering_distance_rows=sampleDists,
         #clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
