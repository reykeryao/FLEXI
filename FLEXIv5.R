rm(list=ls())
library(tidyverse)
library(ggplotify)
library(matrixStats)
library(tidyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggpubr)
library(grid)
library(Rtsne)
library(zinbwave)
library(gridExtra)
library(cowplot)
library(ggfortify)
library(gridGraphics)
library(VennDiagram)
library(ComplexHeatmap)
library(UpSetR)
library(readxl)
library(cluster)
library(dendextend)
library(reshape2)
library(purrr)
library(dplyr)
library(fpc)
library(circlize)
library(cocor)
library(openxlsx)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
set.seed(740714)

### read FLEXIs from all Bio1 dataset
dat<-read.delim("all.FLEXI")

## mapped reads, from TableS1, columns 4; and Summary_dataset.xlsx
mapped_reads<-c(207.491024,692.091831,666.341854,713.775291,715.241521,768.433748,71.116246)
names(mapped_reads)<-colnames(dat[,71:77])
GRCh38<-read.delim("GRCh38.93.intron_deduped.tsv")

#Fig1A, upset plot of FLEXIs, FLEXI host genes from cellular RNA and plasma
#FLEXIs
cut_off=1
pdf("Figures/Fig1A_1.pdf",height=4,width=8)
set_1 <- as.character(dat$ID[dat$UHRR>=cut_off])
set_2 <- as.character(dat$ID[dat$K562>=cut_off])
set_3 <- as.character(dat$ID[dat$HEK>=cut_off])
set_4 <- as.character(dat$ID[dat$Hela>=cut_off])
set_5 <- as.character(dat$ID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK-293T"=set_3,
             "HeLa S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
ss<-set_size(m)
cs=comb_size(m)
UpSet(m,set_order=c("K-562","HEK-293T","HeLa S3","UHRR","Plasma"),
          comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
          column_title="FLEXI RNAs",
          top_annotation = HeatmapAnnotation(
            "Counts" = anno_barplot(cs, 
                                    ylim = c(0, max(cs)*1.1),
                                    border = F,
                                    gp = gpar(border =NA,lty=0,
                                              fill =c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)]), 
                                    height = unit(5, "cm")), 
            annotation_name_side = "left", 
            annotation_name_rot = 90),
          right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {grid.text(cs[1:2], x = 1:2, y = unit(cs[1:2], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6), rot = 45)})
dev.off()
#FLEXI host genes
pdf("Figures/Fig1A_2.pdf",height=4,width=8)
set_1 <- unique(dat$GID[dat$UHRR>=cut_off])
set_2 <- unique(dat$GID[dat$K562>=cut_off])
set_3 <- unique(dat$GID[dat$HEK>=cut_off])
set_4 <- unique(dat$GID[dat$Hela>=cut_off])
set_5 <- unique(dat$GID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK-293T"=set_3,
             "HeLa S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
ss<-set_size(m)
cs=comb_size(m)
UpSet(m,set_order=c("K-562","HEK-293T","HeLa S3","UHRR","Plasma"),
          comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
          column_title="FLEXI host genes",
          top_annotation = HeatmapAnnotation(
            "Counts" = anno_barplot(cs, ylim = c(0, max(cs)*1.1),border = F,
                                    gp = gpar(border =NA,lty=0,
                                              fill =c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)]), 
                                    height = unit(5, "cm")), 
            annotation_name_side = "left", 
            annotation_name_rot = 90),
          right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {
  grid.text(cs[1:2], x = 1:2, y = unit(cs[1:2], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6), rot = 45)})
dev.off()
#remove temporary objects
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","m","cut_off","ss","cs"))
# finish of Fig1A

# Fig 1B: table of mirtron, agotron, etc
# generate the table for Fig1B

'''
Fig1B<-rbind(
  "K-562"=c(length(dat$ID[dat$K562>0]),
            length(dat$ID[dat$K562>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[dat$K562>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[dat$K562>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[dat$K562>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[dat$K562>0 & dat$Has_snoRNA!="."]),
            length(grep ("SCARNA",dat$Has_snoRNA[dat$K562>0 &dat$Has_snoRNA!="."]))),
  "HEK-293T"=c(length(dat$ID[dat$HEK>0]),
            length(dat$ID[dat$HEK>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[dat$HEK>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[dat$HEK>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[dat$HEK>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[dat$HEK>0 & dat$Has_snoRNA!="."]),
            length(grep ("SCARNA",dat$Has_snoRNA[dat$HEK>0 &dat$Has_snoRNA!="."]))),
  "HeLa S3"=c(length(dat$ID[dat$Hela>0]),
            length(dat$ID[dat$Hela>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[dat$Hela>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[dat$Hela>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[dat$Hela>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[dat$Hela>0 & dat$Has_snoRNA!="."]),
            length(grep ("SCARNA",dat$Has_snoRNA[dat$Hela>0 &dat$Has_snoRNA!="."]))),
  "UHRR"=c(length(dat$ID[dat$UHRR>0]),
            length(dat$ID[dat$UHRR>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[dat$UHRR>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[dat$UHRR>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[dat$UHRR>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[dat$UHRR>0 & dat$Has_snoRNA!="."]),
           length(grep ("SCARNA",dat$Has_snoRNA[dat$UHRR>0 &dat$Has_snoRNA!="."]))),
  "Plasma"=c(length(dat$ID[dat$Plasma>0]),
            length(dat$ID[dat$Plasma>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[dat$Plasma>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[dat$Plasma>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[dat$Plasma>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[dat$Plasma>0 & dat$Has_snoRNA!="."]),
            length(grep ("SCARNA",dat$Has_snoRNA[dat$Plasma>0 &dat$Has_snoRNA!="."]))),
  "Total"=c(length(dat$ID[rowSums(dat[,73:77])>0]),
            length(dat$ID[rowSums(dat[,73:77])>0 & dat$Is_agotron!="." & dat$Is_mirtron=="."]),
            length(dat$ID[rowSums(dat[,73:77])>0 & dat$Is_mirtron!="." & dat$Is_agotron=="."]),
            length(dat$ID[rowSums(dat[,73:77])>0 & dat$Is_agotron!="." & dat$Is_mirtron!="."]),
            length(dat$ID[rowSums(dat[,73:77])>0 & (dat$Is_agotron!="." | dat$Is_mirtron!=".")]),
            length(dat$ID[rowSums(dat[,73:77])>0 & dat$Has_snoRNA!="."]),
            length(grep ("SCARNA",dat$Has_snoRNA[rowSums(dat[,73:77])>0 &dat$Has_snoRNA!="."]))),
  "GRCh38"=c(length(GRCh38$ID),
            length(GRCh38$ID[GRCh38$Is_agotron!="." & GRCh38$Is_mirtron=="."]),
            length(GRCh38$ID[GRCh38$Is_mirtron!="." & GRCh38$Is_agotron=="."]),
            length(GRCh38$ID[GRCh38$Is_agotron!="." & GRCh38$Is_mirtron!="."]),
            length(GRCh38$ID[GRCh38$Is_agotron!="." | GRCh38$Is_mirtron!="."]),
            length(GRCh38$ID[GRCh38$Has_snoRNA!="."]),
            length(grep ("SCARNA",GRCh38$Has_snoRNA)))
)
colnames(Fig1B)<-c("FLEXI RNA","Agotron only","Mirtron only","Agotron and Mirtron","Agotron or Mirtron","Embedded snoRNA/scaRNA","scaRNA")
Fig1B<-data.frame(Fig1B)
Fig1B$Dataset<-rownames(Fig1B)
Fig1B$Embedded.snoRNA.scaRNA<-paste0(Fig1B$Embedded.snoRNA.scaRNA,"/",Fig1B$scaRNA," (",
                                     formatC(100*Fig1B$Embedded.snoRNA.scaRNA/Fig1B$FLEXI.RNA,1,format="f"),"%)")
Fig1B$Agotron.or.Mirtron<-paste0(Fig1B$Agotron.or.Mirtron," (",
                                     formatC(100*Fig1B$Agotron.or.Mirtron/Fig1B$FLEXI.RNA,1,format="f"),"%)")
Fig1B$FLEXI.RNA<-formatC(Fig1B$FLEXI.RNA,big.mark = ",")
Fig1B<-Fig1B[,c(8,1:6)]
write.table(Fig1B, "Fig1B_table.tsv",quote=F,sep="\t",row.names=F)
'''
#make Fig 1B
Fig1B<-read.delim(gzfile("Fig1B_table.tsv.gz"))
pdf("Figures/Fig1B.pdf")
tt <- ttheme_default(base_size=8,colhead=list(fg_params = list(parse=TRUE,fontface="plain")))
colnames(Fig1B)<-c("Dataset","FLEXI \nRNAs","Agotron\nonly","Mirtron\nonly","Agotron\nand\nMirtron",
                   "Total\nAgotron/\nMirtron","Embedded\nsnoRNA/\nscaRNA")
grid.table(Fig1B,theme=tt,rows=NULL)
dev.off()
#remove temporary objects
rm("Fig1B","tt")


#Fig2D, weblogo of 5'SS, 3'SS, and BPA consensus
#in Fig2D folder

# Processing circular intron consensus
# from https://elifesciences.org/articles/47809, Fig6
# and from RNAS https://www.pnas.org/doi/abs/10.1073/pnas.1808816115
CI<-read_xlsx("combined_sisRNA.xlsx")
# get 10 nt frequencies at 5SS
for (i in 1:10){
  if (i==1){
    Frq_table<-data.frame(table(substr(CI$`5'ss`,i,i)))
  } else {
    Frq_table<-merge(Frq_table,data.frame(table(substr(CI$`5'ss`,i,i))),by=1,all=T)
  }
}
Frq_table[is.na(Frq_table)]<-0
Frq_table<-t(Frq_table[,-1])
## results copied to txt format (CHEN format) then convert to MEME motif format
## using chen2meme < txt file > motif file
## for weblogo EPS use: ceqlogo -i motif file -m motif name -o motif eps file
## 3SS
for (i in 1:10){
  if (i==1){
    Frq_table<-data.frame(table(substr(CI$`3'ss`,i,i)))
  } else {
    Frq_table<-merge(Frq_table,data.frame(table(substr(CI$`3'ss`,i,i))),by=1,all=T)
  }
}
Frq_table[is.na(Frq_table)]<-0
Frq_table<-t(Frq_table[,-1])
## BPA
## 5 nt up/3 nt down of the BPA
for (i in 1:10){
  if (i==1){
    Frq_table<-data.frame(table(substr(CI$BP_10,i,i)))
  } else {
    Frq_table<-merge(Frq_table,data.frame(table(substr(CI$BP_10,i,i))),
                     by=1,all=T)
  }
}
Frq_table[is.na(Frq_table)]<-0
Frq_table<-t(Frq_table[,-1])
rm ("Frq_table","CI","i")
### saved as motif and then convert to seqlogo

#Fig1C-E, density plot of length, GC%, MFE (kcal/mol)
### with cutoff of 0.002 RPM above and under
pdf("Figures/Fig1CtoE.pdf",width=12,height=12)
col=c("black","red","skyblue")
#length
pdf(NULL)
dev.control(displaylist="enable")
dat1<-dat[rowSums(dat[,73:77])>0,]
dat1$RPM<-(rowSums(dat1[,73:76])/sum(mapped_reads[3:6])>=0.002)
plot(density(dat1[dat1$RPM,"Len"]),bty="n",xlim=c(0,350),lwd=1.5,
     ylim=c(0,0.04),main=NA,xlab="Intron length (bp)",col=col[1])
lines(density(dat1[!dat1$RPM,"Len"]),xlim=c(0,350),lwd=1.5,lty=2,col=col[1])
lines(density(dat1[dat1$Plasma>=1,"Len"]),xlim=c(0,350),lwd=1.5,col=col[2])
lines(density(GRCh38$Len[!GRCh38$ID%in%dat1$ID]),xlim=c(0,350),lwd=1.5,col=col[3])
#lines(density(GRCh38$Len),xlim=c(0,350),lwd=1.5,col=col[6])
legend(150,0.03,lty=c(1,2,1,1),lwd=1.5,col=col[c(1,1,2,3)],
       legend = c("FLEXIs (Cells,>=0.002 RPM)","FLEXIs (Cells,<0.002 RPM)","FLEXIs (Plasma)","Other short introns"),bty="n")
length.line <- recordPlot()
invisible(dev.off())

#GC
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat1[dat1$RPM,"GC"]),bty="n",xlim=c(0,100),lwd=1.5,
     ylim=c(0,0.08),main=NA,xlab="GC (%)",col=col[1])
lines(density(dat1[!dat1$RPM,"GC"]),xlim=c(0,100),lwd=1.5,lty=2,col=col[1])
lines(density(dat1[dat1$Plasma>=1,"GC"]),xlim=c(0,100),lwd=1.5,col=col[2])
lines(density(GRCh38$GC[!GRCh38$ID%in%dat1$ID]),xlim=c(0,100),lwd=1.5,col=col[3])
GC.line <- recordPlot()
invisible(dev.off())

#MEF
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat1[dat1$RPM,"MFE"]),bty="n",xlim=c(-150,0),lwd=1.5,
     ylim=c(0,0.05),main=NA,xlab="MFE (kcal/mol)",col=col[1])
lines(density(dat1[!dat1$RPM,"MFE"]),xlim=c(-150,0),lwd=1.5,col=col[1],lty=2)
lines(density(dat1[dat1$Plasma>=1,"MFE"]),xlim=c(-150,0),lwd=1.5,col=col[2])
lines(density(GRCh38$MFE[!GRCh38$ID%in%dat1$ID]),xlim=c(-150,0),lwd=1.5,col=col[3])
MFE.line <- recordPlot()
invisible(dev.off())

ggarrange(ggarrange(length.line,GC.line,MFE.line,ncol=2,nrow = 2))
dev.off()
#remove objects
rm(list=c("GC.line","length.line","MFE.line","col","dat1"))

#Fig3A density plot of RPM with sncRNA markers
##New ID with scale of sncRNA
##Fig3A, RPM density of FLEXIs, different group
gene_counts<-read.delim(gzfile("combined_counts.tsv.gz"))
gene_counts[gene_counts$Type=="scaRNA",3]<-"snoRNA"
snoRNA<-gene_counts[gene_counts$Type=="snoRNA",c(2,6:10)]
snoRNA<-separate(snoRNA,col = "Name",into = "Name",sep = "-",remove = T)
snoRNA$Name<-gsub(snoRNA$Name,pattern = "[A-Z]$|P[1-9]$",replacement="")
snoRNA$Name[snoRNA$Name=="U3"]<-"SNORD3"
snoRNA$Name[snoRNA$Name=="U8"]<-"SNORD118"
snoRNA$Name[snoRNA$Name=="snoU13"]<-"SNORD13"
snoRNA$Name[snoRNA$Name=="snoU2_19"]<-"snoU2"
snoRNA<-aggregate(.~Name,data=snoRNA,FUN = sum)
snRNA<-gene_counts[gene_counts$Type=="snRNA",c(2,6:10)]
U7<-log10(colSums(snRNA[grep("U7\\b|U7[A-Z]",snRNA$Name),2:6])/mapped_reads[3:7])
U11<-log10(colSums(snRNA[grep("U11\\b|U11[A-Z]",snRNA$Name),2:6])/mapped_reads[3:7])
SNORD74<-log10(colSums(snoRNA[snoRNA$Name=="SNORD74",2:6])/mapped_reads[3:7])
SNORD78<-log10(colSums(snoRNA[snoRNA$Name=="SNORD78",2:6])/mapped_reads[3:7])
YRNA<-log10(colSums(gene_counts[gene_counts$Type=="YRNA",6:10])/mapped_reads[3:7])
VTRNA<-log10(colSums(gene_counts[gene_counts$Type=="VTRNA",6:10])/mapped_reads[3:7])
U6<-log10(colSums(snRNA[grep("U6\\b|U6[A-Z]",snRNA$Name),2:6])/mapped_reads[3:7])
U12<-log10(colSums(snRNA[grep("U12\\b|U12[A-Z]",snRNA$Name),2:6])/mapped_reads[3:7])
U1<-log10(colSums(snRNA[grep("U1\\b|U1[A-Z]",snRNA$Name),2:6])/mapped_reads[3:7])
snc<-data.frame(rbind(SNORD74,SNORD78,U11,U7,U12,U6,YRNA,U1,VTRNA))
title_name<-c("UHRR","K-562","HEK-293T","HeLa S3")
pdf("Figures//Fig3A.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(2,2,2,2,2,2)
for (i in c(74:76,73)){
  plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-70])),bty="n",xlab="RPM",
       xlim=c(-4,6),ylim=c(0,D_height[i-72]),main=title_name[i-72],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-70])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-70])),col="black")
  lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-70])),col="goldenrod")
  lines(density(log10(snoRNA[snoRNA[,i-71]>0,i-71]/mapped_reads[i-70])),col="goldenrod",lty=4)
  abline(v=log10(0.002),lty=2,col="red")
  snc<-snc[order(snc[,i-72],decreasing = F),]
  for (j in (1:dim(snc)[1])){
    if (j==1){
      segments(x0=snc[j,i-72],y0=1.5,y1=1.8,lty=1)
      text (snc[j,i-72],1.85,rownames(snc)[j],cex=0.7,adj=0)
    } else {
      segments(x0=snc[j,i-72],y0=0.5,y1=1.8-j*0.12,lty=1)
      text (snc[j,i-72],1.85-j*0.12,rownames(snc)[j],cex=0.7,adj=0)
    }
  }
  if (i==74){
    legend(2,2,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5,cex=0.7)
  }
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-2),1,bquote(10^2),bquote(10^4),
                  bquote(10^6)),at=seq(-4,6,2))
}
dev.off()

#remove objects
rm(list=c("snoRNA","snRNA","i","U7","U11","SNORD74","SNORD78","U1","U6","U12",
          "YRNA","VTRNA","snc","j","D_height","title_name"))

#Fig3C, distribution of phastCon in FLEXI subsets.
pdf("Figures/Fig3C.pdf")
scol=c("deepskyblue2","firebrick2","goldenrod","black")
plot(density(dat$PhastCon30[dat$Is_mirtron!="."]),ylim=c(0,7),
     col=scol[1],bty="n",xlab="phastCons",main=NA,xlim=c(0,1))
lines(density(dat$PhastCon30[dat$Is_agotron!="."]),col=scol[2])
lines(density(dat$PhastCon30[dat$Has_snoRNA!="."]),col=scol[3])
lines(density(dat$PhastCon30[dat$Is_agotron=="." & dat$Has_snoRNA=="." & dat$Is_mirtron=="."]),
      col=scol[4])
lines(density(GRCh38$PhastCon30[!GRCh38$ID%in%dat$ID]),col="gray50",lty=2)
legend(0.6,5,col=c(scol[c(4,2,1,3)],"gray50"),
       legend = c("All other FLEXIs","Agotron FLEXIs","Mirtron FLEXIs","SnoRNA FLEXI",
                  "Other short introns"),lty=c(1,1,1,1,2),bty="n")
dev.off()
#remove objects
rm(scol)

#Fig2A IGV screen shots in folder
#FigS8 IGV screen shots in folder

#Fig2B histogram of relative length of the intronic reads
file_name<-c("K562","HEK","Hela","UHRR")
for (i in 1:4){
  per_dat<-read.delim(gzfile(paste0(file_name[i],".per.info.gz")))
  if (i==1){
   combined<-per_dat
  } else {
    combined<-rbind(combined,per_dat)
  }
}

p1 <- hist(combined$Per[combined$nonFLEXI!=1]) 
p2 <- hist(combined$Per[combined$nonFLEXI==1])   
pdf("Figures/Fig2B.pdf")
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,100),main=NA,xlab="Fragment length (%)",ylab="Reads")
plot( p2, col=rgb(0,0,1,1/4), xlim=c(0,100), add=T)
legend(0,140000,legend = c("FLEXIs","Other short introns"),bty="n",
       fill =c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),pch=15))
dev.off()
#remove objects
rm(list=c("combined","per_dat","i","name","file_name","p1","p2","col","y_max"))

#Fig2C and S3B are excel spread sheet and their figs
# Figures/Fig2C_S3B.xlsx

#Fig3B
FLEXI_CPM<-dat[,c(1,8,74:76)]
rownames(FLEXI_CPM)<-FLEXI_CPM$ID
Unfrag_total<-mapped_reads[4:6]
FLEXI_CPM[,3:5]<-t(t(FLEXI_CPM[,3:5])/Unfrag_total)
FLEXI_CPM$K562_repo<-apply(dat[,45:52],1,FUN=function(x){sum(x>0)})
FLEXI_CPM$HEK_repo<-apply(dat[,53:60],1,FUN=function(x){sum(x>0)})
FLEXI_CPM$Hela_repo<-apply(dat[,61:70],1,FUN=function(x){sum(x>0)})
FLEXI_CPM<-FLEXI_CPM[rowSums(FLEXI_CPM[,3:5])>0,]
temp<-FLEXI_CPM[,3:5]
temp[temp==0]<-2^-10
FLEXI_CPM[,3:5]<-temp
FLEXI_CPM<-FLEXI_CPM[,2:8]
Cell_counts<-read.delim("combined_run.counts")
Cell_counts$K562<-rowSums(Cell_counts[,c(22:29)])
Cell_counts$HEK<-rowSums(Cell_counts[,c(30:37)])
Cell_counts$Hela<-rowSums(Cell_counts[,c(38:47)])
rownames(Cell_counts)<-Cell_counts$ID
Cell_counts$ID<-as.character(Cell_counts$ID)
#Subset only FLEXI host genes
GID_list<-unique(FLEXI_CPM$GID)
Cell_counts<-Cell_counts[Cell_counts$ID%in%GID_list,c(1,48:50)]
Cell_counts[,2:4]<-t(t(Cell_counts[,2:4])/Unfrag_total)
#Assign log2CPM of -10 for those with 0 count
Cell_counts[Cell_counts==0]<-2^-10

### make downsampled reads
'''
Cell_unfrag_counts<-read.delim("combined_run.counts")
Cell_unfrag_counts$K562<-rowSums(Cell_unfrag_counts[,c(22:29)])
Cell_unfrag_counts$HEK<-rowSums(Cell_unfrag_counts[,c(30:37)])
Cell_unfrag_counts$Hela<-rowSums(Cell_unfrag_counts[,c(38:47)])
rownames(Cell_unfrag_counts)<-Cell_unfrag_counts$ID
Cell_unfrag_counts$ID<-as.character(Cell_unfrag_counts$ID)
Unfrag_total<-colSums(Cell_unfrag_counts[,48:50])
#Subset only FLEXI host genes
GID_list<-unique(FLEXI_CPM$GID)
Cell_counts<-Cell_unfrag_counts[Cell_unfrag_counts$ID%in%GID_list,c(1,48:50)]
FLEXI_host_total<-colSums(Cell_counts[,2:4])
Sampling_ratio<-FLEXI_host_total/FLEXI_total
Sampling_size<-Unfrag_total/Sampling_ratio
#sampling based on ratio of FLEXI counts to FLEXI host gene counts on the whole dataset
for (i in 1:3){
  tmp<-rep(Cell_unfrag_counts$ID,Cell_unfrag_counts[,47+i])
  tmp<-sample(tmp)
  tmp<-sample(tmp,Sampling_size[i],replace = T)
  tmp<-data.frame(table(tmp))
  Cell_unfrag_counts<-merge(Cell_unfrag_counts,tmp,by=1,all=T)
}
Cell_unfrag_counts[is.na(Cell_unfrag_counts)]<-0
GID_list<-unique(FLEXI_CPM$GID)
Cell_counts<-Cell_unfrag_counts[Cell_unfrag_counts$ID%in%GID_list,c(1,51:53)]
Cell_counts[,2:4]<-1e6*t(t(Cell_counts[,2:4])/Unfrag_total)
#Assign log2CPM of -10 for those with 0 count
Cell_counts[Cell_counts==0]<-2^-10
colnames(Cell_counts)[2:4]<-colnames(Cell_unfrag_counts)[48:50]
write.table(Cell_counts,"Fig3A_downsampled.counts",quote=F,sep="\t",row.names=F)
'''
Cell_counts_down<-read.delim(gzfile("Fig3A_downsampled.counts.gz"))

### add FLEXI marker for label
#pick<-c("23I_DOCK6___ENSG00000130158___ENST00000587656___protein_coding___protein_coding",
#        "11I_SGK1___ENSG00000118515___ENST00000237305___protein_coding___protein_coding",
#        "1I_FTH1___ENSG00000167996___ENST00000532601___protein_coding___protein_coding",
#        "11I_THBS1___ENSG00000137801___ENST00000260356___protein_coding___protein_coding",
#        "1I_RPS2___ENSG00000140988___ENST00000526586___protein_coding___protein_coding",
#        "2I_ACTG1___ENSG00000184009___ENST00000575842___protein_coding___protein_coding")

pdf("Figures/Fig3B.pdf",width = 9,height=9)
par(pch=16,mfcol=c(3,3),pty="s")
#FLEXIs scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(2,4)]),col="#FF000050")
points(log2(sig.x[,c(2,4)]),col="#0000FF50")
### identify the most extreme FLEXIs near the axis
top2y<-sig.y %>%
  filter(K562==2^-10) %>%
  arrange(Hela) 
top2y<-rev(rownames(top2y))[1:2]
top2x<-sig.x %>%
  filter(Hela==2^-10) %>%
  arrange(K562) 
top2x<-rev(rownames(top2x))[1:2]
lab_tmp<-unique(rbind(sig.x[rownames(sig.x)%in%top2x,],sig.y[rownames(sig.y)%in%top2y,]))
if (dim(lab_tmp)[1]>0){
  pick_GID<-data.frame("ID"=unique(str_split_i(rownames(lab_tmp),"___",2)),
                       "Name"=unique(str_split_i(rownames(lab_tmp),"_",2)))
  lab_tmp$IID<-rownames(lab_tmp)
  lab_tmp<-separate(lab_tmp,"IID",into="IID",sep="___",extra="drop")
  lab_tmp<-separate(lab_tmp,"IID",into=c("II","Gname"),sep="_",remove = F)
  text(log2(lab_tmp[,c(2,4)]),labels = lab_tmp$IID,cex=0.25)
}
sig<-merge(sig.x,sig.y,by=1:7,all=T)
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,4],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#ccor list for test
cocor_dat<-list(F=data.frame("F1"=FLEXI_by_GID$K562,"F2"=FLEXI_by_GID$Hela))
#Unfrag scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,4)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(2,4)]),labels = lab_tmp$Name,cex=0.25)
}
sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%(union(sig.x,sig.y)),]
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(10,-6,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(10,-7.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,4],method = "spearman"),digits=2, format="f")
text(10,-7.5,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(10,-10.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$K562,"G2"=FLEXI_by_GID$Hela)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#p<0.0001

### downsampled
#Unfrag scatter
FLEXI_by_GID<-Cell_counts_down[!(Cell_counts_down[,2]==2^-10 & Cell_counts_down[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,4)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(2,4)]),labels = lab_tmp$Name,cex=0.25)
}
sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%union(sig.x,sig.y),]
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,4],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$K562,"G2"=FLEXI_by_GID$Hela)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#p<0.0001

#FLEXIs scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,3]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
sig.x<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(3,4)]),col="#FF000080")
points(log2(sig.x[,c(3,4)]),col="#0000FF80")
top2y<-sig.y %>%
  filter(HEK==2^-10) %>%
  arrange(Hela) 
top2y<-rev(rownames(top2y))[1:2]
top2x<-sig.x %>%
  filter(Hela==2^-10) %>%
  arrange(HEK) 
top2x<-rev(rownames(top2x))[1:2]
lab_tmp<-unique(rbind(sig.x[rownames(sig.x)%in%top2x,],sig.y[rownames(sig.y)%in%top2y,]))
if (dim(lab_tmp)[1]>0){
  pick_GID<-data.frame("ID"=unique(str_split_i(rownames(lab_tmp),"___",2)),
                       "Name"=unique(str_split_i(rownames(lab_tmp),"_",2)))
  lab_tmp$IID<-rownames(lab_tmp)
  lab_tmp<-separate(lab_tmp,"IID",into="IID",sep="___",extra="drop")
  lab_tmp<-separate(lab_tmp,"IID",into=c("II","Gname"),sep="_",remove = F)
  text(log2(lab_tmp[,c(3,4)]),labels = lab_tmp$IID,cex=0.25)
}

sig<-merge(sig.x,sig.y,by=1:7,all=T)
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,3],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,3],sig[,4],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#ccor list for test
cocor_dat<-list(F=data.frame("F1"=FLEXI_by_GID$HEK,"F2"=FLEXI_by_GID$Hela))

#Unfrag scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,3]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(3,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(3,4)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(3,4)]),labels = lab_tmp$Name,cex=0.25)
}
sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%(union(sig.x,sig.y)),]
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(10,-6,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(10,-7.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,3],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,3],sig[,4],method = "spearman"),digits=2, format="f")
text(10,-7.5,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(10,-10.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$HEK,"G2"=FLEXI_by_GID$Hela)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#p<0.0001

##UF downsampled
#Unfrag scatter
FLEXI_by_GID<-Cell_counts_down[!(Cell_counts_down[,3]==2^-10 & Cell_counts_down[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(3,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(3,4)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(3,4)]),labels = lab_tmp$Name,cex=0.25)
}
sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%union(sig.x,sig.y),]
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,3],sig[,4],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,3],sig[,4],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$HEK,"G2"=FLEXI_by_GID$Hela)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#p<0.0001

#FELXIs between K562 and HEK
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
points(log2(sig.y[,c(2,3)]),col="#FF000080")
points(log2(sig.x[,c(2,3)]),col="#0000FF80")
top2y<-sig.y
top2y$ratio<-top2y$HEK/top2y$K562
top2y<-top2y[order(top2y$ratio,decreasing = T),]
top2y<-rownames(top2y)[1:2]
top2x<-sig.x
top2x$ratio<-top2x$K562/top2x$HEK
top2x<-top2x[order(top2x$ratio,decreasing = T),]
top2x<-rownames(top2x)[1:2]
lab_tmp<-unique(rbind(sig.x[rownames(sig.x)%in%top2x,],sig.y[rownames(sig.y)%in%top2y,]))
if (dim(lab_tmp)[1]>0){
  pick_GID<-data.frame("ID"=unique(str_split_i(rownames(lab_tmp),"___",2)),
                       "Name"=unique(str_split_i(rownames(lab_tmp),"_",2)))
  lab_tmp$IID<-rownames(lab_tmp)
  lab_tmp<-separate(lab_tmp,"IID",into="IID",sep="___",extra="drop")
  lab_tmp<-separate(lab_tmp,"IID",into=c("II","Gname"),sep="_",remove = F)
  text(log2(lab_tmp[,c(2,3)]),labels = lab_tmp$IID,cex=0.25)
}
sig<-merge(sig.x,sig.y,by=1:7,all=T)
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,3],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#ccor list for test
cocor_dat<-list(F=data.frame("F1"=FLEXI_by_GID$K562,"F2"=FLEXI_by_GID$HEK))
#Unfrag scatter
#FELXIs between K562 and HEK
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,3)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,3)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(2,3)]),labels = lab_tmp$Name,cex=0.25)
}
sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%(union(sig.x,sig.y)),]
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
text(10,-6,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(10,-7.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,3],method = "spearman"),digits=2, format="f")
text(10,-7.5,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(10,-10.5,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$K562,"G2"=FLEXI_by_GID$HEK)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#P=0.002, p<0.01

#Unfrag scatter
FLEXI_by_GID<-Cell_counts_down[!(Cell_counts_down[,2]==2^-10 & Cell_counts_down[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,3)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,3)]),col="#0000FF80")
##label FLEXI host gens
if (dim(pick_GID)[1]>0){
  lab_tmp<-merge(FLEXI_by_GID,pick_GID,by=1)
  text(log2(lab_tmp[,c(2,3)]),labels = lab_tmp$Name,cex=0.25)
}

sig<-FLEXI_by_GID[FLEXI_by_GID$ID%in%union(sig.x,sig.y),]
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
text(2,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
#text(2,-8,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
### reproduced FELXIs
cor_p=formatC(cor(sig[,2],sig[,3],method = "pearson"),digits=2, format="f")
#cor_s=formatC(cor(sig[,2],sig[,3],method = "spearman"),digits=2, format="f")
text(2,-8,bquote(atop(italic(r)== .(cor_p)~phantom())),col="red")
#text(2,-10,bquote(atop(italic(r[s])== .(cor_s)~phantom())),col="red")
#cor-test
cocor_datUF<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$K562,"G2"=FLEXI_by_GID$HEK)))
cocor(~F1 + F2 | G1 + G2, cocor_datUF)
#p<0.0001
dev.off()

#remove objects
rm(list=c("temp","cor_p","GID_list","sig.x","sig.y","sig","cor_s","y_max","cocor_dat","Cell_counts_down",
          "Unfrag_total","Cell_counts","FLEXI_by_GID","FLEXI_CPM"))

#Fig3D table, a excel table with command line to generate it and the final form
# Fig3D.xlsx

#Fig4A bar plot and FigS9, only Cellular FLEXIs in HEk,Hela,K562 and UHRR
RBP_4cell<-read.table(gzfile("4cell_plasma_combined_RBP.info.gz"),col.names=c("ID","RBP"))
FLEXI<-dat[rowSums(dat[,73:76])>0,]
RBP_4cell<-RBP_4cell[RBP_4cell$ID%in%FLEXI$ID,]
RBP_4cell<-unique(RBP_4cell)
RBP_fre<-data.frame(table(RBP_4cell$RBP))
RBP<-read.delim("150_RBP_AGO_DICER.info")
all_intron_RBP<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron_RBP<-unique(all_intron_RBP)
all_intron_RBP_short<-all_intron_RBP[all_intron_RBP$Len<=300,c(1,3)]
all_intron_RBP<-all_intron_RBP[all_intron_RBP$Len>300,c(1,3)]
long_fre<-data.frame(table(all_intron_RBP$RBP))
#make Other short introns
all_intron_RBP_short<-all_intron_RBP_short[!all_intron_RBP_short$FLEXI%in%FLEXI$ID,]
OtherShort_fre<-data.frame(table(all_intron_RBP_short$RBP))
RBP_fre<-merge(RBP_fre,OtherShort_fre,by=1,all=T)
RBP_fre<-merge(RBP_fre,long_fre,by=1,all=T)
RBP_fre[is.na(RBP_fre)]<-0
colnames(RBP_fre)<-c("RBP","Cell","Other Short introns","Long introns")
RBP_fre<-merge(RBP_fre,RBP[,c(1,4,5,11)],by=1)
RBP_fre$col<-(RBP_fre$Splicing.regulation+RBP_fre$Spliceosome)/3+RBP_fre$microRNA.processing
RBP_fre<-RBP_fre[,c(1:4,8)]
RBP_fre[RBP_fre$col>1,5]<-4
RBP_fre[RBP_fre$col==1,5]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,5]<-2
RBP_fre[RBP_fre$col==0,5]<-1


B_col=c("black","red","orange","skyblue","gray50")
RBP_fre$label=""
RBP_fre$label[RBP_fre$Cell>=30]<-"*"

### get numbers of each subset
'''
for (i in c(73:76)){
  names<-as.character(RBP_fre$RBP[RBP_fre$Cell>=30])
  tmp<-dat$ID[dat[,i]>0]
  tmp<-RBP_4cell[RBP_4cell$ID%in%tmp,]
  if (i==73){
    tmp1<-data.frame(table(tmp$RBP))
  } else {
    tmp1<-merge(tmp1,data.frame(table(tmp$RBP)),by=1,all=T)
  }
}
tmp1[is.na(tmp1)]<-0
tmp1<-tmp1[tmp1$Var1%in%names,]
range(tmp1[,-1])
'''
### 6-1086

pdf("Figures/Fig4A_1.pdf",height=8,width=12)
tmp<-RBP_fre[RBP_fre$Cell>=30,]
tmp<-tmp[order(tmp$Cell,decreasing = T),]
### the 6 core spliceosomal proteins are color = 5
tmp$col[tmp$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5
mp<-barplot(tmp$Cell,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[tmp$col],ylab="FLEXI RNAs")
text(mp,7000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)
dev.off()

pdf("Figures/FigS9.pdf",height=24,width=12)
par(mfrow=c(3,1))
tmp<-RBP_fre[RBP_fre$Cell>0,]
tmp<-tmp[order(tmp$Cell,decreasing = T),]
### the 6 core spliceosomal proteins are color = 5
tmp$col[tmp$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5
mp<-barplot(tmp$Cell,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[tmp$col],ylab="FLEXI RNAs")
text(mp,10000,labels=tmp$label,cex=0.5)
text(mp,7000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)

#other short intron
tmp<-RBP_fre[RBP_fre$`Other Short introns`>0,]
tmp<-tmp[order(tmp$`Other Short introns`,decreasing = T),]
### the 6 core spliceosomal proteins are color = 5
tmp$col[tmp$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5
mp<-barplot(tmp$`Other Short introns`,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[tmp$col],
            ylab="Other short introns")
text(mp,10000,labels=tmp$label,cex=0.5)
text(mp,7000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)
#long intron
tmp<-RBP_fre[RBP_fre$`Long introns`>0,]
tmp<-tmp[order(tmp$`Long introns`,decreasing = T),]
### the 6 core spliceosomal proteins are color = 5
tmp$col[tmp$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5
mp<-barplot(tmp$`Long introns`,ylim=c(1,100000),log="y",cex.names=0.5,col=B_col[tmp$col],
            ylab="Long Introns")
text(mp,100000,labels=tmp$label,cex=0.5)
text(mp,70000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)
dev.off()
rm(list=c("tmp","long_fre","mp","OtherShort_fre","B_col","all_intron_RBP","all_intron_RBP_short"))
# Fig4A_2 grids
### new 4A_2
RBP53<-read.delim("53_RBP_info_fig4_7_new.txt")
##mRNA level: if a gene contain binding site affect its expression when KD of that RBP?
## 0 non-sig, 2,B:bidirectional (gray), 3 (sig bias towards LFC<0,down),red, 4 (sig bias towards LFC>0,up), skyblue

pdf("Figures/Fig4A_2.pdf",width=12,height=8)
par(mfrow=c(3,1),mar = c(5,2,2,20))
image(1:53,1:16,as.matrix(RBP53[,18:3]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP53)[18:3],tick = FALSE)
image(1:53,1:12,as.matrix(RBP53[,30:19]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:12,labels = colnames(RBP53)[30:19],tick = FALSE)
image(1:53,1:8,as.matrix(RBP53[,38:31]),col=c("white","blue","gray80","red","skyblue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:8,labels = colnames(RBP53)[38:31],tick = FALSE)
axis(1,las=2,at = 1:53,labels = RBP53$RBP.name,tick = FALSE,cex=0.5)
dev.off()

#Fig5
# Cell: 4cell lines FLEXI
# Other short introns: all other non_FLEXI short introns
# Long introns: all long introns
R_sum<-colSums(RBP_fre[,2:4])
RBP_fre<-RBP_fre[,c(2:5,1)]
#exact fisher test pvalue
# LvF long intron vs FLEXI, 3 vs 1
# SvF other short introns vs FLEXI, 2 vs 1
RBP_fre$LvFpvlue<-1
RBP_fre$SvFpvlue<-1
for (i in 1:152){
  RBP_fre[i,6]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(3,1)],R_sum[c(3,1)]),2,2))$p.value
  RBP_fre[i,7]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,1)],R_sum[c(2,1)]),2,2))$p.value
}
RBP_fre[,1:3]<-data.frame(prop.table(as.matrix(RBP_fre[,1:3]),margin = 2)*100)
percent_cutoff<-1
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue","gray75")
RBP_fre$col[RBP_fre$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5
pdf("Figures/Fig5.pdf",height=8,width=8)
par(mfrow=c(2,2))
par(pch=16,pty="s")
plot(RBP_fre[RBP_fre$col==1,c(3,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,1)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,1)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,1)],col="skyblue",cex=1.5)
points(RBP_fre[RBP_fre$col==5,c(3,1)],col="gray75",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,3]>=4)),c(3,1)],
     labels = RBP_fre$RBP[RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,3]>=4)],
     col=col[RBP_fre$col[RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=4| RBP_fre[,3]>=4)]])
#subpanel
plot(RBP_fre[RBP_fre$col==1,c(3,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,1)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,1)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,1)],col="skyblue",cex=1.5)
points(RBP_fre[RBP_fre$col==5,c(3,1)],col="gray75",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,3]>=1)),c(3,1)],
     labels = RBP_fre$RBP[RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,3]>=1)],
     col=col[RBP_fre$col[RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,3]>=1)]])


plot(RBP_fre[RBP_fre$col==1,c(2,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(2,1)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(2,1)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(2,1)],col="skyblue",cex=1.5)
points(RBP_fre[RBP_fre$col==5,c(2,1)],col="gray75",cex=1.5)
abline(0,1,col="red")
  text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,2]>=4)),c(2,1)],
     labels = RBP_fre$RBP[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,2]>=4)],
     col=col[RBP_fre$col[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,2]>=4)]])

plot(RBP_fre[RBP_fre$col==1,c(2,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(2,1)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(2,1)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(2,1)],col="skyblue",cex=1.5)
points(RBP_fre[RBP_fre$col==5,c(2,1)],col="gray75",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,2]>=1)),c(2,1)],
     labels = RBP_fre$RBP[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,2]>=1)],
     col=col[RBP_fre$col[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,2]>=1)]])
dev.off()
rm(col,i,percent_cutoff,pvalue_cutoff,R_sum)

#Fig6
# make new table of FLEXI only in 4 cells
group_list<-list(group1="agotron")
group_list$group2<-c("mitron")
group_list$group3<-c("sno")
group_list$group4<-c("highP")
group_list$group5<-c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2")
group_list$group6<-c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24","FXR2")
group_list$group7<-c("TIA1","TIAL1")
group_list$group8<-c("U2AF1","U2AF2","KHSRP")
group_list$group9<-c("AATF","NOLC1","DKC1","SMNDC1")
group_list$group10<-c("AGO","DICER")
Name_list<-c("Agotron","Mirtron","snoRNA FLEXIs","Conserved FLEXIs",paste0("Cluster ",1:6))
RBP_fre<-data.frame(table(RBP_4cell$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,46)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]
RBP_fre$col[RBP_fre$RBP%in%c("PRPF8","PPIG","BUD13","EFTUD2","AQR","SF3B4")]<-5

percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue","gray75")
pdf("Figures/Fig6.pdf",width=8,height=20,onefile = T)
par(mfrow=c(10,4),pch=16,pty="s",mar=c(2,2,2,2))
for (i in 1:10){
  if (i==3){
    FLEXI_dat<-FLEXI[FLEXI$Has_snoRNA!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FLEXI[FLEXI$Has_snoRNA==".",]
  } else if (i==4){
    FLEXI_dat<-FLEXI[FLEXI$PhastCon30>=0.75,]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FLEXI[FLEXI$PhastCon30<0.75,]
  } else if (i==1){
    FLEXI_dat<-FLEXI[FLEXI$Is_agotron!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FLEXI[FLEXI$Is_agotron==".",]
  } else if (i==2){
    FLEXI_dat<-FLEXI[FLEXI$Is_mirtron!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FLEXI[FLEXI$Is_mirtron==".",]
  } else {
    RBP_list_sig<-group_list[[i]]
    FLEXI_list<-unique(RBP_4cell$ID[RBP_4cell$RBP%in%RBP_list_sig])
    FLEXI_dat<-FLEXI[FLEXI$ID%in%FLEXI_list,]
    temp2<-FLEXI[!FLEXI$ID%in%FLEXI_list,]
  }
  FLEXI_RBP<-data.frame(table(RBP_4cell$RBP[RBP_4cell$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:dim(RBP_plot)[1]){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=percent_cutoff,6]<-"D"
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="bonferroni")
  axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
  plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,
       bty="n",col=col[RBP_plot$col],ylab=NA,xlab=NA)
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff & (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
         col=col[RBP_plot$col[sig_cutoff]],
         labels = RBP_plot$RBP.name[sig_cutoff])
  }
  Len_t<-0
  GC_t<-0
  MFE_t<-0
  #  PhastCon30_t<-0
  rep_times<-1000
  for (inter in 1:rep_times){
    sample_size<-dim(FLEXI_dat)[1]
    test_set<-FLEXI[sample(1:dim(FLEXI)[1],sample_size,replace = F),]
    if (wilcox.test(FLEXI_dat$Len,test_set$Len,exact = F)$p.value>0.05) {Len_t=Len_t+1}
    if (wilcox.test(FLEXI_dat$GC,test_set$GC,exact = F)$p.value>0.05) {GC_t=GC_t+1}
    if (wilcox.test(FLEXI_dat$MFE,test_set$MFE,exact = F)$p.value>0.05) {MFE_t=MFE_t+1}
    #    if (wilcox.test(FLEXI_dat$PhastCon30,test_set$PhastCon30,exact = F)$p.value>0.05) {PhastCon30_t=PhastCon30_t+1}
  }
  len_max=c(0.04,rep(0.02,9))
  GC_max=c(0.08,0.06,rep(0.05,8))
  MFE_max=c(0.06,rep(0.03,9))
  #Length
  plot(density(FLEXI_dat$Len),bty="n",xlim=c(0,350),lwd=1.5,
       ylim=c(0,len_max[i]),main=NA,xlab="Intron length (bp)",col="red")
  lines(density(temp2$Len),xlim=c(0,350),lwd=1.5,col="black")
  legend(120,0.018,lty=c(1,1),lwd=1.5,col=c("red","black"),
         legend = c(paste0(Name_list[i]," (",dim(FLEXI_dat)[1],")"),"Others"),bty="n")
  pvalue<-wilcox.test(FLEXI_dat$Len,temp2$Len,exact = F)$p.value
  if ((Len_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
  #GC
  plot(density(FLEXI_dat$GC),bty="n",xlim=c(0,100),lwd=1.5,
       ylim=c(0,GC_max[i]),main=NA,xlab="GC (%)",col="red")
  lines(density(temp2$GC),xlim=c(0,100),lwd=1.5,col="black")
  pvalue<-wilcox.test(FLEXI_dat$GC,temp2$GC,exact = F)$p.value
  if ((GC_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
  #MEF
  plot(density(FLEXI_dat$MFE),bty="n",xlim=c(-150,0),lwd=1.5,
       ylim=c(0,MFE_max[i]),main=NA,xlab="Minimal free energy (MFE; kcal/mol)",col="red")
  lines(density(temp2$MFE),xlim=c(-150,0),lwd=1.5,col="black")
  pvalue<-wilcox.test(FLEXI_dat$MFE,temp2$MFE,exact = F)$p.value
  if ((MFE_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
}
dev.off()
rm(axis_max,col,FLEXI_list,GC_max,GC_t,i,inter,j,
   len_max,Len_t,MFE_max,MFE_t,Name_list,RBP_fre,
   percent_cutoff,pvalue,R_sum,pvalue_cutoff,RBP_list_sig,
   rep_times,sample_size,sig_cutoff,group_list,
   RBP_plot,FLEXI_dat,test_set,FLEXI_RBP,temp2)

### Fig7A & 8B
## new 7A
RBP53<-read.delim("53_RBP_info_fig4_7_new.txt")
RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 1))
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,-1]
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
temp<-RBP53[match(RBP_clus_image$RBP.name,RBP53$RBP.name),]
RBP_clus_image<-cbind(RBP_clus_image,temp)
icol<-c("white",brewer.pal(12,"Paired"))

pdf("Figures/Fig7A.pdf",width=8,height=8)
par(mfrow=c(4,1),mar = c(5,2,2,20))
image(1:28,1:4,as.matrix(RBP_clus_image[,4:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:4,labels = colnames(RBP_clus_image)[4:1],tick = FALSE,cex=0.5)
image(1:28,1:15,as.matrix(RBP_clus_image[,c(23:13,11:8)]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:15,labels = colnames(RBP_clus_image)[c(23:13,11:8)],tick = FALSE)
image(1:28,1:12,as.matrix(RBP_clus_image[,35:24]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:12,labels = colnames(RBP_clus_image)[35:24],tick = FALSE)
image(1:28,1:8,as.matrix(RBP_clus_image[,43:36]),col=c("white","blue","gray80","red","skyblue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:8,labels = colnames(RBP_clus_image)[43:36],tick = FALSE)
axis(1,las=2,at = 1:28,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)
dev.off()

#sheet 2 Fig8B
RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 2))
colnames(RBP_clus_image)[c(4,9,30,35)]<-c("SUB1","METAP2","SUB1-1","METAP2-1")
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:40]
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
temp<-RBP53[match(RBP_clus_image$RBP.name,RBP53$RBP.name),]
temp[29,]<-temp[3,]
temp[34,]<-temp[8,]
RBP_clus_image<-cbind(RBP_clus_image,temp)
icol<-c("white",brewer.pal(8,"Set1"),brewer.pal(12,"Paired"))
pdf("Figures/Fig8B.pdf",width=8,height=15)
par(mfrow=c(4,1),mar = c(5,2,2,20))
image(1:39,1:24,as.matrix(RBP_clus_image[,24:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:24,labels = colnames(RBP_clus_image)[24:1],tick = FALSE,cex=0.5)
axis(1,las=2,at = 1:39,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)
dev.off()
rm(icol,temp,RBP_clus_image)

#FigS1 barplot, profile of TGIRT-seq
#FigS1
#Barplot of profile in all sample (combiend)
Cell_counts<-read.delim(gzfile("combined_run.counts.gz"))
Cell_counts$K562<-rowSums(Cell_counts[,c(22:29)])
Cell_counts$HEK<-rowSums(Cell_counts[,c(30:37)])
Cell_counts$Hela<-rowSums(Cell_counts[,c(38:47)])
Cell_counts$UHRR<-rowSums(Cell_counts[,c(14:21)])
Cell_counts$MCF7<-rowSums(Cell_counts[,c(6:13)])
Cell_counts$MDA<-rowSums(Cell_counts[,c(4:5)])
Cell_counts<-Cell_counts[,c(1:3,48:53)]
Cell_counts<-Cell_counts[rowSums(Cell_counts[,4:9])>0,]
Cell_counts$Type<-as.factor(Cell_counts$Type)
Cell_counts<-Cell_counts[Cell_counts$Type!="ERCC",]
Cell_counts<-Cell_counts[Cell_counts$Type!="SP",]
levels(Cell_counts$Type)[1]<-levels(Cell_counts$Type)[2]<-levels(Cell_counts$Type)[3]<-levels(Cell_counts$Type)[4]<-levels(Cell_counts$Type)[15]<-"rRNA"
levels(Cell_counts$Type)[14]<-"Pseudogene"
levels(Cell_counts$Type)[10]<-"misc RNA"
levels(Cell_counts$Type)[8]<-"Other lncRNA"
levels(Cell_counts$Type)[21]<-"VT RNA"
levels(Cell_counts$Type)[22]<-"Y RNA"
levels(Cell_counts$Type)[6]<-levels(Cell_counts$Type)[11]<-levels(Cell_counts$Type)[13]<-levels(Cell_counts$Type)[19]<-"Protein coding"
levels(Cell_counts$Type)[13]<-"snoRNA"
levels(Cell_counts$Type)[11]<-"MT tRNA"
Cell_counts$Type<-droplevels(Cell_counts$Type)
Cell_counts$Type2<-Cell_counts$Type
levels(Cell_counts$Type)[2]<-levels(Cell_counts$Type)[3]<-levels(Cell_counts$Type)[8]<-"sncRNA"
levels(Cell_counts$Type)[7]<-levels(Cell_counts$Type)[10]<-levels(Cell_counts$Type)[11]<-"sncRNA"
levels(Cell_counts$Type)[10]<-levels(Cell_counts$Type)[11]<-"sncRNA"
rownames(Cell_counts)<-Cell_counts$ID

agg<-aggregate(.~Type,data=Cell_counts[,c(3:9)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:7)]),2))
agg.prob$Type=agg$Type
agg.prob<-agg.prob[c(4,8,3,5,6,2,9,7,1),]
agg.prob$Type<-as.character(agg.prob$Type)

pdf("Figures/FigS1_A.pdf")
scol <- brewer.pal(9, "Set3")[c(6,2:5,1,7:9)]
mp<-barplot(as.matrix(agg.prob[,1:6]*100),col=scol,axes=F,
            names.arg=rep(NA,6),width=0.05,space=0.3,
            legend.text=agg.prob$Type,
            adj=0.12,args.legend=list(x=1,y=75,bty="n"),xlim=c(0,1.1))
axis(1,labels=NA,at=c(0,0.4),lwd.ticks = 0)
axis(1,labels=NA,at=mp)
text(mp+0.03, par("usr")[3]-3, labels=colnames(agg.prob[,1:6]),
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()

agg<-aggregate(.~Type2,data=Cell_counts[Cell_counts$Type=="sncRNA",c(4:10)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:7)]),2))
agg.prob$Type=agg$Type2
agg.prob$Type<-as.character(agg.prob$Type)

pdf("Figures/FigS1_B.pdf")
scol <- brewer.pal(8, "Set1")[c(6,2:5,1,7:9)]
mp<-barplot(as.matrix(agg.prob[,1:6]*100),col=scol,axes=F,
            names.arg=rep(NA,6),width=0.05,space=0.3,
            legend.text=agg.prob$Type,
            adj=0.12,args.legend=list(x=1,y=75,bty="n"),xlim=c(0,1.1))
axis(1,labels=NA,at=c(0,0.4),lwd.ticks = 0)
axis(1,labels=NA,at=mp)
text(mp+0.03, par("usr")[3]-3, labels=colnames(agg.prob[,1:6]),
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()


dat1<-read.delim("Protein.info.alldataset")
agg<-t(dat1[,2:5])
colnames(agg)<-dat1$Name
agg.prob<-data.frame(prop.table(as.matrix(agg),2))

pdf("Figures/FigS1_C.pdf")
scol <- rev(c("#6699cc","#999999","#ff9966","#ff6633"))
mp<-barplot(as.matrix(agg.prob[,1:6]*100),col=scol,axes=F,
            names.arg=rep(NA,6),width=0.05,space=0.3,
            legend.text=rownames(agg.prob),
            adj=0.12,args.legend=list(x=1,y=75,bty="n"),xlim=c(0,1.1))
axis(1,labels=NA,at=c(0,0.4),lwd.ticks = 0)
axis(1,labels=NA,at=mp)
text(mp+0.03, par("usr")[3]-3, labels=colnames(agg.prob[,1:6]),
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()
rm(list=c("Cell_counts","agg","agg.prob","scol","dat1","mp"))

#FigS2 IGV screen shot in folder 

#FigS3A FLEXI and long FLEXI distribution 
# len density by all FLEXIs all length
longFLEXI<-list(K562=read.table("K562.longFLEXI.len",col.names=c("Len")))
longFLEXI[[2]]<-read.table("HEK.longFLEXI.len",col.names=c("Len"))
longFLEXI[[3]]<-read.table("Hela.longFLEXI.len",col.names=c("Len"))
longFLEXI[[4]]<-read.table("UHRR.longFLEXI.len",col.names=c("Len"))
tmp<-dat[dat$K562>0,]
tmp<-data.frame("Len"=c(rep(tmp$Len,tmp$K562)))
longFLEXI[[1]]<-rbind(longFLEXI[[1]],tmp)
tmp<-dat[dat$HEK>0,]
tmp<-data.frame("Len"=c(rep(tmp$Len,tmp$HEK)))
longFLEXI[[2]]<-rbind(longFLEXI[[2]],tmp)
tmp<-dat[dat$Hela>0,]
tmp<-data.frame("Len"=c(rep(tmp$Len,tmp$Hela)))
longFLEXI[[3]]<-rbind(longFLEXI[[3]],tmp)
tmp<-dat[dat$UHRR>0,]
tmp<-data.frame("Len"=c(rep(tmp$Len,tmp$UHRR)))
longFLEXI[[4]]<-rbind(longFLEXI[[4]],tmp)
names(longFLEXI)<-c("K562","HEK","Hela","UHRR")

pdf("Figures/FigS3A.pdf",width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
for (i in 1:4){
  d<-density(log10(longFLEXI[[i]]$Len))
  div<-max(which(d$x <= log10(300)))
  plot(density(log10(longFLEXI[[i]]$Len)),bty="n",xlab="Intron length (bp)",
       xlim=c(1,6),main=names(longFLEXI)[i],col="red",ylim=c(0,4),axes=F)
  polygon(c(d$x[1],d$x[1:div], d$x[div]), c(0, d$y[1:div],0),col="red")
  polygon(c(d$x[div], d$x[div:length(d$x)], d$x[length(d$x)]), c(0, d$y[div:length(d$x)], 0),col="blue")
  axis(2,labels=seq(0,4,2),las=1,at=seq(0,4,2),las=2)
  axis(1,labels=c(10,parse(text='10^2'),300,bquote(10^4),bquote(10^6)),
       at=c(1,2,log10(300),4,6))
  if (i==1){
    legend(3,3,pch=16,legend = c("FLEXIs (<=300 nt)","Long FLEXIs (>300 nt)"),col=c("red","blue"),bty="n")
  }
}
dev.off()

## plot insert RPM of short/long FLEXIs
longFLEXI<-list(read.table(gzfile("K562.longFLEXI.counts.gz"),col.names=c("ID","Counts")))
longFLEXI[[2]]<-read.table(gzfile("HEK.longFLEXI.counts.gz"),col.names=c("ID","Counts"))
longFLEXI[[3]]<-read.table(gzfile("Hela.longFLEXI.counts.gz"),col.names=c("ID","Counts"))
longFLEXI[[4]]<-read.table(gzfile("UHRR.longFLEXI.counts.gz"),col.names=c("ID","Counts"))
names(longFLEXI)<-c("K562","HEK","Hela","UHRR")
shortFLEXI<-dat[,c(1,74:76,73)]
tmp_total<-mapped_reads[c(4,5,6,3)]

pdf("Figures/FigS3A_insert.pdf",width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
for (i in 1:4){
  tmp<-shortFLEXI[shortFLEXI[,i+1]>0,i+1]
  tmp1<-longFLEXI[[i]]$Counts
  plot(density(log10(tmp/tmp_total[i])),bty="n",xlab="RPM",
       xlim=c(-4,2),main=names(longFLEXI)[i],col="red",ylim=c(0,2),axes=F)
  lines(density(log10(tmp1/tmp_total[i])),col="blue")
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),parse(text='10^-2'),1,bquote(10^2)),
       at=seq(-4,2,2))
  if (i==1){
    legend(-2,2,pch=16,legend = c("FLEXIs (<=300 nt)","Long FLEXIs (>300 nt)"),col=c("red","blue"),bty="n")
  }
}
dev.off()
rm(i,div,tmp,tmp1,tmp_total,longFLEXI,shortFLEXI,d)

#FigS6 lrm
gene_counts<-read.delim("Combined.corrected.sncRNA.counts")
temp<-gene_counts[gene_counts$Type=="snRNA",c(4,9:16)]
gene_for_correlation<-data.frame("U1"=colSums(temp[grep("U1\\b|U1[A-Z]",temp$Name),2:9]))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U2"=colSums(temp[grep("U2\\b|U2[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4"=colSums(temp[grep("U4\\b|U4[B-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U5"=colSums(temp[grep("U5\\b|U5[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6"=colSums(temp[grep("U6\\b|U6[B-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U7"=colSums(temp[grep("U7\\b|U7[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4ATAC"=colSums(temp[grep("U4[ATAC|atac]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6ATAC"=colSums(temp[grep("U6[ATAC|atac]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SK"=colSums(gene_counts[gene_counts$Type=="7SK",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SL"=colSums(gene_counts[gene_counts$Type=="7SL",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RPPH1"=colSums(gene_counts[gene_counts$Name=="RPPH1",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RMRP"=colSums(gene_counts[gene_counts$Name=="RMRP",9:16])))


temp1<-cbind(gene_for_correlation,
             data.frame("U11"=colSums(temp[grep("U11\\b|U11[A-Z]",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("U12"=colSums(temp[grep("U12\\b|U12[A-Z]",temp$Name),2:9])))
temp<-gene_counts[gene_counts$Type=="snoRNA",c(4,9:16)]
temp1<-cbind(temp1,
             data.frame("SNORD3"=colSums(temp[grep("SNORD3\\b|SNORD3[A-Z]|U3",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("SNORD118"=colSums(temp[grep("SNORD118\\b|SNORD118[A-Z]|U8",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("SNORD13"=colSums(temp[grep("SNORD13\\b|SNORD13[A-Z]|U13",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("SNORD14"=colSums(temp[grep("SNORD14\\b|SNORD14[A-Z]",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("SNORD22"=colSums(temp[grep("SNORD22\\b|SNORD22[A-Z]",temp$Name),2:9])))
temp1<-cbind(temp1,
             data.frame("SNORD44"=colSums(temp[grep("SNORD44\\b|SNORD44[A-Z]",temp$Name),2:9])))

#mappedreads and corrected mapped reads
correct_mapped_reads<-c(715.24152,277.8381,768.43375,368.4548,666.34185,362.6804,713.77529,230.0604)
names(correct_mapped_reads)<-rownames(gene_for_correlation)
gene_for_correlation<-gene_for_correlation/correct_mapped_reads
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",2),"7SK","7SL","RNase P","MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)

temp1<-temp1/correct_mapped_reads
temp1<-data.frame(t(temp1))
temp1<-round(temp1,0)

pdf("Figures/FigS6A_1.pdf",height=6,width=8)
par(bty="n",mfrow=c(2,4),mar=c(2,3,2,1))
scol <- brewer.pal(8, "Set1")
y<-log10(gene_for_correlation[,9])
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i-1])
  lm.out <- lm(y ~ x)
  newx = seq(min(x),max(x),by = 0.05)
  conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95)
  plot(x,y,pch=16,xlim=c(0,5),ylim=c(0,8),
       col=scol[gene_for_correlation$Type],ylab=bquote(log[10]~(molecule/cell)),
       xlab=bquote(log[10] (RPM)),
       main=colnames(gene_for_correlation)[2*i-1])
  abline(lm.out, col="lightblue")
  matlines(newx, conf_interval[,2:3], col = "blue", lty=2)
  if (i==1){
    legend(3,4,legend = levels(gene_for_correlation$Type),col=scol,bty="n",pch=16,cex=0.7)
  }
  cor_p=formatC(cor(x,y,method = "pearson"),digits=2, format="f")
  text(1,7,bquote(atop(italic(r)== .(cor_p)~phantom())))
}
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i])
  lm.out <- lm(y ~ x)
  newx = seq(min(x),max(x),by = 0.05)
  conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95)
  plot(x,y,pch=16,xlim=c(0,5),ylim=c(0,8),
       col=scol[gene_for_correlation$Type],ylab=bquote(log[10]~(molecule/cell)),
       xlab=bquote(log[10] (RPM)),
       main=colnames(gene_for_correlation)[2*i])
  abline(lm.out, col="lightblue")
  matlines(newx, conf_interval[,2:3], col = "blue", lty=2)
  if (i==1){
    legend(3,4,legend = levels(gene_for_correlation$Type),col=scol,bty="n",pch=16,cex=0.7)
  }
  cor_p=formatC(cor(x,y,method = "pearson"),digits=2, format="f")
  text(1,7,bquote(atop(italic(r)== .(cor_p)~phantom())))
}
dev.off()
'''
#dPCR 3 standard, FigS6A part2 table
temp<-gene_counts[gene_counts$Type=="snRNA",c(4,9:16)]
dPCR_standards<-data.frame("U7"=colSums(temp[grep("U7\\b|U7[A-Z]",temp$Name),2:9]))
temp<-gene_counts[gene_counts$Type=="snoRNA",c(4,9:16)]
dPCR_standards<-cbind(dPCR_standards,
                      data.frame("SNORD14B"=colSums(temp[temp$Name=="SNORD14B",2:9])))
dPCR_standards<-cbind(dPCR_standards,
                      data.frame("SNORD44"=colSums(temp[temp$Name=="SNORD44",2:9])))
dPCR_standards<-dPCR_standards/correct_mapped_reads
y<-log10(gene_for_correlation[,9])
for (i in 1:6){
  x<-log10(gene_for_correlation[,i])
  lm.out <- lm(y ~ x)
  newx = log10(unlist(dPCR_standards[i,]))
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx)))
  if (i==1){
    dPCR<-data.frame(round(10^conf_interval,0))
  } else {
    dPCR<-cbind(dPCR,round(10^conf_interval,0))
  }
}
colnames(dPCR)<-colnames(gene_for_correlation)[1:6]
dPCR$Literature_value<-c(4000,10000,10000)
dPCR$Name<-rownames(dPCR)
dPCR<-dPCR[,c(8,7,1,3,5,2,4,6)]
write.table(dPCR,"Figures/FigS6A_table.tsv",quote=F,sep="\t",row.names=F)
'''
dPCR<-read.delim("FigS6A_2_table.tsv")
pdf("Figures/FigS6A_2.pdf")
tt <- ttheme_default(base_size=8,colhead=list(fg_params = list(parse=TRUE,fontface="plain")))
colnames(dPCR)<-c("Copies/cell","Literature\nvalue","HEK-293T\nbefore corr","HeLa S3\nbefore corr","UHRR\nbefore corr",
                  "HEK-293T\nafter corr","HeLa S3\nafter corr","UHRR\nafter corr")
grid.table(dPCR,theme=tt,rows=NULL)
dev.off()

'''
FLEXI<-dat[rowSums(dat[,73:76])>0,]
check_list<-c("1I_H3F3B___ENSG00000132475___ENST00000587560___protein_coding___protein_coding",
              "3I_RAN___ENSG00000132341___ENST00000543796___protein_coding___protein_coding",
              "7I_EIF4A1___ENSG00000161960___ENST00000578495___protein_coding___protein_coding",
              "4I_SEMA4C___ENSG00000168758___ENST00000305476___protein_coding___protein_coding",
              "4I_JUP___ENSG00000173801___ENST00000393930___protein_coding___protein_coding",
              "1I_FTH1___ENSG00000167996___ENST00000532601___protein_coding___protein_coding",
             "11I_THBS1___ENSG00000137801___ENST00000260356___protein_coding___protein_coding",
              "16I_POLG___ENSG00000140521___ENST00000637264___protein_coding___protein_coding",
              "2I_ACTB___ENSG00000075624___ENST00000464611___protein_coding___protein_coding",
              "1I_RPS2___ENSG00000140988___ENST00000526586___protein_coding___protein_coding")
RPM_eve<-FLEXI[match(check_list,FLEXI$ID),c(1,75:76,73)]
RPM_eve$HEK2<-c(10.3392,33.7969,378.067,63.2418,7.00995,0,4.10477,71.1403,557.579,150.335859)
RPM_eve$Hela2<-c(0,40.5708,150.994,35.9247,6458.47,3791.92,0,139.617,1366.71,113.1693)
RPM_eve$UHRR2<-c(17.1409,10.0909,147.737,382.365,20.4748,0.81306,215.188,85.8948,3582.42,59.849117)
RPM_eve<-RPM_eve[,c(1,2,5,3,6,4,7)]
rownames(RPM_eve)<-RPM_eve$ID

RPM_eve<-t(RPM_eve[,2:7])/correct_mapped_reads[1:6]
y<-log10(gene_for_correlation[,9])
for (i in 1:6){
  x<-log10(gene_for_correlation[,i])
  lm.out <- lm(y ~ x)
  newx = log10(unlist(RPM_eve[i,]))
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx)))
  if (i==1){
    dPCR<-data.frame(round(10^conf_interval,0))
  } else {
    dPCR<-cbind(dPCR,round(10^conf_interval,0))
  }
}
colnames(dPCR)<-colnames(gene_for_correlation)[1:6]
write.table(dPCR,"Figures/FigS6B_table_byRPM.tsv",quote=F,sep="\t",row.names=T)
'''
# FigS6B table is in Digital PCR summary 081621.xlsx
rm(dPCR,dPCR_standards,lm.out,RPM_eve,temp,tt,temp1,conf_interval,gene_for_correlation,check_list,cor_p,i,
   newx,PRM_standard,scol,x,y,correct_mapped_reads)


#FigS7A-B cluster PCA, tSNE, ZINB-WAVE
#renumbering based on cutoff, cutoff is based on combined date (rowSums), 
# but applied to individual replicates.
# MDA 2:3
# MCF 4:11
# UHRR 12:19
# Hela 20:29
# K562 batch 2 (newer) 30:37;  batch 1 46:52
# HEK batch 2 (newer) 38:45; batch 1 53:54
#cut off by RPM, 0.01,and 0

# zinb_wave objects are made using make_zinb_object.R

cutoff=c(0,0.01)
P_pch=c(16,17)
dat<-read.delim("all.FLEXI")
dat1<-read.delim("Bio2.FLEXI.counts")
dat<-dat[rowSums(dat[,27:70])>0,c(1,27:70)]
colnames(dat)<-c("ID",paste0("MDA_1_",1:2),paste0("MCF_",1:8),
                 paste0("UHRR_1_",1:8),paste0("K561_1_",1:8),
                 paste0("HEK_1_",1:8),paste0("Hela_",1:10))
dat1<-dat1[,c(1:15,19:28)]
colnames(dat1)<-c("ID",paste0("HEK_2_",1:2),paste0("K562_2_",1:7),
                  paste0("MDA_2_",1:5),paste0("UHRR_2_",1:10))
dat<-merge(dat,dat1,by="ID",all=T)
dat[is.na(dat)]<-0
dat<-dat[rowSums(dat[,2:69])>0,]
dat<-dat[,c(1,4:11,36:45,2:3,12:35,55:69,48:54,46:47)]
sub_mapped_reads<-c(70563696,83971287,93530133,96188889,84195614,95958499,90154958,77528755,
                    91268315,83392799,82750963,62981188,60057226,67704164,86607803,72132500,96955981,64582809,
                    109704771,97786253,
                    87438356,80725874,84193148,84695000,79516751,94126936,80171312,75474477,
                    87388516,82600414,84957488,88141802,95927677,113807174,87425522,73526698,
                    96566478,71023706,90569418,84985720,107210447,92506910,88853873,83524969,
                    37067082,55224062,52426381,46851911,67555020,
                    21744662,38212559,35376687,35131499,40724358,38938071,37683086,38933929,38267462,34349980,
                    5700094,7258578,6925758,5408190,6342889,6941957,9282711,
                    99564097,81230668)
sub_mapped_reads<-sub_mapped_reads/1e6
dat_clus<-dat
dat<-dat_clus
dat1<-data.frame(t(dat[,c(20:69)]))
dat1$label<-c(rep("MDA",2),rep("UHRR",8),rep("K562",8),rep("HEK",8),
              rep("MDA",5),rep("UHRR",10),rep("K562",7),rep("HEK",2))
dat1$batch<-c(rep("1",26),rep("2",24))
dat1$label<-factor(dat1$label)
dat1$batch<-factor(dat1$batch)

col=c("black","royalblue1","goldenrod","orchid","greenyellow","tomato")
col<-col2hex(col)
col<-paste0(col,"A0")
name<-c("MCF-7","HeLa S3","MDA-MB-231","UHRR","K-562","HEK 293T")
file_name<-c("A","B")
for (i in 1:2){
  pdf_name<-paste0("Figures/FigS7",file_name[i],".pdf")
  pdf(pdf_name,height=9,width=9)
  par(pty="s",mfrow=c(3,3))
  #44 dataset
  obj_name<-paste0("44_zinbwave",cutoff[i],"RPM.pbj")
  dat<-dat_clus
  dat<-dat[rowMaxs(t(t(as.matrix(dat[,2:45]))/sub_mapped_reads[1:44]))>=cutoff[i],2:45]
  dat<-dat[rowSums(dat)>0,]
  #PCA
  pca_dat<-data.frame(t(dat/sub_mapped_reads[1:44]))
  pca_dat$label<-c(rep(1,8),rep(2,10),rep(3,2),rep(4:6,each=8))
  pca_dat$label<-as.factor(pca_dat$label)
  pca<-prcomp(pca_dat[,1:(dim(pca_dat)[2]-1)])
  plot(pca$x[,1:2],col=col[pca_dat$label],main="PCA",
       xlab="PC1",ylab="PC2",pch=16,cex=1.5)
  legend(0,4,legend = name,col=col,pch=16,bty="n")
  #t-SNE
  tsne <- Rtsne(pca_dat[,1:dim(pca_dat)[2]-1], dims = 2, perplexity=4, theta = 0.5,normalize=F,
                num_threads=4,verbose=TRUE, max_iter = 5000)
  plot(tsne$Y,main="t-SNE ",col=col[pca_dat$label],
       xlab="t-SNE1",ylab="t-SNE2",pch=16,cex=1.5)
  #ZINB-WaVE
  zinb<-readRDS(obj_name)
  W <- data.frame(reducedDim(zinb))
  W$label<-pca_dat$label
  plot(W$W1,W$W2,main="ZINB-WaVE",col=col[W$label],
       xlab="W1",ylab="W2",pch=16,cex=1.5)
  #biological replicate
  dat<-dat_clus
  dat<-dat[rowMaxs(t(t(as.matrix(dat[,20:69]))/sub_mapped_reads[19:68]))>=cutoff,20:69]
  dat<-as.matrix(dat[rowSums(dat)>0,])
  
  #without batch correction
  BR_col<-col[c(6,5,3,4)]
  obj_name<-paste0("Bio_rep_zinbwave",cutoff[i],"RPM.pbj")
  #PCA
  pca_dat<-data.frame(t(dat/sub_mapped_reads[19:68]))
  pca_dat$label<-dat1$label
  pca_dat$batch<-dat1$batch
  pca_dat$label<-factor(pca_dat$label)
  pca_dat$batch<-factor(pca_dat$batch)
  pca<-prcomp(pca_dat[,1:(dim(pca_dat)[2]-2)])
  plot(pca$x[,1:2],col=BR_col[pca_dat$label],main="PCA",
       xlab="PC1",ylab="PC2",pch=P_pch[pca_dat$batch],cex=1.5)
  
  #t-SNE
  tsne <- Rtsne(dat1[,1:(dim(dat1)[2]-2)], dims = 2,
                perplexity=8, theta = 0.5,normalize=F,
                num_threads=4,verbose=TRUE, max_iter = 5000)
  plot(tsne$Y,  main="t-SNE",col=BR_col[pca_dat$label],pch=P_pch[pca_dat$batch],
       xlab="t-SNE1",ylab="t-SNE2",cex=1.5)
  #ZINB-Wave w/o batch correction
  zinb<-readRDS(obj_name)
  W <- data.frame(reducedDim(zinb))
  W$label<-pca_dat$label
  W$batch<-pca_dat$batch
  plot(W$W1,W$W2,main="ZINB-WaVE (batch effect corrected)",col=BR_col[W$label],pch=P_pch[W$batch],
       xlab="W1",ylab="W2",cex=1.5)
  #ZINB-WaVE-batch corrected
  obj_name<-paste0("Bio_rep_zinbwave_Batch_cor_",cutoff[i],"RPM.pbj")
  zinb<-readRDS(obj_name)
  W <- data.frame(reducedDim(zinb))
  W$label<-pca_dat$label
  W$batch<-pca_dat$batch
  #PCA using ZINB-wave batch corrected normalized counts
  pca_dat<-data.frame(t(assays(zinb)$normalizedValues))
  pca_dat$label<-dat1$label
  pca_dat$batch<-dat1$batch
  pca_dat$label<-factor(pca_dat$label)
  pca_dat$batch<-factor(pca_dat$batch)
  pca<-prcomp(pca_dat[,1:(dim(pca_dat)[2]-2)])
  plot(pca$x[,1:2],col=BR_col[pca_dat$label],main="PCA",
       xlab="PC1",ylab="PC2",pch=P_pch[pca_dat$batch],cex=1.5)
  ##t-SNE after batch correction
  tsne <- Rtsne(W[,1:2], dims = 2,pca=F,
                perplexity=8, theta = 0.5,
                num_threads=4,verbose=TRUE, max_iter = 5000)
  plot(tsne$Y,  main="t-SNE (batch effect corrected)",col=BR_col[W$label],pch=P_pch[W$batch],
       xlab="t-SNE1",ylab="t-SNE2",cex=1.5)
  #zinb-wave batch corrrection
  plot(W$W1,W$W2,main="ZINB-WaVE (batch effect corrected)",col=BR_col[W$label],pch=P_pch[W$batch],
       xlab="W1",ylab="W2",cex=1.5)
  legend(-2,1,legend = levels(W$label),col=BR_col,pch=16,bty="n")
  dev.off()
}
rm(tsne,W,zinb,col,BR_col,i,file_name,name,obj_name,cutoff,P_pch,pdf_name,pca,pca_dat,dat1)

#FigS9A
### FLEXI biomarker
dat<-read.delim("all.FLEXI")
dat1<-read.delim("Bio2.FLEXI.counts")
dat<-dat[rowSums(dat[,27:70])>0,c(1,27:70)]
colnames(dat)<-c("ID",paste0("MDA_1_",1:2),paste0("MCF_",1:8),
                 paste0("UHRR_1_",1:8),paste0("K561_1_",1:8),
                 paste0("HEK_1_",1:8),paste0("Hela_",1:10))
dat1<-dat1[,c(1:15,19:28)]
colnames(dat1)<-c("ID",paste0("HEK_2_",1:2),paste0("K562_2_",1:7),
                  paste0("MDA_2_",1:5),paste0("UHRR_2_",1:10))
dat<-merge(dat,dat1,by="ID",all=T)
dat[is.na(dat)]<-0
dat<-dat[rowSums(dat[,2:69])>0,]
dat<-dat[,c(1,4:11,36:45,2:3,12:35,55:69,48:54,46:47)]
sub_mapped_reads<-c(70563696,83971287,93530133,96188889,84195614,95958499,90154958,77528755,
                    91268315,83392799,82750963,62981188,60057226,67704164,86607803,72132500,96955981,64582809,
                    109704771,97786253,
                    87438356,80725874,84193148,84695000,79516751,94126936,80171312,75474477,
                    87388516,82600414,84957488,88141802,95927677,113807174,87425522,73526698,
                    96566478,71023706,90569418,84985720,107210447,92506910,88853873,83524969,
                    37067082,55224062,52426381,46851911,67555020,
                    21744662,38212559,35376687,35131499,40724358,38938071,37683086,38933929,38267462,34349980,
                    5700094,7258578,6925758,5408190,6342889,6941957,9282711,
                    99564097,81230668)
sub_mapped_reads<-sub_mapped_reads/1e6
dat[,-1]<-t(t(dat[,-1])/sub_mapped_reads)
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)

dat1<-dat
rownames(dat1)<-dat1$ID
dat1<-dat1[,-1]
dat1<-dat1[,order(colnames(dat1))]
dat1<-dat1[rowSums(dat1)>0,]
dat1[dat1==0]<-2^-7
dat1<-log2(dat1)
#ht<-heatmap.2(as.matrix(dat1))
#saveRDS(ht,"FLEXI_cell_all_rep_heatmap.obj")
ht<-readRDS("FLEXI_cell_all_rep_heatmap.obj")
pick<-rev(ht$rowInd)[1:100]
tmp1<-dat1[pick,]
rowLab<-str_split_i(rownames(tmp1),pattern = "___",i=1)
pdf("Figures/FigS9A.pdf",height=11,width=8.5)
heatmap.2(as.matrix(tmp1),labRow = rowLab,
          dendrogram = "none",labCol=colnames(tmp1),
          scale="none",breaks=seq(-7,1,0.08),cexRow = 0.1,
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,Colv=F,Rowv=T,cexCol=0.25,
          col =heatPalette)
dev.off()

### FigS9B
### examples of cell-specific FLEXIs, also marked in Fig3B
pick<-c("2I_ACTG1___ENSG00000184009___ENST00000575842___protein_coding___protein_coding",
        "1I_RPS2___ENSG00000140988___ENST00000526586___protein_coding___protein_coding",
        "1I_FTH1___ENSG00000167996___ENST00000532601___protein_coding___protein_coding",
        "11I_SGK1___ENSG00000118515___ENST00000237305___protein_coding___protein_coding",
        "4I_SESN2___ENSG00000130766___ENST00000253063___protein_coding___protein_coding",
        "3I_TAOK2___ENSG00000149930___ENST00000308893___protein_coding___protein_coding",
        "23I_DOCK6___ENSG00000130158___ENST00000587656___protein_coding___protein_coding",
        "2I_H3F3B___ENSG00000132475___ENST00000587560___protein_coding___protein_coding",
        "1I_ERRFI1___ENSG00000116285___ENST00000467067___protein_coding___protein_coding",
        "11I_THBS1___ENSG00000137801___ENST00000260356___protein_coding___protein_coding")
dat1<-data.frame(label=c(rep("MCF",8),rep("Hela",10),
                         rep("MDA_Bio1",2),rep("UHRR_Bio1",8),rep("K562_Bio1",8),rep("HEK_Bio1",8),
                         rep("MDA_Bio2",5),rep("UHRR_Bio2",10),rep("K562_Bio2",7),rep("HEK_Bio2",2)))
dat2<-dat[dat$ID%in%pick,]
dat2<-dat2[match(pick,dat2$ID),]
pdf("Figures/FigS9B.pdf",height=12,width=8)
par(mfrow=c(5,2),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=2,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

#FigS11 created by Shelby, RBP binding site location
#FigS12 made by Shelby, FLEXIreads in eCLIP dataset
# raw data in S5 table
#FigS13 plotted by ENCODE_DE_ver2.R
#FigS14 plotted by ENCODE_rMAT.R
#FigS15 created by RBP_clus_byCellType.R


#FigS16 individual RBP scatter and density plots
#individual 4 panel for each 47 rbPs
FLEXI<-dat$ID[rowSums(dat[,73:76])>0]
FLEXI<-dat[dat$ID%in%FLEXI,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXIFLEXI,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,46)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]

percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")
RBP_list<-sort(RBP53$RBP.name)
RBP_list<-RBP_list[c(7,8,19,22,23,25,37,43,50,
                     4,6,12,16,32,46,49,52,53,
                     44,45,
                     18,47,48,
                     1,10,24,40,
                     2,9,
                     3,5,11,29,31,39,
                     13,14,15,17,20,21,26,27,28,30,33,34,35,36,38,41,42,51)]

pdf("Figures/FigS16.pdf",width=8,height=20,onefile = T)
par(mfrow=c(10,4),pch=16,pty="s",mar=c(2,2,2,2))
for (i in 1:53) {
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list[i]])
  FLEXI_dat<-FLEXI[FLEXI$ID%in%FLEXI_list,]
  temp2<-FLEXI[!FLEXI$ID%in%FLEXI_list,]
  FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:dim(RBP_plot)[1]){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=percent_cutoff,6]<-"D"
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="bonferroni")
  axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
  plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=paste0(RBP_list[i]," FLEXIs"),
       bty="n",col=col[RBP_plot$col],ylab=NA,xlab=NA)
  axis(1,at=seq(0,axis_max,5),label=F)
  axis(2,at=seq(0,axis_max,5),label=F)
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff & (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
         col=col[RBP_plot$col[sig_cutoff]],
         labels = RBP_plot$RBP.name[sig_cutoff])
  }
  Len_t<-0
  GC_t<-0
  MFE_t<-0
  #  PhastCon30_t<-0
  rep_times<-1000
  for (inter in 1:rep_times){
    sample_size<-dim(FLEXI_dat)[1]
    test_set<-FLEXI[sample(1:dim(FLEXI)[1],sample_size,replace = F),]
    if (wilcox.test(FLEXI_dat$Len,test_set$Len,exact = F)$p.value>0.05) {Len_t=Len_t+1}
    if (wilcox.test(FLEXI_dat$GC,test_set$GC,exact = F)$p.value>0.05) {GC_t=GC_t+1}
    if (wilcox.test(FLEXI_dat$MFE,test_set$MFE,exact = F)$p.value>0.05) {MFE_t=MFE_t+1}
    #    if (wilcox.test(FLEXI_dat$PhastCon30,test_set$PhastCon30,exact = F)$p.value>0.05) {PhastCon30_t=PhastCon30_t+1}
  }
  len_max=0.02
  GC_max=0.08
  MFE_max=0.03
  
  #Length
  plot(density(FLEXI_dat$Len),bty="n",xlim=c(0,350),lwd=1.5,xlab=NA,ylab=NA,
       ylim=c(0,len_max),main=NA,axes=F,col="red")
  lines(density(temp2$Len),xlim=c(0,350),lwd=1.5,col="black")
  legend(120,0.018,lty=c(1,1),lwd=1.5,col=c("red","black"),
         legend = c(paste0(RBP_list[i]," (",dim(FLEXI_dat)[1],")"),"Others"),bty="n")
  axis(1,at=seq(0,350,50),labels = F)
  axis(2,at=seq(0,len_max,0.01),labels = F)
  pvalue<-wilcox.test(FLEXI_dat$Len,temp2$Len,exact = F)$p.value
  if ((Len_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
  #GC
  plot(density(FLEXI_dat$GC),bty="n",xlim=c(0,100),lwd=1.5,
       ylim=c(0,GC_max),main=NA,xlab=NA,ylab=NA,axes=F,col="red")
  lines(density(temp2$GC),xlim=c(0,100),lwd=1.5,col="black")
  axis(1,at=seq(0,100,25),labels = F)
  axis(2,at=seq(0,GC_max,0.02),labels = F)
  pvalue<-wilcox.test(FLEXI_dat$GC,temp2$GC,exact = F)$p.value
  if ((GC_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
  #MEF
  plot(density(FLEXI_dat$MFE),bty="n",xlim=c(-150,0),lwd=1.5,
       ylim=c(0,MFE_max),main=NA,xlab=NA,ylab=NA,axes=F,col="red")
  lines(density(temp2$MFE),xlim=c(-150,0),lwd=1.5,col="black")
  axis(1,at=seq(-150,0,50),labels = F)
  axis(2,at=seq(0,MFE_max,0.01),labels = F)
  pvalue<-wilcox.test(FLEXI_dat$MFE,temp2$MFE,exact = F)$p.value
  if ((MFE_t/rep_times<=0.05) & (pvalue<=0.05)) {
    if (pvalue<0.01) {
      legend("topleft",legend=paste0("p<0.01"),bty="n")
    } else {
      legend("topleft",legend=paste0("p = ",format(pvalue,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  } 
  
}
dev.off()

## FigS5A from Supplementary file
## FigS5B raw image as FigS5B.tiff
## sanger sequencing results file (ab1) in FigS5_sanger_sequencing folder

## Gower distance object and heatmap in make_RBP_clusTable.R

### figure S17
### RBP sites overlap
RBP_over<-read.table("FLEXI_RBP.overlapMax1Overlap.info",
                     col.names=c("ID","Overlap"))
RBP_name<-sort(unique(RBP_over$ID))
### make tables
for (i in 1:146){
  c_name<-RBP_name[i]
  if (i==1){
    RBP_matrix<-data.frame(table(RBP_over$Overlap[RBP_over$ID==c_name]))
    colnames(RBP_matrix)<-c("ID",c_name)
  } else {
    tmp<-data.frame(table(RBP_over$Overlap[RBP_over$ID==c_name]))
    colnames(tmp)<-c("ID",c_name)
    RBP_matrix<-merge(RBP_matrix,tmp,by=1,all=T)
  }
}
rownames(RBP_matrix)<-RBP_matrix$ID
RBP_matrix<-RBP_matrix[,-1]
RBP_matrix[is.na(RBP_matrix)]<-0
RBP_matrix<-log10(RBP_matrix+1)
### only 47 RBPs
RBP53<-read.delim("53_RBP_info_fig4_7.txt")
Core_splicesome<-c("PRPF8","SF3B4","AQR","EFTUD2","BUD13","PPIG")
RBP47<-RBP53$RBP.name[!RBP53$RBP.name%in%Core_splicesome]
RBP_matrix<-RBP_matrix[RBP47,RBP47]
#RBP_matrix<-apply(RBP_matrix,2,function(x){x/max(x)})
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)
L_col<-c("black","red","orange","skyblue","orchid","blue1","goldenrod4")
Clusters<-list(Cluster1=c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"),
               Cluster2=c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24","FXR2"),
               Cluster3=c("TIA1","TIAL1"),
               Cluster4=c("U2AF1","U2AF2","KHSRP"),
               Cluster5=c("AATF","DKC1","NOLC1","SMNDC1"),
               Cluster6=c("AGO","DICER"))

pdf("Figures/FigS17.pdf",height=11,width=8)
for (i in 1:6){
  tmp<-RBP_matrix[,Clusters[[i]]]
  y_max<-2
  tmp$Col<-1
  tmp$Col[rownames(tmp)%in%Clusters[[1]]]<-2
  tmp$Col[rownames(tmp)%in%Clusters[[2]]]<-3
  tmp$Col[rownames(tmp)%in%Clusters[[3]]]<-4
  tmp$Col[rownames(tmp)%in%Clusters[[4]]]<-5
  tmp$Col[rownames(tmp)%in%Clusters[[5]]]<-6
  tmp$Col[rownames(tmp)%in%Clusters[[6]]]<-7
  heatmap.2(as.matrix(tmp[,1:(dim(tmp)[2]-1)]),labRow =rownames(tmp),
            dendrogram = "none",
            scale="none",breaks=seq(0,y_max,y_max/100),margins = c(10, 10),
            main=paste0("Cluster ",i),colRow=L_col[tmp[,dim(tmp)[2]]],
            density.info = "none",trace="none",symm=F,key=T,
            key.title = "log10(binding sites)",key.xlab = NA,key.ylab=NA,
            keysize = 1,colCol = rep(L_col[i+1],length(Clusters[[i]])),
            symbreaks=F,labCol =colnames(tmp),cexRow=0.5,
            col =heatPalette)
}
dev.off()

## Fig10,11,12A plotted by Cellular_cluster.r
## Fig12B from String database results
## FigS18, S19 in Cellular_cluster.r
## FigS20 plotted by Peter. in Cell_frac_peter.R

plot_upset <- function(m, i) {
  cs = comb_size(m)  
  ss<-set_size(m)
  ht <- UpSet(m, set_order=c("FLEXI","Other short introns","Long introns"),
              comb_col = c("black","skyblue","orchid")[comb_degree(m)],
              #comb_order=od,
              column_title=paste0("Cluster ",i," RBP"),
              top_annotation = HeatmapAnnotation(
                "Host genes" = anno_barplot(cs, 
                                            ylim = c(0, max(cs)*1.1),
                                            border = F,
                                            gp = gpar(border =NA,lty=0,
                                                      fill =c("black","skyblue","orchid")[comb_degree(m)]), 
                                            height = unit(6, "cm")), 
                annotation_name_side = "left", 
                annotation_name_rot = 90),
              right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
  ht = draw(ht)
  decorate_annotation("Host genes", {grid.text(formatC(cs,big.mark = ","), x = 1:length(cs), 
                                               y = unit(cs, "native") + unit(2, "pt"), 
                                               default.units = "native", just = c("left", "bottom"), 
                                               gp = gpar(fontsize = 6), rot = 45)})
}
### FLEXI: combined 4 cells (HEK, K-562, Hela nad UHRR)
dat<-read.delim("all.FLEXI")
FLEXI<-dat[rowSums(dat[,73:76])>0,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
FLEXI<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
FLEXI<-unique(FLEXI)
FLEXI<-separate(FLEXI,ID,into=c("IID","GID","TID","GType","TType"),sep="___",remove = F)
### long|short
all_intron_RBP<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron_RBP<-separate(all_intron_RBP,col="FLEXI",into=c("IID","GID","TID","GType","TType"),sep="___",remove = F)
all_intron_RBP_short<-all_intron_RBP[all_intron_RBP$Len<=300,]
all_intron_RBP<-all_intron_RBP[all_intron_RBP$Len>300,]
#make Other short introns
all_intron_RBP_short<-all_intron_RBP_short[!all_intron_RBP_short$FLEXI%in%FLEXI$ID,]

pdf("Figures/Fig8A.pdf",height=11,width=8)
ht_list =vector(6, mode="list")
for (i in 1: length(Clusters)){
  tmp<-list(unique(FLEXI$GID[FLEXI$RBP%in%Clusters[[i]]]),
            unique(all_intron_RBP_short$GID[all_intron_RBP_short$RBP%in%Clusters[[i]]]),
            unique(all_intron_RBP$GID[all_intron_RBP$RBP%in%Clusters[[i]]]))
  names(tmp)<-c("FLEXI","Other short introns","Long introns")
  m = make_comb_mat(tmp)
  dht <- as.ggplot( ~ plot_upset(m, i))
  ht_list[[i]] <- dht
}
patchwork::wrap_plots(ht_list, ncol=2)
dev.off()

## FigS21
## A/B numbers
mitocarat<-read.table("mitocarta.genes")
NR<-read.table("RIbosome.gene")
MR<-read.table("mito.ribosome.gene")
###
##
for (i in 1:length(Clusters[[1]])){
  RBP_name<-Clusters[[1]][i]
  tmp2<-unique(tmp1$GID[grep(RBP_name,tmp1$RBP)])
  print(c(RBP_name,length(tmp2),sum(tmp2%in%mitocarat$V1),
          sum(tmp2%in%NR$V1),sum(tmp2%in%MR$V1)))
  ### exon site host gene %
  #tmp<-read.table(paste0("/stor/work/Lambowitz/ref/TGSEQ/RBP/RBP_cov/",RBP_name,"_exon_UTR.hotsgene"))
  #tmp<-unique(tmp$V1)
  #print(c(RBP_name,length(tmp),sum(tmp%in%mitocarat$V1)))
}
### FLEXIs
#[1] "LARP4" "45"    "10"    "15"    "1"    
#[1] "PABPC4" "54"     "8"      "23"     "1"     
#[1] "SUB1" "26"   "3"    "8"    "0"   
#[1] "DDX3X" "137"   "35"    "14"    "9"    
#[1] "RPS3" "99"   "14"   "10"   "2"   
#[1] "NCBP2" "36"    "3"     "3"     "0"    
#[1] "DDX55" "28"    "4"     "4"     "0"    
#[1] "METAP2" "48"     "4"      "7"      "0" 
tmp1<-unique(tmp$GID)
print(c("FLEXI host genes",length(tmp1),sum(tmp1%in%mitocarat$V1),
        sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
#[1] "FLEXI host genes" "3729"             "280"              "60"               "23"       
### other short intron/long intron
all_intron<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron<-unique(all_intron)
### remove all FLEXI introns
all_intron<-all_intron[!all_intron$FLEXI%in%unique(tmp$ID),]
### seperate GID
all_intron<-separate(all_intron,"FLEXI",into=c("IID","GID","TID","GType","TType"),sep="___",remove = F)
### make long introns
long_intron<-all_intron[all_intron$Len>300,c(1,3,8)]
#make Other short introns
other_short_intron<-all_intron[all_intron$Len<=300,c(1,3,8)]

for (i in 1:length(Clusters[[1]])){
  RBP_name<-Clusters[[1]][i]
  tmp1<-unique(other_short_intron$GID[other_short_intron$RBP==RBP_name])
  print(c(RBP_name,length(tmp1),sum(tmp1%in%mitocarat$V1),
          sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
}
## OSI
#[1] "LARP4" "132"   "25"    "18"    "4"    
#[1] "PABPC4" "112"    "15"     "25"     "0"     
#[1] "SUB1" "57"   "6"    "9"    "0"   
#[1] "DDX3X" "510"   "83"    "23"    "11"   
#[1] "RPS3" "223"  "28"   "16"   "3"   
#[1] "NCBP2" "170"   "19"    "6"     "0"    
#[1] "DDX55" "87"    "15"    "3"     "0"    
#[1] "METAP2" "96"     "7"      "8"      "0" 
tmp1<-unique(other_short_intron$GID)
print(c("OSI host genes",length(tmp1),sum(tmp1%in%mitocarat$V1),
        sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
#[1] "OSI host genes" "5051"           "399"            "74"             "29"     

for (i in 1:length(Clusters[[1]])){
  RBP_name<-Clusters[[1]][i]
  tmp1<-unique(long_intron$GID[long_intron$RBP==RBP_name])
  print(c(RBP_name,length(tmp1),sum(tmp1%in%mitocarat$V1),
          sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
}
## LI
#[1] "LARP4" "801"   "118"   "49"    "14"   
#[1] "PABPC4" "661"    "86"     "68"     "9"     
#[1] "SUB1" "464"  "62"   "31"   "5"   
#[1] "DDX3X" "3107"  "354"   "47"    "34"   
#[1] "RPS3" "872"  "112"  "35"   "4"   
#[1] "NCBP2" "1554"  "133"   "20"    "3"    
#[1] "DDX55" "410"   "44"    "11"    "3"    
#[1] "METAP2" "445"    "39"     "28"     "2"   
 
### as a cluster
FLEXI<-read.delim("all.FLEXI")
tmp<-FLEXI[rowSums(FLEXI[,73:76])>0,]

for (i in c(1,5)){
  RBP_name<-Clusters[[i]]
  for (j in 1:length(RBP_name)){
    name_tmp<-RBP_name[j]
    if (j==1){
      tmp2<-unique(tmp$GID[grep(name_tmp,tmp$RBP)])
    } else {
      tmp2<-unique(c(tmp2,tmp$GID[grep(name_tmp,tmp$RBP)]))
    }
  }
  print(c(names(Clusters)[i],length(tmp2),sum(tmp2%in%mitocarat$V1),
          sum(tmp2%in%NR$V1),sum(tmp2%in%MR$V1)))
}
#[1] "Cluster1" "324"      "59"       "31"       "10"      
#[1] "Cluster5" "116"      "8"        "13"       "0"       
c(length(unique(tmp$GID)),length(unique(tmp$GID[tmp$GID%in%mitocarat$V1])),
  length(unique(tmp$GID[tmp$GID%in%NR$V1])),length(unique(tmp$GID[tmp$GID%in%MR$V1])))
## all FLEIX host gene
# [1] 3729  280   60   23
tmp1<-unique(other_short_intron$GID)
print(c("OSI host genes",length(tmp1),sum(tmp1%in%mitocarat$V1),
        sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
#[1] "OSI host genes" "5051"           "399"            "74"             "29"  
tmp1<-unique(long_intron$GID)
print(c("LI host genes",length(tmp1),sum(tmp1%in%mitocarat$V1),
        sum(tmp1%in%NR$V1),sum(tmp1%in%MR$V1)))
#[1] "LI host genes" "13166"         "970"           "95"            "73"  
all_protein<-read.table("protein.bed")
c(length(unique(all_protein$V8)),
  length(unique(all_protein$V8[all_protein$V8%in%mitocarat$V1])),
  length(unique(all_protein$V8[all_protein$V8%in%NR$V1])),
  length(unique(all_protein$V8[all_protein$V8%in%MR$V1])))
#[1] 19780  1122   101    78

