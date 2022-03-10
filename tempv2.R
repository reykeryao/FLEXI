library(matrixStats)
library(tidyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggfortify)
library(gridGraphics)
library(VennDiagram)
library(UpSetR)
library(ComplexHeatmap)
library(DESeq2)
library("Rtsne")
library("umap")
library(zinbwave)
library(cluster)
library("dendextend")
library("reshape2")
library("purrr")
library("dplyr")
library(fpc)
library(circlize)
setwd("Documents/NGS/full_length_intron/Paper/FLEXI_git/")
load("../FLEXI.image")


dat<-read.delim("all.FLEXI")
mapped_reads<-c(305.069837,251.558067,268.210336,477.543790,207.491024,
                692.091831,666.341854,713.775291,715.241521,768.433748,71.116246)
names(mapped_reads)<-colnames(dat[,82:92])
GRCh38<-read.delim("GRCh38.93.intron_deduped.tsv")
RBP<-read.delim("150_RBP_AGO_DICER.info")
Cell<-read.table("4cell_plasma_RBP.counts",col.names=c("ID","Cells"))
FLEXI<-read.table("all_FLEXI_RBP.counts",col.names=c("ID","All_FLEXI"))
Intron<-read.table("all_intron_RBP.counts",col.names=c("ID","All_Intron"))
Genome<-read.table("GRCh38_RBP.counts",col.names=c("ID","GRCh38"))
RBP<-merge(RBP,Cell,by=1,all=T)
RBP<-merge(RBP,FLEXI,by=1,all=T)
RBP<-merge(RBP,Intron,by=1,all=T)
RBP<-merge(RBP,Genome,by=1,all=T)
RBP[is.na(RBP)]<-0
RBP_info<-read.delim(pipe("cut -f 1,6 all_FLEXI_RBP_intersect.info"),col.names=c("ID","RBP"))
RBP_info<-RBP_info[RBP_info$RBP!=".",]
RBP_info<-separate(RBP_info,col="ID",into=c("IID","GID","TID","GType","TType"),
                   sep="___",remove=F)
RBP_info<-separate(RBP_info,col="IID",into=c("IID","GName"),sep="_")
#Scatter plots (Fig1B and 7B), FLEXIs and host gene counts, not RBP scatter plots
#subset to get K562, HEK293T, and Hela FLEXI counts
#Calculate CPM of FLEXIs, by all mapped reads
FLEXI_CPM<-dat[,c(8,89:91)]
FLEXI_CPM<-FLEXI_CPM[rowSums(FLEXI_CPM[,2:4])>0,]
Unfrag_total<-mapped_reads[8:10]
FLEXI_CPM[,2:4]<-t(t(FLEXI_CPM[,2:4])/Unfrag_total)
#Assign log2CPM of -10 for those with 0 count
FLEXI_CPM[FLEXI_CPM==0]<-2^-10
#Calculate CPM of Unfragmented cellular RNA
Cell_counts<-read.delim("combined_run.counts")
Cell_counts$K562<-rowSums(Cell_counts[,c(34:41)])
Cell_counts$HEK<-rowSums(Cell_counts[,c(42:49)])
Cell_counts$Hela<-rowSums(Cell_counts[,c(50:59)])
rownames(Cell_counts)<-Cell_counts$ID
Cell_counts$ID<-as.character(Cell_counts$ID)
#Subset only FLEXI host genes
GID_list<-unique(FLEXI_CPM$GID)
Cell_counts<-Cell_counts[Cell_counts$ID%in%GID_list,c(1,60:62)]
Cell_counts[,2:4]<-t(t(Cell_counts[,2:4])/Unfrag_total)
#Assign log2CPM of -10 for those with 0 count
Cell_counts[Cell_counts==0]<-2^-10

postscript("Figures/Fig1B.eps",width = 15,height=10,horizontal = F)
par(pch=16,mfrow=c(2,3),pty="s")
#FLEXIs scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#FELXIs between HEK and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,3]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#FELXIs between K562 and HEK
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#Unfrag scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#FELXIs between HEK and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,3]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#FELXIs between K562 and HEK
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
dev.off()
#remove the temporay object created for this plot
rm(list=c("Cell_counts","FLEXI_by_GID","FLEXI_CPM","cor_p","cor_s","Unfrag_total","GID_list"))

#Scatter plots in Fig7B
#Make reproducibility marker (detected in ≥ 50% of dataset)
#Make combined healthy column  
Repo<-dat[,c(1,32:34,29:31,26:28,35:37,38:47)]
Repo<-Repo[rowSums(Repo[,2:23])>0,]
Repo$BCH3_repro<-rowMedians(as.matrix(Repo[,2:4]))>0
Repo$BCH4_repro<-rowMedians(as.matrix(Repo[,5:7]))>0
Repo$BC3_repro<-rowMedians(as.matrix(Repo[,8:10]))>0
Repo$BC4_repro<-rowMedians(as.matrix(Repo[,11:13]))>0
Repo$MDA_repro<-rowMins(as.matrix(Repo[,14:15]))>0
Repo$MCF_repro<-apply(Repo[,16:23],MARGIN = 1,FUN=function(x){sort(x)[5]>0})
Repo$PatientA_H<-rowSums(Repo[,c(2:4)])
Repo$PatientB_H<-rowSums(Repo[,c(5:7)])
Repo$PatientA_C<-rowSums(Repo[,c(8:10)])
Repo$PatientB_C<-rowSums(Repo[,c(11:13)])
Repo$MDA<-rowSums(Repo[,c(14:15)])
Repo$MCF<-rowSums(Repo[,c(16:23)])
Unfrag_total<-mapped_reads[1:6]
Repo[,30:35]<-t(t(Repo[,30:35])/Unfrag_total)
Repo<-Repo[,c(1,24:35)]
dat1<-Repo[,8:13]
dat1[dat1==0]<-2^-10
Repo[,8:13]<-dat1
rm(dat1)
Repo$Combined_H<-rowMeans(Repo[,8:9])
Repo$Combined_H_repo<-Repo$BCH3_repro | Repo$BCH4_repro
Repo<-Repo[,c(1,15,2:7,14,8:13)]
Repo$ID<-as.character(Repo$ID)
Repo<-separate(Repo,ID,into=c("IID","GID"),sep = "___",remove = F,extra="drop")
Repo<-separate(Repo,IID,into=c("Intron","GName"),sep = "_",remove = T,extra="drop")
colnames(Repo)[6:9]<-c("PatientA_H_repro","PatientB_H_repro","PatientA_C_repro","PatientB_C_repro")

#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("IBC_frag.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_total<-c(50.48402,33.084889,56.983289,54.093371)
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(1,4:7)]
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-colnames(FLEXI_CPM)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                                                "PatientA_Cancer","PatientB_Cancer")


postscript("Figures/Fig7B.eps",width = 10,height=10,horizontal = F)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
FLEXI_by_GID<-Repo[!(Repo[,13]==2^-10 & Repo[,15]==2^-10),]
plot(log2(FLEXI_by_GID[,c(13,15)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientA_H==2^-10 & FLEXI_by_GID$PatientA_C>=0.05 & 
                           FLEXI_by_GID$PatientA_C_repro,c(13,15)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientA_H==2^-10 & FLEXI_by_GID$PatientA_C>=0.05 & 
                                    FLEXI_by_GID$PatientA_C_repro])
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,15],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>=0.05 & 
                           FLEXI_by_GID$PatientB_C_repro,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>=0.05 & 
                                    FLEXI_by_GID$PatientB_C_repro])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()
#remove temperary object
rm(list=c("FLEXI_by_GID","FLEXI_CPM","Frag_by_FLEXI","Repo","cor_p","cor_s","Frag_total","Unfrag_total","GID_list"))

#Fig7C hallmark geneset enrichment heatmap
#GID of detected FLEXIs from each sample used for ShinyGO
Hall<-read.delim("hallmark.txt")
Hall<-Hall[,c(1,7,2:6)]
pdf("Figures/Fig7C.pdf",height=15,width=6)
scol<-colorRampPalette(c(rev(brewer.pal(11, "RdYlBu")[3:8])))(50)
heatmap.2(as.matrix(Hall[,2:7]),labRow = Hall$Hallmark,dendrogram = "none",scale="none",margins = c(10, 20),
          density.info = "none",trace="none",symm=F,symbreaks=F,keysize=1,Colv = F,symkey=F,
          col = rev(scol))
dev.off()
#remove temperary object
rm(list=c("Hall","scol"))

#FigS1 density plot 
col=c("tomato","royalblue1","greenyellow","goldenrod","orchid","black")
name<-c("HEK","Hela","K562","UHRR")
pdf("Figures/FigS1.pdf",width=10,height=10)
par(mfrow=c(2,2))
for (i in 1:4){
  per_dat<-read.delim(gzfile(paste0(name[i],".per.info.gz")))
  if (i==1){
    plot(density(per_dat$Per[per_dat$AGO==1]),bty="n",lwd=1, ylim=c(0,0.3),xlim=c(0,100),main=name[i],
         xlab="Length (%)",col=col[1])
    legend(20,0.28,lty=c(1,1,1,1,1,1),lwd=1,col=col,
           legend = c("AGO","DICER","Core spliceosome RBPs","Other RBPs",
                      "FLEXIs w/o RBP","Other short introns"),bty="n")
    combined<-per_dat
  } else {
    plot(density(per_dat$Per[per_dat$AGO==1]),bty="n",lwd=1, ylim=c(0,0.15),xlim=c(0,100),main=name[i],
         xlab="Length (%)",col=col[1])
    combined<-rbind(combined,per_dat)
  }
  lines(density(per_dat$Per[per_dat$DICER==1]),lwd=1,col=col[2])
  lines(density(per_dat$Per[per_dat$SPLICE==1]),lwd=1,col=col[3])
  lines(density(per_dat$Per[per_dat$RBP_other==1]),lwd=1,col=col[4])
  lines(density(per_dat$Per[per_dat$nonRBP==1]),lwd=1,col=col[5])
  lines(density(per_dat$Per[per_dat$nonFLEXI==1]),lwd=1,col=col[6])
}
dev.off()
#histogram version
# temp, stll working
pdf("Figures/FigS4B_his.pdf",width=10,height=10)
par(mfrow=c(2,2))
for (i in 1:4){
  per_dat<-read.delim(gzfile(paste0(name[i],".per.info.gz")))
  if (i==1){
    p1 <- hist(per_dat$Per[per_dat$nonFLEXI!=1])
    p2 <- hist(per_dat$Per[per_dat$nonFLEXI==1])   
    plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,100),main=NA,xlab="Length (% of intron)",ylab="Reads")
    plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,100), add=T)
    legend(20,0.28,lty=c(1,1,1,1,1,1),lwd=1,col=col,
           legend = c("AGO","DICER","Core spliceosome RBPs","Other RBPs",
                      "FLEXIs w/o RBP","Other short introns"),bty="n")
    combined<-per_dat
  } else {
    plot(density(per_dat$Per[per_dat$AGO==1]),bty="n",lwd=1, ylim=c(0,0.15),xlim=c(0,100),main=name[i],
         xlab="Length (%)",col=col[1])
    combined<-rbind(combined,per_dat)
  }
  lines(density(per_dat$Per[per_dat$DICER==1]),lwd=1,col=col[2])
  lines(density(per_dat$Per[per_dat$SPLICE==1]),lwd=1,col=col[3])
  lines(density(per_dat$Per[per_dat$RBP_other==1]),lwd=1,col=col[4])
  lines(density(per_dat$Per[per_dat$nonRBP==1]),lwd=1,col=col[5])
  lines(density(per_dat$Per[per_dat$nonFLEXI==1]),lwd=1,col=col[6])
}
dev.off()

#combined cellular RNA, density plot of FLEXIs vs other short intron, Fig1C left panel
pdf("Figures/Fig1C_1.pdf")
plot(density(combined$Per[combined$nonFLEXI==1]),bty="n",lwd=1, ylim=c(0,0.1),xlim=c(0,100),main=NA,
         xlab="Length (%)",col="black")
lines(density(combined$Per[combined$nonFLEXI!=1]),lwd=1,col="red")
legend(20,0.08,lty=c(1,1),lwd=1,col=c("red","black"),
           legend = c("Other short introns","FLEXIs"),bty="n")
dev.off()
pdf("Figures/Fig1C_1h.pdf")
p1 <- hist(combined$Per[combined$nonFLEXI!=1]) 
p2 <- hist(combined$Per[combined$nonFLEXI==1])   
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,100),main=NA,xlab="Length (% of intron)",ylab="Reads")
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,100), add=T)
legend(0,140000,legend = c("FLEXIs","Other short introns"),bty="n",
       fill =c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),pch=15))
dev.off()
#remove objects
rm(list=c("combined","per_dat","i","name"))
#Density plot,Fig1C right 4 panels
pdf("Figures/Fig1C_2.pdf",width=12,height=12)
#length
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat[dat$UHRR>=1,22]),bty="n",xlim=c(0,350),lwd=1.5,
     ylim=c(0,0.04),main=NA,xlab="Intron length (bp)",col=col[1])
lines(density(dat[dat$K562>=1,22]),xlim=c(0,350),lwd=1.5,col=col[2])
lines(density(dat[dat$HEK>=1,22]),xlim=c(0,350),lwd=1.5,col=col[3])
lines(density(dat[dat$Hela>=1,22]),xlim=c(0,350),lwd=1.5,col=col[4])
lines(density(dat[dat$Plasma>=1,22]),xlim=c(0,350),lwd=1.5,col=col[5])
lines(density(GRCh38$Len),xlim=c(0,350),lwd=1.5,col=col[6])
legend(200,0.03,lty=c(1,1,1,1,1,1),lwd=1.5,col=col,legend = c("UHRR","K-562","HEK-293T","Hela S3","Plasma","GRCh38"),bty="n")
length.line <- recordPlot()
invisible(dev.off())

#GC
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat[dat$UHRR>=1,21]),bty="n",xlim=c(0,100),lwd=1.5,
     ylim=c(0,0.08),main=NA,xlab="GC (%)",col=col[1])
lines(density(dat[dat$K562>=1,21]),xlim=c(0,350),lwd=1.5,col=col[2])
lines(density(dat[dat$HEK>=1,21]),xlim=c(0,350),lwd=1.5,col=col[3])
lines(density(dat[dat$Hela>=1,21]),xlim=c(0,350),lwd=1.5,col=col[4])
lines(density(dat[dat$Plasma>=1,21]),xlim=c(0,350),lwd=1.5,col=col[5])
lines(density(GRCh38$GC),xlim=c(0,350),lwd=1.5,col=col[6])
GC.line <- recordPlot()
invisible(dev.off())

#MEF
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat[dat$UHRR>=1,20]),bty="n",xlim=c(-150,0),lwd=1.5,
     ylim=c(0,0.05),main=NA,xlab="MFE (kcal/mol)",col=col[1])
lines(density(dat[dat$K562>=1,20]),xlim=c(-150,0),lwd=1.5,col=col[2])
lines(density(dat[dat$HEK>=1,20]),xlim=c(-150,0),lwd=1.5,col=col[3])
lines(density(dat[dat$Hela>=1,20]),xlim=c(-150,0),lwd=1.5,col=col[4])
lines(density(dat[dat$Plasma>=1,20]),xlim=c(-150,0),lwd=1.5,col=col[5])
lines(density(GRCh38$MFE),xlim=c(-150,0),lwd=1.5,col=col[6])
MFE.line <- recordPlot()
invisible(dev.off())

#PhastCons
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat[dat$UHRR>=1,17]),bty="n",xlim=c(0,1),lwd=1.5,
     ylim=c(0,7),main=NA,xlab="PhastCons score",col=col[1])
lines(density(dat[dat$K562>=1,17]),lwd=1.5,col=col[2])
lines(density(dat[dat$HEK>=1,17]),lwd=1.5,col=col[3])
lines(density(dat[dat$Hela>=1,17]),lwd=1.5,col=col[4])
lines(density(dat[dat$Plasma>=1,17]),lwd=1.5,col=col[5])
lines(density(GRCh38$PhastCon30),lwd=1.5,col=col[6])
lines(density(dat[rowMaxs(as.matrix(dat[,88:92]))>0 & dat$Has_snoRNA!=".",17]),lwd=1.5,col="gray50",lty=2)
PhastCons30.line <- recordPlot()
invisible(dev.off())
ggarrange(ggarrange(length.line,GC.line,ncol=2),
          ggarrange(MFE.line,PhastCons30.line,ncol = 2),nrow = 2)
dev.off()
#remove objects
rm(list=c("GC.line","length.line",'MFE.line',"PhastCons30.line"))

#Fig1D, RPM density of FLEXIs, different group
gene_counts<-read.delim("combined_counts.tsv")
snoRNA_list<-unique(dat$Has_snoRNA[dat$Has_snoRNA!="."])
snoRNA<-gene_counts[gene_counts$Name%in%snoRNA_list,c(2,4:14)]
pdf("Figures/Fig1D.pdf",onefile = T,width=8,height=12)
par(mfrow=c(3,2),lwd=1.5)
D_height<-c(1.5,1.5,2,2,2)
for (i in c(88:92)){
    plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
         xlim=c(-4,4),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
    lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
    lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                                i]/mapped_reads[i-81])),col="black")
    if (i<92){
      lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
    }
    lines(density(log10(snoRNA[snoRNA[,i-80]>0,i-80]/mapped_reads[i-81])),col="goldenrod",lty=4)
    if (i==88){
      legend(0,1.5,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
             col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
             lty=c(1,1,1,1,4),lwd=1.5)
    }
    if (D_height[i-87]==2){
      axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
    } else {
      axis(2,labels=seq(0,1.5,0.5),las=1,at=seq(0,1.5,0.5),las=2)
    }
    
    axis(1,labels=c(parse(text='10^-4'),bquote(10^-3),bquote(10^-2),bquote(10^-1),1,10,bquote(10^2),bquote(10^3),
                    bquote(10^4)),
         at=seq(-4,4,1))
}
dev.off()
rm(list=c("snoRNA","D_height","i","snoRNA_list"))

#Fig3B Venn of mirtrons and agotrons
pdf(paste0("Figures/Fig3B.pdf"),width=12,height=6)
#venn diagram of agotrons (including mirtron)
set_1 <- as.character(dat$ID[dat$UHRR>0 & dat$Is_agotron=="+"])
set_2 <- as.character(dat$ID[dat$K562>0 & dat$Is_agotron=="+"])
set_3 <- as.character(dat$ID[dat$HEK>0 & dat$Is_agotron=="+"])
set_4 <- as.character(dat$ID[dat$Hela>0 & dat$Is_agotron=="+"])
set_5 <- as.character(dat$ID[dat$Plasma>0 & dat$Is_agotron=="+"])
set <- list ("UHRR"=set_1,"K-562"=set_2,"HEK 293T"=set_3,"Hela S3"=set_4,"Plasma"=set_5)
vennplot1 <- venn.diagram (set, filename=NULL,category.names=names(set),
                           cat.col = col[1:5], main="Agotrons",
                           fill = col[1:5],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(0,0,180,180,0),
                           main.cex=1, cat.cex = 1) 
#venn diagram of mirtrons (including agotron)
set_1 <- as.character(dat$ID[dat$UHRR>0 & dat$Is_mirtron=="+"])
set_2 <- as.character(dat$ID[dat$K562>0 & dat$Is_mirtron=="+"])
set_3 <- as.character(dat$ID[dat$HEK>0 & dat$Is_mirtron=="+"])
set_4 <- as.character(dat$ID[dat$Hela>0 & dat$Is_mirtron=="+"])
set_5 <- as.character(dat$ID[dat$Plasma>0 & dat$Is_mirtron=="+"])
set <- list ("UHRR"=set_1,"K-562"=set_2,"HEK 293T"=set_3,"Hela S3"=set_4,"Plasma"=set_5)
vennplot2 <- venn.diagram (set, filename=NULL,category.names=names(set),
                           cat.col = col[1:5], main="Mirtrons",
                           fill = col[1:5],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(0,0,180,180,0),
                           main.cex=1, cat.cex = 1) 
ggarrange(vennplot1, vennplot2,nrow = 1,ncol=2)
dev.off()
#cleanup Veen log file
unlink("*.log")
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","vennplot1","vennplot2"))

#Fig1A, upset plot of FLEXIs, FLEXI host genes from cellular RNA and plasma
#FLEXIs
cutoff<-c(0,0.01,0.05,0.1,0.5)
dat<-read.delim("all.FLEXI")
mapped_reads<-c(305.069837,251.558067,268.210336,477.543790,207.491024,
                692.091831,666.341854,713.775291,715.241521,768.433748,71.116246)
names(mapped_reads)<-colnames(dat[,82:92])
for (i in 1:5){
  set_1 <- as.character(dat$ID[dat$UHRR/mapped_reads[7]>=cutoff[i] & dat$UHRR>0])
  set_2 <- as.character(dat$ID[dat$K562/mapped_reads[8]>=cutoff[i] & dat$K562>0])
  set_3 <- as.character(dat$ID[dat$HEK/mapped_reads[9]>=cutoff[i] & dat$HEK>0])
  set_4 <- as.character(dat$ID[dat$Hela/mapped_reads[10]>=cutoff[i] & dat$Hela>0])
  set_5 <- as.character(dat$ID[dat$Plasma/mapped_reads[11]>=cutoff[i] & dat$Plasma>0])
  set <- list ("UHRR"=set_1,
               "K-562"=set_2,"HEK 293T"=set_3,
               "Hela S3"=set_4,"Plasma"=set_5)
  m = make_comb_mat(set)
  fil_name<-paste0("temp_fig/Fig1A_",cutoff[i],"RPM.eps")
  postscript(fil_name,height=4,width=8)
  UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
        comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
  dev.off()
  #FLEXI host genes
  set_1 <- unique(dat$GID[dat$UHRR/mapped_reads[7]>=cutoff[i] & dat$UHRR>0])
  set_2 <- unique(dat$GID[dat$K562/mapped_reads[8]>=cutoff[i] & dat$K562>0])
  set_3 <- unique(dat$GID[dat$HEK/mapped_reads[9]>=cutoff[i] & dat$HEK>0])
  set_4 <- unique(dat$GID[dat$Hela/mapped_reads[10]>=cutoff[i] & dat$Hela>0])
  set_5 <- unique(dat$GID[dat$Plasma/mapped_reads[11]>=cutoff[i] & dat$Plasma>0])
  set <- list ("UHRR"=set_1,
               "K-562"=set_2,"HEK 293T"=set_3,
               "Hela S3"=set_4,"Plasma"=set_5)
  m = make_comb_mat(set)
  fil_name<-paste0("temp_fig/Fig1A_",cutoff[i],"RPM_G.eps")
  postscript(fil_name,height=4,width=8)
  UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
        comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
  dev.off()
}
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","m"))

#Fig7A and FigS7, FLEXI by 0 cutoff and ≥ 0.01 RPM FLEXIs
cut_off=0.01
FLEXI_CPM<-dat[,c(1,8,82:87)]
FLEXI_CPM<-FLEXI_CPM[rowSums(FLEXI_CPM[,3:8])>0,]
FLEXI_CPM[,3:8]<-t(t(FLEXI_CPM[,3:8])/mapped_reads[1:6])
#Fig7A, cutoff 0.01 RPM, by FLEXI
postscript("Figures/Fig7A_1.eps",height=4,width=8)
set_1 <- as.character(FLEXI_CPM$ID[(FLEXI_CPM$BCH3+FLEXI_CPM$BCH4)>=cut_off])
set_2 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$BC3>=cut_off])
set_3 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$BC4>=cut_off])
set_4 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$MDA>=cut_off])
set_5 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$MCF7>=cut_off])
set <- list ("Healthy (Combined)"=set_1,
             "Patient A (Cancer)"=set_2,"Patient B (Cancer)"=set_3,
             "MDA-MB-231"=set_4,"MCF7"=set_5)
m = make_comb_mat(set)
UpSet(m,set_order=c("MDA-MB-231","MCF7","Patient A (Cancer)",
                    "Patient B (Cancer)","Healthy (Combined)"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
      comb_order=c(6,15,16,13,14,26:21,31:30,29:28,1,5:2,12:7,20:17,27))
dev.off()
postscript(paste0("Figures/Fig7A_2.eps"),height=4,width=8)
set_1 <- unique(FLEXI_CPM$GID[(FLEXI_CPM$BCH3+FLEXI_CPM$BCH4)>=cut_off])
set_2 <- unique(FLEXI_CPM$GID[FLEXI_CPM$BC3>=cut_off])
set_3 <- unique(FLEXI_CPM$GID[FLEXI_CPM$BC4>=cut_off])
set_4 <- unique(FLEXI_CPM$GID[FLEXI_CPM$MDA>=cut_off])
set_5 <- unique(FLEXI_CPM$GID[FLEXI_CPM$MCF7>=cut_off])
set <- list ("Healthy (Combined)"=set_1,
             "Patient A (Cancer)"=set_2,"Patient B (Cancer)"=set_3,
             "MDA-MB-231"=set_4,"MCF7"=set_5)
m = make_comb_mat(set)
UpSet(m,set_order=c("MDA-MB-231","MCF7","Patient A (Cancer)",
                    "Patient B (Cancer)","Healthy (Combined)"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
      comb_order=c(6,15,16,13,14,26:21,31:30,29:28,1,5:2,12:7,20:17,27))
dev.off()

#FigS7, cutoff 0 RPM, by FLEXI
cut_off=0
postscript("Figures/FigS7_1.eps",height=4,width=8)
set_1 <- as.character(FLEXI_CPM$ID[(FLEXI_CPM$BCH3+FLEXI_CPM$BCH4)>cut_off])
set_2 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$BC3>cut_off])
set_3 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$BC4>cut_off])
set_4 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$MDA>cut_off])
set_5 <- as.character(FLEXI_CPM$ID[FLEXI_CPM$MCF7>cut_off])
set <- list ("Healthy (Combined)"=set_1,
             "Patient A (Cancer)"=set_2,"Patient B (Cancer)"=set_3,
             "MDA-MB-231"=set_4,"MCF7"=set_5)
m = make_comb_mat(set)
UpSet(m,set_order=c("MDA-MB-231","MCF7","Patient A (Cancer)",
                    "Patient B (Cancer)","Healthy (Combined)"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
      comb_order=c(6,15,16,13,14,26:21,31:30,29:28,1,5:2,12:7,20:17,27))
dev.off()
postscript(paste0("Figures/FigS7_2.eps"),height=4,width=8)
set_1 <- unique(FLEXI_CPM$GID[(FLEXI_CPM$BCH3+FLEXI_CPM$BCH4)>cut_off])
set_2 <- unique(FLEXI_CPM$GID[FLEXI_CPM$BC3>cut_off])
set_3 <- unique(FLEXI_CPM$GID[FLEXI_CPM$BC4>cut_off])
set_4 <- unique(FLEXI_CPM$GID[FLEXI_CPM$MDA>cut_off])
set_5 <- unique(FLEXI_CPM$GID[FLEXI_CPM$MCF7>cut_off])
set <- list ("Healthy (Combined)"=set_1,
             "Patient A (Cancer)"=set_2,"Patient B (Cancer)"=set_3,
             "MDA-MB-231"=set_4,"MCF7"=set_5)
m = make_comb_mat(set)
UpSet(m,set_order=c("MDA-MB-231","MCF7","Patient A (Cancer)",
                    "Patient B (Cancer)","Healthy (Combined)"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
      comb_order=c(6,15,16,13,14,26:21,31:30,29:28,1,5:2,12:7,20:17,27))
dev.off()
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","m","cut_off","FLEXI_CPM"))

#Fig8A-D, upset plot of oncogene TSG in BC dataset
BC_FLEXI<-dat[,c(1:25,93,82:92)]
Onco<-read.delim("OncoGenes.table")
Onco$Onco<-1
Onco<-Onco[,c(2,8)]
colnames(Onco)[1]<-"GName"
TSG<-read.delim("TSGs.table")
TSG$TSG<-1
TSG<-TSG[,c(2,9)]
colnames(TSG)[1]<-"GName"
BC_FLEXI[,27:37]<-t(t(BC_FLEXI[,27:37])/mapped_reads)
BC_FLEXI$Healthy<-(BC_FLEXI$BCH3+BC_FLEXI$BCH4)/2
temp<-BC_FLEXI[,27:38]
temp[temp==0]<-2^-10
BC_FLEXI[,27:38]<-log2(temp)
BC_FLEXI<-merge(BC_FLEXI,Onco,by="GName",all.x=T)
BC_FLEXI<-merge(BC_FLEXI,TSG,by="GName",all.x=T)
BC_FLEXI[is.na(BC_FLEXI)]<-0
pdf("Figures/Fig8A_D.pdf",height=8,width=11)
set1<-BC_FLEXI[((BC_FLEXI$BC3-BC_FLEXI$BCH3)>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
up1 <- upset(fromList(set),number.angles = 45,mainbar.y.label = "Oncogene FLEXIs",
             sets.x.label = "Oncogene FLEXIs (Up)",keep.order = T,sets=rev(names(set)),
             main.bar.color = c(rep("tomato",4),rep("royalblue1",6),rep("goldenrod",3),"orchid"))
up1 <- cowplot::plot_grid(NULL, up1$Main_bar, up1$Sizes, up1$Matrix,
                          nrow=2, align='hv', rel_heights = c(3,1),
                          rel_widths = c(2,3))
set1<-BC_FLEXI[(BC_FLEXI$BCH3-BC_FLEXI$BC3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BCH4-BC_FLEXI$BC4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MDA>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MCF7>=1) & BC_FLEXI$Onco==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
up2 <- upset(fromList(set),number.angles = 45,mainbar.y.label = "Oncogene FLEXIs",
             sets.x.label = "Oncogene FLEXIs (Down)",keep.order = T,sets=rev(names(set)),
             main.bar.color = c(rep("tomato",4),rep("royalblue1",4),rep("goldenrod",3),"orchid"))
up2 <- cowplot::plot_grid(NULL, up2$Main_bar, up2$Sizes, up2$Matrix,
                          nrow=2, align='hv', rel_heights = c(3,1),
                          rel_widths = c(2,3))
set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$TSG==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$TSG==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$TSG==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$TSG==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
up3 <-upset(fromList(set),number.angles = 45,mainbar.y.label = "TSG FLEXIs",
            sets.x.label = "TSG FLEXIs (Up)",keep.order = T,sets=rev(names(set)),
            main.bar.color = c(rep("tomato",4),rep("royalblue1",6),rep("goldenrod",3),"orchid"))
up3 <- cowplot::plot_grid(NULL, up3$Main_bar, up3$Sizes, up3$Matrix,
                          nrow=2, align='hv', rel_heights = c(3,1),
                          rel_widths = c(2,3))
set1<-BC_FLEXI[(BC_FLEXI$BCH3-BC_FLEXI$BC3>=1) & BC_FLEXI$TSG==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BCH4-BC_FLEXI$BC4>=1) & BC_FLEXI$TSG==1,2]
set3<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MDA>=1) & BC_FLEXI$TSG==1,2]
set4<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MCF7>=1) & BC_FLEXI$TSG==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
up4 <- upset(fromList(set),number.angles = 45,mainbar.y.label = "TSG FLEXIs",
             sets.x.label = "TSG FLEXIs (Down)",keep.order = T,sets=rev(names(set)),
             main.bar.color = c(rep("tomato",3),rep("royalblue1",5),rep("goldenrod",3),"orchid"))
up4 <- cowplot::plot_grid(NULL, up4$Main_bar, up4$Sizes, up4$Matrix,
                          nrow=2, align='hv', rel_heights = c(3,1),
                          rel_widths = c(2,3))
grid.arrange(up1, up2, up3,up4,nrow = 2,ncol=2)
dev.off()
rm(list=c("temp","dat1","up1","up2","up3","up4","set1","set2","set3","set4","set"))

#FIg3D, RBP scatter plot in conserved FELXIs, phastCons ≥ 0.99, various cutoff 0.5, 0.75, 0.8,0.9,0.95,0.99
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>0,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP","Cell")
cutoff<-c(0.5,0.75,0.8,0.9,0.95,0.99)
pdf("temp_fig//Fig3D2.pdf",height=12,width=8)
par(mfrow=c(3,2),pty="s",pch=16)
for (i in 1:6){
  Phast99<-dat$ID[dat$PhastCon30>=cutoff[i] & rowSums(dat[,88:91])>0
                  & dat$Is_agotron=="." & dat$Is_mirtron=="."
                  & dat$Has_snoRNA=="."]
  Phast99_fre<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%Phast99,]
  Phast99_fre<-data.frame(table(Phast99_fre$RBP))
  colnames(Phast99_fre)<-c("RBP","Phast99")
  Phast<-merge(RBP_fre,Phast99_fre,by=1,all=T)
  Phast[is.na(Phast)]<-0
  rownames(Phast)<-Phast$ID
  Phast$pvalue<-1
  R_sum<-colSums(Phast[,2:3])
  for (j in 1:126){
    Phast[j,4]<-fisher.test(as.matrix(rbind(Phast[j,2:3],R_sum)))$p.value
  }
  Phast$padj<-p.adjust(Phast$pvalue,method="fdr")
  Phast[,2:3]<-data.frame(prop.table(as.matrix(Phast[,2:3]),margin = 2)*100)
  Phast<-merge(Phast,RBP[,c(1,4,5,11)],by=1)
  Phast$col<-(Phast$Splicing.regulation+Phast$Spliceosome)/3+Phast$microRNA.processing
  plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
       main=paste0("phastCons>=",cutoff[i],", n=",length(Phast99)),
       xlab="All FELXIs (% RBP sites)",ylab=paste0("Conserved FELXIs (% RBP sites)"))
  points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
  points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
  points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  text(Phast[Phast$padj<=0.05 & Phast$Phast99>=2,c(2,3)],pos = 3,cex=0.5,
       labels = Phast$RBP[Phast$padj<=0.05 & Phast$Phast99>=2])
  
#  plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
#       xlab="All FELXIs (% RBP sites)",ylab="Conserved FELXIs (% RBP sites)")
#  points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
#  points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
#  points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
#  abline(0,1,col="red")
#  text(Phast[Phast$padj<=0.05 & Phast$Phast99>=2,c(2,3)],pos = 3,cex=0.5,
#       labels = Phast$RBP[Phast$padj<=0.05 & Phast$Phast99>=2])
}
dev.off()

rm(list=c("Phast","Phast99"))
#Fig4B
pdf("Figures/Fig4B.pdf",height=6,width=12)
par(mfrow=c(1,2),pty="s",pch=16)
RBP_fre<-RBP[,c(46,48,49)]
RBP_fre$All_Intron<-RBP$All_Intron-RBP$All_FLEXI
RBP_fre<-data.frame(prop.table(as.matrix(RBP_fre),margin = 2)*100)
RBP_fre$Name<-RBP$RBP.name
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-RBP_fre[,c(2,1,3,4,5)]
plot(RBP_fre[RBP_fre$col==0,1:2],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,1:2],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,1:2],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,1:2],col="skyblue",cex=1.5)
text(RBP_fre[RBP_fre$Cells>4,1:2],labels = RBP_fre$Name[RBP_fre$Cells>4])
abline(0,1,col="red")

plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="All RBP binding sites (%)",ylab="FLEXIs (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[RBP_fre$Cells>4,c(3,2)],labels = RBP_fre$Name[RBP_fre$Cells>4])
dev.off()
rm(RBP_fre)
#Fig4A and S3
#All RBPs by occurrence in FLEXIs (log scale)
pdf("Figures/FigS3.pdf",height=12,width=12)
Fun<-read.table("4cell_plasma_RBP_by_FLEXI.counts",col.names=c("Name","RBP_by_FLEXI"))
Fun<-merge(Fun,RBP[,c(1,4,5,11)],by=1,all.x=T)
Fun<-Fun[order(Fun$RBP_by_FLEXI,decreasing = T),]
Fun$Color<-(Fun$Splicing.regulation+Fun$Spliceosome)/3+Fun$microRNA.processing
Fun[Fun$Color==1,6]<-2
Fun[Fun$Color==4/3,6]<-3
Fun[Fun$Color==1/3 | Fun$Color==2/3,6]<-1
Fun$Color<-factor(Fun$Color)
B_col=c("black","red","orange","skyblue")
mp<-barplot(Fun$RBP_by_FLEXI,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[Fun$Color],
            names.arg =Fun$Name,las=2)
legend("top",pch=16,col=B_col[c(2,3,4,1)],legend = c("RNA splicing","miRNA related","Both","Other"),bty="n")
par(xpd=T)
legend("topright",legend = Fun$Name,text.col=B_col[Fun$Color],cex=0.3,bty="n")
par(xpd=F)
dev.off()
#RBPs by occurrence in FLEXIs (≥30) log scale
pdf("Figures/Fig4A.pdf",height=12,width=12)
mp<-barplot(Fun$RBP_by_FLEXI[Fun$RBP_by_FLEXI>=30],ylim=c(1,10000),log="y",
            cex.names=0.5,col=B_col[Fun$Color[Fun$RBP_by_FLEXI>=30]],
            names.arg =Fun$Name[Fun$RBP_by_FLEXI>=30],las=2)
legend("top",pch=16,col=B_col[c(2,3,4,1)],legend = c("RNA splicing","miRNA related","Both","Other"),bty="n")
par(xpd=T)
legend("topright",legend = Fun$Name[Fun$RBP_by_FLEXI>=30],text.col=B_col[Fun$Color],cex=0.3)
par(xpd=F)
dev.off()
rm(list=c("mp","B_col"))

#Fig6 bottom panel
pick<-c(113,12,32,9,94,5,117,111,68,151,96,131,102,124,49,69,85,30,93,87,74,81,105,20,145,101,24,
        2,76,146,29,150,136,47,134,130,132,119,91,19,135,11,53,86,63,118,72,88,103,100,51,71,48)
#redefine SG protein from MSGP database
RBP_53<-RBP[pick,1:43]
RBP_53$StressGranule<-0
SG<-read.table("SG_MSGP_RGD.list",col.names="ID")
RBP_53[RBP_53$RBP.name%in%SG$ID,44]<-1
RBP_53<-RBP_53[,c(1:31,44,32:43)]
RBP_53[is.na(RBP_53)]<-0
RBP_53<-RBP_53[,c(TRUE,TRUE,colSums(RBP_53[,3:44])>0)]
RBP_53[6,27]<-1
pdf("Figures/Fig6.pdf",width=12,height=8)
par(mfrow=c(3,1),mar = c(5,2,2,20))
image(1:53,1:16,as.matrix(RBP_53[,18:3]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP_53)[18:3],tick = FALSE)
image(1:53,1:14,as.matrix(RBP_53[,32:19]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP_53)[32:19],tick = FALSE)
image(1:53,1:6,as.matrix(RBP_53[38:33]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:6,labels = colnames(RBP_53)[38:33],tick = FALSE)
axis(1,las=2,at = 1:53,labels = RBP_53$RBP.name,tick = FALSE,cex=0.5)
dev.off()
rm(RBP_53)
rm(pick)
rm(SG)
#FigS5
FLEXI_ID<-dat$ID[rowSums(dat[,88:92])>0]
RBP_info_4cell<-RBP_info[RBP_info$ID%in%FLEXI_ID,]
FourCellPlasma<-read.delim("4cell_plasma_FLEXI.tsv")
pdf("Figures/FigS5.pdf",width=8,height=5)
Splicesome<-c("SF3B4","PRPF8","EFTUD2","BUD13","AQR")
set_1 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[1]])
set_2 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[2]])
set_3 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[3]])
set_4 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[4]])
set_5 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[5]])
RBP_3<-c("AATF","DKC1","NOLC1")
set_6 <- union(set_1,union(set_2,union(set_3,union(set_4,set_5))))
set_7 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[1]])
set_8 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[2]])
set_9 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[3]])
set <- list (set_6,set_7,set_8,set_9)
names(set)<-c("Splicesomal proteins",RBP_3)
upset(fromList(set),keep.order=T,set_size.show = T,
      mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",
      sets=names(set),nintersects = 1000,nsets = 4)
dev.off()
#Fig5
RBP_30<-c("AGO1-4","DICER","BCLAF1","DDX3X","G3BP1","PABPN1","YBX3","ZNF622")
Title_name<-c("(51%)","(44%)","(43%)","(40%)","(41%)","(43)","(52%)","(33%)")
Title_name<-paste(RBP_30,Title_name)
plot_l <- list()
for (i in 1:8){
  if (i==2){
    set_6 <- unique(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
    set <- list (set_1,set_2,set_3,set_4,set_5,set_6)
    names(set)<-c(Splicesome,RBP_30[i])
    vp <-upset(fromList(set),keep.order=T,set_size.show = T,
               mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",mainbar.y.max=600,
               sets=names(set),nintersects = 1000,nsets = 6)
    vp <- plot_grid(NULL, vp$Main_bar, vp$Sizes, vp$Matrix,
                    nrow=2, align='hv', rel_heights = c(3,1),
                    rel_widths = c(1,3))
    plot_l[[RBP_30[i]]] <-vp
  } else if (i==1){
    set_6 <- unique(FourCellPlasma$ID[FourCellPlasma$AGO_CCR!="."])
    set <- list (set_1,set_2,set_3,set_4,set_5,set_6)
    names(set)<-c(Splicesome,RBP_30[i])
    vp <-upset(fromList(set),keep.order=T,set_size.show = T,
               mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",mainbar.y.max=600,
               sets=names(set),nintersects = 1000,nsets = 6)
    vp <- plot_grid(NULL, vp$Main_bar, vp$Sizes, vp$Matrix,
                    nrow=2, align='hv', rel_heights = c(3,1),
                    rel_widths = c(1,3))
    plot_l[[RBP_30[i]]] <-vp
  } else {
    set_6 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_30[i]])
    set <- list (set_1,set_2,set_3,set_4,set_5,set_6)
    names(set)<-c(Splicesome,RBP_30[i])
    vp <-upset(fromList(set),keep.order=T,set_size.show = T,
               mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",mainbar.y.max=600,
               sets=names(set),nintersects = 1000,nsets = 6)
    vp <- plot_grid(NULL, vp$Main_bar, vp$Sizes, vp$Matrix,
                    nrow=2, align='hv', rel_heights = c(3,1),
                    rel_widths = c(1,3))
    plot_l[[RBP_30[i]]] <-vp
  }
}
pdf("Figures/Fig5.pdf",width=16,height=22)
grid.arrange(arrangeGrob(plot_l[[RBP_30[1]]],top=Title_name[1]),
             arrangeGrob(plot_l[[RBP_30[2]]],top=Title_name[2]),
             arrangeGrob(plot_l[[RBP_30[3]]],top=Title_name[3]),
             arrangeGrob(plot_l[[RBP_30[4]]],top=Title_name[4]),
             arrangeGrob(plot_l[[RBP_30[5]]],top=Title_name[5]),
             arrangeGrob(plot_l[[RBP_30[6]]],top=Title_name[6]),
             arrangeGrob(plot_l[[RBP_30[7]]],top=Title_name[7]),
             arrangeGrob(plot_l[[RBP_30[8]]],top=Title_name[8]),
             ncol = 2,nrow=4)
dev.off()
rm(list=c("set_1","set_2","set_3","set_4","set_5","set_6","set_7","set_8","set_9",
          "RBP_3","RBP_30","Splicesome","Title_name","vp","set","i","plot_l","FLEXI_ID"))
#Fig8E
#RBP of BC(onco)
set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
RBP_info<-RBP_info[,c(1,8)]
RBP_info<-unique(RBP_info)
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))
AGO_sites<-c(dim(dat[dat$ID%in%MCF_unique & dat$AGO_CCR!=".",])[1],
             dim(dat[dat$ID%in%MDA_unique & dat$AGO_CCR!=".",])[1],
             dim(dat[dat$ID%in%All_cancer & dat$AGO_CCR!=".",])[1])
DICER_sites<-c(dim(dat[dat$ID%in%MCF_unique & dat$DICER_CCR!=".",])[1],
               dim(dat[dat$ID%in%MDA_unique & dat$DICER_CCR!=".",])[1],
               dim(dat[dat$ID%in%All_cancer & dat$DICER_CCR!=".",])[1])

MDA_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MDA_unique]))
colnames(MDA_unique)<-c("RBP","MDA_unique")
MCF_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MCF_unique]))
colnames(MCF_unique)<-c("RBP","MCF_unique")
All_cancer<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%All_cancer]))
colnames(All_cancer)<-c("RBP","All_cancer")
Cell<-RBP_fre
colnames(Cell)<-c("RBP","Cell")
Four_cancer_RBP<-merge(MCF_unique,MDA_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,All_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Cell,by=1,all=T)
Four_cancer_RBP[is.na(Four_cancer_RBP)]<-0
Four_cancer_RBP$RBP<-as.character(Four_cancer_RBP$RBP)

Four_cancer_RBP[Four_cancer_RBP$RBP=="AGO",2:4]<-AGO_sites
Four_cancer_RBP[Four_cancer_RBP$RBP=="DICER",2:4]<-DICER_sites

pdf("Figures/Fig7D.pdf",height=6,width=4)
par(mfrow=c(3,2),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
Name=c("MCF7","MDA-MB-231","All cancer")
for (i in 1:3){
  dat_set<-Four_cancer_RBP[,c(1,i+1,5)]
  RBP_sum<-colSums(dat_set[,2:3])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)])))$p.value})
  dat_set$C_sig<-p.adjust(dat_set$C_sig,method = "fdr")
  dat_set[,2:3]<-100*prop.table(as.matrix(dat_set[,2:3]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name[i],axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  sig_set<-dat_set[dat_set$C_sig<=0.05 & (dat_set[,3]>=2 |dat_set[,2]>=2),]
  if (dim(sig_set)[1]>0){
    points(sig_set[,c(3,2)],col="red")
    text(sig_set[,c(3,2)],
         labels = sig_set[,1],pos = 3,cex=0.7)
  }
    if (i<3) {
    axis(1,at=seq(0,25,5),labels = NA)
  } else {
    axis(1,at=seq(0,25,5),labels = seq(0,25,5))
    title(xlab="All FLEXIs (% RBP sites)")
  }
  if (i%%2==1){
    axis(2,at=seq(0,25,5),labels = seq(0,25,5))
    title(ylab="Oncogene FLEXIs (% RBP sites)")
  } else {
    axis(2,at=seq(0,25,5),labels = NA)
  }
}
dev.off()
#FigS6
FourCellFLEXI<-dat$ID[rowSums(dat[,88:91])>0]
FourCell<-dat[dat$ID%in%FourCellFLEXI,]

pdf("temp_fig/Fig5.pdf",width=11,height=14,onefile = T)
RBP_list_sig<-c("DKC1","NOLC1","AATF","AGO","DICER","highP")
par(mfrow=c(6,4))
for (i in 1:length(RBP_list_sig)){
  # subset
  if (RBP_list_sig[i] =="AGO") {
    temp1<-FourCell[FourCell$AGO_CCR!=".",]
    temp2<-FourCell[FourCell$AGO_CCR==".",]
  } else if (RBP_list_sig[i] =="DICER") {
    temp1<-FourCell[FourCell$DICER_CCR!=".",]
    temp2<-FourCell[FourCell$DICER_CCR==".",]
  } else if (RBP_list_sig[i] =="highP") {
    temp1<-FourCell[FourCell$PhastCon30>=0.8 & FourCell$Has_snoRNA=="." &
                      FourCell$Is_agotron=="." & FourCell$Is_mirtron==".",]
    temp2<-FourCell[!FourCell$ID%in%temp1$ID,]
  } else {
    temp1<-FourCell[grep(RBP_list_sig[i],FourCell$RBP),]
    temp2<-FourCell[grep(RBP_list_sig[i],FourCell$RBP,invert=T),]
  }
  Len_t<-0
  GC_t<-0
  MFE_t<-0
  PhastCon30_t<-0
  rep_times<-100
  for (inter in 1:rep_times){
    sample_size<-dim(temp1)[1]
    test_set<-FourCell[sample(1:8098,sample_size,replace = F),]
    if (wilcox.test(temp1$Len,test_set$Len,exact = F)$p.value>0.05) {Len_t=Len_t+1}
    if (wilcox.test(temp1$GC,test_set$GC,exact = F)$p.value>0.05) {GC_t=GC_t+1}
    if (wilcox.test(temp1$MFE,test_set$MFE,exact = F)$p.value>0.05) {MFE_t=MFE_t+1}
    if (wilcox.test(temp1$PhastCon30,test_set$PhastCon30,exact = F)$p.value>0.05) {PhastCon30_t=PhastCon30_t+1}
  }
  if (Len_t>0 | GC_t >0 | MFE_t >0 | PhastCon30_t >0 ) {
    #Length
    plot(density(temp1$Len),bty="n",xlim=c(0,350),lwd=1.5,
         ylim=c(0,0.015),main=NA,xlab="Intron length (bp)",col="red")
    lines(density(temp2$Len),xlim=c(0,350),lwd=1.5,col="black")
    legend(120,0.012,lty=c(1,1),lwd=1.5,col=c("red","black"),
           legend = c(paste0(RBP_list_sig[i]," (",dim(temp1)[1],")"),"Others"),bty="n")
    if (Len_t==0) {
      legend("topleft",legend="FDR < 0.01",bty="n") 
    } else {
      legend("topleft",legend=paste0("FDR = ",format(Len_t/rep_times,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
    #GC
    plot(density(temp1$GC),bty="n",xlim=c(0,100),lwd=1.5,
         ylim=c(0,0.05),main=NA,xlab="GC (%)",col="red")
    lines(density(temp2$GC),xlim=c(0,100),lwd=1.5,col="black")
    if (GC_t==0) {
      legend("topleft",legend="FDR < 0.01",bty="n") 
    } else {
      legend("topleft",legend=paste0("FDR = ",format(GC_t/rep_times,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
    #MEF
    plot(density(temp1$MFE),bty="n",xlim=c(-150,0),lwd=1.5,
         ylim=c(0,0.03),main=NA,xlab="Minimal free energy (MFE; kcal/mol)",col="red")
    lines(density(temp2$MFE),xlim=c(-150,0),lwd=1.5,col="black")
    if (MFE_t==0) {
      legend("topleft",legend="FDR < 0.01",bty="n") 
    } else {
      legend("topleft",legend=paste0("FDR = ",format(MFE_t/rep_times,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
    #PhastCons
    plot(density(temp1$PhastCon30),bty="n",xlim=c(0,1),lwd=1.5,
         ylim=c(0,7),main=NA,xlab="PhastCons score",col="red")
    lines(density(temp2$PhastCon30),lwd=1.5,col="black")
    if (PhastCon30_t==0) {
      legend("topleft",legend="FDR < 0.01",bty="n") 
    } else {
      legend("topleft",legend=paste0("FDR = ",format(PhastCon30_t/rep_times,digits=2,nsmall=2,big.mark = ".")),bty="n") 
    }
  }
}
dev.off()



#Barplot of profile in all sample (combiend)
Cell_counts<-read.delim("combined_run.counts")
Cell_counts$K562<-rowSums(Cell_counts[,c(34:41)])
Cell_counts$HEK<-rowSums(Cell_counts[,c(42:49)])
Cell_counts$Hela<-rowSums(Cell_counts[,c(50:59)])
Cell_counts$UHRR<-rowSums(Cell_counts[,c(26:33)])
Cell_counts$MCF7<-rowSums(Cell_counts[,c(18:25)])
Cell_counts$MDA<-rowSums(Cell_counts[,c(16:17)])
Cell_counts$BC3H<-rowSums(Cell_counts[,c(10:12)])
Cell_counts$BC4H<-rowSums(Cell_counts[,c(7:9)])
Cell_counts$BC3<-rowSums(Cell_counts[,c(4:6)])
Cell_counts$BC4<-rowSums(Cell_counts[,c(13:15)])
Frag_by_FLEXI<-read.delim("IBC_frag.counts")
dat<-Cell_counts[,c(1:3,60:69)]
dat<-merge(dat,Frag_by_FLEXI,by=1:3,all=T)
dat[is.na(dat)]<-0
dat<-dat[rowSums(dat[,4:17])>0,]
dat$Type<-as.factor(dat$Type)
dat<-dat[dat$Type!="ERCC",]
dat<-dat[dat$Type!="SP",]
levels(dat$Type)[1]<-levels(dat$Type)[2]<-levels(dat$Type)[3]<-levels(dat$Type)[4]<-levels(dat$Type)[15]<-"rRNA"
levels(dat$Type)[14]<-"Pseudogene"
levels(dat$Type)[10]<-"misc RNA"
levels(dat$Type)[8]<-"Other lncRNA"
levels(dat$Type)[21]<-"VT RNA"
levels(dat$Type)[22]<-"Y RNA"
levels(dat$Type)[6]<-levels(dat$Type)[11]<-levels(dat$Type)[13]<-levels(dat$Type)[19]<-"Protein coding"
levels(dat$Type)[13]<-"snoRNA"
levels(dat$Type)[10]<-"MT tRNA"
dat$Type<-droplevels(dat$Type)
dat$Type2<-dat$Type
levels(dat$Type)[2]<-levels(dat$Type)[3]<-levels(dat$Type)[8]<-"sncRNA"
levels(dat$Type)[7]<-levels(dat$Type)[10]<-levels(dat$Type)[11]<-"sncRNA"
levels(dat$Type)[10]<-levels(dat$Type)[11]<-"sncRNA"
rownames(dat)<-dat$ID
agg<-aggregate(.~Type,data=dat[,c(3:17)],FUN=sum)
agg.prob<-data.frame(prop.table(as.matrix(agg[,c(2:15)]),2))
agg.prob$Type=agg$Type
agg.prob<-agg.prob[c(4,8,3,5,6,2,9,7,1),]
agg.prob$Type<-as.character(agg.prob$Type)
pdf("../Pass1_IGV/bar_all_genes.pdf")
scol <- brewer.pal(9, "Set3")
mp<-barplot(as.matrix(agg.prob[,1:14]*100),col=scol,axes=F,
            names.arg=rep(NA,14),width=0.05,space=0.3,
            legend.text=agg.prob$Type,
            adj=0.12,args.legend=list(x=1.3,y=75,bty="n"),xlim=c(0,1.4))
axis(1,labels=NA,at=c(0,mp,0.95))
text(mp+0.05, par("usr")[3]-3, labels=colnames(agg.prob[,1:14]),
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()

dat<-read.delim("Protein.info.alldataset")
agg<-dat[,6:7]
agg.prob<-data.frame(prop.table(as.matrix(agg),1))
pdf("bar_protein_sense_anti.pdf")
mp<-barplot(t(agg.prob*100),col=scol,axes=F,
            names.arg=rep(NA,14),width=0.05,space=0.3,
            legend.text=c("Sense","Antisense"),
            adj=0.12,args.legend=list(x=1.3,y=75,bty="n"),xlim=c(0,1.4))
axis(1,labels=NA,at=c(0,mp,0.95))
text(mp+0.05, par("usr")[3]-3, labels=dat$Name,
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()

agg<-dat[,2:5]
agg.prob<-data.frame(prop.table(as.matrix(agg),1))
pdf("bar_protein_CDS.pdf")
mp<-barplot(t(agg.prob*100),col=rev(scol),axes=F,
            names.arg=rep(NA,14),width=0.05,space=0.3,
            legend.text=c("CDS","UTR","Intron","Intergenic"),
            adj=0.12,args.legend=list(x=1.3,y=75,bty="n"),xlim=c(0,1.4))
axis(1,labels=NA,at=c(0,mp,0.95))
text(mp+0.05, par("usr")[3]-3, labels=dat$Name,
     srt=45, pos=2,xpd=TRUE,cex=0.7)
axis(2,labels=c(0,25,50,75,100),las=1,at=c(0,25,50,75,100),pos=0,cex.axis=0.7)
dev.off()

# genebody cov
pdf("GenebodyCov.pdf")
par(mfrow=c(2,1))
dat<-read.delim("combined.genebodycov.per",row.names=1)
scol <- brewer.pal(9, "Set1")
plot(dat$UHRR,xlim=c(1,100),ylim=c(0,1),col=scol[1],type = "l",bty="n",xlab="Gene body percentile (5'->3')",ylab="Coverage")
lines(dat$K562,col=scol[2])
lines(dat$HEK,col=scol[3])
lines(dat$Hela,col=scol[4])
lines(dat$MDA,col=scol[5])
lines(dat$MCF7,col=scol[7])
legend(40,0.9,col=scol[c(1:5,7)],lty = 1,
       legend = c("UHRR","K 562","HEK 293T","HeLa S3","MCF7","MDA-MB-231"),bty="n")

plot(dat$BCH3,xlim=c(1,100),ylim=c(0,1),col=scol[1],type = "l",bty="n",xlab="Gene body percentile (5'->3')",ylab="Coverage")
lines(dat$BCH4,col=scol[2])
lines(dat$BC3,col=scol[3])
lines(dat$BC4,col=scol[4])
lines(dat$BCH3F,col=scol[1],lty=2,lwd=1.2)
lines(dat$BCH4F,col=scol[2],lty=2,lwd=1.2)
lines(dat$BC3F,col=scol[3],lty=2,lwd=1.2)
lines(dat$BC4F,col=scol[4],lty=2,lwd=1.2)
legend(40,0.8,col=scol[c(1:4,1:4)],lty = rep(c(1,2),each=4),
       legend = c("Patient A (H)","Patient B (H)","Patient A (C)","Patient B (C)",
                  "Patient A (H)","Patient B (H)","Patient A (C)","Patient B (C)"),bty="n")
dev.off()

dat<-read.delim("../Pass1_IGV/polyA.len")
pdf("../Pass1_IGV/polyA.pdf")
plot(density(rep(dat$length,dat$UHRR)),bty="n",ylim=c(0,1),
     xlim=c(0,100),col=scol[1],main=NA,xlab="Length of polyA tail (nt)")
lines(density(rep(dat$length,dat$K562)),col=scol[2])     
lines(density(rep(dat$length,dat$HEK)),col=scol[3])  
lines(density(rep(dat$length,dat$Hela)),col=scol[4])  
legend(70,0.4,lty=1,col=scol[1:4],legend = c("UHRR","K-562","HEK 293T","HeLa S3"),bty="n")
dev.off()



dat<-read.delim("all.FLEXI",row.names = 1)
dat<-dat[,25:80]
dat$Type="FELXI"
Cell_counts<-read.delim("combined_run.counts")
Cell_counts<-Cell_counts[,4:59]
Cell_counts$Type="Other"
dat<-rbind(Cell_counts,dat)
coldata <- data.frame(rep(c("BC3","BCH4","BCH4","BC4","MDA","MCF","UHRR","K562","HEK","Hela"),
                            times=c(3,3,3,3,2,8,8,8,8,10)),
                      colnames(dat[,1:56]),
                      row.names=colnames(dat[,1:56]))
colnames(coldata) <- c("Dataset","Name")
#FLEXI DE
FLEXI_dds <- DESeqDataSetFromMatrix(countData = dat[,1:56],
                                   colData = coldata,design = ~Dataset)
FLEXI_dds <- DESeq(FLEXI_dds,parallel = T)
scol=c("red","black")

Compare=c("MCF","Hela")
FLEXI_res <- data.frame(results(FLEXI_dds,contrast = c("Dataset",Compare), alpha = 0.05))
FLEXI_res$Type=dat$Type
FLEXI_res$Type=factor(FLEXI_res$Type)
FLEXI_res<-FLEXI_res[complete.cases(FLEXI_res),]
plot(FLEXI_res$log2FoldChange, -log10(FLEXI_res$padj),col=scol[FLEXI_res$Type],
     main=paste0(Compare[1]," vs ",Compare[2]),
     pch=16,bty="n",xlab="log2(FC)",ylab="-log10 (adjusted p-value)")

#Fig1E
pdf("temp_fig/Phast.pdf")
plot(density(dat$PhastCon30[dat$Is_mirtron!="."]),ylim=c(0,7),
     col=scol[1],bty="n",xlab="phastCons",main=NA,xlim=c(0,1))
lines(density(dat$PhastCon30[dat$Is_agotron!="."]),col=scol[2])
lines(density(dat$PhastCon30[dat$Has_snoRNA!="."]),col=scol[3])
lines(density(dat$PhastCon30[dat$Is_agotron=="." & dat$Has_snoRNA=="." & dat$Is_mirtron=="."]),
      col=scol[4])
lines(density(GRCh38$PhastCon30[!GRCh38$ID%in%dat$ID]),col="black")
legend(0.6,5,col=c(scol[1:4],"black"),
       legend = c("Mirtron","Agotron","SnoRNA FLEXI","Other FLEXI RNAs",
                  "Other short introns"),lty=1,bty="n")
dev.off()


#new fig 1C

pdf("Figures/Fig1C_2.pdf",width=12,height=12)
col=c("red","orchid","black")
#length
pdf(NULL)
dev.control(displaylist="enable")
dat1<-dat[rowSums(dat[,88:92])>0,]
plot(density(dat1[,22]),bty="n",xlim=c(0,350),lwd=1.5,
     ylim=c(0,0.04),main=NA,xlab="Intron length (bp)",col=col[1])
lines(density(dat1[dat1$Plasma>=1,22]),xlim=c(0,350),lwd=1.5,col=col[2])
lines(density(GRCh38$Len[!GRCh38$ID%in%dat1$ID]),xlim=c(0,350),lwd=1.5,col=col[3],lty=2)
#lines(density(GRCh38$Len),xlim=c(0,350),lwd=1.5,col=col[6])
legend(200,0.03,lty=c(1,1,2),lwd=1.5,col=col,
       legend = c("FLEXIs (Cellular)","FLEXIs (Plasma)","Other short introns"),bty="n")
length.line <- recordPlot()
invisible(dev.off())

#GC
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat1[,21]),bty="n",xlim=c(0,100),lwd=1.5,
     ylim=c(0,0.08),main=NA,xlab="GC (%)",col=col[1])
lines(density(dat1[dat1$Plasma>=1,21]),xlim=c(0,350),lwd=1.5,col=col[2])
lines(density(GRCh38$GC[!GRCh38$ID%in%dat1$ID]),xlim=c(0,350),lwd=1.5,col=col[3],lty=2)
GC.line <- recordPlot()
invisible(dev.off())

#MEF
pdf(NULL)
dev.control(displaylist="enable")
plot(density(dat1[,20]),bty="n",xlim=c(-150,0),lwd=1.5,
     ylim=c(0,0.05),main=NA,xlab="MFE (kcal/mol)",col=col[1])
lines(density(dat1[dat1$Plasma>=1,20]),xlim=c(-150,0),lwd=1.5,col=col[2])
lines(density(GRCh38$MFE[!GRCh38$ID%in%dat1$ID]),xlim=c(-150,0),lwd=1.5,col=col[3],lty=2)
MFE.line <- recordPlot()
invisible(dev.off())

ggarrange(ggarrange(length.line,GC.line,ncol=2),
          ggarrange(MFE.line),nrow = 2)
dev.off()
#remove objects
rm(list=c("GC.line","length.line",'MFE.line'))

#RBP comparason between FLEXI and non-FLEXI
dat1<-RBP
dat1$All_FLEXI<-dat1$All_FLEXI-dat1$Cells
dat1$Color<-(dat1$Splicing.regulation+dat1$Spliceosome)/3+dat1$microRNA.processing
dat1<-dat1[,c(1,2,50,46,47)]
dat1[dat1$Color==1,3]<-2
dat1[dat1$Color==4/3,3]<-3
dat1[dat1$Color==1/3 | dat1$Color==2/3,3]<-1
dat1$Color<-factor(dat1$Color)
B_col=c("black","red","orange","skyblue")
dat1<-dat1[rowSums(dat1[,4:5])>0,]
dat1[,4:5]<-prop.table(as.matrix(dat1[,4:5]),margin = 2)*100
pdf("../Pass1_IGV/RBP_new.pdf")
par(pty="s")
plot(dat1[,4:5],pch=16,xlim=c(0,15),ylim=c(0,15),bty="n",
     xlab="FLEXIs (Cellular)",ylab="Other short introns")
abline(0,1,col="red")
text(5,12,"r = 0.98")
text(5,11,"rs = 0.91")
dev.off


#try t-sne

dat1<- dat[,26:81]
sub_mapped_reads<-c(89264248,83532383,95413705,79734689,81242936,90580442,113477444,
                    106080387,85512006,282545532,93703609,101294649,109704771,97786253,
                    70563696,83971287,93530133,96188889,84195614,95958499,90154958,
                    77528755,87438356,80725874,84193148,84695000,79516751,94126936,
                    80171312,75474477,87388516,82600414,84957488,88141802,95927677,
                    113807174,87425522,73526698,96566478,71023706,90569418,84985720,
                    107210447,92506910,88853873,83524969,91268315,83392799,82750963,
                    62981188,60057226,67704164,86607803,72132500,96955981,64582809)
dat1<-1e6*t(t(dat1/sub_mapped_reads))
rownames(dat1)<-dat$ID
dat1<-data.frame(t(dat1))
dat1$label<-c(rep(1:4,each=3),rep(5,2),rep(6:9,each=8),rep(10,10))
dat1$label<-as.factor(dat1$label)
colors = rainbow(length(unique(dat1$label)))
names(colors) = unique(dat1$label)
pick=c(13:56)
dat2<-dat1[pick,1:8687]
cut_off_title<-c(0,0.05,0.25,0.5,0.75,0.95,0.99,0.999)
cut_off<-quantile(colSums(dat2),cut_off_title)
pdf("../Pass1_IGV/t-SNE_cutoff.pdf",height=12,width=12)
par(pty="s",mfrow=c(4,4))
col=c("black","orchid","tomato","royalblue1","greenyellow","goldenrod")
col<-col2hex(col)
col<-paste0(col,"A0")
name<-c("MDA-MB-231","MCF-7","UHRR","K-562","HEK 293T","HeLa S3")
for (i in 1:length(cut_off)){
  dat2<-dat1[pick,1:8687]
  dat2<-dat2[,colSums(dat2)>=cut_off[i]]
  dat2$label<-dat1$label[pick]
  dat2$label<-droplevels(dat2$label)
  tsne <- Rtsne(dat2[,1:dim(dat2)[2]-1], dims = 2, perplexity=3, theta = 0.5,normalize=F,
                num_threads=4,verbose=TRUE, max_iter = 5000)
  plot(tsne$Y,main=paste0("t-SNE ",cut_off_title[i]),col=col[dat2$label],
       xlab=NA,ylab=NA,pch=16)
  if (i==1){
    legend(round(quantile(tsne$Y[,1],0.1)/10)*10,
           round(quantile(tsne$Y[,2],0.8)/10)*10,
           legend = name, col=col,pch=16,bty="n")
  }
  pca<-prcomp(dat2[,1:dim(dat2)[2]-1])
  plot(pca$x[,1:2],main=paste0("PCA ",cut_off_title[i]),col=col[dat2$label],
       xlab=NA,ylab=NA,pch=16)
}
dev.off()

pick=c(13:56)
dat2<-dat1[pick,1:8687]
dat2<-dat2[,colSums(dat2)>0]
dat2$label<-dat1$label[pick]
dat2$label<-droplevels(dat2$label)
tsne <- Rtsne(dat2[,1:dim(dat2)[2]-1], dims = 2,
              perplexity=3, theta = 0.5,normalize=F,
              num_threads=4,verbose=TRUE, max_iter = 5000)
pdf("../Pass1_IGV/t-SNE_RPM.pdf",height=4,width=4)
par(pty="s")

col=c("black","orchid","tomato","royalblue1","greenyellow","goldenrod")
col<-col2hex(col)
col<-paste0(col,"A0")
name<-c("MDA-MB-231","MCF-7","UHRR","K-562","HEK 293T","HeLa S3")
plot(tsne$Y,  main="t-SNE",col=col[dat2$label],
     xlab="t-SNE1",ylab="t-SNE2",pch=16)
legend(round(quantile(tsne$Y[,1],0.1)/10)*10,
       round(quantile(tsne$Y[,2],0.8)/10)*10,
       legend = name, col=col,pch=16,bty="n")
dev.off()

#PCA
pdf("../Pass1_IGV/PCA.pdf",height=4,width=4)
par(pty="s")
pca<-prcomp(dat2[,1:dim(dat2)[2]-1])
plot(pca$x[,1:2],main="Cell lines",col=col[dat2$label],
     xlab="PC1",ylab="PC2",pch=16)
legend(0,6,
       legend = name, col=col,pch=16,bty="n")
dev.off()


pdf("../Pass1_IGV/Umap.pdf",height=4,width=4)
par(pty="s")
custom.settings = umap.defaults
custom.settings$n_neighbors<-6
custom.settings$min_dist<-0.01
FLEXIumap<-umap(dat2[,1:8391],config=custom.settings)
FLEXIumap<-data.frame(FLEXIumap$layout)
FLEXIumap$label<-dat2$label
plot(FLEXIumap[,1:2],main="UMAP",col=col[FLEXIumap$label],
     xlab="W1",ylab="W2",pch=16)
dev.off()

pick=c(13:56)
dat2<-dat[,38:81]
dat2<-dat2[rowSums(dat2)>0,]
Fc<-SummarizedExperiment(assays = list(counts = as.matrix(dat2)),
                     colData = data.frame(label=c(rep(5,2),rep(6:9,each=8),rep(10,10))))
zinb<-zinbwave(Fc)
W <- data.frame(reducedDim(zinb))
col=c("black","orchid","tomato","royalblue1","greenyellow","goldenrod")
col<-col2hex(col)
col<-paste0(col,"A0")
name<-c("MDA-MB-231","MCF-7","UHRR","K-562","HEK 293T","HeLa S3")
W$label<-as.factor(Fc$label)
pdf("../Pass1_IGV/ZINB-WaVE.pdf",width=4,height=4)
par(pty="s")
plot(W$W1,W$W2,main="ZINB-WaVE",col=col[W$label],
     xlab="W1",ylab="W2",pch=16)
#legend(1,0.5,legend = name, col=col,pch=16,bty="n")
dev.off()

'''
pick=c(1:22)
dat2<-dat1[pick,]
dat2$label<-droplevels(dat2$label)
tsne <- Rtsne(dat2[,1:8687], dims = 2, perplexity=5, theta = 0.5,normalize=F,
              num_threads=4,
              verbose=TRUE, max_iter = 5000)
plot(tsne$Y, t='n', main="All breast cancer + Healthy")
text(tsne$Y, labels=dat2$label, col=colors[dat2$label])

pick=c(1:12)
dat2<-dat1[pick,]
dat2$label<-droplevels(dat2$label)
tsne <- Rtsne(dat2[,1:8687], dims = 2, perplexity=2, theta = 0.5,normalize=F,
              num_threads=4,
              verbose=TRUE, max_iter = 5000)
plot(tsne$Y, t='n', main="Patients")
text(tsne$Y, labels=dat2$label, col=colors[dat2$label])

pick=c(1:3,10:22)
dat2<-dat1[pick,]
dat2$label<-droplevels(dat2$label)
tsne <- Rtsne(dat2[,1:8687], dims = 2, perplexity=2, theta = 0.5,normalize=F,
              num_threads=4,
              verbose=TRUE, max_iter = 5000)
plot(tsne$Y, t='n', main="Breast cancers only")
text(tsne$Y, labels=dat2$label, col=colors[dat2$label])

dev.off()
'''
#density check for FLEXI with or w/o RBP sites
dat1<-dat[,c(1:25,93,88:91)]
dat1[,27:30]<-t(t(dat1[,27:30]/mapped_reads[7:10]))
dat2<-dat1[,27:30]
dat2[dat2==0]<-2^-10
dat1[,27:30]<-dat2
pdf("../Pass1_IGV/RPM_RBP.pdf")
par(bty="n",mfrow=c(2,2))
dat2<-log10(dat1$UHRR[dat1$RBP=="."&dat1$UHRR>2^-10])
plot(density(dat2),type="l",xlim=c(-4,1),ylim=c(0,1.5),main="UHRR",xlab="log10 (RPM)")
dat2<-log10(dat1$UHRR[!dat1$RBP=="."&dat1$UHRR>2^-10])
lines(density(dat2),col="red")
legend(-3,1.4,legend = c("without RBP site","with RBP site"),col=c("black","red"),lty=1,bty="n")

dat2<-log10(dat1$K562[dat1$RBP=="."&dat1$K562>2^-10])
plot(density(dat2),type="l",xlim=c(-4,1),ylim=c(0,1.5),main="K-562",xlab="log10 (RPM)")
dat2<-log10(dat1$K562[!dat1$RBP=="."&dat1$K562>2^-10])
lines(density(dat2),col="red")

dat2<-log10(dat1$HEK[dat1$RBP=="."&dat1$HEK>2^-10])
plot(density(dat2),type="l",xlim=c(-4,1),ylim=c(0,1.5),main="HEk 293T",xlab="log10 (RPM)")
dat2<-log10(dat1$HEK[!dat1$RBP=="."&dat1$HEK>2^-10])
lines(density(dat2),col="red")

dat2<-log10(dat1$Hela[dat1$RBP=="."&dat1$Hela>2^-10])
plot(density(dat2),type="l",xlim=c(-4,1),ylim=c(0,1.5),main="HeLa S3",xlab="log10 (RPM)")
dat2<-log10(dat1$Hela[!dat1$RBP=="."&dat1$Hela>2^-10])
lines(density(dat2),col="red")
dev.off()

#5percentile vs 95 percentile RBP distribution
dat1<-dat[,c(1,88:91)]
pdf("../Pass1_IGV/RBP_quantile.pdf")
par(mfrow=c(2,2),pty="s")
for (i in 2:5){
  perc<-quantile(dat1[dat1[,i]>0,i],0.95)
  RBP95<-dat1$ID[dat1[,i]>perc]
  RBP95<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%RBP95]))
  RBP05<-data.frame(table(RBP_info$RBP[!RBP_info$ID%in%RBP95]))
  RBP95<-merge(RBP95,RBP05,by=1,all=T)
  RBP95[is.na(RBP95)]<-0
  RBP95[,2:3]<-100*prop.table(as.matrix(RBP95[,2:3]),margin = 2)
  plot(RBP95[,2:3],pch=16,xlim=c(0,20),ylim=c(0,20),xlab="5% quantile",ylab="Others",
       main=colnames(dat1)[i])
  abline(0,1,col="red")
}
dev.off()

#BPA position
pdf("../Pass1_IGV/BPA.pdf")
#par(mfrow=c(2,1))
dat1<-read.table("../../JA20150_combined/Branch_point/meme_out40/fimo_bpa/bpa.pos")
dat2<-read.table("../../JA20150_combined/Branch_point/meme_out40/non_FLEXI_bpa/bpa.pos")
plot(density(dat1$V1),xlim=c(0,40),ylim=c(0,0.15),bty="n",xlab="Distance of BPA to 3'SS (nt)",main=NA,col="red")
lines(density(dat2$V1),col="black")
legend(0,0.1,legend = c("FLEXIs","Other short introns"),col=c("red","black"),lty=1,bty="n")
#dat1<-read.table("../../JA20150_combined/Branch_point/meme_out_nonFELXI40/fimo_bpa/bpa.pos")
#dat2<-read.table("../../JA20150_combined/Branch_point/meme_out_nonFELXI40/non_FLEXI_bpa/bpa.pos")
#plot(density(dat1$V1),xlim=c(0,40),ylim=c(0,0.1),bty="n",xlab="Distance of BPA to 3'SS (nt)",main=NA,col="red")
#lines(density(dat2$V1),col="black")
dev.off()

##New ID with scale of sncRNA
##Fig1D, RPM density of FLEXIs, different group
gene_counts<-read.delim("combined_counts.tsv")
gene_counts[gene_counts$Type=="scaRNA",3]<-"snoRNA"
snoRNA<-gene_counts[gene_counts$Type=="snoRNA",c(2,8:14)]
snoRNA<-separate(snoRNA,col = "Name",into = "Name",sep = "-",remove = T)
snoRNA$Name<-gsub(snoRNA$Name,pattern = "[A-Z]$|P[1-9]$",replacement="")
snoRNA$Name[snoRNA$Name=="U3"]<-"SNORD3"
snoRNA$Name[snoRNA$Name=="U8"]<-"SNORD118"
snoRNA$Name[snoRNA$Name=="snoU13"]<-"SNORD13"
snoRNA$Name[snoRNA$Name=="snoU2_19"]<-"snoU2"
snoRNA<-aggregate(.~Name,data=snoRNA,FUN = sum)


snRNA<-gene_counts[gene_counts$Type=="snRNA",c(2,8:13)]
U7<-log10(colSums(snRNA[grep("U7\\b|U7[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
U11<-log10(colSums(snRNA[grep("U11\\b|U11[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
SNORD74<-log10(colSums(snoRNA[snoRNA$Name=="SNORD74",2:7])/mapped_reads[5:10])
SNORD78<-log10(colSums(snoRNA[snoRNA$Name=="SNORD78",2:7])/mapped_reads[5:10])
RN7SK<-log10(colSums(gene_counts[gene_counts$Type=="7SK",8:13])/mapped_reads[5:10])
RN7SL<-log10(colSums(gene_counts[gene_counts$Type=="7SL",8:13])/mapped_reads[5:10])
YRNA<-log10(colSums(gene_counts[gene_counts$Type=="YRNA",8:13])/mapped_reads[5:10])
VTRNA<-log10(colSums(gene_counts[gene_counts$Type=="VTRNA",8:13])/mapped_reads[5:10])
RMRP=log10(colSums(gene_counts[gene_counts$Name=="RMRP",8:13])/mapped_reads[5:10])
RPPH1=log10(colSums(gene_counts[gene_counts$Name=="RPPH1",8:13])/mapped_reads[5:10])
snc<-data.frame(rbind(SNORD74,SNORD78,U7,U11,YRNA,VTRNA,RN7SK,RN7SL,RMRP,RPPH1))

pdf("Figures/Fig1D.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(2,2,2,2,2)
for (i in c(88:91)){
  plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,6),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-81])),col="black")
  if (i<92){
    lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
  }
  lines(density(log10(snoRNA[snoRNA[,i-86]>0,i-86]/mapped_reads[i-81])),col="goldenrod",lty=4)
  snc<-snc[order(snc[,i-87],decreasing = F),]
  for (j in (1:dim(snc)[1])){
    if (j==1){
      segments(x0=snc[j,i-87],y0=1.5,y1=1.8,lty=1)
      text (snc[j,i-87],1.85,rownames(snc)[j],cex=0.7,adj=0)
    } else {
      segments(x0=snc[j,i-87],y0=0.5,y1=1.8-j*0.12,lty=1)
      text (snc[j,i-87],1.85-j*0.12,rownames(snc)[j],cex=0.7,adj=0)
    }
  }
  if (i==88){
    legend(0.5,2,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5,cex=0.7)
  }
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-2),1,bquote(10^2),bquote(10^4),
                  bquote(10^6)),at=seq(-4,6,2))
}
dev.off()

pdf("temp_fig/Fig1F_MDA_MCF.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(3,3)
for (i in c(86:87)){
  plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,6),ylim=c(0,D_height[i-85]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-81])),col="black")
  if (i<92){
    lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
  }
  lines(density(log10(snoRNA[snoRNA[,i-84]>0,i-84]/mapped_reads[i-81])),col="goldenrod",lty=4)
  snc<-snc[order(snc[,i-85],decreasing = F),]
  for (j in (1:dim(snc)[1])){
    if (j==1){
      segments(x0=snc[j,i-85],y0=1.5,y1=1.8,lty=1)
      text (snc[j,i-85],1.85,rownames(snc)[j],cex=0.7,adj=0)
    } else {
      segments(x0=snc[j,i-85],y0=0.5,y1=1.8-j*0.12,lty=1)
      text (snc[j,i-85],1.85-j*0.12,rownames(snc)[j],cex=0.7,adj=0)
    }
  }
  if (i==86){
    legend(0.5,2,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5,cex=0.7)
  }
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-2),1,bquote(10^2),bquote(10^4),
                  bquote(10^6)),at=seq(-4,6,2))
}
dev.off()


gene_counts<-read.delim("combined_counts.tsv")
temp<-gene_counts[gene_counts$Type=="snRNA",c(2,10:14)]
gene_for_correlation<-data.frame("U1"=colSums(temp[grep("U1\\b|U1[A-Z]",temp$Name),2:6]))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U2"=colSums(temp[grep("U2\\b|U2[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4"=colSums(temp[grep("U4\\b|U4[B-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U5"=colSums(temp[grep("U5\\b|U5[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6"=colSums(temp[grep("U6\\b|U6[B-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U7"=colSums(temp[grep("U7\\b|U7[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U11"=colSums(temp[grep("U11\\b|U11[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U12"=colSums(temp[grep("U12\\b|U12[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4ATAC"=colSums(temp[grep("U4[ATAC|atac]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6ATAC"=colSums(temp[grep("U6[ATAC|atac]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SK"=colSums(gene_counts[gene_counts$Type=="7SK",10:14])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SL"=colSums(gene_counts[gene_counts$Type=="7SL",10:14])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RPPH1"=colSums(gene_counts[gene_counts$Name=="RPPH1",10:14])))
temp<-gene_counts[gene_counts$Type=="snoRNA",c(2,10:14)]
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD3"=colSums(temp[grep("SNORD3\\b|SNORD3[A-Z]|U3",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD118"=colSums(temp[grep("SNORD118\\b|SNORD118[A-Z]|U8",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD13"=colSums(temp[grep("SNORD13\\b|SNORD13[A-Z]|U13",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD14"=colSums(temp[grep("SNORD14\\b|SNORD14[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD22"=colSums(temp[grep("SNORD22\\b|SNORD22[A-Z]",temp$Name),2:6])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RMRP"=colSums(gene_counts[gene_counts$Name=="RMRP",10:14])))
gene_for_correlation<-gene_for_correlation/mapped_reads[7:11]
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,1e4,5e3,2e3,2e3,2e5,5e5,2e5,
                                                   2e5,4e4,1e4,1e4,1e4,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",4),"7SK","7SL","RNase P",
                            rep("C/D box snoRNA",5),"MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)
pdf("../Pass1_IGV/copy_cor.pdf",height=8,width=8)
par(bty="n",pty="s",mfrow=c(2,2))
scol <- brewer.pal(8, "Set1")
y<-log10(gene_for_correlation[,6])
for (i in 1:4){
  x<-log10(gene_for_correlation[,i])
  lm.out <- lm(y ~ x)
  newx = seq(min(x),max(x),by = 0.05)
  conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95)
  plot(x,y,pch=16,xlim=c(0,10),ylim=c(0,10),
       col=scol[gene_for_correlation$Type],ylab=bquote(log[10]~(molecule/cell)),
       xlab=bquote(log[10] (RPM)),
       main=colnames(gene_for_correlation)[i])
  abline(lm.out, col="lightblue")
  matlines(newx, conf_interval[,2:3], col = "blue", lty=2)
  if (i==1){
    legend(5,4,legend = levels(gene_for_correlation$Type),col=scol,bty="n",pch=16,cex=0.7)
  }
  cor_s=formatC(cor(x,y,method = "spearman"),digits=2, format="f")
  cor_p=formatC(cor(x,y,method = "pearson"),digits=2, format="f")
  text(1,8,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
}
dev.off()

# calculate expected RPM for FLEXIs used in digital PCR
check_list<-c(2232,5094,7563,6003,5726,2183,549,1489,3377)
RPM_eve<-FourCell[rownames(FourCell)%in%check_list,c(1,88:91)]
RPM_eve<-RPM_eve[c(4,6,9,8,7,3,1,2,5),]
y<-log10(gene_for_correlation[,6])
for (i in 1:4){
  x<-log10(gene_for_correlation[,i])
  xlm.out <- lm(y ~ x)
  newx = RPM_eve[,i+1]/mapped_reads[i+6]
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95))
  t(t(round(10^conf_interval$fit,0)))
}
# Venn for Alan

pdf(paste0("Figures/Fig3B.pdf"),width=12,height=6)
#venn diagram of agotrons (including mirtron)
set_1 <- as.character(dat$ID[dat$UHRR>0])
set_2 <- as.character(dat$ID[dat$K562>0])
set_3 <- as.character(dat$ID[dat$HEK>0])
set_4 <- as.character(dat$ID[dat$Hela>0])
set_5 <- as.character(dat$ID[dat$Plasma>0])
set <- list ("UHRR"=set_1,"K-562"=set_2,"HEK 293T"=set_3,"Hela S3"=set_4,"Plasma"=set_5)
vennplot1 <- venn.diagram (set, filename=NULL,category.names=names(set),
                           cat.col = col[1:5],
                           fill = col[1:5],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(0,0,180,180,0),
                           main.cex=1, cat.cex = 1) 
#venn diagram of mirtrons (including agotron)
set_1 <- as.character(dat$ID[dat$MDA>0])
set_2 <- as.character(dat$ID[dat$MCF7>0])
set_3 <- as.character(dat$ID[dat$BC3>0])
set_4 <- as.character(dat$ID[dat$BC4>0])
set_5 <- as.character(dat$ID[dat$Plasma>0])
set <- list ("MDA-MB-231"=set_1,"MCF7"=set_2,"PatientA"=set_3,"PatientB"=set_4,"Plasma"=set_5)
vennplot2 <- venn.diagram (set, filename=NULL,category.names=names(set),
                           cat.col = col[1:5],
                           fill = col[1:5],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(0,0,180,180,0),
                           main.cex=1, cat.cex = 1) 
ggarrange(vennplot1, vennplot2,nrow = 1,ncol=2)
dev.off()
#cleanup Veen log file
unlink("*.log")

all_CPM<-dat[,c(1:25,93,82:92)]
all_CPM[,27:37]<-t(t(all_CPM[,27:37])/mapped_reads)
common<-intersect(set_5,intersect(set_1,intersect(set_2,intersect(set_3,set_4))))
non_common<-setdiff(set_5,union(set_1,union(set_2,union(set_3,set_4))))
write.table(rbind(all_CPM[all_CPM$ID%in%common,c(1,7,10,11,13,14,17,20:26,29:32,37)],
                  all_CPM[all_CPM$ID%in%non_common,c(1,7,10,11,13,14,17,20:26,29:32,37)]),
            "../Old/Table_Plasma_commonCancer.tsv",quote=F,row.names=F,sep="\t")
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","vennplot1","vennplot2"))

#reproduciable
dat<-read.delim("all.FLEXI")
dat1<-read.delim("../../Old_UF_fragmented/full_length_intron.counts")
dat<-dat[rowSums(dat[,38:81])>0,c(1,38:81)]
colnames(dat)<-c("ID",paste0("MDA_",1:2),paste0("MCF_",1:8),
                 paste0("UHRR_",1:8),paste0("K561_2_",1:8),
                 paste0("HEK_2_",1:8),paste0("Hela",1:10))
colnames(dat1)<-c("ID",paste0("HEK_1_",1:2),paste0("K562_1_",1:7))
dat<-merge(dat,dat1,by="ID",all=T)
dat<-dat[,c(1:19,36:45,20:35,48:54,46:47)]
dat[is.na(dat)]<-0
sub_mapped_reads<-c(109704771,97786253,
                    70563696,83971287,93530133,96188889,84195614,95958499,90154958,77528755,
                    87438356,80725874,84193148,84695000,79516751,94126936,80171312,75474477,
                    91268315,83392799,82750963,62981188,60057226,67704164,86607803,72132500,96955981,64582809,
                    87388516,82600414,84957488,88141802,95927677,113807174,87425522,73526698,
                    96566478,71023706,90569418,84985720,107210447,92506910,88853873,83524969,
                    5700094,7258578,6925758,5408190,6342889,6941957,9282711,
                    99564097,81230668)
sub_mapped_reads<-sub_mapped_reads/1e6
'''
dat$K562B2_repro<-dat$K562B1_repro<-dat$HEKB2_repro<-dat$HEKB1_repro<-0
dat$K562B2_over<-dat$K562B1_over<-dat$HEKB2_over<-dat$HEKB1_over<-0
cutoff=3
#more than half of the dataset has FLEXIs  greater/equal than cutoff
dat[rowSums(dat[,2:3]>=cutoff)>1,27]<-1
dat[rowSums(dat[,4:11]>=cutoff)>3,28]<-1
dat[rowSums(dat[,12:18]>=cutoff)>3,29]<-1
dat[rowSums(dat[,19:26]>=cutoff)>3,30]<-1
#total number of reads greater/equal than cutoff
dat[rowSums(dat[,2:3])>=cutoff,31]<-1
dat[rowSums(dat[,4:11])>=cutoff,32]<-1
dat[rowSums(dat[,12:18])>=cutoff,33]<-1
dat[rowSums(dat[,19:26])>=cutoff,34]<-1
col=c("tomato","royalblue1","greenyellow","goldenrod","orchid","black")
pdf("../../Old_UF_fragmented/Venn_FLEXI.pdf")
set_1 <- as.character(dat$ID[dat$HEKB1_repro>0])
set_2 <- as.character(dat$ID[dat$HEKB2_repro>0])
set_3 <- as.character(dat$ID[dat$K562B1_repro>0])
set_4 <- as.character(dat$ID[dat$K562B2_repro>0])
set_5 <- as.character(dat$ID[dat$HEKB1_over>0])
set_6 <- as.character(dat$ID[dat$HEKB2_over>0])
set_7 <- as.character(dat$ID[dat$K562B1_over>0])
set_8 <- as.character(dat$ID[dat$K562B2_over>0])

set_repo <- list ("HEK 293T (Batch1)"=set_1,"HEK 293T (Batch2)"=set_2,
             "K-562 (Batch1)"=set_3,"K-562 (Batch2)"=set_4)
set_over <- list ("HEK 293T (Batch1)"=set_5,"HEK 293T (Batch2)"=set_6,
                  "K-562 (Batch1)"=set_7,"K-562 (Batch2)"=set_8)
vennplot1 <- venn.diagram (set_repo[1:2], filename=NULL,category.names=names(set_repo)[1:2],
                           cat.col = col[1:2],
                           fill = col[1:2],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(180,180),
                           main.cex=1, cat.cex = 1) 
vennplot2 <- venn.diagram (set_repo[3:4], filename=NULL,category.names=names(set_repo)[3:4],
                           cat.col = col[1:2],
                           fill = col[1:2],cat.pos=c(180,180),
                           height = 300, width = 300, units = "px",
                           cex = 1,
                           main.cex=1, cat.cex = 1) 
vennplot3 <- venn.diagram (set_over[1:2], filename=NULL,category.names=names(set_repo)[1:2],
                           cat.col = col[1:2],
                           fill = col[1:2],cat.pos=c(180,180),
                           height = 300, width = 300, units = "px",
                           cex = 1,
                           main.cex=1, cat.cex = 1) 
vennplot4 <- venn.diagram (set_over[3:4], filename=NULL,category.names=names(set_repo)[3:4],
                           cat.col = col[1:2],
                           fill = col[1:2],cat.pos=c(180,180),
                           height = 300, width = 300, units = "px",
                           cex = 1,
                           main.cex=1, cat.cex = 1) 

ggarrange(vennplot1, vennplot2,vennplot3, vennplot4,nrow = 2,ncol=2)

dev.off()
#cleanup Veen log file
unlink("*.log")
'''

#renumbering based on cutoff, cutoff is based on combined date (rowSums), 
# but applied to individual replicates.
# MDA 2:3
# MCF 4:11
# UHRR 12:19
# Hela 20:29
# K562 batch 2 (newer) 30:37;  batch 1 46:52
# HEK batch 2 (newer) 38:45; batch 1 53:54
#cut off by RPM, 0.01,0.05,0.1,and 0.5
cutoff=c(0,0.01)
P_pch=c(16,17)
dat<-read.delim("all.FLEXI")
dat1<-read.delim("../../Old_UF_fragmented/full_length_intron.counts")
dat<-dat[rowSums(dat[,38:81])>0,c(1,38:81)]
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

for (i in 1:2){
  pdf_name<-paste0("temp_fig/cluster",cutoff[i],"RPM.pdf")
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

#alternate Fig1D without 1 read FLEXI

pdf("temp_fig/Fig1EV2.pdf",onefile = T,width=8,height=12)
par(mfrow=c(3,2),lwd=1.5)
D_height<-c(2,2,2,2,2)
for (i in c(88:92)){
  plot(density(log10(dat[dat[,i]>1 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,4),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>1 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>1 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-81])),col="black")
  if (i<92){
    lines(density(log10(dat[dat[,i]>1 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
  }
  lines(density(log10(snoRNA[snoRNA[,i-80]>1,i-80]/mapped_reads[i-81])),col="goldenrod",lty=4)
  if (i==88){
    legend(0,1.5,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5)
  }
  if (D_height[i-87]==2){
    axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  } else {
    axis(2,labels=seq(0,1.5,0.5),las=1,at=seq(0,1.5,0.5),las=2)
  }
  
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-3),bquote(10^-2),bquote(10^-1),1,10,bquote(10^2),bquote(10^3),
                  bquote(10^4)),
       at=seq(-4,4,1))
}
dev.off()


pdf("temp_fig/Fig1EV1.pdf",onefile = T,width=8,height=12)
par(mfrow=c(3,2),lwd=1.5)
D_height<-c(2,2,2,2,2)
for (i in c(88:92)){
  plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,4),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-81])),col="black")
  if (i<92){
    lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
  }
  lines(density(log10(snoRNA[snoRNA[,i-80]>0,i-80]/mapped_reads[i-81])),col="goldenrod",lty=4)
  if (i==88){
    legend(0,1.5,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5)
  }
  if (D_height[i-87]==2){
    axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  } else {
    axis(2,labels=seq(0,1.5,0.5),las=1,at=seq(0,1.5,0.5),las=2)
  }
  
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-3),bquote(10^-2),bquote(10^-1),1,10,bquote(10^2),bquote(10^3),
                  bquote(10^4)),
       at=seq(-4,4,1))
}
dev.off()



# new fig4B/C/D

# cell: 4cell lines+plasma FLEXI
# all FLEXI: all short introns
# all intron: all short and long introns
# GRCh38: all positions
RBP_fre<-RBP[,c(46:49)]
# now all intron means all long introns
RBP_fre$All_Intron<-RBP_fre$All_Intron-RBP_fre$All_FLEXI
# now all FLEXI are all other short introns (not in 4cell + plasma)
RBP_fre$All_FLEXI<-RBP$All_FLEXI-RBP$Cells
R_sum<-colSums(RBP_fre)
#exact fisher test pvalue
# SvF short intron vs FLEXI, 2 vs 1
# LvF long intron vs FLEXI, 3 vs 1
# AvF all RBP sites  vs FLEXI, 4 vs 1
# AvL all RBP sites  vs long intron, 4 vs 3
RBP_fre$SvFpvlue<-1
RBP_fre$LvFpvlue<-1
RBP_fre$AvFpvlue<-1
RBP_fre$AvLpvlue<-1
for (i in 1:152){
  RBP_fre[i,5]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(1,2)],R_sum[c(1,2)]),2,2))$p.value
  RBP_fre[i,6]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(1,3)],R_sum[c(1,3)]),2,2))$p.value
  RBP_fre[i,7]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(1,4)],R_sum[c(1,4)]),2,2))$p.value
  RBP_fre[i,8]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(3,4)],R_sum[c(3,4)]),2,2))$p.value
}
RBP_fre[,1:4]<-data.frame(prop.table(as.matrix(RBP_fre[,1:4]),margin = 2)*100)
RBP_fre$Name<-RBP$RBP.name
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre[RBP_fre$col>1,10]<-4
RBP_fre[RBP_fre$col==1,10]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,10]<-2
RBP_fre[RBP_fre$col==0,10]<-1
Splicesome<-c("SF3B4","PRPF8","EFTUD2","BUD13","AQR","PPIG")
col<-c("black","red","orange","skyblue")

pdf("temp_fig/Fig4BCD.pdf",height=6,width=12)
par(mfcol=c(2,4))
par(pch=16,pty="s")
plot(RBP_fre[,c(2,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,1]>=4)),c(2,1)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,2]>=4)])
abline(0,1,col="red")
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(2,1)],col="red",pos = 2,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])
     
plot(RBP_fre[,c(2,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,1]>=1) ),c(2,1)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$SvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,2]>=1) ]],
     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,2]>=1) ])
abline(0,1,col="red")

plot(RBP_fre[,c(3,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,3]>=4 | RBP_fre[,1]>=4)),c(3,1)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,3]>=4)])
abline(0,1,col="red")
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(3,1)],col="red",pos = 4,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])

plot(RBP_fre[,c(3,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
text(RBP_fre[(RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,3]>=1 | RBP_fre[,1]>=1) ),c(3,1)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$LvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,3]>=1) ]],
     labels = RBP_fre$Name[RBP_fre$LvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,3]>=1) ])
abline(0,1,col="red")

plot(RBP_fre[,c(4,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="All RBP binding sites (%)",ylab="FLEXIs (% RBP sites)")
abline(0,1,col="red")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,1]>=4)),c(4,1)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,4]>=4)])
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(4,1)],col="red",pos = 4,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])


plot(RBP_fre[,c(4,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="All RBP binding sites (%)",ylab="FLEXIs (% RBP sites)")
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$AvFpvlue<=0.05 & (RBP_fre[,4]>=1 | RBP_fre[,1]>=1)),c(4,1)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$AvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,4]>=1)]],
     labels = RBP_fre$Name[RBP_fre$AvFpvlue<=0.05 & (RBP_fre[,1]>=1 | RBP_fre[,4]>=1)])

plot(RBP_fre[,c(4,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="All RBP binding sites (%)",ylab="Long introns (% RBP sites)")
abline(0,1,col="red")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,3]>=4)),c(4,3)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,3]>=4 | RBP_fre[,4]>=4)])
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(4,3)],col="red",pos = 4,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])

plot(RBP_fre[,c(4,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="All RBP binding sites (%)",ylab="Long introns (% RBP sites)")
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$AvLpvlue<=0.05 & (RBP_fre[,4]>=1 | RBP_fre[,3]>=1)),c(4,3)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$AvLpvlue<=0.05 & (RBP_fre[,3]>=1 | RBP_fre[,4]>=1)]],
     labels = RBP_fre$Name[RBP_fre$AvLpvlue<=0.05 & (RBP_fre[,3]>=1 | RBP_fre[,4]>=1)])
dev.off()
#rm(RBP_fre)

# new fig4B/C/D_ver2

# cell: 4cell lines+plasma FLEXI
# all FLEXI: all short introns
# all intron: all short and long introns
# GRCh38: all positions
RBP_fre<-RBP[,c(46:49)]
# now all intron means all long introns
RBP_fre$All_Intron<-RBP_fre$All_Intron-RBP_fre$All_FLEXI
# now all FLEXI are all other short introns (not in 4cell + plasma)
RBP_fre$All_FLEXI<-RBP$All_FLEXI-RBP$Cells
R_sum<-colSums(RBP_fre)
#exact fisher test pvalue
# LvF long intron vs FLEXI, 3 vs 1
# LvS long intron vs other short introns, 3 vs 2
RBP_fre$LvFpvlue<-1
RBP_fre$LvSpvlue<-1
for (i in 1:152){
  RBP_fre[i,5]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(3,1)],R_sum[c(3,1)]),2,2))$p.value
  RBP_fre[i,6]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,3)],R_sum[c(2,3)]),2,2))$p.value
}
RBP_fre[,1:4]<-data.frame(prop.table(as.matrix(RBP_fre[,1:4]),margin = 2)*100)
RBP_fre$Name<-RBP$RBP.name
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre[RBP_fre$col>1,8]<-4
RBP_fre[RBP_fre$col==1,8]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,8]<-2
RBP_fre[RBP_fre$col==0,8]<-1
Splicesome<-c("SF3B4","PRPF8","EFTUD2","BUD13","AQR","PPIG")
col<-c("black","red","orange","skyblue")

pdf("temp_fig/Fig4BV2.pdf",height=3,width=12)
par(mfrow=c(1,4))
par(pch=16,pty="s")
plot(RBP_fre[,c(3,1)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,3]>=4 | RBP_fre[,1]>=4)),c(3,1)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,3]>=4)])
abline(0,1,col="red")
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(3,1)],col="red",pos = 4,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])

plot(RBP_fre[,c(3,1)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
text(RBP_fre[(RBP_fre$LvFpvlue<=0.05 & (RBP_fre[,3]>=1 | RBP_fre[,1]>=1) ),c(3,1)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$LvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,3]>=1) ]],
     labels = RBP_fre$Name[RBP_fre$LvFpvlue<=0.05& (RBP_fre[,1]>=1 | RBP_fre[,3]>=1) ])
abline(0,1,col="red")

plot(RBP_fre[,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="Other short introns (% RBP sites)",ylab="Long introns (% RBP sites)")
abline(0,1,col="red")
#text(RBP_fre[(RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,1]>=4)),c(4,1)],
#     labels = RBP_fre$Name[RBP_fre$SvFpvlue<=0.05 & (RBP_fre[,1]>=4 | RBP_fre[,4]>=4)])
text(RBP_fre[RBP_fre$Name%in%Splicesome,c(3,2)],col="red",pos = 4,cex=0.5,
     labels = RBP_fre$Name[RBP_fre$Name%in%Splicesome])


plot(RBP_fre[,c(3,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_fre$col],
     xlab="Other short introns (% RBP sites)",ylab="Long introns (% RBP sites)")
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$LvSpvlue<=0.05 & (RBP_fre[,3]>=1 | RBP_fre[,2]>=1)),c(3,2)],cex=0.5,
     col=col[RBP_fre$col[RBP_fre$LvSpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,3]>=1)]],
     labels = RBP_fre$Name[RBP_fre$LvSpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,3]>=1)])
dev.off()
#rm(RBP_fre)


#FigS10
dat<-read.delim("all.FLEXI")
FLEXI_ID<-dat$ID[rowSums(dat[,88:92])>0]
RBP_list<-c("BCLAF1","DX3X","RPS3","ZNF800","XRN2","DKC1","BCLAF")
#pdf("Figures/FigS9.pdf",width=8,height=5)
Splicesome<-c("SF3B4","PRPF8","EFTUD2","BUD13","AQR")
set_1 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[1]])
set_2 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[2]])
set_3 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[3]])
set_4 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[4]])
set_5 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[5]])
set_6 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[6]])
set_7 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_list[7]])
set <- list (set_1,set_2,set_3,set_4,set_5,set_6,set_7)
names(set)<-RBP_list
upset(fromList(set),keep.order=T,set_size.show = T,
      mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",
      sets=names(set),nintersects = 1000,nsets = 4)
dev.off()

#RBP-specific FLEXI (one by one) to all FLEXI
RBP_fre<-RBP[,c(1,46)]
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-RBP_fre[RBP_fre$Cells>0,]
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_list<-sort(unique(RBP_4cell_plasma$RBP))

pdf("temp_fig/RBP_boundFLEXI_RBP.pdf",width=12,height=18)
col<-c("black","red","orange","skyblue")
par(mfrow=c(6,4))
par(pch=16,pty="s")
for (i in 1:126) {
  RBP_name<-RBP_list[i]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name])
  FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Padj_withSign<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:126){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum,2,2)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=1,7]<-"D"
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=1,7]<-"U"
  RBP_plot[RBP_plot$Cells<=RBP_plot$Freq,6]<- -log10(RBP_plot[RBP_plot$Cells<=RBP_plot$Freq,5])
  RBP_plot[RBP_plot$Cells>RBP_plot$Freq,6]<- log10(RBP_plot[RBP_plot$Cells>RBP_plot$Freq,5])
  if (i==1){
    RBP_clus<-RBP_plot[,c(1,2,7)]
    colnames(RBP_clus)[3]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-RBP_plot[,c(1,2,6)]
    colnames(RBP_clus_log10)[3]<-paste0(RBP_name,"BF")
    
  } else {
    RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,2,7)],by=1:2)
    colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,2,6)],by=1:2)
    colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
  }
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="bonferroni")
  axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
  plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,bty="n",col=col[RBP_plot$col],
       ylab=paste0(RBP_name," bound FLEXIs (% RBP sites)"),xlab="ALl FLEXIs (% RBP sites)")
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=0.01 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
       col=col[RBP_plot$col[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
  abline(0,1,col="red")
  legend(axis_max/2,axis_max,cex=0.5,
         text.col=col[RBP_plot$col[sig_cutoff& (RBP_plot$Cells<RBP_plot[,4])]],
         legend = RBP_plot$RBP.name[sig_cutoff & (RBP_plot$Cells<RBP_plot[,4])])
  legend(axis_max-5,axis_max-5,cex=0.5,
         text.col=col[RBP_plot$col[sig_cutoff& (RBP_plot$Cells>RBP_plot[,4])]],
         legend = RBP_plot$RBP.name[sig_cutoff & (RBP_plot$Cells>RBP_plot[,4])])
}
dev.off()

rownames(RBP_clus)<-RBP_clus$RBP.name
RBP_clus<-RBP_clus[,c(3:128)]
RBP_clus<-data.frame(t(RBP_clus))
for (i in 1:126){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}
RBP_col<-RBP_fre[,c(1,2)]
rownames(RBP_col)<-RBP_col$RBP.name

rownames(RBP_clus_log10)<-RBP_clus_log10$RBP.name
RBP_clus_col<-RBP_clus_log10[,1:2]
RBP_clus_col$BP.name<-paste0(RBP_clus_col$RBP.name,"BF")
RBP_clus_log10<-RBP_clus_log10[,3:128]
RBP_clus_log10.t<-data.frame(t(RBP_clus_log10))

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")

#gower.dist.r <- daisy(RBP_clus_log10.t, metric = c("euclidean"))
gower.dist.c <- daisy(RBP_clus_log10, metric = c("euclidean"))
#aggl.clust.r <- hclust(gower.dist.r, method = "complete")
aggl.clust.c <- hclust(gower.dist.c, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10.t
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (i in 1:126){
  colname<-paste0(colnames(RBP_clus_log10_convertedvalue)[i],"BF")
  RBP_clus_log10_convertedvalue[colname,i]<-NA
}

pdf("temp_fig/heatmap_RBP.pdf",height=10,width=10)
labels.c<-colnames(RBP_clus_log10_convertedvalue)
col.c<-RBP_clus_col$col[RBP_clus_col$RBP.name%in%labels.c]
labels.r<-labels(aggl.clust.r)
labels.r<-substr(labels.r,1,nchar(labels.r)-2)
col.r<-RBP_clus_col$col[RBP_clus_col$RBP.name%in%labels.r]
Heatmap(RBP_clus_log10_convertedvalue, 
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        na_col="black",
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        column_names_gp = gpar(fontsize = 7,col=col[col.c]),
        cluster_rows =aggl.clust.r )
dev.off()

#53 RBP
RBP_53<-paste0(Fun$Name[Fun$RBP_by_FLEXI>=30],"BF")
RBP_clus_log10.53<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_53]
RBP_clus_log10.53.t<-data.frame(t(RBP_clus_log10.53))
RBP_clus.53<-RBP_clus[rownames(RBP_clus)%in%RBP_53,]
gower.dist <- daisy(RBP_clus.53, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")

#gower.dist.r <- daisy(RBP_clus_log10.53.t, metric = c("euclidean"))
gower.dist.c <- daisy(RBP_clus_log10.53, metric = c("euclidean"))
#aggl.clust.r <- hclust(gower.dist.r, method = "complete")
aggl.clust.c <- hclust(gower.dist.c, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10.53.t
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (i in 1:53){
  colname<-rownames(RBP_clus_log10_convertedvalue)[i]
  colname<-substr(colname,1,nchar(colname)-2)
  RBP_clus_log10_convertedvalue[i,colnames(RBP_clus_log10_convertedvalue)==colname]<-NA
}

pdf("temp_fig/53RBP_heatmap_cluster.pdf",height=6,width=12)
labels.c<-colnames(RBP_clus_log10_convertedvalue)
col.c<-RBP_clus_col$col[RBP_clus_col$RBP.name%in%labels.c]
labels.r<-labels(aggl.clust.r)
labels.r<-substr(labels.r,1,nchar(labels.r)-2)
col.r<-RBP_clus_col$col[RBP_clus_col$RBP.name%in%labels.r]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        column_names_gp = gpar(fontsize = 7,col=col[col.c]),
        cluster_rows = color_branches(aggl.clust.r, k =7),
        #cluster_columns = aggl.clust.c
        )
dev.off()

#53 RBP -no core splice RBP
# make RBP_clus and RBP_clus_log10 data frame
# table is made in RBP_clus_byCEllType.R

RBP_fre<-RBP[,c(1,46)]
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-RBP_fre[RBP_fre$Cells>0,]
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]

dat<-read.delim("4cell_plasma_FLEXI.tsv")
cutoff<-0
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>cutoff,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_list<-sort(unique(RBP_4cell_plasma$RBP))

RBP_clus<-read.delim("RBP_clus_froGower.txt")
marker_text<-RBP_clus
marker_text[marker_text=="U"]<-"X"
marker_text[marker_text=="N"]<-""
marker_text[marker_text=="D"]<-"X"
RBP_clus<-data.frame(t(RBP_clus))
for (i in 1:dim(RBP_clus)[2]){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}

RBP_clus_log10<-read.delim("RBP_clus_froHeatmap.txt")

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}
RBP_col<-RBP[,c(1,46)]
RBP_col$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_col[RBP_col$col>1,3]<-4
RBP_col[RBP_col$col==1,3]<-3
RBP_col[RBP_col$col<1 & RBP_col$col>0,3]<-2
RBP_col[RBP_col$col==0,3]<-1
RBP_col<-RBP_col[,c(1,3)]
rownames(RBP_col)<-RBP_col$RBP.name

col<-c("black","red","orange","skyblue")
pdf("temp_fig/126RBP_heatmap_cluster3cutoff.pdf",height=8,width=8)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 4,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        show_column_names = F,
        cluster_rows = F, cluster_columns =aggl.clust.r)
dev.off()

RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_47<-data.frame(table(RBP_4cell_plasma$RBP))
RBP_47<-RBP_47[order(RBP_47$Freq,decreasing = T),]
RBP_47_name<-as.character(RBP_47$Var1[RBP_47$Freq>=30])
RBP_47_name<-RBP_47_name[7:length(RBP_47_name)]
RBP_47<-paste0(RBP_47_name,"BF")

RBP_clus_log10.47<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
RBP_clus_log10.47<-RBP_clus_log10.47[rownames(RBP_clus_log10.47)%in%RBP_47_name,]
marker_text.47<-marker_text[,colnames(marker_text)%in%RBP_47]
marker_text.47<-marker_text.47[rownames(marker_text.47)%in%RBP_47_name,]
RBP_clus.47<-RBP_clus[rownames(RBP_clus)%in%RBP_47,]
RBP_clus.47<-RBP_clus.47[,colnames(RBP_clus.47)%in%RBP_47_name]

gower.dist <- daisy(RBP_clus.47, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10.47
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

pdf("temp_fig/47RBP_heatmap_cluster3cutoff.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
marker_text.47<-marker_text.47[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        show_column_names = F,
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        cluster_rows = F, cluster_columns =aggl.clust.r,
        cell_fun = function(i, j, x, y,w, h, col) {
          grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))})
dev.off()

labels.r<-labels(aggl.clust.r)
labels.r<-substr(labels.r,1,nchar(labels.r)-2)
RBP_47_info<-RBP[match(RBP_fre$RBP.name[order(RBP_fre$Cells,decreasing = T)][1:53],RBP$RBP.name),1:43]

RBP_47_info$StressGranule<-0
SG<-read.table("SG_MSGP_RGD.list",col.names="ID")
RBP_47_info[RBP_47_info$RBP.name%in%SG$ID,44]<-1
RBP_47_info<-RBP_47_info[,c(1:31,44,32:43)]
RBP_47_info[is.na(RBP_47_info)]<-0
RBP_47_info<-RBP_47_info[,c(TRUE,TRUE,colSums(RBP_47_info[,3:44])>0)]
#AGO is in P body, stress granule, nuclei and cytoplasm
RBP_47_info[RBP_47_info$RBP.name=="AGO",c(19,23,26,27)]<-1
#DICER is in P body, ER, nuclei and cytoplasm
RBP_47_info[RBP_47_info$RBP.name=="DICER",c(19,23,26,28)]<-1
#BCLAF1 has spicing regulation function
RBP_47_info[RBP_47_info$RBP.name=="BCLAF1",4]<-1
#additonl infor bout binding site location
RBP_47_info$I<-0 
RBP_47_info$SS_I<-0 
RBP_47_info$SS_E<-0 
name_list<-c("AATF", "DICER", "DKC1" , "EFTUD2", "KHSRP", "NOLC1", "PCBP1", 
             "PCBP2", "PRPF4", "RBFOX2", "RBM5", "RBM15", "RBM22", "SMNDC1", "TIA1",  "XRN2")
RBP_47_info$I[RBP_47_info$RBP.name%in%name_list]<-1
name_list<-c("PRPF8", "SF3B4", "AQR" , "BUD13", "AGO", "SF3A3", "U2AF1", 
             "GPKOW", "LSM11", "TIAL1", "GEMIN5", "NCBP2", "IGF2BP1")
RBP_47_info$SS_I[RBP_47_info$RBP.name%in%name_list]<-1
RBP_47_info$SS_E<-1-(RBP_47_info$I+RBP_47_info$SS_I)

#additional infor on mRNA  level and alternative expression
RBP_47_info$mRNA<-0 
RBP_47_info$AltSPl_RI<-0
RBP_47_info$AltSPl_SE<-0
RBP_47_info$LowGC<-0
RBP_47_info$phastCons<-0
name_list<-c("DDX55", "GPKOW", "LARP4", "NCBP2", "PCBP2", "RBM15", "RBM22", "RPS3", 
             "SF3A3", "SRSF1", "TIAL1", "U2AF2", "UCHL5", "ZNF622","SF3B4","AQR","U2AF1")
RBP_47_info$mRNA[RBP_47_info$RBP.name%in%name_list]<-1

name_list<-c("AQR", "RBM22", "SF3B4",
             "AATF", "BUD13", "DKC1", "EFTUD2", 
            "METAP2", "NOLC1", "PABPC4", "TIA1")
RBP_47_info$AltSPl_RI[RBP_47_info$RBP.name%in%name_list]<-1

name_list<-c("AATF","AQR","EFTUD2", "GEMIN5", "GRWD1", "IGF2BP1","LARP4","LIN28B", "LSM11", "METAP2",
             "NCBP2","NOLC1","PABPC4", "PRPF8", "RBM15", "RBM22", "RPS3", "SF3B4", "SUB1", "TIA1",
             "TIAL1", "U2AF1","UCHL5", "XRN2","YBX3")
RBP_47_info$AltSPl_SE[RBP_47_info$RBP.name%in%name_list]<-1


name_list<-c("BCLAF1","TRA2A","ZNF800","TIA1","TIAL1","KHSRP","U2AF2","U2AF1","GPKOW",
             "NOLC1","AATF","DKC1")
RBP_47_info$LowGC[RBP_47_info$RBP.name%in%name_list]<-1
name_list<-c("ZNF622","BCLAF1","U2AF2","U2AF1","TIAL1")
RBP_47_info$phastCons[RBP_47_info$RBP.name%in%name_list]<-1

write.table(RBP_47_info,"53_RBP_info_fig4_7_S16.txt",sep="\t",quote=F,row.names=F)
### need modify
library(readxl)
RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 1))
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:28]
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
temp<-RBP_47_info[match(RBP_clus_image$RBP.name,RBP_47_info$RBP.name),]
RBP_clus_image<-cbind(RBP_clus_image,temp)


icol<-c(brewer.pal(8,"Set1"),brewer.pal(12,"Paired"))
pdf("temp_fig/RBP_clus_image.pdf",width=8,height=8)
par(mfrow=c(4,1),mar = c(5,2,2,20))
image(1:27,1:4,as.matrix(RBP_clus_image[,4:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:4,labels = colnames(RBP_clus_image)[4:1],tick = FALSE,cex=0.5)

image(1:27,1:8,as.matrix(RBP_clus_image[,51:44]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:8,labels = colnames(RBP_clus_image)[51:44],tick = FALSE)

image(1:27,1:16,as.matrix(RBP_clus_image[,23:8]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP_clus_image)[23:8],tick = FALSE)

image(1:27,1:14,as.matrix(RBP_clus_image[,37:24]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP_clus_image)[37:24],tick = FALSE)
axis(1,las=2,at = 1:27,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)
dev.off()

#sheet 2
RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 2))
colnames(RBP_clus_image)[c(4,9,30,35)]<-c("SUB1","METAP2","SUB1-1","METAP2-1")
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:40]
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
temp<-RBP_47_info[match(RBP_clus_image$RBP.name,RBP_47_info$RBP.name),]
temp[29,]<-temp[3,]
temp[34,]<-temp[8,]
RBP_clus_image<-cbind(RBP_clus_image,temp)
icol<-c("white",brewer.pal(8,"Set1"),brewer.pal(12,"Paired"))

pdf("temp_fig/RBP_clus_image_p005.pdf",width=8,height=8)
par(mfrow=c(4,1),mar = c(5,2,2,20))
image(1:39,1:24,as.matrix(RBP_clus_image[,24:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:24,labels = colnames(RBP_clus_image)[24:1],tick = FALSE,cex=0.5)

image(1:39,1:4,as.matrix(RBP_clus_image[,67:64]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:4,labels = colnames(RBP_clus_image)[67:64],tick = FALSE)

image(1:39,1:16,as.matrix(RBP_clus_image[,43:28]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP_clus_image)[43:28],tick = FALSE)

image(1:39,1:14,as.matrix(RBP_clus_image[,57:44]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP_clus_image)[57:44],tick = FALSE)
axis(1,las=2,at = 1:39,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)
dev.off()

'''
#sheet 5
RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 5))
colnames(RBP_clus_image)[c(4,9,20,30,35,41)]<-c("SUB1","METAP2","FXR2","SUB1-1","METAP2-1","FXR2-1")
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:42]
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
temp<-RBP_47_info[match(RBP_clus_image$RBP.name,RBP_47_info$RBP.name),]
temp[29,]<-temp[3,]
temp[34,]<-temp[8,]
temp[40,]<-temp[19,]
RBP_clus_image<-cbind(RBP_clus_image,temp)
icol<-c("white",brewer.pal(8,"Set1"),brewer.pal(12,"Paired"))

pdf("temp_fig/RBP_clus_image_exon.pdf",width=8,height=8)
par(mfrow=c(4,1),mar = c(5,2,2,20))
image(1:41,1:34,as.matrix(RBP_clus_image[,34:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:34,labels = colnames(RBP_clus_image)[34:1],tick = FALSE,cex=0.5)

image(1:41,1:4,as.matrix(RBP_clus_image[,77:74]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:4,labels = colnames(RBP_clus_image)[77:74],tick = FALSE)

image(1:41,1:16,as.matrix(RBP_clus_image[,53:38]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP_clus_image)[53:38],tick = FALSE)

image(1:41,1:14,as.matrix(RBP_clus_image[,67:54]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP_clus_image)[67:54],tick = FALSE)
axis(1,las=2,at = 1:41,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)
dev.off()
'''
# new fig4A
RBP53<-read.delim("53_RBP_info_fig4_7_S16.txt")

pdf("temp_fig/RBP53_image.pdf",width=8,height=8)
par(mfrow=c(3,1),mar = c(5,2,2,20))
image(1:53,1:8,as.matrix(RBP53[,46:39]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:8,labels = colnames(RBP53)[46:39],tick = FALSE)

image(1:53,1:16,as.matrix(RBP53[,18:3]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP53)[18:3],tick = FALSE)

image(1:53,1:14,as.matrix(RBP53[,32:19]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP53)[32:19],tick = FALSE)
axis(1,las=2,at = 1:53,labels = RBP53$RBP.name,tick = FALSE,cex=0.5)
dev.off()

#53 RBP-by group
group_list<-list(group1=c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"))
group_list$group2<-c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24")
group_list$group3<-c("TIA1","TIAL1")
group_list$group4<-c("U2AF1","U2AF2","KHSRP")
group_list$group5<-c("AATF","NOLC1","DKC1","SMNDC1")
group_list$group6<-c("AGO","DICER")


RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,50)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]

percent_cutoff<-2
pvalue_cutoff<-0.05

pdf("temp_fig/RBP_boundFLEXI_RBP_by_group.pdf",width=12,height=9)
col<-c("black","red","orange","skyblue")
par(mfrow=c(4,3),mar=c(3,3,1,1))
par(pch=16,pty="s")
for (i in 1:6) {
  RBP_name<-group_list[[i]]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_name])
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
  plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,
       bty="n",col=col[RBP_plot$col],
       ylab=paste0("Group ",i,"RBP bound FLEXIs (% RBP sites)"),xlab="ALl FLEXIs (% RBP sites)")
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff & (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
       col=col[RBP_plot$col[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
  }
}
dev.off()
#Fig1G
pdf("temp_fig/Fig1G.pdf",height=8,width=8)

group_CPM<-dat_CPM[,c(1,33:36)]
group_FLEXI<-list(goup1=unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[1]]]))
group_FLEXI$group2<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[2]]])
group_FLEXI$group3<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[3]]])
group_FLEXI$group4<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[4]]])
group_FLEXI$group5<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[5]]])
group_FLEXI$group6<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[6]]])
group_FLEXI$group7<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[7]]])
group_FLEXI$group8<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[8]]])
group_FLEXI$group9<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[9]]])
group_FLEXI$group10<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[10]]])
group_FLEXI$group11<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[11]]])
par(mfrow=c(3,2))
for (i in 1:11){
  plot(density(log10(rowMeans(group_CPM[group_CPM$ID%in%group_FLEXI[[i]],2:5]))),col="red",
       bty="n",xlim=c(-4,1),ylim=c(0,2),xlab="RPM",main=paste0("Complex ",i))
  lines(density(log10(rowMeans(group_CPM[!group_CPM$ID%in%group_FLEXI[[i]],2:5]))))
}
dev.off()

#individual 47 RBP scatter
dat<-read.delim("4cell_plasma_FLEXI.tsv")
#this can be chanhged
cutoff<-0
percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")
#simple cutoff
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>cutoff,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)

RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,50)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]


pdf("temp_fig/47RBP_boundFLEXI_RBP.pdf",width=9,height=15)
col<-c("black","red","orange","skyblue")
par(mfrow=c(5,3))
par(pch=16,pty="s")
for (i in 1:length(labels.r)) {
    Name<-labels.r[i]
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==Name])
    FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    R_sum<-colSums(RBP_plot[,3:4])
    for (j in 1:dim(RBP_plot)[1]){
      RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
    }
    RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
    plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,bty="n",col=col[RBP_plot$col],
         ylab=paste0(Name," bound FLEXIs (% RBP sites)"),xlab="ALl FLEXIs (% RBP sites)")
    #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
    sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff & (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
    abline(0,1,col="red")
    text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
         col=col[RBP_plot$col[sig_cutoff]],
         labels = RBP_plot$RBP.name[sig_cutoff])
}
dev.off()

#GO enrich via ShinyGO
library(readxl)
GO<-read_xlsx("Gene_by_clusterofRBPs.xlsx",sheet = 2)
colnames(GO)[1]<-"C1"
GO<-GO[GO$C1<0.001,c(4,1)]
for (i in 1:5){
  temp<-read_xlsx("Gene_by_clusterofRBPs.xlsx",sheet = i+2)
  colnames(temp)[1]<-paste0("C",i+1)
  temp<-temp[,c(4,1)]
  temp<-temp[temp[,2]<0.001,]
  GO<-merge(GO,temp,by=1,all=T)
}
GO[is.na(GO)]<-0.001
GO[,2:7]<- -log10(GO[,2:7])
rownames(GO)<-GO$`Functional Category`
pdf("temp_fig/GO_RBP_by_cluster.pdf")
Heatmap(GO[,2:7], 
        row_names_gp = gpar(fontsize = 4),
        col = colorRamp2(seq(3,9,length = 2),c("white", "red")),
        column_names_gp = gpar(fontsize = 4),
        cluster_rows = T, cluster_columns =F)
dev.off()
## Hallmark
GO<-read_xlsx("hallmark_RBP_by_cluster.xlsx",sheet = 1)
colnames(GO)[1]<-"C1"
GO<-GO[,c(4,1)]
for (i in 1:5){
  temp<-read_xlsx("hallmark_RBP_by_cluster.xlsx",sheet = i+1)
  colnames(temp)[1]<-paste0("C",i+1)
  temp<-temp[,c(4,1)]
  GO<-merge(GO,temp,by=1,all=T)
}
GO[is.na(GO)]<-1
rownames(GO)<-GO$`Functional Category`
GO[,2:7]<- -log10(GO[,2:7])
pdf("temp_fig/Hallmark_RBP_by_cluster.pdf")
Heatmap(GO[,2:7], 
        row_names_gp = gpar(fontsize = 4),
        col = colorRamp2(seq(0,9,length = 2),c("white", "red")),
        column_names_gp = gpar(fontsize = 4),
        cluster_rows = T, cluster_columns =F)
dev.off()

#snoRNA RBP scatter

pdf("temp_fig/snoRNAFLEXI_RBP.pdf",width=10,height=5)
col<-c("black","red","orange","skyblue")
par(mfrow=c(1,2),pch=16,pty="s")
snoRNA_FLEXI<-unique(dat$ID[dat$Has_snoRNA!="."])
FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%snoRNA_FLEXI]))
RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
RBP_plot[is.na(RBP_plot)]<-0
RBP_plot$Padj<-1
R_sum<-colSums(RBP_plot[,3:4])
for (j in 1:126){
  RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum,2,2)))$p.value
}
RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,bty="n",col=col[RBP_plot$col],
         ylab="snoRNA FLEXIs (% RBP sites)",xlab="ALl FLEXIs (% RBP sites)")
sig_cutoff<-(RBP_plot$Padj<=0.01 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
abline(0,1,col="red")
text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
    col=col[RBP_plot$col[sig_cutoff]],
    labels = RBP_plot$RBP.name[sig_cutoff])

plot(RBP_plot[,c(3,4)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",col=col[RBP_plot$col],
     ylab="snoRNA FLEXIs (% RBP sites)",xlab="ALl FLEXIs (% RBP sites)")
abline(0,1,col="red")
sig_cutoff<-(RBP_plot$Padj<=0.01 & (RBP_plot[,3]>=1 | RBP_plot[,4]>=1))

text(RBP_plot[sig_cutoff,3:4],cex=0.5,pos = 4,
     col=col[RBP_plot$col[sig_cutoff]],
     labels = RBP_plot$RBP.name[sig_cutoff])
dev.off()


####other chunk
RBP_53<-Fun$Name[Fun$RBP_by_FLEXI>=30]
RBP_53_col<-RBP_col[RBP_53,]
RBP_53<-RBP_clus[rownames(RBP_clus)%in%RBP_53,]

RBP_23<-RBP_53_col$RBP.name[RBP_53_col$col==1]
RBP_23<-RBP_clus[RBP_23,]

gower.dist <- daisy(RBP_23, metric = c("gower"))
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = TRUE, keep.diss = TRUE)
#check silhouette plot
aggl.clust.c <- hclust(gower.dist, method = "complete")
ggplot(data = data.frame(t(cstats.table(gower.dist, divisive.clust, 15))), 
       aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("Divisive clustering") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(data = data.frame(t(cstats.table(gower.dist, divisive.clust, 15))), 
       aes(x=cluster.number, y=avg.silwidth)) + 
  geom_point()+
  geom_line()+
  ggtitle("Divisive clustering") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))
ggplot(data = data.frame(t(cstats.table(gower.dist, aggl.clust.c, 15))), 
       aes(x=cluster.number, y=avg.silwidth)) + 
  geom_point()+
  geom_line()+
  ggtitle("Agglomerative clustering") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))

tsne<-Rtsne(gower.dist, dims = 2,
                 perplexity=3, theta = 0.5,is_distance = TRUE)

plot(tsne$Y,  main="t-SNE",
     xlab="t-SNE1",ylab="t-SNE2",cex=1.5)

dendro <- as.dendrogram(aggl.clust.c)
labels<-labels(dendro.col)
labels_col<-rep(0,23)
labels_hei<-rep(0,23)
labels_fre<-rep(0,23)
for (i in 1:23){
  labels_hei[i]<-length(unique(RBP_info_4cell$ID[RBP_info_4cell$RBP==labels[i]]))
  labels_fre[i]<-RBP_fre$Cells[RBP_fre$RBP.name==labels[i]]
  labels_col[i]<-col[RBP_col$col[rownames(RBP_col)==labels[i]]]
}
labels_hei<-labels_hei/max(labels_hei)
labels_fre<-labels_fre/max(labels_fre)
dendro.col <- dendro %>%
  set("branches_k_color", 
      k = 3)%>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", value = labels_col) %>% 
  set("labels_cex", 0.5)
pdf("temp_fig/Polar_dendrogram_plot53.pdf")
circlize_dendrogram(dendro.col)
dev.off()

pdf("temp_fig/verticalden_53.pdf",width=11,height=6)
ggplot(dendro.col, horiz = TRUE, theme = NULL)
RBP_53_info<-RBP
rownames(RBP_53_info)<-RBP_53_info$RBP.name
RBP_53_info<-RBP_53_info[labels,4:22]
par(mar = c(5,2,2,10))
image(1:53,1:19,as.matrix(RBP_53_info),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:19,labels = colnames(RBP_53_info),cex=0.5,tick = FALSE)
axis(1,las=2,at = 1:53,labels = rownames(RBP_53_info),cex=0.5,tick = FALSE)
dev.off()

pdf("temp_fig/verticalden_23.pdf",width=11,height=6)
ggplot(dendro.col, horiz = TRUE, theme = NULL)
RBP_23_info<-RBP
rownames(RBP_23_info)<-RBP_23_info$RBP.name
RBP_23_info<-RBP_23_info[labels,4:22]
par(mar = c(5,2,2,10))
image(1:23,1:19,as.matrix(RBP_23_info),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:19,labels = colnames(RBP_23_info),cex=0.5,tick = FALSE)
axis(1,las=2,at = 1:23,labels = rownames(RBP_23_info),cex=0.5,tick = FALSE)
dev.off()


gower.dist.r <- daisy(RBP_53, metric = c("gower"))
RBP_53.c<-data.frame(t(RBP_53))
for (i in 1:53){
  RBP_53.c[,i]<-factor(RBP_53.c[,i],levels = c("D","N","U"))
}
gower.dist.c <- daisy(RBP_53.c, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist.r, method = "complete")
aggl.clust.c <- hclust(gower.dist.c, method = "complete")
pdf("temp_fig/RBP_cluster.pdf",height=6,width=12)
Heatmap(RBP_53, 
        row_names_gp = gpar(fontsize = 7),col=c("royalblue","gray75","red"),
        column_names_gp = gpar(fontsize = 7),
        cluster_rows = color_branches(aggl.clust.r, k =3),
        cluster_columns = color_branches(aggl.clust.c, k = 3))
dev.off()


#cluster by RBP and FLEXIs contain them, work not well,
for (i in 1:126) {
  RBP_name<-RBP_list[i]
  if (i==1){
    RBP_by_FLEXI<-data.frame("ID"=unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name]))
    RBP_by_FLEXI$Counts<-1
    colnames(RBP_by_FLEXI)[2]<-RBP_name
  } else {
    temp_table<-data.frame("ID"=unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name]))
    temp_table$Counts<-1
    colnames(temp_table)[2]<-RBP_name
    RBP_by_FLEXI<-merge(RBP_by_FLEXI,temp_table,by=1,all=T)
  }
}
RBP_by_FLEXI[is.na(RBP_by_FLEXI)]<-0
rownames(RBP_by_FLEXI)<-RBP_by_FLEXI$ID
RBP_by_FLEXI.t<-data.frame(t(RBP_by_FLEXI[,2:127]))
RBP_by_FLEXI<-RBP_by_FLEXI[,2:127]

for (i in 1:126){
  RBP_by_FLEXI[,i]<-factor(RBP_by_FLEXI[,i],levels = c(0,1))
}
for (i in 1:4505){
  RBP_by_FLEXI.t[,i]<-factor(RBP_by_FLEXI.t[,i],levels = c(0,1))
}
gower.dist.r <- daisy(RBP_by_FLEXI.t, metric = c("gower"))
gower.dist.c <- daisy(RBP_by_FLEXI, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist.r, method = "complete")
aggl.clust.c <- hclust(gower.dist.c, method = "complete")
pdf("temp_fig/RBP_cluster.pdf",height=6,width=12)
Heatmap(RBP_by_FLEXI.t, 
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = color_branches(aggl.clust.r, k =3),
        cluster_columns = color_branches(aggl.clust.c, k = 3))
dev.off()


## density for long FLEXIs

##Fig1D, RPM density of FLEXIs, different group
longFLEXI<-list(longHEK=read.table("UHRR.longFLEXI.counts",col.names=c("ID","counts")))
longFLEXI[[2]]<-read.table("K562.longFLEXI.counts",col.names=c("ID","counts"))
longFLEXI[[3]]<-read.table("HEK.longFLEXI.counts",col.names=c("ID","counts"))
longFLEXI[[4]]<-read.table("Hela.longFLEXI.counts",col.names=c("ID","counts"))
dat<-read.delim("all.FLEXI")
pdf("temp_fig/longFLEXIDensity.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(2,2,2,2,2)
for (i in c(88:91)){
  plot(density(log10(longFLEXI[[i-87]]$counts/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,2),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="red",axes=F)
  lines(density(log10(dat[dat[,i]>0,i]/mapped_reads[i-81])),col="black")
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-2),1,bquote(10^2)),
       at=seq(-4,2,2))
  legend("topleft",lty=1,
         bty="n",col=c("black","red"),legend = c("FLEXIs","Long FLEXIs"))
}
dev.off()

pdf("temp_fig/longFLEXIDensity1.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(200,200,200,50)
for (i in c(88:91)){
  plot(density((longFLEXI[[i-87]]$counts/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(0,0.1),ylim=c(0,D_height[i-87]),main=colnames(dat)[i],col="red",axes=F)
  lines(density((dat[dat[,i]>0,i]/mapped_reads[i-81])),col="black")
  axis(2,labels=seq(0,200,50),las=1,at=seq(0,200,50),las=2)
  axis(1,labels=c(0,0.01,0.05,0.1),
       at=c(0,0.01,0.05,0.1))
  legend("topright",lty=1,
         bty="n",col=c("black","red"),legend = c("FLEXIs","Long FLEXIs"))
}
dev.off()
# len density by all FLEXIs all length

longFLEXI<-list(K562=read.table("K562.longFLEXI.len",col.names=c("Len")))
longFLEXI[[2]]<-read.table("HEK.longFLEXI.len",col.names=c("Len"))
longFLEXI[[3]]<-read.table("Hela.longFLEXI.len",col.names=c("Len"))
longFLEXI[[4]]<-read.table("UHRR.longFLEXI.len",col.names=c("Len"))
dat<-read.delim("all.FLEXI")
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

pdf("temp_fig/longFLEXIDensity_len.pdf",onefile = T,width=8,height=8)
par(mfrow=c(2,2),lwd=1.5)
D_height<-c(2,2,2,2,2)
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
    legend(4,3,pch=16,legend = c("FLEXIs (≤300 nt)","Long FLEXIs (>300 nt)"),col=c("red","blue"),bty="n")
  }
}
dev.off()


##new Fig5
# make new table of FLEXI only in 4 cells
FourCellFLEXI<-dat$ID[rowSums(dat[,88:91])>0]
FourCell<-dat[dat$ID%in%FourCellFLEXI,]

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

RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FourCellFLEXI,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,50)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]

percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")
pdf("temp_fig/Fig5-1.pdf",width=8,height=20,onefile = T)
par(mfrow=c(10,4),pch=16,pty="s",mar=c(2,2,2,2))
for (i in 1:10){
  if (i==3){
    FLEXI_dat<-FourCell[FourCell$Has_snoRNA!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FourCell[FourCell$Has_snoRNA==".",]
  } else if (i==4){
    FLEXI_dat<-FourCell[FourCell$PhastCon30>=0.75,]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FourCell[FourCell$PhastCon30<0.75,]
  } else if (i==1){
    FLEXI_dat<-FourCell[FourCell$Is_agotron!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FourCell[FourCell$Is_agotron==".",]
  } else if (i==2){
    FLEXI_dat<-FourCell[FourCell$Is_mirtron!=".",]
    FLEXI_list<-unique(FLEXI_dat$ID)
    temp2<-FourCell[FourCell$Is_mirtron==".",]
  } else {
    RBP_list_sig<-group_list[[i]]
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig])
    FLEXI_dat<-FourCell[FourCell$ID%in%FLEXI_list,]
    temp2<-FourCell[!FourCell$ID%in%FLEXI_list,]
  }
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
    test_set<-FourCell[sample(1:dim(FourCell)[1],sample_size,replace = F),]
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

#RPM of each clusters
#1. by cell type
pdf("temp_fig/CLuster_RPM1.pdf",width=8,height=8)
par(mfrow=c(2,2),pch=16,pty="s",mar=c(2,2,2,2))
col=brewer.pal(6,"Paired")
y_height<-c(1.5,1.5,1.5,1.5)
for (i in 1:4){
  map_cell<-mapped_reads[i+6]
  FLEXI_dat<-FourCell[FourCell[87+i]>0,]
  for (j in 1:6) {
    RBP_list_sig<-group_list[[j+4]]
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig])
    if(j==1) {
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      temp2<-FLEXI_dat[!FLEXI_dat$ID%in%FLEXI_list,]
      plot(density(log10(temp1[,i+87]/map_cell)),bty="n",xlab="RPM",ylim=c(0,y_height[i]),
           xlim=c(-4,1),main=colnames(FourCell)[87+i],col=col[j],axes=T)
    } else if (j==6) {
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      temp2<-temp2[!temp2$ID%in%FLEXI_list,]
      lines(density(log10(temp1[,i+87]/map_cell)),col=col[j])
      lines(density(log10(temp2[,i+87]/map_cell)),col="black")
    } else {
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      temp2<-temp2[!temp2$ID%in%FLEXI_list,]
      lines(density(log10(temp1[,i+87]/map_cell)),col=col[j])
    }
    if (i==1){
      legend(0,0.8,lty=1,col=c(col,"black"),legend = c(paste0("Cluster ",1:6),"Other"),bty="n")
    }
  }
}
dev.off()
  
pdf("temp_fig/CLuster_RPM2.pdf",width=8,height=12)
par(mfrow=c(2,3),pch=16,pty="s",mar=c(2,2,2,2))
col=brewer.pal(6,"Dark2")
y_height<-c(1.5,1.5,1.5,1.5,1.5,1.5)
for (i in 1:6){
  RBP_list_sig<-group_list[[i+4]]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig])
  for (j in 1:4){
    map_cell<-mapped_reads[j+6]
    FLEXI_dat<-FourCell[FourCell[87+j]>0,]
    if (j==1) {
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      plot(density(log10(temp1[,j+87]/map_cell)),bty="n",xlab="RPM",ylim=c(0,y_height[i]),
           xlim=c(-4,2),main=paste0("Cluster ",i),col=col[j],axes=T)
    } else {
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      lines(density(log10(temp1[,j+87]/map_cell)),col=col[j])
    }
  }
  if (i==1){
      legend(0,1.2,lty=1,col=c(col[1:4]),legend = names(FLEXI_dat)[88:91],bty="n")
    }
}
dev.off()

#individual RBP in clusters I/II
pdf("temp_fig/cluster1.pdf",width=6,height=12,onefile = T)
par(mfrow=c(4,2),pch=16,pty="s",mar=c(2,2,2,2))
col=brewer.pal(6,"Paired")
y_height<-1.5
RBP_list_sig<-sort(group_list[[5]])
for (i in 1:4){
  map_cell<-mapped_reads[i+6]
  FLEXI_dat<-FourCell[FourCell[87+i]>0,]
  for (j in 1:8) {
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig[j]])
    temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
    temp2<-FLEXI_dat[!FLEXI_dat$ID%in%FLEXI_list,]
    plot(density(log10(temp1[,i+87]/map_cell)),bty="n",xlab="RPM",ylim=c(0,2),
         xlim=c(-4,1),main=paste(colnames(FourCell)[87+i],"-",RBP_list_sig[j]),col="red",axes=T)
    lines(density(log10(temp2[,i+87]/map_cell)),col="black")
    legend(-2,1.8,lty=1,col=c("red","black"),
           legend = c(paste0(RBP_list_sig[j]," (",dim(temp1)[1],")")
                      ,paste0("Other (",dim(temp2)[1],")")),
           bty="n")
  }
}
dev.off()

pdf("temp_fig/cluster2.pdf",width=6,height=12,onefile = T)
par(mfrow=c(4,2),pch=16,pty="s",mar=c(2,2,2,2))
col=brewer.pal(6,"Paired")
y_height<-1.5
RBP_list_sig<-sort(group_list[[6]])
for (i in 1:4){
  map_cell<-mapped_reads[i+6]
  FLEXI_dat<-FourCell[FourCell[87+i]>0,]
  for (j in 1:8) {
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig[j]])
    temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
    temp2<-FLEXI_dat[!FLEXI_dat$ID%in%FLEXI_list,]
    plot(density(log10(temp1[,i+87]/map_cell)),bty="n",xlab="RPM",ylim=c(0,2),
         xlim=c(-4,1),main=paste(colnames(FourCell)[87+i],"-",RBP_list_sig[j]),col="red",axes=T)
    lines(density(log10(temp2[,i+87]/map_cell)),col="black")
    legend(-2,1.8,lty=1,col=c("red","black"),
           legend = c(paste0(RBP_list_sig[j]," (",dim(temp1)[1],")")
                      ,paste0("Other (",dim(temp2)[1],")")),
           bty="n")
  }
}
dev.off()

##various density plot related to RBP
pdf("temp_fig/RBP_group_density.pdf",width=6,height=6,onefile = T)
par(mfrow=c(2,2),pch=16,pty="s",mar=c(2,2,2,2))
col=brewer.pal(6,"Paired")
group_name<-c("FLEXIs w/o RBP sites","FLEXIs w/RBP sites (without core RBP)","FLEXIs w/RBP sites (47 RBPs)","All FLEXIs")
tmp<-data.frame(table(RBP_4cell_plasma$RBP))
tmp<-tmp[order(tmp$Freq,decreasing = T),]
tmp<-tmp[7:126,]
RBP_list1<-unique(as.character(tmp$Var1))
RBP_list2<-unique(as.character(tmp$Var1[tmp$Freq>=30]))
for (i in 1:4){
  map_cell<-mapped_reads[i+6]
  FLEXI_dat<-FourCell[FourCell[87+i]>0,]
  for (j in 1:4) {
    if (j==1){
      FLEXI_list<-unique(RBP_4cell_plasma$ID)
      temp1<-FLEXI_dat[!FLEXI_dat$ID%in%FLEXI_list,]
      n1<-dim(temp1)[1]
      plot(density(log10(temp1[,i+87]/map_cell)),bty="n",xlab="RPM",ylim=c(0,2),
           xlim=c(-4,1),main=paste(colnames(FourCell)[87+i]),col=col[1],axes=T)
    } else if (j==2){
      FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list1])
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      n2<-dim(temp1)[1]
      lines(density(log10(temp1[,i+87]/map_cell)),col=col[4])
    } else if (j==3){
      FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list2])
      temp1<-FLEXI_dat[FLEXI_dat$ID%in%FLEXI_list,]
      n3<-dim(temp1)[1]
      lines(density(log10(temp1[,i+87]/map_cell)),col=col[5])
    } else {
      n4<-dim(FLEXI_dat)[1]
      lines(density(log10(FLEXI_dat[,i+87]/map_cell)),col="black")
    }
  }
  legend(-3,2,lty=1,col=c(col[c(1,4,5)],"black"),cex=0.5,
           legend = paste0(group_name,": ",c(n1,n2,n3,n4)),
           bty="n")
}
dev.off()


#individual 4 panel for each 47 rbPs
FourCellFLEXI<-dat$ID[rowSums(dat[,88:91])>0]
FourCell<-dat[dat$ID%in%FourCellFLEXI,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FourCellFLEXI,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,50)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]

percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")
RBP_list<-RBP_47_info$RBP.name
RBP_list<-c(sort(RBP_list[c(6:8,15,17,25,27,31,35)]),
            sort(RBP_list[c(1:5,13:14,22,47)]),
            sort(RBP_list[c(18:19)]),
            sort(RBP_list[c(42,20:21)]),
            sort(RBP_list[c(9:12)]),
            sort(RBP_list[c(45:46)]),
            sort(RBP_list[c(16,23,24,26,28:30,32:34,36:41,43:44)]))
pdf("temp_fig/Fig5_sup.pdf",width=8,height=20,onefile = T)
par(mfrow=c(10,4),pch=16,pty="s",mar=c(2,2,2,2))
for (i in 1:47) {
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list[i]])
  FLEXI_dat<-FourCell[FourCell$ID%in%FLEXI_list,]
  temp2<-FourCell[!FourCell$ID%in%FLEXI_list,]
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
    test_set<-FourCell[sample(1:dim(FourCell)[1],sample_size,replace = F),]
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

# plot 6 core splicesome
RBP_list<-c("AQR","BUD13","EFTUD2","PPIG","PRPF8","SF3")
pdf("temp_fig/Fig5_sup2.pdf",width=8,height=20,onefile = T)
par(mfrow=c(10,4),pch=16,pty="s",mar=c(2,2,2,2))
for (i in 1:6) {
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list[i]])
  FLEXI_dat<-FourCell[FourCell$ID%in%FLEXI_list,]
  temp2<-FourCell[!FourCell$ID%in%FLEXI_list,]
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
    test_set<-FourCell[sample(1:dim(FourCell)[1],sample_size,replace = F),]
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
# upset of clusters using humap2
#group_list$group5<-c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2")
#group_list$group6<-c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24")
#group_list$group7<-c("TIA1","TIAL1")
#group_list$group8<-c("U2AF1","U2AF2","KHSRP")
#group_list$group9<-c("AATF","NOLC1","DKC1","SMNDC1")
#group_list$group10<-c("AGO","DICER")
humap<-read.table("humap2_complexes_20200809.txt",sep=",",header=T)
i=9

set_1 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][1],humap$genenames)])
set_2 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][2],humap$genenames)])
set_3 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][3],humap$genenames)])
set_4 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][4],humap$genenames)])
set_5 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][5],humap$genenames)])
set_6 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][6],humap$genenames)])
set_7 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][7],humap$genenames)])
set_8 <- as.character(humap$HuMAP2_ID[grep(group_list[[i]][8],humap$genenames)])

set <- list (set_1,set_2,set_3,set_4,set_5,set_6,set_7,set_8)
set<-list(set_1,set_2,set_3,set_4)
names(set)<-group_list[[i]]
m = make_comb_mat(set)
UpSet(m,
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black","tomato","royalblue1","goldenrod")[comb_degree(m)])

#new Fig1F
##New ID with scale of sncRNA
##Fig1D, RPM density of FLEXIs, different group
dat<-read.delim("all.FLEXI")

gene_counts<-read.delim("combined_counts.tsv")
gene_counts[gene_counts$Type=="scaRNA",3]<-"snoRNA"
snoRNA<-gene_counts[gene_counts$Type=="snoRNA",c(2,8:14)]
snoRNA<-separate(snoRNA,col = "Name",into = "Name",sep = "-",remove = T)
snoRNA$Name<-gsub(snoRNA$Name,pattern = "[A-Z]$|P[1-9]$",replacement="")
snoRNA$Name[snoRNA$Name=="U3"]<-"SNORD3"
snoRNA$Name[snoRNA$Name=="U8"]<-"SNORD118"
snoRNA$Name[snoRNA$Name=="snoU13"]<-"SNORD13"
snoRNA$Name[snoRNA$Name=="snoU2_19"]<-"snoU2"
snoRNA<-aggregate(.~Name,data=snoRNA,FUN = sum)


snRNA<-gene_counts[gene_counts$Type=="snRNA",c(2,8:13)]
U7<-log10(colSums(snRNA[grep("U7\\b|U7[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
U11<-log10(colSums(snRNA[grep("U11\\b|U11[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
SNORD74<-log10(colSums(snoRNA[snoRNA$Name=="SNORD74",2:7])/mapped_reads[5:10])
SNORD78<-log10(colSums(snoRNA[snoRNA$Name=="SNORD78",2:7])/mapped_reads[5:10])
RN7SK<-log10(colSums(gene_counts[gene_counts$Type=="7SK",8:13])/mapped_reads[5:10])
RN7SL<-log10(colSums(gene_counts[gene_counts$Type=="7SL",8:13])/mapped_reads[5:10])
YRNA<-log10(colSums(gene_counts[gene_counts$Type=="YRNA",8:13])/mapped_reads[5:10])
VTRNA<-log10(colSums(gene_counts[gene_counts$Type=="VTRNA",8:13])/mapped_reads[5:10])
RMRP=log10(colSums(gene_counts[gene_counts$Name=="RMRP",8:13])/mapped_reads[5:10])
RPPH1=log10(colSums(gene_counts[gene_counts$Name=="RPPH1",8:13])/mapped_reads[5:10])
U6<-log10(colSums(snRNA[grep("U6\\b|U6[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
U12<-log10(colSums(snRNA[grep("U12\\b|U12[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
SNORD14<-log10(colSums(gene_counts[gene_counts$Name=="SNORD14B",8:13])/mapped_reads[5:10])
U1<-log10(colSums(snRNA[grep("U1\\b|U1[A-Z]",snRNA$Name),2:7])/mapped_reads[5:10])
#SNORD44<-log10(c(84550,40187,29285,565526)/mapped_reads[7:10])
snc<-data.frame(rbind(SNORD74,SNORD78,U11,U7,U12,U6,YRNA,U1,SNORD14,VTRNA,RN7SK,RN7SL))
snc<-snc[c(1:6,8:9,12),]
pdf("temp_fig/Fig1F.pdf",onefile = T,width=12,height=8)
par(mfrow=c(3,2),lwd=1.5)
D_height<-c(2,2,2,2,2,2)
for (i in c(86:91)){
  plot(density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81])),bty="n",xlab="RPM",
       xlim=c(-4,6),ylim=c(0,D_height[i-85]),main=colnames(dat)[i],col="deepskyblue2",axes=F)
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81])),col="firebrick2")
  lines(density(log10(dat[dat[,i]>0 & dat$Is_mirtron=="." & dat$Is_agotron=="." & dat$Has_snoRNA==".",
                          i]/mapped_reads[i-81])),col="black")
  lines(density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81])),col="goldenrod")
  lines(density(log10(snoRNA[snoRNA[,i-84]>0,i-84]/mapped_reads[i-81])),col="goldenrod",lty=4)
  snc<-snc[order(snc[,i-85],decreasing = F),]
  for (j in (1:dim(snc)[1])){
    if (j==1){
      segments(x0=snc[j,i-85],y0=1.5,y1=1.8,lty=1)
      text (snc[j,i-85],1.85,rownames(snc)[j],cex=0.7,adj=0)
    } else {
      segments(x0=snc[j,i-85],y0=0.5,y1=1.8-j*0.12,lty=1)
      text (snc[j,i-85],1.85-j*0.12,rownames(snc)[j],cex=0.7,adj=0)
    }
  }
  if (i==86){
    legend(0.5,2,bty="n",legend = c("Other FLEXIs", "Agotron","Mirtron","snoRNA FLEXIs","snoRNAs"),
           col=c("black","firebrick2","deepskyblue2","goldenrod","goldenrod"),
           lty=c(1,1,1,1,4),lwd=1.5,cex=0.7)
  }
  axis(2,labels=seq(0,2,1),las=1,at=seq(0,2,1),las=2)
  axis(1,labels=c(parse(text='10^-4'),bquote(10^-2),1,bquote(10^2),bquote(10^4),
                  bquote(10^6)),at=seq(-4,6,2))
}
dev.off()

for (i in 86:90){
#  a<-density(log10(dat[dat[,i]>0 & dat$Is_agotron!=".",i]/mapped_reads[i-81]))
#  a<-density(log10(dat[dat[,i]>0 & dat$Is_mirtron!=".",i]/mapped_reads[i-81]))
  a<-density(log10(dat[dat[,i]>0 & dat$Has_snoRNA!=".",i]/mapped_reads[i-81]))
  print (10^a$x[a$y==max(a$y)])
}

#new fig S10
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>0,]
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))

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


B_col=c("black","red","orange","skyblue")
RBP_fre$label=""
RBP_fre$label[RBP_fre$Cell>=30]<-"*"

pdf("temp_fig//FigS9.pdf",height=24,width=12)
par(mfrow=c(3,1))
tmp<-RBP_fre[RBP_fre$Cell>0,]
tmp<-tmp[order(tmp$Cell,decreasing = T),]
mp<-barplot(tmp$Cell,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[tmp$col],ylab="FLEXI RNAs")
text(mp,10000,labels=tmp$label,cex=0.5)
text(mp,9000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)

#other short intron
tmp<-RBP_fre[RBP_fre$`Other Short introns`>0,]
tmp<-tmp[order(tmp$`Other Short introns`,decreasing = T),]
mp<-barplot(tmp$`Other Short introns`,ylim=c(1,10000),log="y",cex.names=0.5,col=B_col[tmp$col],
            ylab="Other short introns")
text(mp,10000,labels=tmp$label,cex=0.5)
text(mp,9000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)
#long intron
tmp<-RBP_fre[RBP_fre$`Long introns`>0,]
tmp<-tmp[order(tmp$`Long introns`,decreasing = T),]
mp<-barplot(tmp$`Long introns`,ylim=c(1,100000),log="y",cex.names=0.5,col=B_col[tmp$col],
            ylab="Long Introns")
text(mp,100000,labels=tmp$label,cex=0.5)
text(mp,80000,labels=tmp$RBP,col=B_col[tmp$col],cex=0.5,srt=90)

dev.off()

'''
#S10C
Fun<-read.table("4cell_plasma_RBP_by_FLEXI.counts",col.names=c("Name","RBP_by_FLEXI"))
topRBP<-Fun$Name[Fun$RBP_by_FLEXI>=30]
Fun<-read.table(gzfile("All_exon_inter_RBP152_collapsed.info.gz")
                ,col.names=c("ID","RBP"))
Fun<-data.frame(table(Fun$RBP))
colnames(Fun)<-c("Name","RBP_byExon")
Fun<-merge(Fun,RBP[,c(1,4,5,11)],by=1,all.x=T)
Fun<-Fun[order(Fun$RBP_byExon,decreasing = T),]
Fun$Color<-(Fun$Splicing.regulation+Fun$Spliceosome)/3+Fun$microRNA.processing
Fun[Fun$Color==1,6]<-2
Fun[Fun$Color==4/3,6]<-3
Fun[Fun$Color==1/3 | Fun$Color==2/3,6]<-1
Fun$Color<-factor(Fun$Color)
B_col=c("black","red","orange","skyblue")
Fun$label=""
Fun$label[Fun$Name%in%topRBP]<-"*"
pdf("temp_fig//FigS10C.pdf",height=12,width=12)
mp<-barplot(Fun$RBP_byExon,ylim=c(1,100000),log="y",cex.names=0.5,col=B_col[Fun$Color],
            names.arg =Fun$Name,las=2)
text(mp,100000,labels=Fun$label,cex=0.5)
par(xpd=T)
legend("topright",legend = Fun$Name,text.col=B_col[Fun$Color],cex=0.3,bty="n")
par(xpd=F)
dev.off()


#S10def
Fun<-read.table("4cell_plasma_RBP_by_FLEXI.counts",col.names=c("Name","RBP_by_FLEXI"))
topRBP<-Fun$Name[Fun$RBP_by_FLEXI>=30]
Fun<-read.table("tRNA_RBP_inter.info"
                ,col.names=c("Name","ID","RBP"))
Fun<-data.frame(table(Fun$RBP))
colnames(Fun)<-c("Name","RBP_byExon")
Fun<-merge(Fun,RBP[,c(1,4,5,11)],by=1,all.x=T)
Fun<-Fun[order(Fun$RBP_byExon,decreasing = T),]
Fun$Color<-(Fun$Splicing.regulation+Fun$Spliceosome)/3+Fun$microRNA.processing
Fun[Fun$Color==1,6]<-2
Fun[Fun$Color==4/3,6]<-3
Fun[Fun$Color==1/3 | Fun$Color==2/3,6]<-1
Fun$Color<-factor(Fun$Color)
B_col=c("black","red","orange","skyblue")
Fun$label=""
Fun$label[Fun$Name%in%topRBP]<-"*"
pdf("temp_fig//FigS10f.pdf",height=12,width=12)
mp<-barplot(Fun$RBP_byExon,ylim=c(1,1000),log="y",cex.names=0.5,col=B_col[Fun$Color],
            names.arg =Fun$Name,las=2)
text(mp,1000,labels=Fun$label,cex=0.5)
par(xpd=T)
legend("topright",legend = Fun$Name,text.col=B_col[Fun$Color],cex=0.3,bty="n")
par(xpd=F)
dev.off()
'''
#linear regression with old and new corrcted counts

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

#mappedreads and corrected mapped reads
correct_mapped_reads<-c(715.24152,277.8381,768.43375,368.4548,666.34185,362.6804,713.77529,230.0604)
names(correct_mapped_reads)<-rownames(gene_for_correlation)
gene_for_correlation<-gene_for_correlation/correct_mapped_reads
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",2),"7SK","7SL","RNase P","MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)

temp1<-temp1/correct_mapped_reads
#temp1<-rbind(temp1,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
temp1<-data.frame(t(temp1))
temp1<-round(temp1,0)
  
pdf("temp_fig/FIgS61.pdf",height=6,width=4)
par(bty="n",mfrow=c(2,2),mar=c(1,1,1,1))
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
  cor_s=formatC(cor(x,y,method = "spearman"),digits=2, format="f")
  cor_p=formatC(cor(x,y,method = "pearson"),digits=2, format="f")
  text(1,7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
}
dev.off()
pdf("temp_fig/FIgS62.pdf",height=6,width=4)
par(bty="n",mfrow=c(2,2),mar=c(1,1,1,1))
scol <- brewer.pal(8, "Set1")
y<-log10(gene_for_correlation[,9])
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
  cor_s=formatC(cor(x,y,method = "spearman"),digits=2, format="f")
  cor_p=formatC(cor(x,y,method = "pearson"),digits=2, format="f")
  text(1,7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
}
dev.off()
# calculate expected RPM for FLEXIs used in digital PCR
check_list<-c(2232,5094,7563,6003,5726,2183,549,1489,3377)
RPM_eve<-FourCell[rownames(FourCell)%in%check_list,c(1,88:91)]
RPM_eve<-RPM_eve[c(4,6,9,8,7,3,1,2,5),]
PRM_standard<-c(0.01,1)
y<-log10(gene_for_correlation[,9])
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i])
  lm.out <- lm(y ~ x)
#  newx = log10(RPM_eve[,i+1]/mapped_reads[i+6])
#  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
#                                      level = 0.95))
#  print(t(t(round(10^conf_interval$fit,0))))
  newx = log10(PRM_standard)
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                                      level = 0.95))
  print(t(t(round(10^conf_interval$fit,0))))
}

#dPCR 4 standard
temp<-gene_counts[gene_counts$Type=="snRNA",c(4,9:16)]
dPCR_standards<-data.frame("U7"=colSums(temp[grep("U7\\b|U7[A-Z]",temp$Name),2:9]))
dPCR_standards<-cbind(dPCR_standards,"U11"=colSums(temp[grep("U11\\b|U11[A-Z]",temp$Name),2:9]))
dPCR_standards<-cbind(dPCR_standards,
                            data.frame("U12"=colSums(temp[grep("U12\\b|U12[A-Z]",temp$Name),2:9])))
temp<-gene_counts[gene_counts$Type=="snoRNA",c(4,9:16)]
dPCR_standards<-cbind(dPCR_standards,
                            data.frame("SNORD14B"=colSums(temp[temp$Name=="SNORD14B",2:9])))
dPCR_standards<-cbind(dPCR_standards,
                      data.frame("SNORD44"=colSums(temp[temp$Name=="SNORD44",2:9])))
dPCR_standards<-dPCR_standards/correct_mapped_reads
y<-log10(gene_for_correlation[,9])
for (i in 1:8){
  x<-log10(gene_for_correlation[,i])
  lm.out <- lm(y ~ x)
  newx = log10(unlist(dPCR_standards[i,]))
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx)))
  print(t(t(round(10^conf_interval,0))))
}

check_list<-c(2232,5094,7563,6003,5726,2183,549,1489,3377,2568)
RPM_eve<-FourCell[rownames(FourCell)%in%check_list,c(1,90:91,88)]
RPM_eve<-RPM_eve[c(4,7,10,9,8,3,1,2,6,5),]
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
  print(t(t(round(10^conf_interval,0))))
}
#FLEIX_corrected
FLEXI_corr<-read.delim("FLEXI_corrected.counts")
ago_list<-FourCell$ID[FourCell$Is_agotron!="."]
tmp<-FLEXI_corr[FLEXI_corr$ID%in%ago_list,]
corrected_range<-data.frame(apply(tmp[,2:4],2,function(x){range(x[x>0])}))
ago_list<-FourCell$ID[FourCell$Is_mirtron!="."]
tmp<-FLEXI_corr[FLEXI_corr$ID%in%ago_list,]
corrected_range<-rbind(corrected_range,data.frame(apply(tmp[,2:4],2,function(x){range(x[x>0])})))
ago_list<-FourCell$ID[FourCell$Has_snoRNA!="."]
tmp<-FLEXI_corr[FLEXI_corr$ID%in%ago_list,]
corrected_range<-rbind(corrected_range,data.frame(apply(tmp[,2:4],2,function(x){range(x[x>0])})))
rownames(corrected_range)<-c("AgoL","AgoH","MirL","MirH","SnoL","SnoH")
corrected_range<-t(t(corrected_range)/correct_mapped_reads[c(2,4,6)])
y<-log10(gene_for_correlation[,9])
for (i in 1:3){
  pick_col<-2*i
  x<-log10(gene_for_correlation[,pick_col])
  lm.out <- lm(y ~ x)
  newx = log10(unlist(corrected_range[,i]))
  conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx)))
  print(t(t(round(10^conf_interval,0))))
}


# compare cellular FLEXI to plasma FLEXI
Plasma_RBP<-read.table("Plasma_FLEXI.frq",col.names="RBP")
Plasma_RBP<-data.frame(table(c(Plasma_RBP$RBP,rep("AGO",48),rep("DICER",59))))
colnames(Plasma_RBP)<-c("RBP","Plasma")
Plasma_RBP<-merge(Plasma_RBP,RBP_fre,by=1,all=T)
Plasma_RBP<-merge(Plasma_RBP,Genome,by=1,all=T)
Plasma_RBP[is.na(Plasma_RBP)]<-0

Plasma_RBP<-merge(Plasma_RBP,RBP[,c(1,4,5,11)],by=1)
Plasma_RBP$col<-(Plasma_RBP$Splicing.regulation+Plasma_RBP$Spliceosome)/3+Plasma_RBP$microRNA.processing
Plasma_RBP<-Plasma_RBP[,c(1,3,2,4,5)]
Plasma_RBP$Bpvlue<-1
Plasma_RBP$Cpvlue<-1
#Dpvalue is for plasma FLEXI vs Other FLEXI
Plasma_RBP$Dpvlue<-1
#set BCLAF1 splicing function
Plasma_RBP[11,2]<-1/3
R_sum<-colSums(Plasma_RBP[,3:5])
for (i in 1:152){
  Plasma_RBP[i,6]<-fisher.test(as.matrix(rbind(Plasma_RBP[i,c(3,5)],R_sum[c(1,3)])))$p.value
  Plasma_RBP[i,7]<-fisher.test(as.matrix(rbind(Plasma_RBP[i,c(4,5)],R_sum[c(2,3)])))$p.value
  Plasma_RBP[i,8]<-fisher.test(as.matrix(rbind(Plasma_RBP[i,c(3,4)],R_sum[c(1,2)])))$p.value
}
Plasma_RBP[,3:5]<-data.frame(prop.table(as.matrix(Plasma_RBP[,3:5]),margin = 2)*100)
Plasma_RBP$Bpvlue<-p.adjust(Plasma_RBP$Bpvlue,method="fdr")
Plasma_RBP$Cpvlue<-p.adjust(Plasma_RBP$Cpvlue,method="fdr")
Plasma_RBP$Dpvlue<-p.adjust(Plasma_RBP$Dpvlue,method="fdr")

pdf("temp_fig/Plasma_vsALL.pdf",height=12,width=8)
par(mfrow=c(3,2),pty="s",pch=16)
plot(Plasma_RBP[Plasma_RBP$col==0,c(5,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="GRCh38 (% RBP sites)",ylab="Plasma FLEXIs (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(5,3)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(5,3)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(5,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Bpvlue<=0.01 & (Plasma_RBP[,3]>=4 | Plasma_RBP[,5]>=4)),c(5,3)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Bpvlue<=0.01 & (Plasma_RBP[,3]>=4 | Plasma_RBP[,5]>=4)],
     cex=0.75)

plot(Plasma_RBP[Plasma_RBP$col==0,c(5,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     xlab="GRCh38 (% RBP sites)",ylab="Plasma FLEXIs (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(5,3)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(5,3)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(5,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Bpvlue<=0.01 & (Plasma_RBP[,3]>=2 | Plasma_RBP[,5]>=2)),c(5,3)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Bpvlue<=0.01 & (Plasma_RBP[,3]>=2 | Plasma_RBP[,5]>=2)],
     cex=0.75)


plot(Plasma_RBP[Plasma_RBP$col==0,c(5,4)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="All FLEXIs (% RBP sites)",xlab="GRCh38 (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(5,4)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(5,4)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(5,4)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Cpvlue<=0.01 & (Plasma_RBP[,4]>=4 | Plasma_RBP[,5]>=4)),c(5,4)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Cpvlue<=0.01 & (Plasma_RBP[,5]>=4 | Plasma_RBP[,4]>=4)],
     cex=0.75)

plot(Plasma_RBP[Plasma_RBP$col==0,c(5,4)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="All FLEXIs (% RBP sites)",xlab="GRCh38 (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(5,4)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(5,4)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(5,4)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Cpvlue<=0.01 & (Plasma_RBP[,4]>=2 | Plasma_RBP[,5]>=2)),c(5,4)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Cpvlue<=0.01 & (Plasma_RBP[,4]>=2 | Plasma_RBP[,5]>=2)],
     cex=0.75)

plot(Plasma_RBP[Plasma_RBP$col==0,c(4,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="Plasma FLEXIs (% RBP sites)",xlab="Celllular FLEXIs (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(4,3)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(4,3)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(4,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Dpvlue<=0.01 & (Plasma_RBP[,4]>=4 | Plasma_RBP[,3]>=4)),c(4,3)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Dpvlue<=0.01 & (Plasma_RBP[,3]>=4 | Plasma_RBP[,4]>=4)],
     cex=0.75)

plot(Plasma_RBP[Plasma_RBP$col==0,c(4,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="Plasma FLEXIs (% RBP sites)",xlab="Celllular FLEXIs (% RBP sites)")
points(Plasma_RBP[Plasma_RBP$col<1 & Plasma_RBP$col>0,c(4,3)],col="red",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col==1,c(4,3)],col="orange",cex=1.5)
points(Plasma_RBP[Plasma_RBP$col>1,c(4,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Plasma_RBP[(Plasma_RBP$Dpvlue<=0.01 & (Plasma_RBP[,4]>=2 | Plasma_RBP[,3]>=2)),c(4,3)],
     labels = Plasma_RBP$RBP[Plasma_RBP$Dpvlue<=0.01 & (Plasma_RBP[,4]>=2 | Plasma_RBP[,3]>=2)],
     cex=0.75)
dev.off()


#RBP scatter for 0.01 RPM
K562001RPM<-dat_CPM$ID[dat_CPM$K562>=0.01]
HEK001RPM<-dat_CPM$ID[dat_CPM$HEK>=0.01]
Hela001RPM<-dat_CPM$ID[dat_CPM$Hela>=0.01]
MDA001RPM<-dat_CPM$ID[dat_CPM$MDA>=0.01]
MCF001RPM<-dat_CPM$ID[dat_CPM$MCF7>=0.01]
All_FLEXIRBPInfo<-read.table("all_FLEXI_152RBP.info",col.names=c("ID","RBP"))

K562001RPM<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%K562001RPM,]
K562001RPM<-data.frame(table(K562001RPM$RBP))

HEK001RPM<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%HEK001RPM,]
HEK001RPM<-data.frame(table(HEK001RPM$RBP))

Hela001RPM<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%Hela001RPM,]
Hela001RPM<-data.frame(table(Hela001RPM$RBP))

MDA001RPM<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%MDA001RPM,]
MDA001RPM<-data.frame(table(MDA001RPM$RBP))

MCF001RPM<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%MCF001RPM,]
MCF001RPM<-data.frame(table(MCF001RPM$RBP))

RPM001_RBP<-merge(K562001RPM,HEK001RPM,by=1,all=T)
RPM001_RBP<-merge(RPM001_RBP,Hela001RPM,by=1,all=T)
RPM001_RBP<-merge(RPM001_RBP,MDA001RPM,by=1,all=T)
RPM001_RBP<-merge(RPM001_RBP,MCF001RPM,by=1,all=T)

RPM001_RBP[is.na(RPM001_RBP)]<-0
colnames(RPM001_RBP)<-c("RBP","K562","HEK","Hela","MDA","MCF")

RPM001_RBP<-merge(RPM001_RBP,RBP[,c(1,4,5,11)],by=1)
RPM001_RBP$col<-(RPM001_RBP$Splicing.regulation+RPM001_RBP$Spliceosome)/3+RPM001_RBP$microRNA.processing
RPM001_RBP<-RPM001_RBP[,c(1:6,10)]
RPM001_RBP$pvlue<-1
#set BCLAF1 splicing function
RPM001_RBP[9,7]<-1/3
R_sum<-colSums(RPM001_RBP[,2:6])

com_list<-list(c(2,3),c(2,4),c(4,3),c(5,6))
pdf("temp_fig/RPM001.pdf",height=16,width=8)
par(mfrow=c(4,2),pty="s",pch=16)
for (i in (1:4)){
  comp<-com_list[[i]]
  tmp<-RPM001_RBP
  tmp<-tmp[rowSums(tmp[,comp])>0,]
  for (j in 1:dim(tmp)[1]){
    tmp[j,8]<-fisher.test(as.matrix(rbind(tmp[j,comp],R_sum[comp-1])))$p.value
  }
  tmp[,2:6]<-data.frame(prop.table(as.matrix(tmp[,2:6]),margin = 2)*100)
  tmp$pvlue<-p.adjust(tmp$pvlue,method="fdr")
  
  plot(tmp[tmp$col==0,comp],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n")
  points(tmp[tmp$col<1 & tmp$col>0,comp],col="red",cex=1.5)
  points(tmp[tmp$col==1,comp],col="orange",cex=1.5)
  points(tmp[tmp$col>1,comp],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)),comp]
  if (dim(sig)[1]>0){
    text(tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)),comp],
         labels = tmp$RBP[tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)],
         cex=0.75)
  }
  
  plot(tmp[tmp$col==0,comp],,xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n")
  points(tmp[tmp$col<1 & tmp$col>0,comp],col="red",cex=1.5)
  points(tmp[tmp$col==1,comp],col="orange",cex=1.5)
  points(tmp[tmp$col>1,comp],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)),comp]
  if (dim(sig)[1]>0){
    text(tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)),comp],
         labels = tmp$RBP[tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)],
         cex=0.75)
  }
}
dev.off()

#Ago
PlasmaAgo<-dat$ID[dat$Plasma>0 &dat$Is_agotron!="."]
CellAgo<-FourCell$ID[FourCell$Is_agotron!="."]
PlasmaAgo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%PlasmaAgo,]
PlasmaAgo<-data.frame(table(PlasmaAgo$RBP))
CellAgo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%CellAgo,]
CellAgo<-data.frame(table(CellAgo$RBP))
AgoRBP<-merge(PlasmaAgo,CellAgo,by=1,all=T)
AgoRBP[is.na(AgoRBP)]<-0
colnames(AgoRBP)<-c("RBP","Plasma","Cellular")

AgoRBP<-merge(AgoRBP,RBP[,c(1,4,5,11)],by=1)
AgoRBP$col<-(AgoRBP$Splicing.regulation+AgoRBP$Spliceosome)/3+AgoRBP$microRNA.processing
AgoRBP<-AgoRBP[,c(1:3,7)]
AgoRBP$pvlue<-1
#set BCLAF1 splicing function
AgoRBP[3,4]<-1/3
R_sum<-colSums(AgoRBP[,2:3])
for (j in 1:dim(AgoRBP)[1]){
  AgoRBP[j,5]<-fisher.test(as.matrix(rbind(AgoRBP[j,2:3],R_sum)))$p.value
}
AgoRBP[,2:3]<-data.frame(prop.table(as.matrix(AgoRBP[,2:3]),margin = 2)*100)
AgoRBP$pvlue<-p.adjust(AgoRBP$pvlue,method="fdr")

pdf("temp_fig/Plasma_cell.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)

plot(AgoRBP[AgoRBP$col==0,2:3],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",main="Agotron FLEXIs")
points(AgoRBP[AgoRBP$col<1 & AgoRBP$col>0,2:3],col="red",cex=1.5)
points(AgoRBP[AgoRBP$col==1,2:3],col="orange",cex=1.5)
points(AgoRBP[AgoRBP$col>1,2:3],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,3]>=4 | AgoRBP[,2]>=4)),2:3]
if (dim(sig)[1]>0){
  text(AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=4 | AgoRBP[,3]>=4)),comp],
       labels = AgoRBP$RBP[AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=4 | AgoRBP[,3]>=4)],
       cex=0.75)
}

plot(AgoRBP[AgoRBP$col==0,2:3],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",main="Agotron FLEXIs")
points(AgoRBP[AgoRBP$col<1 & AgoRBP$col>0,2:3],col="red",cex=1.5)
points(AgoRBP[AgoRBP$col==1,2:3],col="orange",cex=1.5)
points(AgoRBP[AgoRBP$col>1,2:3],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,3]>=2 | AgoRBP[,2]>=2)),2:3]
if (dim(sig)[1]>0){
  text(AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=2 | AgoRBP[,3]>=2)),comp],
       labels = AgoRBP$RBP[AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=2 | AgoRBP[,3]>=2)],
       cex=0.75)
}

PlasmaAgo<-dat$ID[dat$Plasma>0 &dat$Is_mirtron!="."]
CellAgo<-FourCell$ID[FourCell$Is_mirtron!="."]
PlasmaAgo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%PlasmaAgo,]
PlasmaAgo<-data.frame(table(PlasmaAgo$RBP))
CellAgo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%CellAgo,]
CellAgo<-data.frame(table(CellAgo$RBP))
AgoRBP<-merge(PlasmaAgo,CellAgo,by=1,all=T)
AgoRBP[is.na(AgoRBP)]<-0
colnames(AgoRBP)<-c("RBP","Plasma","Cellular")

AgoRBP<-merge(AgoRBP,RBP[,c(1,4,5,11)],by=1)
AgoRBP$col<-(AgoRBP$Splicing.regulation+AgoRBP$Spliceosome)/3+AgoRBP$microRNA.processing
AgoRBP<-AgoRBP[,c(1:3,7)]
AgoRBP$pvlue<-1
#set BCLAF1 splicing function
AgoRBP[3,4]<-1/3
R_sum<-colSums(AgoRBP[,2:3])
for (j in 1:dim(AgoRBP)[1]){
  AgoRBP[j,5]<-fisher.test(as.matrix(rbind(AgoRBP[j,2:3],R_sum)))$p.value
}
AgoRBP[,2:3]<-data.frame(prop.table(as.matrix(AgoRBP[,2:3]),margin = 2)*100)
AgoRBP$pvlue<-p.adjust(AgoRBP$pvlue,method="fdr")
plot(AgoRBP[AgoRBP$col==0,2:3],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",main="Mirtron FLEXIs")
points(AgoRBP[AgoRBP$col<1 & AgoRBP$col>0,2:3],col="red",cex=1.5)
points(AgoRBP[AgoRBP$col==1,2:3],col="orange",cex=1.5)
points(AgoRBP[AgoRBP$col>1,2:3],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,3]>=4 | AgoRBP[,2]>=4)),2:3]
if (dim(sig)[1]>0){
  text(AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=4 | AgoRBP[,3]>=4)),comp],
       labels = AgoRBP$RBP[AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=4 | AgoRBP[,3]>=4)],
       cex=0.75)
}

plot(AgoRBP[AgoRBP$col==0,2:3],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",main="Mirtron FLEXIs")
points(AgoRBP[AgoRBP$col<1 & AgoRBP$col>0,2:3],col="red",cex=1.5)
points(AgoRBP[AgoRBP$col==1,2:3],col="orange",cex=1.5)
points(AgoRBP[AgoRBP$col>1,2:3],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,3]>=2 | AgoRBP[,2]>=2)),2:3]
if (dim(sig)[1]>0){
  text(AgoRBP[(AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=2 | AgoRBP[,3]>=2)),comp],
       labels = AgoRBP$RBP[AgoRBP$pvlue<=0.01 & (AgoRBP[,2]>=2 | AgoRBP[,3]>=2)],
       cex=0.75)
}
dev.off()

#reproducible
K562repo<-rownames(FLEXI_by_GID)[FLEXI_by_GID$K562_repo>=8]
Helarepo<-rownames(FLEXI_by_GID)[FLEXI_by_GID$Hela_repo>=10]
HEKrepo<-rownames(FLEXI_by_GID)[FLEXI_by_GID$HEK_repo>=8]

All_FLEXIRBPInfo<-read.table("all_FLEXI_152RBP.info",col.names=c("ID","RBP"))

K562repo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%K562repo,]
K562repo<-data.frame(table(K562repo$RBP))

Helarepo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%Helarepo,]
Helarepo<-data.frame(table(Helarepo$RBP))

HEKrepo<-All_FLEXIRBPInfo[All_FLEXIRBPInfo$ID%in%HEKrepo,]
HEKrepo<-data.frame(table(HEKrepo$RBP))

RBPrepo<-merge(K562repo,Helarepo,by=1,all=T)
RBPrepo<-merge(RBPrepo,HEKrepo,by=1,all=T)
RBPrepo[is.na(RBPrepo)]<-0


colnames(RBPrepo)<-c("RBP","K562","Hela","HEK")

RBPrepo<-merge(RBPrepo,RBP[,c(1,4,5,11)],by=1)
RBPrepo$col<-(RBPrepo$Splicing.regulation+RBPrepo$Spliceosome)/3+RBPrepo$microRNA.processing
RBPrepo<-RBPrepo[,c(1:4,8)]
RBPrepo$pvlue<-1
#set BCLAF1 splicing function
RBPrepo[5,5]<-1/3
R_sum<-colSums(RBPrepo[,2:4])

com_list<-list(c(2,3),c(2,4),c(4,3))
pdf("temp_fig/reproducible.pdf",height=12,width=8)
par(mfrow=c(3,2),pty="s",pch=16)
for (i in (1:3)){
  comp<-com_list[[i]]
  tmp<-RBPrepo
  tmp<-tmp[rowSums(tmp[,comp])>0,]
  for (j in 1:dim(tmp)[1]){
    tmp[j,6]<-fisher.test(as.matrix(rbind(tmp[j,comp],R_sum[comp-1])))$p.value
  }
  tmp[,2:4]<-data.frame(prop.table(as.matrix(tmp[,2:4]),margin = 2)*100)
  tmp$pvlue<-p.adjust(tmp$pvlue,method="fdr")
  
  plot(tmp[tmp$col==0,comp],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n")
  points(tmp[tmp$col<1 & tmp$col>0,comp],col="red",cex=1.5)
  points(tmp[tmp$col==1,comp],col="orange",cex=1.5)
  points(tmp[tmp$col>1,comp],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)),comp]
  if (dim(sig)[1]>0){
    text(tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)),comp],
         labels = tmp$RBP[tmp$pvlue<=0.01 & (tmp[,comp[1]]>=4 | tmp[,comp[2]]>=4)],
         cex=0.75)
  }
  
  plot(tmp[tmp$col==0,comp],,xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n")
  points(tmp[tmp$col<1 & tmp$col>0,comp],col="red",cex=1.5)
  points(tmp[tmp$col==1,comp],col="orange",cex=1.5)
  points(tmp[tmp$col>1,comp],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)),comp]
  if (dim(sig)[1]>0){
    text(tmp[(tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)),comp],
         labels = tmp$RBP[tmp$pvlue<=0.01 & (tmp[,comp[1]]>=2 | tmp[,comp[2]]>=2)],
         cex=0.75)
  }
}
dev.off()

#2-D countour
conm<-read.delim("4cell.2D.txt")
conm<-rbind(conm,c(-3,3,0))
pdf("temp_fig/2d-density.pdf")
ggplot(conm, aes(x = SS5, y = SS3,z=Cellular)) + 
  geom_contour_filled(color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
dev.off()

## ohastcon 0.75 scatter plot
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>0,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP","Cell")
Phast99<-dat$ID[dat$PhastCon30>=0.75 & rowSums(dat[,88:91])>0]
Phast99_fre<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%Phast99,]
Phast99_fre<-data.frame(table(Phast99_fre$RBP))
colnames(Phast99_fre)<-c("RBP","Phast99")
Phast<-merge(RBP_fre,Phast99_fre,by=1,all=T)
Phast[is.na(Phast)]<-0
rownames(Phast)<-Phast$ID
Phast$pvalue<-1
R_sum<-colSums(Phast[,2:3])
for (i in 1:126){
  Phast[i,4]<-fisher.test(as.matrix(rbind(Phast[i,2:3],R_sum)))$p.value
}
Phast$padj<-p.adjust(Phast$pvalue,method="fdr")
Phast[,2:3]<-data.frame(prop.table(as.matrix(Phast[,2:3]),margin = 2)*100)
Phast<-merge(Phast,RBP[,c(1,4,5,11)],by=1)
Phast$col<-(Phast$Splicing.regulation+Phast$Spliceosome)/3+Phast$microRNA.processing

pdf("Figures/Fig3D.pdf",height=6,width=12)
par(mfrow=c(1,2),pty="s",pch=16)
plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="All FELXIs (% RBP sites)",ylab="Conserved FELXIs (% RBP sites)")
points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-Phast$padj<=0.05 & Phast$Phast99>=4
text(Phast[sig,c(2,3)],pos = 3,cex=0.5,
     labels = Phast$RBP[sig])

plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     xlab="All FELXIs (% RBP sites)",ylab="Conserved FELXIs (% RBP sites)")
points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
sig<-Phast$padj<=0.05 & Phast$Phast99>=2 & Phast$Phast99<4
text(Phast[sig,c(2,3)],pos = 3,cex=0.5,
     labels = Phast$RBP[sig])
dev.off()