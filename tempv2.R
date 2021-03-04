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

#combined cellular RNA, density plot of FLEXIs vs other short intron, Fig1C left panel
pdf("Figures/Fig1C_1.pdf")
plot(density(combined$Per[combined$nonFLEXI==1]),bty="n",lwd=1, ylim=c(0,0.1),xlim=c(0,100),main=NA,
         xlab="Length (%)",col="black")
lines(density(combined$Per[combined$nonFLEXI!=1]),lwd=1,col="red")
legend(20,0.08,lty=c(1,1),lwd=1,col=c("red","black"),
           legend = c("Other short introns","FLEXIs"),bty="n")
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
cut_off=1
set_1 <- as.character(dat$ID[dat$UHRR>=cut_off])
set_2 <- as.character(dat$ID[dat$K562>=cut_off])
set_3 <- as.character(dat$ID[dat$HEK>=cut_off])
set_4 <- as.character(dat$ID[dat$Hela>=cut_off])
set_5 <- as.character(dat$ID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK 293T"=set_3,
             "Hela S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
postscript("Figures/Fig1A_1.eps",height=4,width=8)
UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
dev.off()
#FLEXI host genes
set_1 <- unique(dat$GID[dat$UHRR>=cut_off])
set_2 <- unique(dat$GID[dat$K562>=cut_off])
set_3 <- unique(dat$GID[dat$HEK>=cut_off])
set_4 <- unique(dat$GID[dat$Hela>=cut_off])
set_5 <- unique(dat$GID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK 293T"=set_3,
             "Hela S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
postscript("Figures/Fig1A_2.eps",height=4,width=8)
UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
dev.off()
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","m","cut_off"))

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

#FIg3D, RBP scatter plot in conserved FELXIs, phastCons ≥ 0.99
Phast99<-read.table("44phastCon99.RBP.counts",col.names=c("ID","P99"))
Phast<-merge(Cell,Phast99,by=1,all=T)
Phast[is.na(Phast)]<-0
rownames(Phast)<-Phast$ID
Phast<-data.frame(prop.table(as.matrix(Phast[,2:3]),margin = 2)*100)
Phast$Name<-rownames(Phast)
Phast<-Phast[,c(3,1:2)]
Phast<-merge(Phast,RBP[,c(1,4,5,11)],by=1)
Phast$col<-(Phast$Splicing.regulation+Phast$Spliceosome)/3+Phast$microRNA.processing
Phast<-Phast[,c(1:3,7)]
pdf("Figures/Fig3D.pdf",height=6,width=12)
par(mfrow=c(1,2),pty="s",pch=16)
plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="All FELXIs (% RBP sites)",ylab="Conserved FELXIs (% RBP sites)")
points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Phast[Phast$P99/Phast$Cells>2 & Phast$P99>3,c(2,3)],pos = 3,cex=0.5,
     labels = Phast$Name[Phast$P99/Phast$Cells>2 & Phast$P99>3])

plot(Phast[Phast$col==0,c(2,3)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     xlab="All FELXIs (% RBP sites)",ylab="Conserved FELXIs (% RBP sites)")
points(Phast[Phast$col<1 & Phast$col>0,c(2,3)],col="red",cex=1.5)
points(Phast[Phast$col==1,c(2,3)],col="orange",cex=1.5)
points(Phast[Phast$col>1,c(2,3)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Phast[Phast$P99/Phast$Cells>2 & Phast$P99>=2,c(2,3)],pos = 3,cex=0.5,
     labels = Phast$Name[Phast$P99/Phast$Cells>2 & Phast$P99>=2])
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

Four_cancer_RBP<-merge(MCF_unique,MDA_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,All_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Cell,by=1,all=T)
Four_cancer_RBP[is.na(Four_cancer_RBP)]<-0
Four_cancer_RBP$RBP<-as.character(Four_cancer_RBP$RBP)

Four_cancer_RBP[Four_cancer_RBP$RBP=="AGO",2:4]<-AGO_sites
Four_cancer_RBP[Four_cancer_RBP$RBP=="DICER",2:4]<-DICER_sites

pdf("Figures/Fig8E.pdf",height=6,width=4)
par(mfrow=c(3,2),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
Name=c("MCF7","MDA-MB-231","All cancer")
for (i in 1:3){
  dat_set<-Four_cancer_RBP[,c(1,i+1,5)]
  RBP_sum<-colSums(dat_set[,2:3])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)]),2,2))$p.value})
  dat_set[,2:3]<-100*prop.table(as.matrix(dat_set[,2:3]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name[i],axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
  text(dat_set[dat_set$C_sig<=0.05 & (dat_set[,3]>=6 |dat_set[,2]>=6),c(3,2)],
       labels = dat_set[dat_set$C_sig<=0.05 & (dat_set[,3]>=6 |dat_set[,2]>=6),1],pos = 3,cex=0.7)
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
  plot(dat_set[,c(3,2)],xlim=c(0,6),ylim=c(0,6),bty="n",axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
  text(dat_set[dat_set$C_sig<=0.05,c(3,2)],
       labels = dat_set[dat_set$C_sig<=0.05,1],pos = 3,cex=0.7)
  if (i<3) {
    axis(1,at=seq(0,6,2),labels = NA)
  } else {
    axis(1,at=seq(0,6,2),labels = seq(0,6,2))
    title(xlab="All FLEXIs (% RBP sites)")
  }
  axis(2,at=seq(0,6,2),labels = NA)
}
dev.off()
#FigS6
pdf("Figures/FigS6.pdf",width=11,height=8,onefile = T)
RBP_list_sig<-c("DKC1","NOLC1","AATF","BCLAF1","GRWD1","SRSF1","TIA1","UCHL5","U2AF1","U2AF2","ZNF622","TRA2A")
par(mfrow=c(3,4))
for (i in 1:length(RBP_list_sig)){
  # subset
  if (RBP_list_sig[i] =="AGO") {
    temp1<-FourCellPlasma[FourCellPlasma$AGO_CCR!=".",]
    temp2<-FourCellPlasma[FourCellPlasma$AGO_CCR==".",]
  } else if (RBP_list_sig[i] =="DICER") {
    temp1<-FourCellPlasma[FourCellPlasma$DICER_CCR!=".",]
    temp2<-FourCellPlasma[FourCellPlasma$DICER_CCR==".",]
  } else {
    temp1<-FourCellPlasma[grep(RBP_list_sig[i],FourCellPlasma$RBP),]
    temp2<-FourCellPlasma[grep(RBP_list_sig[i],FourCellPlasma$RBP,invert=T),]
  }
  Len_t<-0
  GC_t<-0
  MFE_t<-0
  PhastCon30_t<-0
  rep_times<-100
  for (inter in 1:rep_times){
    sample_size<-Fun$RBP_by_FLEXI[Fun$Name==RBP_list_sig[i]]
    test_set<-FourCellPlasma[sample(1:8144,sample_size,replace = F),]
    if (wilcox.test(temp1$Len,test_set$Len,exact = F)$p.value>0.05) {Len_t=Len_t+1}
    if (wilcox.test(temp1$GC,test_set$GC,exact = F)$p.value>0.05) {GC_t=GC_t+1}
    if (wilcox.test(temp1$MFE,test_set$MFE,exact = F)$p.value>0.05) {MFE_t=MFE_t+1}
    if (wilcox.test(temp1$PhastCon30,test_set$PhastCon30,exact = F)$p.value>0.05) {PhastCon30_t=PhastCon30_t+1}
  }
  if (Len_t<=5 | GC_t <=5 | MFE_t <=5 | PhastCon30_t <=5 ) {
    #Length
    plot(density(temp1$Len),bty="n",xlim=c(0,350),lwd=1.5,
         ylim=c(0,0.015),main=NA,xlab="Intron length (bp)",col="red")
    lines(density(temp2$Len),xlim=c(0,350),lwd=1.5,col="black")
    legend(120,0.012,lty=c(1,1),lwd=1.5,col=c("red","black"),
           legend = c(paste0(RBP_list_sig[i]," (",Fun$RBP_by_FLEXI[Fun$Name==RBP_list_sig[i]],")"),"Others"),bty="n")
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


library(DESeq2)
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
pdf("../Pass1_IGV/Phast.pdf")
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
library("Rtsne")
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

library("umap")
custom.settings = umap.defaults
custom.settings$n_neighbors<-6
custom.settings$min_dist<-0.01
FLEXIumap<-umap(dat2[,1:8391],config=custom.settings)
FLEXIumap<-data.frame(FLEXIumap$layout)
FLEXIumap$label<-dat2$label
plot(FLEXIumap[,1:2],main="UMAP",col=col[FLEXIumap$label],
     xlab="W1",ylab="W2",pch=16)

library(zinbwave)
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
snoRNA<-gene_counts[gene_counts$Type=="snoRNA",c(2,10:14)]
snoRNA<-separate(snoRNA,col = "Name",into = "Name",sep = "-",remove = T)
snoRNA$Name<-gsub(snoRNA$Name,pattern = "[A-Z]$|P[1-9]$",replacement="")
snoRNA$Name[snoRNA$Name=="U3"]<-"SNORD3"
snoRNA$Name[snoRNA$Name=="U8"]<-"SNORD118"
snoRNA$Name[snoRNA$Name=="snoU13"]<-"SNORD13"
snoRNA$Name[snoRNA$Name=="snoU2_19"]<-"snoU2"
snoRNA<-aggregate(.~Name,data=snoRNA,FUN = sum)


snRNA<-gene_counts[gene_counts$Type=="snRNA",c(2,10:13)]
U7<-log10(colSums(snRNA[grep("U7\\b|U7[A-Z]",snRNA$Name),2:5])/mapped_reads[7:10])
U11<-log10(colSums(snRNA[grep("U11\\b|U11[A-Z]",snRNA$Name),2:5])/mapped_reads[7:10])
SNORD74<-log10(colSums(snoRNA[snoRNA$Name=="SNORD74",2:5])/mapped_reads[7:10])
SNORD78<-log10(colSums(snoRNA[snoRNA$Name=="SNORD78",2:5])/mapped_reads[7:10])
RN7SK<-log10(colSums(gene_counts[gene_counts$Type=="7SK",10:13])/mapped_reads[7:10])
RN7SL<-log10(colSums(gene_counts[gene_counts$Type=="7SL",10:13])/mapped_reads[7:10])
YRNA<-log10(colSums(gene_counts[gene_counts$Type=="YRNA",10:13])/mapped_reads[7:10])
VTRNA<-log10(colSums(gene_counts[gene_counts$Type=="VTRNA",10:13])/mapped_reads[7:10])
RMRP=log10(colSums(gene_counts[gene_counts$Name=="RMRP",10:13])/mapped_reads[7:10])
RPPH1=log10(colSums(gene_counts[gene_counts$Name=="RPPH1",10:13])/mapped_reads[7:10])
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





