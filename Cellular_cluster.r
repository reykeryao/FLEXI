rm(list=ls())
library(tidyverse)
#library(DESeq2)
library(ggalluvial)
library(ggplotify)
library(ComplexHeatmap)
set.seed(740714)
dat<-read.csv("Cellular_frac/raw_counts.csv",row.names=1)
'''
### pre-process

coldata<-data.frame("Cell"=rep(c("Hela","K562","MDA","MCF"),each=6),
                    "Frac"=rep(c("Total","Nuc","Cyto"),8),
                    "Sample"=colnames(dat)[5:28])
dds<-DESeqDataSetFromMatrix(countData = dat[,5:28],colData = coldata,design=~Frac+Cell)
dds<-estimateSizeFactors(dds)
size<-sizeFactors(dds)

Hela<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,5:10])>0,5:10],
                             colData = coldata[coldata$Cell=="Hela",],design=~Frac)
sizeFactors(Hela)<-size[1:6]
Hela<-DESeq(Hela,parallel = T)
saveRDS(Hela,"Cellular_frac/Hela.deseq")

K562<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,11:16])>0,11:16],
                             colData = coldata[coldata$Cell=="K562",],design=~Frac)
sizeFactors(K562)<-size[7:12]
K562<-DESeq(K562,parallel = T)
saveRDS(K562,"K562.deseq")

MDA<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,17:22])>0,17:22],
                             colData = coldata[coldata$Cell=="MDA",],design=~Frac)
sizeFactors(MDA)<-size[13:18]
MDA<-DESeq(MDA,parallel = T)
saveRDS(MDA,"Cellular_frac/MDA.deseq")

MCF<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,23:28])>0,23:28],
                             colData = coldata[coldata$Cell=="MCF",],design=~Frac)
sizeFactors(MCF)<-size[19:24]
MCF<-DESeq(MCF,parallel = T)
saveRDS(MCF,"Cellular_frac/MCF.deseq")

All<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,5:28])>0,5:28],
                             colData = coldata,design=~Frac)
sizeFactors(All)<-size
All<-DESeq(All,parallel = T)
saveRDS(All,"Cellular_frac/All.deseq")

Res<-c(list(data.frame(results(K562,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(Hela,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(MDA,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(MCF,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(All,contrast=c("Frac","Cyto","Nuc")))))
names(Res)<-c("K562","Hela","MDA","MCF","All")
saveRDS(Res,"Cellular_frac/Res_CytovsNuc")
'''
RBP<-read.delim("150_RBP_AGO_DICER.info")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP[RBP$col>1,46]<-4
RBP[RBP$col==1,46]<-3
RBP[RBP$col<1 & RBP$col>0,46]<-2
RBP[RBP$col==0,46]<-1

FLEXI_RBP<-read.table("all_short_intron_RBP.bed12",
                      col.names=c("FLEXI_chr","FLEXI_st","FLEXI_ed","ID","Score",
                                  "Strand","RBP_chr","RBP_st","RBP_ed","RBP",
                                  "Score1","Strand1","ENCODE_ID","Overlaps"))
Res<-readRDS("Cellular_frac/Res_CytovsNuc")
percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("gray80","red","orange","skyblue","orchid","blue1","goldenrod4")
L_col<-c("black","red","orange","skyblue","orchid","blue1","goldenrod4")
Clusters<-list(Cluster1=c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"),
               Cluster2=c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24","FXR2"),
               Cluster3=c("TIA1","TIAL1"),
               Cluster4=c("U2AF1","U2AF2","KHSRP"),
               Cluster5=c("AATF","DKC1","NOLC1","SMNDC1"),
               Cluster6=c("AGO","DICER"))

Core_splicesome<-c("PRPF8","SF3B4","AQR","EFTUD2","BUD13","PPIG")

### LFC no number plot
pdf("Figures/Fig10A_B.pdf",
    width=6,height=12)
par(mfrow=c(4,2),pty="s",mar=c(5,5,5,5))
## def 1, fraction specific FLEXI is defined as FC>1.5
LFC_cutoff<-log2(1.5)
for (i in 1:4){
  tmp<-Res[[i]]
  Up<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>LFC_cutoff]),c(4,10)]
  Up<-unique(Up)
  Down<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange<LFC_cutoff]),c(4,10)]
  Down<-unique(Down)
  
  RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
  colnames(RBP_plot)<-c("RBP.name","Up","Down")
  RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot<-merge(RBP_plot,RBP[,c(1,13,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,5,2:4)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP.name%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  Up_num<-length(unique(Up$ID))
  Down_num<-length(unique(Down$ID))
  RBP_plot$Clu<-0
  ### plot
  for (j in 1:dim(RBP_plot)[1]){
    R_tmp<-as.character(RBP_plot$RBP.name[j])
    RBP_plot$Padj[j]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
    Clu_tmp<-which(sapply(Clusters, function(x) (RBP_plot$RBP.name[j] %in% x)))
    if (length(Clu_tmp)>0){
      RBP_plot$Clu[j]<-Clu_tmp 
    }
  }
  RBP_plot$Clu<-RBP_plot$Clu+1
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  axis_max<-8
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  tmp<-RBP_plot[RBP_plot$Clu>1,]
  points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
  points(RBP_plot[RBP_plot$RNA.export==1,c(4,3)],pch=1,cex=1.5,lwd=1.5,col="black")
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  rect(0,0,2,2,lwd=0.75)
  segments(2,2,10,8,lty=2,lwd=0.75,xpd=T)
  segments(2,0,10,0,lty=2,lwd=0.75,xpd=T)
  par(xpd=F)
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 | RBP_plot$RNA.export==1)&
                 (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff | RBP_plot$RNA.export==1))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = RBP_plot$RBP.name[sig_cutoff])
  }
  if(i==1){
    legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other","RNA export"),xpd=TRUE,
           pch=c(rep(16,7),1),col = c(col[c(2:7,1)],"black"))
  }
  axis_max<-2
  RBP_plot<-RBP_plot[RBP_plot$Up<=axis_max & RBP_plot$Down<=axis_max,]
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  tmp<-RBP_plot[RBP_plot$Clu>1,]
  points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
  points(RBP_plot[RBP_plot$RNA.export==1,c(4,3)],pch=1,cex=1.5,lwd=1.5,col="black")
  axis(1,at=seq(0,axis_max,1),
       label=seq(0,axis_max,1))
  axis(2,at=seq(0,axis_max,1),las=2,
       label=seq(0,axis_max,1))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 | RBP_plot$RNA.export==1)
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = RBP_plot$RBP.name[sig_cutoff])
  }
}
dev.off()


#### short/long intron RBPs
long<-read.delim("Cellular_frac/longIF.RBP")
short<-read.delim("Cellular_frac/shortIF.RBP")

### can not apply cutoff
pdf("Figures/Fig10D.pdf",
    width=6,height=12)
par(mfrow=c(4,1),pty="s",mar=c(5,5,5,5))
for (i in 2:5){
  tmp<-long[,c(1,i,i+5)]
  RBP_plot<-merge(tmp,RBP[,c(1,13,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,5,2:4)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  RBP_plot$Clu<-0
  ### plot
  for (j in 1:dim(RBP_plot)[1]){
    R_tmp<-as.character(RBP_plot$RBP[j])
    RBP_plot$Padj[j]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
    Clu_tmp<-which(sapply(Clusters, function(x) (RBP_plot$RBP[j] %in% x)))
    if (length(Clu_tmp)>0){
      RBP_plot$Clu[j]<-Clu_tmp 
    }
  }
  RBP_plot$Clu<-RBP_plot$Clu+1
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  print(max(RBP_plot[,3:4]))
  print(RBP_plot$Padj[RBP_plot$RBP=="NONO"])
  axis_max<-6
  plot(RBP_plot[RBP_plot$Clu==1,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=strsplit(colnames(tmp)[3],"_")[[1]][1],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm (RBP sites, %)"),
       xlab=paste0("Nuclear (RBP sites, %)"))
  tmp<-RBP_plot[RBP_plot$Clu>1,]
  points(tmp[,c(3,4)],pch=16,cex=1.5,col=col[tmp$Clu])
  points(RBP_plot[RBP_plot$RNA.export==1,c(3,4)],pch=1,cex=1.5,lwd=1.5,col="black")
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  rect(0,0,2,2,lwd=0.75)
  segments(2,2,8,6,lty=2,lwd=0.75,xpd=T)
  segments(2,0,8,0,lty=2,lwd=0.75,xpd=T)
  par(xpd=F)
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Clu==2 | RBP_plot$Clu==3 | RBP_plot$RNA.export==1) & 
                 (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(3:4)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = RBP_plot$RBP[sig_cutoff])
  }
  if(i==2){
    legend(x=0,y=6,bty="n",legend = c(paste0("Cluster ",as.roman(1:6)),"Other","RNA export"),xpd=TRUE,
           pch=c(rep(16,7),1),col = c(col[c(2:7,1)],"black"))
  }
}
dev.off()

### short
pdf("Figures/Fig10C.pdf",
    width=6,height=12)
par(mfrow=c(4,1),pty="s",mar=c(5,5,5,5))
for (i in 2:5){
  tmp<-short[,c(1,i,i+5)]
  RBP_plot<-merge(tmp,RBP[,c(1,13,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,5,2:4)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  RBP_plot$Clu<-0
  ### plot
  for (j in 1:dim(RBP_plot)[1]){
    R_tmp<-as.character(RBP_plot$RBP[j])
    RBP_plot$Padj[j]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
    Clu_tmp<-which(sapply(Clusters, function(x) (RBP_plot$RBP[j] %in% x)))
    if (length(Clu_tmp)>0){
      RBP_plot$Clu[j]<-Clu_tmp 
    }
  }
  RBP_plot$Clu<-RBP_plot$Clu+1
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  print(max(RBP_plot[,3:4]))
  print(RBP_plot$Padj[RBP_plot$RBP=="DDX3X"])
  axis_max<-6
  plot(RBP_plot[RBP_plot$Clu==1,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=strsplit(colnames(tmp)[3],"_")[[1]][1],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm (RBP sites, %)"),
       xlab=paste0("Nuclear (RBP sites, %)"))
  tmp<-RBP_plot[RBP_plot$Clu>1,]
  points(tmp[,c(3,4)],pch=16,cex=1.5,col=col[tmp$Clu])
  points(RBP_plot[RBP_plot$RNA.export==1,c(3,4)],pch=1,cex=1.5,lwd=1.5,col="black")
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  rect(0,0,2,2,lwd=0.75)
  segments(2,2,8,6,lty=2,lwd=0.75,xpd=T)
  segments(2,0,8,0,lty=2,lwd=0.75,xpd=T)
  par(xpd=F)
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Clu==2 | RBP_plot$Clu==3 | RBP_plot$RNA.export==1) & 
                 (RBP_plot[,3]>=percent_cutoff | RBP_plot[,4]>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(3:4)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = RBP_plot$RBP[sig_cutoff])
  }
  if(i==2){
    legend(x=0,y=6,bty="n",legend = c(paste0("Cluster ",as.roman(1:6)),"Other","RNA export"),xpd=TRUE,
           pch=c(rep(16,7),1),col = c(col[c(2:7,1)],"black"))
  }
}
dev.off()


### Fig11
dat<-read.csv("Cellular_frac/raw_counts.csv",row.names = 1)
dat<-dat[dat$Type=="FLEXI",]
dat$K562_cyto<-dat$K562_Cyto_1+dat$K562_Cyto_2
dat$K562_nuc<-dat$K562_Nuc_1+dat$K562_Nuc_2
dat$Hela_cyto<-dat$HeLa_Cyto_1+dat$HeLa_Cyto_2
dat$Hela_nuc<-dat$HeLa_Nuc_1+dat$HeLa_Nuc_2
dat$MDA_cyto<-dat$MDA_Cyto_1+dat$MDA_Cyto_2
dat$MDA_nuc<-dat$MDA_Nuc_1+dat$MDA_Nuc_2
dat$MCF7_cyto<-dat$MCF7_Cyto_1+dat$MCF7_Cyto_2
dat$MCF7_nuc<-dat$MCF7_Nuc_1+dat$MCF7_Nuc_2
dat$cyto<-dat$K562_cyto+dat$Hela_cyto+dat$MDA_cyto+dat$MCF7_cyto
dat$nuc<-dat$K562_nuc+dat$Hela_nuc+dat$MDA_nuc+dat$MCF7_nuc
dat<-dat[,c(1,37,38,31,32,29,30,35,36,33,34)]

clu_RBP<-unlist(Clusters)
for (i in 1: length(clu_RBP)){
  name_tmp<-clu_RBP[i]
  Fid_tmp<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%name_tmp])
  dat_tmp<-dat[dat$Ensbl%in%Fid_tmp,]
  for (j in c(2,4,6,8,10)){
    nUm<-c(sum((dat_tmp[,j]>0) & (dat_tmp[,j+1]==0)),
           sum((dat_tmp[,j]==0) & (dat_tmp[,j+1]>0)),
           sum((dat_tmp[,j]>0) & (dat_tmp[,j+1]>0)))
    if (i==1 & j==2){
      res_tmp<-nUm
    } else {
      res_tmp<-rbind(res_tmp,nUm)
    }
  }
}
res_tmp<-data.frame(res_tmp)
colnames(res_tmp)<-c("Cytoplasm only","Nucleus only","Both")
res_tmp$RBP<-rep(clu_RBP,each=5)
res_tmp$Cell<-rep(c("Combined","HeLa S3","K-562","MCF7","MDA-MB-231"),28)
res_tmp1<-res_tmp
res_tmp<-gather(res_tmp, Type, FLEXIs, "Cytoplasm only":Both, factor_key=TRUE)
res_tmp$Type<-factor(res_tmp$Type,levels=c("Nucleus only","Both","Cytoplasm only"))
res_tmp$Cluster<-1
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[2]]]<-2
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[3]]]<-3
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[4]]]<-4
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[5]]]<-5
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[6]]]<-6
res_tmp$RBP<-factor(res_tmp$RBP,levels=clu_RBP)
lab_tmp<-res_tmp
res_tmp<-res_tmp1
res_tmp[,1:3]<-100*prop.table(as.matrix(res_tmp[,1:3]),1)
res_tmp<-gather(res_tmp, Type, FLEXIs, "Cytoplasm only":Both, factor_key=TRUE)
res_tmp$Type<-factor(res_tmp$Type,levels=c("Nucleus only","Both","Cytoplasm only"))
res_tmp$Cluster<-1
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[2]]]<-2
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[3]]]<-3
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[4]]]<-4
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[5]]]<-5
res_tmp$Cluster[res_tmp$RBP%in%Clusters[[6]]]<-6
res_tmp$RBP<-factor(res_tmp$RBP,levels=clu_RBP)
res_tmp$lab<-lab_tmp$FLEXIs
pdf("Figures/Fig11.pdf",height=4,width=30)
ggplot(res_tmp,aes(x = Cell, y = FLEXIs, fill = Type,label=lab)) +
  geom_bar(
    position = "stack",
    stat = "identity") +
  scale_fill_manual(values = c("skyblue","goldenrod2","red")) +
  geom_text(stat = "stratum", aes(stratum = Type),size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("FLEXI (%)") +
  facet_grid(~ Cluster+RBP)
dev.off()







### get a few numbers
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
dat$HeLa_Nuc<-dat$HeLa_Nuc_1+dat$HeLa_Nuc_2
dat$HeLa_Cyto<-dat$HeLa_Cyto_1+dat$HeLa_Cyto_2
dat$K562_Nuc<-dat$K562_Nuc_1+dat$K562_Nuc_2
dat$K562_Cyto<-dat$K562_Cyto_1+dat$K562_Cyto_2
dat$MDA_Nuc<-dat$MDA_Nuc_1+dat$MDA_Nuc_2
dat$MDA_Cyto<-dat$MDA_Cyto_1+dat$MDA_Cyto_2
dat$MCF7_Nuc<-dat$MCF7_Nuc_1+dat$MCF7_Nuc_2
dat$MCF7_Cyto<-dat$MCF7_Cyto_1+dat$MCF7_Cyto_2
### print some numbers
tmp<-dat[dat$Type=="FLEXI",c(1:4,29:36)]
for (i in c(7,5,9,11)){
  nUm<-c(sum((tmp[,i]>0 )& (tmp[,i+1]==0)),
         sum((tmp[,i]==0 )& (tmp[,i+1]>0)),
         sum((tmp[,i]>0 )& (tmp[,i+1]>0)))
  print(nUm)
}

for (i in c(7,5,9,11)){
  nUm<-c(sum(tmp[,i][(tmp[,i]>0) & (tmp[,i+1]==0)]),
         sum(tmp[,i+1][(tmp[,i]==0) & (tmp[,i+1]>0)]),
         sum(tmp[,i][(tmp[,i]>0) & (tmp[,i+1]>0)]),
         sum(tmp[,i+1][(tmp[,i]>0) & (tmp[,i+1]>0)]))
  print(nUm)
}
###


dat<-dat[,c(1:4,29:36)]
Sums<-colSums(dat[,5:12])/1e6
dat[,5:12]<-t(t(dat[,5:12])/Sums)
dat<-dat[dat$Type=="FLEXI",]
### get some num
tmp<-dat
for (i in c(7,5,9,11)){
  nUm<-c(sum(tmp[,i][(tmp[,i]>0) & (tmp[,i+1]==0)]),
         sum(tmp[,i+1][(tmp[,i]==0) & (tmp[,i+1]>0)]),
         sum(tmp[,i][(tmp[,i]>0) & (tmp[,i+1]>0)]),
         sum(tmp[,i+1][(tmp[,i]>0) & (tmp[,i+1]>0)]))
  print(nUm)
}
###
# number of FLEXIs
x<-apply(dat[,5:12],2,function(x){sum(x>0)})
#nuclear FLEXI #
range(x[seq(1,8,2)])
#[1]  226 1092
#cytoplasma FLEXI #
range(x[seq(0,8,2)])
#[1]  82 237

# total RPM of FLEXIs
x<-apply(dat[,5:12],2,function(x){sum(x)})
#nuclear FLEXI #
range(x[seq(1,8,2)])
#[1] 34.88362 39.43325
#cytoplasma FLEXI #
range(x[seq(0,8,2)])
#[1] 12.55991 24.82834

# average individual RPM of FLEXIs
x<-apply(dat[,5:12],2,function(x){mean(x[x>0])})
#nuclear FLEXI #
range(x[seq(1,8,2)])
#[1] 0.03611104 0.15435230
#cytoplasma FLEXI #
range(x[seq(0,8,2)])
#[1] 0.07397975 0.3027846

### overlaps
length(intersect(dat$Ensbl[dat$HeLa_Nuc>0],dat$Ensbl[dat$HeLa_Cyto>0]))/length(union(dat$Ensbl[dat$HeLa_Nuc>0],dat$Ensbl[dat$HeLa_Cyto>0]))
length(intersect(dat$Ensbl[dat$K562_Nuc>0],dat$Ensbl[dat$K562_Cyto>0]))/length(union(dat$Ensbl[dat$K562_Nuc>0],dat$Ensbl[dat$K562_Cyto>0]))
length(intersect(dat$Ensbl[dat$MDA_Nuc>0],dat$Ensbl[dat$MDA_Cyto>0]))/length(union(dat$Ensbl[dat$MDA_Nuc>0],dat$Ensbl[dat$MDA_Cyto>0]))
length(intersect(dat$Ensbl[dat$MCF7_Nuc>0],dat$Ensbl[dat$MCF7_Cyto>0]))/length(union(dat$Ensbl[dat$MCF7_Nuc>0],dat$Ensbl[dat$MCF7_Cyto>0]))
#[1] 0.06944444
#[1] 0.07104413
#[1] 0.04901118
#[1] 0.08603239

### FigS17
dat<-read.delim("all.FLEXI")
Fourcell<-dat[rowSums(dat[,73:76])>0,c(1,7,8)]
dat<-read.csv("Cellular_frac/raw_counts.csv",row.names = 1)
dat<-dat[dat$Type=="FLEXI",]
MDA<-dat[rowSums(dat[,17:22])>0,1:2]
MCF<-dat[rowSums(dat[,23:28])>0,1:2]
MDA<-separate(data = MDA,col="Gene_name",into=c("IID","GName","a","b","GID"),sep="_",remove = F,extra = "drop")
MDA<-MDA[,c(2,4,7)]
MCF<-separate(data = MCF,col="Gene_name",into=c("IID","GName","a","b","GID"),sep="_",remove = F,extra = "drop")
MCF<-MCF[,c(2,4,7)]

Onco<-read.delim("OncoGenes.table")
Onco$Onco<-1
Onco<-Onco[,c(2,8)]
colnames(Onco)[1]<-"GName"
TSG<-read.delim("TSGs.table")
TSG$TSG<-1
TSG<-TSG[,c(2,9)]
colnames(TSG)[1]<-"GName"
###
plot_upset <- function(m, i) {
  cs = comb_size(m)  
  ss<-set_size(m)
  ht <- UpSet(m, set_order=c("Other cell lines","MCF7","MDA-MB-231"),
              comb_col = c("black","skyblue","orchid")[comb_degree(m)],
              #comb_order=od,
              column_title=i,
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
### FLEIXs by RNA
set1<-Fourcell$ID
set2<-MCF$Gene_name
set3<-MDA$Gene_name
US_list <- list("FLEXI"=list ("Other cell lines"=set1,"MCF7"=set2,"MDA-MB-231"=set3))
### FLEIXs by gene
set1<-unique(Fourcell$GID)
set2<-unique(MCF$GID)
set3<-unique(MDA$GID)
US_list <- c(US_list,list("Gene"=list ("Other cell lines"=set1,"MCF7"=set2,"MDA-MB-231"=set3)))
### FLEIXs by gene, oncogene
set1<-unique(Fourcell$GName)
set1<-set1[set1%in%Onco$GName]
set2<-unique(MCF$GName)
set2<-set2[set2%in%Onco$GName]
set3<-unique(MDA$GName)
set3<-set3[set3%in%Onco$GName]
US_list <- c(US_list,list("Onco"=list ("Other cell lines"=set1,"MCF7"=set2,"MDA-MB-231"=set3)))
### FLEIXs by gene
set1<-unique(Fourcell$GName)
set1<-set1[set1%in%TSG$GName]
set2<-unique(MCF$GName)
set2<-set2[set2%in%TSG$GName]
set3<-unique(MDA$GName)
set3<-set3[set3%in%TSG$GName]
US_list <- c(US_list,list("TSG"=list ("Other cell lines"=set1,"MCF7"=set2,"MDA-MB-231"=set3)))

pdf("Figures/FigS17.pdf",height=8,width=8)
ht_list =vector(4, mode="list")
for (i in 1: length(US_list)){
  tmp<-US_list[[i]]
  m = make_comb_mat(tmp)
  dht <- as.ggplot( ~ plot_upset(m, names(US_list)[i]))
  ht_list[[i]] <- dht
}
patchwork::wrap_plots(ht_list, ncol=2,nrow=2)
dev.off()
#get gene names listed below the corresponding plot
###Oncogene
tmp<-US_list[[3]]
sort(t(t(setdiff(intersect(tmp$MCF7,tmp$`MDA-MB-231`),tmp$`Other cell lines`))))
sort(t(t(setdiff(setdiff(tmp$MCF7,tmp$`MDA-MB-231`),tmp$`Other cell lines`))))
sort(sort(t(t(setdiff(setdiff(tmp$`MDA-MB-231`,tmp$MCF7),tmp$`Other cell lines`)))))
# TSGs
tmp<-US_list[[4]]
sort(t(t(setdiff(intersect(tmp$MCF7,tmp$`MDA-MB-231`),tmp$`Other cell lines`))))
sort(t(t(setdiff(setdiff(tmp$MCF7,tmp$`MDA-MB-231`),tmp$`Other cell lines`))))
sort(sort(t(t(setdiff(setdiff(tmp$`MDA-MB-231`,tmp$MCF7),tmp$`Other cell lines`)))))

##FigS18
### scatter plots , CPM of frac
'''
dat<-read.csv("Cellular_frac/raw_counts.csv",row.names=1)
dat<-dat[dat$Type!="ERCC",]
dat<-dat[dat$Type!="18S_rRNA",]
dat<-dat[dat$Type!="28S_rRNA",]
dat<-dat[dat$Type!="5.8S_rRNA",]
dat<-dat[dat$Type!="5S_rRNA",]
dat<-dat[dat$Type!="FLEXI",]
dat<-dat[dat$Type!="Mt_rRNA",]
dat$Type2<-dat$Type
dat$Type2[dat$Type=="protein_coding"]<-"Protein coding"
dat$Type2[dat$Type=="TR"]<-"Protein coding"
dat$Type2[dat$Type=="IG"]<-"Protein coding"
dat$Type2[dat$Type=="pseudogene"]<-"lncRNA"
dat$Type2[dat$Type=="lincRNA"]<-"lncRNA"
dat$Type2[dat$Type=="Antisense"]<-"lncRNA"
dat$Type2[dat$Type=="7SK"]<-"sncRNA"
dat$Type2[dat$Type=="7SL"]<-"sncRNA"
dat$Type2[dat$Type=="miRNA"]<-"sncRNA"
dat$Type2[dat$Type=="miscRNA"]<-"sncRNA"
dat$Type2[dat$Type=="scaRNA"]<-"sncRNA"
dat$Type2[dat$Type=="snRNA"]<-"sncRNA"
dat$Type2[dat$Type=="snoRNA"]<-"sncRNA"
dat$Type2[dat$Type=="VTRNA"]<-"sncRNA"
dat$Type2[dat$Type=="YRNA"]<-"sncRNA"
dat$Type2[dat$Type=="DNA"]<-"Repeat element"
dat$Type2[dat$Type=="LINE"]<-"Repeat element"
dat$Type2[dat$Type=="Low_complexity"]<-"Repeat element"
dat$Type2[dat$Type=="LTR"]<-"Repeat element"
dat$Type2[dat$Type=="RC"]<-"Repeat element"
dat$Type2[dat$Type=="Retroposon"]<-"Repeat element"
dat$Type2[dat$Type=="Satellite"]<-"Repeat element"
dat$Type2[dat$Type=="Simple_repeat"]<-"Repeat element"
dat$Type2[dat$Type=="SINE"]<-"Repeat element"
dat$Type2[dat$Type=="Unknown"]<-"Repeat element"
dat$Type2[dat$Type=="Mt_tRNA"]<-"MT"
dat$Type2[dat$Type=="Mt_protein_coding"]<-"MT"

dat$HeLa_Nuc<-dat$HeLa_Nuc_1+dat$HeLa_Nuc_2
dat$HeLa_Cyto<-dat$HeLa_Cyto_1+dat$HeLa_Cyto_2
dat$K562_Nuc<-dat$K562_Nuc_1+dat$K562_Nuc_2
dat$K562_Cyto<-dat$K562_Cyto_1+dat$K562_Cyto_2
dat$MDA_Nuc<-dat$MDA_Nuc_1+dat$MDA_Nuc_2
dat$MDA_Cyto<-dat$MDA_Cyto_1+dat$MDA_Cyto_2
dat$MCF7_Nuc<-dat$MCF7_Nuc_1+dat$MCF7_Nuc_2
dat$MCF7_Cyto<-dat$MCF7_Cyto_1+dat$MCF7_Cyto_2
dat<-dat[,c(1:4,29:37)]
dat<-dat[rowSums(dat[,6:13])>0,]
Sums<-colSums(dat[,6:13])/1e6
dat[,6:13]<-t(t(dat[,6:13])/Sums)
write.table(dat,"Cell_frac_cleaned_combined.CPM",quote=F,sep="\t")
'''
dat<-read.delim("Cell_frac_cleaned_combined.CPM")
dat[,6:13]<-log2(dat[,6:13]+1)
dat$Type2<-factor(dat$Type2,levels = c("Protein coding","lncRNA","sncRNA","tRNA","MT","Repeat element"))
dat<-dat[order(dat$Type2),]

### scatter plots
pcol<-c("gray80","red","skyblue","orange","brown","blue1")
pcol<-col2hex(pcol)
pcol<-paste0(pcol,"80")
### 
pdf("Figures/FigS18.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)
plot(dat$HeLa_Nuc~dat$HeLa_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (all RNAs)")
abline(0,1)
legend("bottomright",legend = levels(dat$Type2),col=pcol,pch=16,bty="n")

plot(dat$K562_Nuc~dat$K562_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (all RNAs)")
abline(0,1)

plot(dat$MDA_Nuc~dat$MDA_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (all RNAs)")
abline(0,1)

plot(dat$MCF7_Nuc~dat$MCF7_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (all RNAs)")
abline(0,1)

### protein
tmp<-dat[dat$Type2=="Protein coding",]
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,,xlim=c(0,20),ylim=c(0,20),col=pcol[1],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (Protein coding)")
abline(0,1)

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[1],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (Protein coding)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (Protein coding)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (Protein coding)")
abline(0,1)

#lncRNA
tmp<-dat[dat$Type2=="lncRNA",]
tmp$Type<-factor(tmp$Type)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[tmp$Type],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (lncRNA)")
abline(0,1)
legend("bottomright",legend = levels(tmp$Type),col=pcol,pch=16,bty="n")

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (lncRNA)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (lncRNA)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (lncRNA)")
abline(0,1)

### sncRNA
tmp<-dat[dat$Type2=="sncRNA",]
tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (sncRNA)")
abline(0,1)
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (sncRNA)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (sncRNA)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (sncRNA)")
abline(0,1)

### repeat elements
tmp<-dat[dat$Type2=="Repeat element",]
tmp$Type[tmp$Type=="Low_complexity"]<-"Simple repeat"
tmp$Type[tmp$Type=="Satellite"]<-"Simple repeat"
tmp$Type[tmp$Type=="Simple_repeat"]<-"Simple repeat"
tmp$Type[tmp$Type=="RC"]<-"Other"
tmp$Type[tmp$Type=="Retroposon"]<-"Other"
tmp$Type[tmp$Type=="Unknown"]<-"Other"

tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (repeat element)")
abline(0,1)
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (repeat element)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (repeat element)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (repeat element)")
abline(0,1)

### tRNA
tmp<-dat[dat$Type2=="tRNA",]
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[4],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (tRNA)")
abline(0,1)

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (tRNA)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (tRNA)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (tRNA)")
abline(0,1)

## MT genes
tmp<-dat[dat$Type2=="MT",]
tmp$Type<-factor(tmp$Type)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[tmp$Type],
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="HeLa S3 (MT)")
abline(0,1)
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="K-562 (MT)")
abline(0,1)

plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MDA-MB-231 (MT)")
abline(0,1)

plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab="Cytoplasm RNA",ylab="Nuclear Rna",main="MCF7 (MT)")
abline(0,1)
dev.off()

# make png
png("Figures/FigS18_1.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
plot(dat$HeLa_Nuc~dat$HeLa_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)

plot(dat$K562_Nuc~dat$K562_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(dat$MDA_Nuc~dat$MDA_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(dat$MCF7_Nuc~dat$MCF7_Cyto,col=pcol[dat$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_2.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
### protein
tmp<-dat[dat$Type2=="Protein coding",]
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,,xlim=c(0,20),ylim=c(0,20),col=pcol[1],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[1],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type2],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_3.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
#lncRNA
tmp<-dat[dat$Type2=="lncRNA",]
tmp$Type<-factor(tmp$Type)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[tmp$Type],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_4.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
### sncRNA
tmp<-dat[dat$Type2=="sncRNA",]
tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)

plot(tmp$K562_Nuc~tmp$K562_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_5.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
### repeat elements
tmp<-dat[dat$Type2=="Repeat element",]
tmp$Type[tmp$Type=="Low_complexity"]<-"Simple repeat"
tmp$Type[tmp$Type=="Satellite"]<-"Simple repeat"
tmp$Type[tmp$Type=="Simple_repeat"]<-"Simple repeat"
tmp$Type[tmp$Type=="RC"]<-"Other"
tmp$Type[tmp$Type=="Retroposon"]<-"Other"
tmp$Type[tmp$Type=="Unknown"]<-"Other"

tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=ncol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_6.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))
### tRNA
tmp<-dat[dat$Type2=="tRNA",]
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[4],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[4],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

png("Figures/FigS18_7.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))

## MT genes
tmp<-dat[dat$Type2=="MT",]
tmp$Type<-factor(tmp$Type)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=pcol[tmp$Type],
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$K562_Nuc~tmp$K562_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MDA_Nuc~tmp$MDA_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
plot(tmp$MCF7_Nuc~tmp$MCF7_Cyto,col=pcol[tmp$Type],xlim=c(0,20),ylim=c(0,20),
     xlab=NA,ylab=NA,main=NA,axes=F)
abline(0,1)
axis(1,at=seq(0,20,5),labels = NA)
axis(2,at=seq(0,20,5),labels = NA)
dev.off()

### make legend
pdf("Figures/FigS18_legend.pdf",height=9,width=6)
pcol<-c("gray80","red","skyblue","orange","brown","blue1")
pcol<-col2hex(pcol)
par(mfrow=c(3,2),pty="s",pch=16)
plot.new()
legend("bottomright",legend = levels(dat$Type2),col=pcol,pch=16,bty="n")
#lncRNA
tmp<-dat[dat$Type2=="lncRNA",]
tmp$Type<-factor(tmp$Type)
plot.new()
legend("bottomright",legend = levels(tmp$Type),col=pcol,pch=16,bty="n")

### sncRNA
tmp<-dat[dat$Type2=="sncRNA",]
tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot.new()
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")

### repeat elements
tmp<-dat[dat$Type2=="Repeat element",]
tmp$Type[tmp$Type=="Low_complexity"]<-"Simple repeat"
tmp$Type[tmp$Type=="Satellite"]<-"Simple repeat"
tmp$Type[tmp$Type=="Simple_repeat"]<-"Simple repeat"
tmp$Type[tmp$Type=="RC"]<-"Other"
tmp$Type[tmp$Type=="Retroposon"]<-"Other"
tmp$Type[tmp$Type=="Unknown"]<-"Other"

tmp$Type<-factor(tmp$Type)
ncol<-brewer.pal(10,"Set3")[-2]
ncol<-col2hex(ncol)
plot.new()
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")
## MT genes
tmp<-dat[dat$Type2=="MT",]
tmp$Type<-factor(tmp$Type)
plot.new()
legend("bottomright",legend = levels(tmp$Type),col=pcol,pch=16,bty="n")
dev.off()


### Fig12A
### all 4 cell cytoplasm FLEXIs
### Res names: "K562" "Hela" "MDA"  "MCF"  "All" 

pdf("Figures/Fig12A_1.pdf",width=11,height=6)
par(mfrow=c(1,4))
for (k in 1:4){
  tmp<-Res[[k]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  # Cluster 1
  for (m in 1:6){
    # Cluster 1
    tmp1<-tmp[tmp$RBP%in%Clusters[[m]],]
    F_id<-unique(tmp1$ID)
    for (i in 1:length(F_id)){
      print (i)
      if (i==1){bar_tmp<-rep(0,5)}
      id<-F_id[i]
      FLEXI_tmp<-tmp1[tmp1$ID%in%id,c(1:4,10,8:9)]
      FLEXI_tmp<-FLEXI_tmp[FLEXI_tmp$RBP%in%Clusters[[m]],]
      RBP_tmp<-unique(FLEXI_tmp$RBP)
      for (j in 1:dim(FLEXI_tmp)[1]){
        FLEXI_tmp$range[j]<-list(c(FLEXI_tmp$RBP_st[j]:FLEXI_tmp$RBP_ed[j]))
      }
      if (length(RBP_tmp)==1){
        bar_tmp[1]<-bar_tmp[1]+1
      } else if (length(RBP_tmp)==2){
        range1<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[1]]))
        range2<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[2]]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          bar_tmp[2]<-bar_tmp[2]+1
        } else {
          bar_tmp[3]<-bar_tmp[3]+1
        }
      } else {
        for (j in 1:length(RBP_tmp)){
          if (j==1){
            range1<-list(unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[j]])))
            over<-0
          } else {
            range2<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[j]]))
            range1<-c(range1,list(range2))
            o_flag<-lapply(range1[1:(j-1)],function(x){intersect(unlist(x),range2)})
            over<-over+sum(lapply(o_flag,length)>0)
          }
        }
        if (over>0){
          bar_tmp[4]<-bar_tmp[4]+1
        } else {
          bar_tmp[5]<-bar_tmp[5]+1
        }
      }
    }
    if (m==1){
      bar_4cell<-bar_tmp
    } else {
      bar_4cell<-cbind(bar_4cell,bar_tmp)
    }
  }
  colnames(bar_4cell)<-names(Clusters)
  barplot(prop.table(bar_4cell,2),beside=F,col=col,main=names(Res)[k])
}
dev.off()

pdf("Figures/Fig12A_2.pdf",width=11,height=6)
par(mfrow=c(1,6))
for (m in 1:6){
  ### overalpping in cluster 1-6
  for (k in 1:4){
    tmp<-Res[[k]]
    tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
    tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
    # Cluster 1
    tmp1<-tmp[tmp$RBP%in%Clusters[[m]],]
    F_id<-unique(tmp1$ID)
    for (i in 1:length(F_id)){
      print (i)
      if (i==1){bar_tmp<-rep(0,5)}
      id<-F_id[i]
      FLEXI_tmp<-tmp1[tmp1$ID%in%id,c(1:4,10,8:9)]
      FLEXI_tmp<-FLEXI_tmp[FLEXI_tmp$RBP%in%Clusters[[m]],]
      RBP_tmp<-unique(FLEXI_tmp$RBP)
      for (j in 1:dim(FLEXI_tmp)[1]){
        FLEXI_tmp$range[j]<-list(c(FLEXI_tmp$RBP_st[j]:FLEXI_tmp$RBP_ed[j]))
      }
      if (length(RBP_tmp)==1){
        bar_tmp[1]<-bar_tmp[1]+1
      } else if (length(RBP_tmp)==2){
        range1<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[1]]))
        range2<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[2]]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          bar_tmp[2]<-bar_tmp[2]+1
        } else {
          bar_tmp[3]<-bar_tmp[3]+1
        }
      } else {
        for (j in 1:length(RBP_tmp)){
          if (j==1){
            range1<-list(unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[j]])))
            over<-0
          } else {
            range2<-unique(unlist(FLEXI_tmp$range[FLEXI_tmp$RBP==RBP_tmp[j]]))
            range1<-c(range1,list(range2))
            o_flag<-lapply(range1[1:(j-1)],function(x){intersect(unlist(x),range2)})
            over<-over+sum(lapply(o_flag,length)>0)
          }
        }
        if (over>0){
          bar_tmp[4]<-bar_tmp[4]+1
        } else {
          bar_tmp[5]<-bar_tmp[5]+1
        }
      }
    }
    if (k==1){
      bar_4cell<-bar_tmp
    } else {
      bar_4cell<-cbind(bar_4cell,bar_tmp)
    }
  }
  colnames(bar_4cell)<-names(Res)[1:4]
  barplot(prop.table(bar_4cell,2),beside=F,col=col,main=paste0("Cluster ",m),cex.names =0.5)
}
dev.off()
