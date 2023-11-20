rm(list=ls())
library(tidyverse)
#library(DESeq2)
set.seed(740714)
setwd("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/")
'''
### pre-process
dat<-read.csv("raw_counts.csv",row.names=1)
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
saveRDS(Hela,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/Hela.deseq")

K562<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,11:16])>0,11:16],
                             colData = coldata[coldata$Cell=="K562",],design=~Frac)
sizeFactors(K562)<-size[7:12]
K562<-DESeq(K562,parallel = T)
saveRDS(K562,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/K562.deseq")

MDA<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,17:22])>0,17:22],
                             colData = coldata[coldata$Cell=="MDA",],design=~Frac)
sizeFactors(MDA)<-size[13:18]
MDA<-DESeq(MDA,parallel = T)
saveRDS(MDA,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/MDA.deseq")

MCF<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,23:28])>0,23:28],
                             colData = coldata[coldata$Cell=="MCF",],design=~Frac)
sizeFactors(MCF)<-size[19:24]
MCF<-DESeq(MCF,parallel = T)
saveRDS(MCF,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/MCF.deseq")

All<-DESeqDataSetFromMatrix(countData = dat[dat$Type=="FLEXI" & rowSums(dat[,5:28])>0,5:28],
                             colData = coldata,design=~Frac)
sizeFactors(All)<-size
All<-DESeq(All,parallel = T)
saveRDS(All,"All.deseq")
K562<-readRDS("K562.deseq")
Hela<-readRDS("Hela.deseq")
MCF<-readRDS("MCF.deseq")
MDA<-readRDS("MDA.deseq")
All<-readRDS("All.deseq")

Res<-c(list(data.frame(results(K562,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(Hela,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(MDA,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(MCF,contrast=c("Frac","Cyto","Nuc")))),
       list(data.frame(results(All,contrast=c("Frac","Cyto","Nuc")))))
names(Res)<-c("K562","Hela","MDA","MCF","All")
saveRDS(Res,"Res_CytovsNuc")
'''
RBP<-read.delim("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Paper/FLEXI_git/150_RBP_AGO_DICER.info")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP[RBP$col>1,46]<-4
RBP[RBP$col==1,46]<-3
RBP[RBP$col<1 & RBP$col>0,46]<-2
RBP[RBP$col==0,46]<-1

FLEXI_RBP<-read.table("/stor/work/Lambowitz/ref/TGSEQ/RBP/FLEXI_RBP.bed12",
                      col.names=c("FLEXI_chr","FLEXI_st","FLEXI_ed","ID","Score",
                                  "Strand","RBP_chr","RBP_st","RBP_ed","RBP",
                                  "Score1","Strand1","ENCODE_ID","Overlaps"))
Res<-readRDS("Res_CytovsNuc")
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

LFC_cutoff<-log2(1.5)
Core_splicesome<-c("PRPF8","SF3B4","AQR","EFTUD2","BUD13","PPIG")

pdf("RBP_various_def_6coreRemoved.pdf",
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
  RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,4,2:3)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP.name%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  Up_num<-length(unique(Up$ID))
  Down_num<-length(unique(Down$ID))
  RBP_plot$Clu<-0
  RBP_plot$x<-lapply(RBP_plot$RBP.name,function(x){length(unique(Up$ID[Up$RBP%in%x]))})
  RBP_plot$y<-lapply(RBP_plot$RBP.name,function(x){length(unique(Down$ID[Down$RBP%in%x]))})
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
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
                 (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
  if(i==1){
    legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
           pch=16,col = col[c(2:7,1)])
  }
  axis_max<-2
  RBP_plot<-RBP_plot[RBP_plot$Up<=axis_max & RBP_plot$Down<=axis_max,]
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
  points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
  axis(1,at=seq(0,axis_max,1),
       label=seq(0,axis_max,1))
  axis(2,at=seq(0,axis_max,1),las=2,
       label=seq(0,axis_max,1))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
}
mtext("Cellular fraction-enriched FLEXI: FC>1.5 in the corresponding fraction than the other",side=3,
      cex=0.7,line = -1.5,outer=T)


### def 2, FC>1
LFC_cutoff<-0
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
  RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,4,2:3)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP.name%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  Up_num<-length(unique(Up$ID))
  Down_num<-length(unique(Down$ID))
  RBP_plot$Clu<-0
  RBP_plot$x<-lapply(RBP_plot$RBP.name,function(x){length(unique(Up$ID[Up$RBP%in%x]))})
  RBP_plot$y<-lapply(RBP_plot$RBP.name,function(x){length(unique(Down$ID[Down$RBP%in%x]))})
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
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
                 (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
  segments(0,5.1,5.1,5.1,lty = 2)
  segments(5.1,0,5.1,5.1,lty = 2)
  if(i==1){
    legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
           pch=16,col = col[c(2:7,1)])
  }
  axis_max<-2
  RBP_plot<-RBP_plot[RBP_plot$Up<=axis_max & RBP_plot$Down<=axis_max,]
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
  points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
  axis(1,at=seq(0,axis_max,1),
       label=seq(0,axis_max,1))
  axis(2,at=seq(0,axis_max,1),las=2,
       label=seq(0,axis_max,1))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
}
mtext("Cellular fraction-enriched FLEXI: normalized counts in the corresponding fraction > that of the other fraction",side=3,
      cex=0.7,line = -1.5,outer=T)

##def 3, detected in only one fraction
dat<-read.csv("raw_counts.csv",row.names = 1)
dat<-dat[dat$Type=="FLEXI",]
dat$K562_Cyto<-dat$K562_Cyto_1+dat$K562_Cyto_2
dat$K562_Nuc<-dat$K562_Nuc_1+dat$K562_Nuc_2
dat$HeLa_Cyto<-dat$HeLa_Cyto_1+dat$HeLa_Cyto_2
dat$HeLa_Nuc<-dat$HeLa_Nuc_1+dat$HeLa_Nuc_2
dat$MDA_Cyto<-dat$MDA_Cyto_1+dat$MDA_Cyto_2
dat$MDA_Nuc<-dat$MDA_Nuc_1+dat$MDA_Nuc_2
dat$MCF7_Cyto<-dat$MCF7_Cyto_1+dat$MCF7_Cyto_2
dat$MCF7_Nuc<-dat$MCF7_Nuc_1+dat$MCF7_Nuc_2
dat$Cyto<-dat$K562_Cyto+dat$HeLa_Cyto+dat$MDA_Cyto+dat$MCF7_Cyto
dat$Nuc<-dat$K562_Nuc+dat$HeLa_Nuc+dat$MDA_Nuc+dat$MCF7_Nuc
dat<-dat[,c(1,29:38)]


for (i in 1:4){
  tmp<-dat[,c(1,(i*2),(i*2+1))]
  Up<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]>0&tmp[,3]==0],c(4,10)]
  Up<-unique(Up)
  Down<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]==0&tmp[,3]>0],c(4,10)]
  Down<-unique(Down)
  RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
  colnames(RBP_plot)<-c("RBP.name","Up","Down")
  RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,4,2:3)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP.name%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  Up_num<-length(unique(Up$ID))
  Down_num<-length(unique(Down$ID))
  RBP_plot$Clu<-0
  RBP_plot$x<-lapply(RBP_plot$RBP.name,function(x){length(unique(Up$ID[Up$RBP%in%x]))})
  RBP_plot$y<-lapply(RBP_plot$RBP.name,function(x){length(unique(Down$ID[Down$RBP%in%x]))})
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
  axis_max<-10
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  tmp<-RBP_plot[RBP_plot$Clu>1,]
  points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
                 (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
  segments(0,5.1,5.1,5.1,lty = 2)
  segments(5.1,0,5.1,5.1,lty = 2)
  if(i==1){
    legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
           pch=16,col = col[c(2:7,1)])
  }
  axis_max<-2
  RBP_plot<-RBP_plot[RBP_plot$Up<=axis_max & RBP_plot$Down<=axis_max,]
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
  points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
  axis(1,at=seq(0,axis_max,1),
       label=seq(0,axis_max,1))
  axis(2,at=seq(0,axis_max,1),las=2,
       label=seq(0,axis_max,1))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
}
mtext("Cellular fraction-only FLEXI: only detected in the corresponding fraction",side=3,
      cex=0.7,line = -1.5,outer=T)

## def 4, detcted in the fraction
for (i in 1:4){
  tmp<-dat[,c(1,(i*2),(i*2+1))]
  Up<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]>0],c(4,10)]
  Up<-unique(Up)
  Down<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,3]>0],c(4,10)]
  Down<-unique(Down)
  RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
  colnames(RBP_plot)<-c("RBP.name","Up","Down")
  RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
  RBP_plot<-RBP_plot[,c(1,4,2:3)]
  RBP_plot$Padj<-1
  RBP_plot<-RBP_plot[!RBP_plot$RBP.name%in%Core_splicesome,]
  R_sum<-colSums(RBP_plot[,3:4])
  Up_num<-length(unique(Up$ID))
  Down_num<-length(unique(Down$ID))
  RBP_plot$Clu<-0
  RBP_plot$x<-lapply(RBP_plot$RBP.name,function(x){length(unique(Up$ID[Up$RBP%in%x]))})
  RBP_plot$y<-lapply(RBP_plot$RBP.name,function(x){length(unique(Down$ID[Down$RBP%in%x]))})
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
  axis(1,at=seq(0,axis_max,2),
       label=seq(0,axis_max,2))
  axis(2,at=seq(0,axis_max,2),las=2,
       label=seq(0,axis_max,2))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
                 (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
  segments(0,5.1,5.1,5.1,lty = 2)
  segments(5.1,0,5.1,5.1,lty = 2)
  if(i==1){
    legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
           pch=16,col = col[c(2:7,1)])
  }
  axis_max<-2
  RBP_plot<-RBP_plot[RBP_plot$Up<=axis_max & RBP_plot$Down<=axis_max,]
  plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
       main=names(Res)[i],
       pch=16,bty="n",col=col[1],
       ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
       xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
  RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
  points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
  axis(1,at=seq(0,axis_max,1),
       label=seq(0,axis_max,1))
  axis(2,at=seq(0,axis_max,1),las=2,
       label=seq(0,axis_max,1))
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
  abline(0,1,col="red")
  if(sum(sig_cutoff)>0){
    text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
         col=L_col[RBP_plot$Clu[sig_cutoff]],
         labels = paste0(RBP_plot$RBP.name[sig_cutoff]," (",RBP_plot$x[sig_cutoff],",",RBP_plot$y[sig_cutoff],")"))
  }
}
mtext("Cellular fraction FLEXI: all FLEXIs detected in the corresponding fraction",side=3,
      cex=0.7,line = -1.5,outer=T)

dev.off()

### LFC no number plot
pdf("RBP_6coreRemoved.pdf",
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


### barplot of CLuster 1-6
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
res_tmp$Cell<-rep(c("K-562","HeLa S3","MDA-MB-231","MCF7","Combined"),28)
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
library(ggalluvial)
pdf("frac_stack_bar.pdf",height=6,width=25)
ggplot(res_tmp,aes(x = Cell, y = FLEXIs, fill = Type,label=FLEXIs)) +
  geom_bar(
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  geom_text(stat = "stratum", aes(stratum = Type),size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(~ Cluster+RBP)
dev.off()


## 100% stacked
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
pdf("frac_stack_bar_percentage.pdf",height=4,width=30)
ggplot(res_tmp,aes(x = Cell, y = FLEXIs, fill = Type,label=lab)) +
  geom_bar(
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  geom_text(stat = "stratum", aes(stratum = Type),size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("FLEXI (%)") +
  facet_grid(~ Cluster+RBP)
dev.off()

### by clusters
### barplot of CLuster 1-6
for (i in 1: length(Clusters)){
  name_tmp<-Clusters[[i]]
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
res_tmp$Cluster<-rep(paste0("Cluster ",as.roman(1:6)),each=5)
res_tmp$Cell<-rep(c("K-562","HeLa S3","MDA-MB-231","MCF7","Combined"),6)
res_tmp1<-res_tmp
res_tmp<-gather(res_tmp, Type, FLEXIs, "Cytoplasm only":Both, factor_key=TRUE)
res_tmp$Type<-factor(res_tmp$Type,levels=c("Nucleus only","Both","Cytoplasm only"))
lab_tmp<-res_tmp
library(ggalluvial)
pdf("frac_stack_bar_byCluster.pdf",height=6,width=25)
ggplot(res_tmp,aes(x = Cell, y = FLEXIs, fill = Type,label=FLEXIs)) +
  geom_bar(
    position = "stack",
    stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  geom_text(stat = "stratum", aes(stratum = Type),size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(~ Cluster)
dev.off()

## 100% stacked
res_tmp<-res_tmp1
res_tmp[,1:3]<-100*prop.table(as.matrix(res_tmp[,1:3]),1)
res_tmp<-gather(res_tmp, Type, FLEXIs, "Cytoplasm only":Both, factor_key=TRUE)
res_tmp$Type<-factor(res_tmp$Type,levels=c("Nucleus only","Both","Cytoplasm only"))
res_tmp$lab<-lab_tmp$FLEXIs
pdf("frac_stack_bar_percentage_by_cluster.pdf",height=6,width=25)
ggplot(res_tmp,aes(x = Cell, y = FLEXIs, fill = Type,label=lab)) +
  geom_bar(
    position = "stack",
    stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  geom_text(stat = "stratum", aes(stratum = Type),size=2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("FLEXI (%)") +
  facet_grid(~ Cluster)
dev.off()



##
### for other Clu only
name_type<-c("K-562","HeLa S3","MDA-MB-231","MCF7","Combined")
pdf("RBP_distribution_in_FLEXIs_with_specific_RBPbinding_site.pdf")
for (i in 1: length(clu_RBP)){
  name_tmp<-clu_RBP[i]
  Fid_tmp<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%name_tmp])
  dat_tmp<-dat[dat$Ensbl%in%Fid_tmp,]
  for (j in c(2,4,6,8,10)){
    Cyto<-unique(dat_tmp$Ensbl[dat_tmp[,j]>0 & dat_tmp[,j+1]==0])
    ### solve 0 issue
    if (length(Cyto)==0) {
      Cyto<-data.frame(matrix(ncol = 2, nrow = 0))
      colnames(Cyto)<-c("RBP",paste0(name_type[j/2],":Cyto"))
    } else {
      Cyto<-FLEXI_RBP[FLEXI_RBP$ID%in%Cyto,c(4,10)]
      Cyto<-unique(Cyto)
      Cyto<-data.frame(table(Cyto$RBP))
      colnames(Cyto)<-c("RBP",paste0(name_type[j/2],":Cyto"))
      Cyto<-Cyto[Cyto$RBP!=name_tmp,]
    }
    Nuc<-unique(dat_tmp$Ensbl[dat_tmp[,j]==0 & dat_tmp[,j+1]>0])
    if (length(Nuc)==0){
      Nuc<-data.frame(matrix(ncol = 2, nrow = 0))
      colnames(Nuc)<-c("RBP",paste0(name_type[j/2],":Nuc"))
    } else {
      Nuc<-FLEXI_RBP[FLEXI_RBP$ID%in%Nuc,c(4,10)]
      Nuc<-unique(Nuc)
      Nuc<-data.frame(table(Nuc$RBP))
      colnames(Nuc)<-c("RBP",paste0(name_type[j/2],":Nuc"))
      Nuc<-Nuc[Nuc$RBP!=name_tmp,]
    }
    
    if (j==2){
      res_tmp<-merge(Cyto,Nuc,by="RBP",all=T)
    } else {
      res_tmp<-merge(res_tmp,Cyto,by="RBP",all=T)
      res_tmp<-merge(res_tmp,Nuc,by="RBP",all=T)
    }
  }
  res_tmp[is.na(res_tmp)]<-0
  res_tmp<-res_tmp[,c(1,2,4,6,8,10,3,5,7,9,11)]
  res_tmp[,2:11]<-100*prop.table(as.matrix(res_tmp[,2:11]),2)
  RBP_name<-res_tmp$RBP
  res_tmp<-gather(res_tmp, Type, Percent, "K-562:Cyto":"Combined:Nuc", factor_key=TRUE)
  res_tmp<-separate(res_tmp,Type,into=c("Cell","Fraction"),sep=":")
  gp<-ggplot(res_tmp) +
    geom_bar(aes(x = Cell, y = Percent, fill = RBP),
             position = "stack",
             stat = "identity") +
    ylab("RBP (%)") +
    ggtitle(paste0(name_tmp," bound FLEXIs")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(~ Fraction)
  print(gp)
}
dev.off()


### combined cell type
for (i in 1: length(clu_RBP)){
  name_tmp<-clu_RBP[i]
  Fid_tmp<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%name_tmp])
  dat_tmp<-dat[dat$Ensbl%in%Fid_tmp,]
  # order of column: cyto:nuc:both
  nUm<-c(sum((dat_tmp$Cyto>0) & (dat_tmp$Nuc==0)),
         sum((dat_tmp$Cyto==0) & (dat_tmp$Nuc>0)),
         sum((dat_tmp$Cyto>0) & (dat_tmp$Nuc>0)))
  if (i==1){
    res_tmp<-nUm
  } else {
    res_tmp<-rbind(res_tmp,nUm)
  }
}
res_tmp<-data.frame(res_tmp)
colnames(res_tmp)<-c("Cytoplasm only","Nucleus only","Both")
res_tmp$RBP<-clu_RBP
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
pdf("frac_stack_bar_combiend.pdf",height=6,width=8)
ggplot(res_tmp) +
  geom_bar(aes(x = RBP, y = FLEXIs, fill = Type),
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
dev.off()
## 100% stacked
res_tmp<-res_tmp1
Number<-rowSums(res_tmp[,1:3])
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
pdf("frac_stack_bar_percentage_combined.pdf",height=6,width=20)
ggplot(res_tmp) +
  geom_bar(aes(x = RBP, y = FLEXIs, fill = Type),
           position = "stack",
           stat = "identity") +
  scale_fill_manual(values = c("skyblue","tomato","pink")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
dev.off()


### for other Clu only
clu_RBP<-unlist(Clusters)
for (i in 1: length(clu_RBP)){
  name_tmp<-clu_RBP[i]
  Fid_tmp<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%name_tmp])
  dat_tmp<-dat[dat$Ensbl%in%Fid_tmp,]
  for (j in c(2,4,6,8)){
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
res_tmp$RBP<-rep(clu_RBP,each=4)
res_tmp$Cell<-rep(c("K-562","HeLa S3","MDA-MB-231","MCF7"),28)


### RBP using combiend various def (6 core removed)
pdf("RBP_various_def_6coreRemoved_combined_cell.pdf",
    width=6,height=12)
par(mfrow=c(4,2),pty="s",mar=c(5,5,5,5))
## def 1, fraction specific FLEXI is defined as FC>1.5
LFC_cutoff<-log2(1.5)
## Res use all FC>1.5
tmp<-Res[[5]]
Up<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>LFC_cutoff]),c(4,10)]
Up<-unique(Up)
Down<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange<LFC_cutoff]),c(4,10)]
Down<-unique(Down)
RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
colnames(RBP_plot)<-c("RBP.name","Up","Down")
RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
RBP_plot[is.na(RBP_plot)]<-0
RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
RBP_plot<-RBP_plot[,c(1,4,2:3)]
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
axis_max<-6
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
     main="Cell line combined: enriched as FC>1.5",
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
tmp<-RBP_plot[RBP_plot$Clu>1,]
points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
axis(1,at=c(seq(0,axis_max,2)),
     label=c(seq(0,axis_max,2)))
axis(2,at=c(seq(0,axis_max,2)),las=2,
     label=c(seq(0,axis_max,2)))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
               (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}
#segments(0,5.1,5.1,5.1,lty = 2)
#segments(5.1,0,5.1,5.1,lty = 2)
#if(i==1){
#  legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
#         pch=16,col = col[c(2:7,1)])
#}

RBP_plot<-RBP_plot[RBP_plot$Up<=2 & RBP_plot$Down<=2,]
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,2),ylim=c(0,2),cex=1.5,axes=F,
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
axis(1,at=c(0,1,2),label=c(0,1,2))
axis(2,at=c(0,1,2),label=c(0,1,2),las=2)
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}



### def 2, FC>1
LFC_cutoff<-0
tmp<-Res[[5]]
Up<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>LFC_cutoff]),c(4,10)]
Up<-unique(Up)
Down<-FLEXI_RBP[FLEXI_RBP$ID%in%sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange<LFC_cutoff]),c(4,10)]
Up<-unique(Up)
RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
colnames(RBP_plot)<-c("RBP.name","Up","Down")
RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
RBP_plot[is.na(RBP_plot)]<-0
RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
RBP_plot<-RBP_plot[,c(1,4,2:3)]
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
axis_max<-10
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
     main="Cell line combined: enriched as FC>1",
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
tmp<-RBP_plot[RBP_plot$Clu>1,]
points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
axis(1,at=c(seq(0,axis_max,2)),label=c(seq(0,axis_max,2)))
axis(2,las=2,at=c(seq(0,axis_max,2)),label=c(seq(0,axis_max,2)))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
               (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}
#segments(0,5.1,5.1,5.1,lty = 2)
#segments(5.1,0,5.1,5.1,lty = 2)
#if(i==1){
#  legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
#         pch=16,col = col[c(2:7,1)])
#}

RBP_plot<-RBP_plot[RBP_plot$Up<=2 & RBP_plot$Down<=2,]
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,2),ylim=c(0,2),cex=1.5,axes=F,
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
axis(1,at=c(0,1,2),label=c(0,1,2))
axis(2,las=2,at=c(0,1,2),label=c(0,1,2))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}

##def 3, detected in only one fraction
dat<-read.csv("raw_counts.csv",row.names = 1)
dat<-dat[dat$Type=="FLEXI",]
dat$Cyto<-rowSums(dat[,c(6,9,12,15,18,21,24,27)+1])
dat$Nuc<-rowSums(dat[,c(6,9,12,15,18,21,24,27)])

tmp<-dat[,c(1,29,30)]
Up<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]>0&tmp[,3]==0],c(4,10)]
Up<-unique(Up)
Down<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]==0&tmp[,3]>0],c(4,10)]
Down<-unique(Down)
RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
colnames(RBP_plot)<-c("RBP.name","Up","Down")
RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
RBP_plot[is.na(RBP_plot)]<-0
RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
RBP_plot<-RBP_plot[,c(1,4,2:3)]
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
     main="Cell line combined: fraction-only FLEXIS",
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
tmp<-RBP_plot[RBP_plot$Clu>1,]
points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
axis(1,at=seq(0,axis_max,2),label=seq(0,axis_max,2))
axis(2,las=2,at=seq(0,axis_max,2),label=seq(0,axis_max,2))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
               (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}
#segments(0,5.1,5.1,5.1,lty = 2)
#segments(5.1,0,5.1,5.1,lty = 2)
##if(i==1){
#  legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
#         pch=16,col = col[c(2:7,1)])
#}

RBP_plot<-RBP_plot[RBP_plot$Up<=2 & RBP_plot$Down<=2,]
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,2),ylim=c(0,2),cex=1.5,axes=F,
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
axis(1,at=c(0,1,2),label=c(0,1,2))
axis(2,at=c(0,1,2),las=2,label=c(0,1,2))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}

## def 4, detcted in the fraction
tmp<-dat[,c(1,29,30)]
Up<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,2]>0],c(4,10)]
Up<-unique(Up)
Down<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp$Ensbl[tmp[,3]>0],c(4,10)]
Down<-unique(Down)
RBP_plot<-merge(data.frame(table(Up$RBP)),data.frame(table(Down$RBP)),by=1,all=T)
colnames(RBP_plot)<-c("RBP.name","Up","Down")
RBP_plot<-RBP_plot[RBP_plot$RBP!=".",]
RBP_plot[is.na(RBP_plot)]<-0
RBP_plot<-merge(RBP_plot,RBP[,c(1,46)],by=1)
RBP_plot<-RBP_plot[,c(1,4,2:3)]
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
axis_max<-6
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,axes=F,
     main="Cell line combined: all FLEXIS",
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
tmp<-RBP_plot[RBP_plot$Clu>1,]
points(tmp[,c(4,3)],pch=16,cex=1.5,col=col[tmp$Clu])
axis(1,at=seq(0,axis_max,2),
     label=seq(0,axis_max,2))
axis(2,at=seq(0,axis_max,2),las=2,
     label=seq(0,axis_max,2))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-((RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3 )&
               (RBP_plot$Up>=percent_cutoff | RBP_plot$Down>=percent_cutoff))
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}
#segments(0,5.1,5.1,5.1,lty = 2)
#segments(5.1,0,5.1,5.1,lty = 2)
#if(i==1){
#  legend(x=6.5,y=6.5,bty="n",legend = c(paste0("Cluster ",1:6),"Other"),xpd=TRUE,
##         pch=16,col = col[c(2:7,1)])
#}

RBP_plot<-RBP_plot[RBP_plot$Up<=2 & RBP_plot$Down<=2,]
plot(RBP_plot[RBP_plot$Clu==1,c(4,3)],xlim=c(0,2),ylim=c(0,2),cex=1.5,axes=F,
     pch=16,bty="n",col=col[1],
     ylab=paste0("Cytoplasm FLEXIs (RBP sites, %)"),
     xlab=paste0("Nuclear FLEXIs (RBP sites, %)"))
RBP_plot<-RBP_plot[RBP_plot$Clu>1,]
points(RBP_plot[,c(4,3)],pch=16,cex=1.5,col=col[RBP_plot$Clu])
axis(1,at=c(0,1,2),
     label=c(0,1,2))
axis(2,at=c(0,1,2),las=2,
     label=c(0,1,2))
#sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
sig_cutoff<-(RBP_plot$Padj<=pvalue_cutoff |  RBP_plot$Clu==2 | RBP_plot$Clu==3)
abline(0,1,col="red")
if(sum(sig_cutoff)>0){
  text(RBP_plot[sig_cutoff,c(4:3)],cex=0.5,pos = 4,
       col=L_col[RBP_plot$Clu[sig_cutoff]],
       labels = RBP_plot$RBP.name[sig_cutoff])
}
dev.off()












### venn

library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
Core_set<-c(list(Clusters[[1]][c(1,2,8,3)]),
            list(Clusters[[1]][c(1,2,8,4)]),
                 list(Clusters[[1]][c(1,2,8,5)]),
                      list(Clusters[[1]][c(1,2,8,6)]),
                           list(Clusters[[1]][c(1,2,8,7)]))

for (i in 1: length(Core_set)){
  tmp<-c(list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][1]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][2]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][3]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][4]])))
  v1<-venn.diagram(tmp,category.names = Core_set[[i]],filename = NULL,
                   fill = col[2:5])
  if (i==1){
    v_list<-list(v1)
  } else {
    v_list<-c(v_list,list(v1))
  }
}

pdf("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/venn1.pdf")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 3, nc = 2,
                                           heights = unit(c(2,2),"null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(v_list[[1]])
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(v_list[[2]])
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(v_list[[3]])
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(v_list[[4]])
upViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(v_list[[5]])
upViewport()

upViewport()
dev.off()


Core_set<-c(list(Clusters[[1]][c(1,2,8,3,4)]),
            list(Clusters[[1]][c(1,2,8,3,5)]),
            list(Clusters[[1]][c(1,2,8,3,6)]),
            list(Clusters[[1]][c(1,2,8,3,7)]),
            list(Clusters[[1]][c(1,2,8,4,5)]),
            list(Clusters[[1]][c(1,2,8,4,6)]),
            list(Clusters[[1]][c(1,2,8,4,7)]),
            list(Clusters[[1]][c(1,2,8,5,6)]),
            list(Clusters[[1]][c(1,2,8,5,7)]),
            list(Clusters[[1]][c(1,2,8,6,7)]))

for (i in 1: length(Core_set)){
  tmp<-c(list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][1]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][2]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][3]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][4]])),
         list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Core_set[[i]][5]])))
  v1<-venn.diagram(tmp,category.names = Core_set[[i]],filename = NULL,
                   fill = col[2:6])
  if (i==1){
    v_list<-list(v1)
  } else {
    v_list<-c(v_list,list(v1))
  }
}

pdf("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/venn2.pdf",height=15,width=6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 5, nc = 2,
                                           heights = unit(c(2,2),"null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(v_list[[1]])
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(v_list[[2]])
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(v_list[[3]])
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(v_list[[4]])
upViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(v_list[[5]])
upViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
grid.draw(v_list[[6]])
upViewport()

pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
grid.draw(v_list[[7]])
upViewport()

pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
grid.draw(v_list[[8]])
upViewport()

pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 1))
grid.draw(v_list[[9]])
upViewport()

pushViewport(viewport(layout.pos.row = 5, layout.pos.col = 2))
grid.draw(v_list[[10]])
upViewport()

upViewport()
dev.off()



'''
  ###inlet bar
  adu<-0.5
  for (j in 1:length(Clusters)){
    clu<-Clusters[[j]]
    clu<-RBP_plot[RBP_plot$RBP.name%in%clu,]
    if (j==1){
      bar_in<-data.frame(c(sum(clu$Up>clu$Down & clu$Up>=adu),
                           length(clu$RBP.name)-sum(clu$Up>clu$Down & clu$Up>=adu)-sum(clu$Up<clu$Down & clu$Down>=adu),
                           sum(clu$Up<clu$Down & clu$Down>=adu)))
      colnames(bar_in)<-"Cluster1"
    } else {
      bar_in<-cbind(bar_in,data.frame(c(sum(clu$Up>clu$Down & clu$Up>=adu),
                                        length(clu$RBP.name)-sum(clu$Up>clu$Down & clu$Up>=adu)-sum(clu$Up<clu$Down & clu$Down>=adu),
                                        sum(clu$Up<clu$Down & clu$Down>=adu))))
      colnames(bar_in)[dim(bar_in)[2]]<-paste0("Cluster",j)
    }
  }
  bar_in[is.na(bar_in)]<-0
  bar_in$Type<-c("Cytoplasma","Not enriched","Nuclear")
  
  barplot(as.matrix(bar_in[,1:6]),ylim=c(0,10),col=c("skyblue","gray80","salmon"),names.arg =colnames(bar_in[1:6]),
          ylab="Number of RBPs",legend.text = bar_in$Type,args.legend = c(bty="n",title="Enriched RBP sites"))
'''

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

###check FLEXI sequence similarity in RBP sites
### first venn of FLEXIs in different fractions
### 
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
dat<-dat[dat$Type=="FLEXI",]
set1<-dat$Ensbl[(dat$HeLa_Nuc_1+dat$HeLa_Nuc_2)>0]
set2<-dat$Ensbl[(dat$HeLa_Cyto_1+dat$HeLa_Cyto_2)>0]
set3<-dat$Ensbl[(dat$K562_Nuc_1+dat$K562_Nuc_2)>0]
set4<-dat$Ensbl[(dat$K562_Cyto_1+dat$K562_Cyto_2)>0]
set5<-dat$Ensbl[(dat$MDA_Nuc_1+dat$MDA_Nuc_2)>0]
set6<-dat$Ensbl[(dat$MDA_Cyto_1+dat$MDA_Cyto_2)>0]
set7<-dat$Ensbl[(dat$MCF7_Nuc_1+dat$MCF7_Nuc_2)>0]
set8<-dat$Ensbl[(dat$MCF7_Cyto_1+dat$MCF7_Cyto_2)>0]
FLEXI_set<-list(set1,set2,set3,set4,set5,set6,set7,set8)

library(VennDiagram)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

v1<-venn.diagram(FLEXI_set[1:2],category.names = c("Hela Nuc","Hela Cyto"),
             fill = c("#E69F00", "#56B4E9"),
             cat.pos = c(0, 0),filename = NULL,
             cat.dist = c(0.05, 0.05))
v2<-venn.diagram(FLEXI_set[3:4],category.names = c("K562 Nuc","K562 Cyto"),
             fill = c("#E69F00", "#56B4E9"),
             cat.pos = c(0, 0),
             cat.dist = c(0.05, 0.05),filename = NULL)
v3<-venn.diagram(FLEXI_set[5:6],category.names = c("MDA Nuc","MDA Cyto"),
             fill = c("#E69F00", "#56B4E9"),
             cat.pos = c(0, 0),
             cat.dist = c(0.05, 0.05),filename = NULL)
v4<-venn.diagram(FLEXI_set[7:8],category.names = c("MCF7 Nuc","MCF7 Cyto"),
             fill = c("#E69F00", "#56B4E9"),
             cat.pos = c(0, 0),
             cat.dist = c(0.05, 0.05),filename = NULL)

pdf("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/venn.pdf")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2,
                                           heights = unit(c(2,2),"null"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(v1)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(v2)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(v3)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(v4)
upViewport()

upViewport()
dev.off()


### combined cell line venn
set1<-dat$Ensbl[rowSums(dat[,c(6,9,12,15,18,21,24,27)])>0]
set2<-dat$Ensbl[rowSums(dat[,c(6,9,12,15,18,21,24,27)+1])>0]
all_set<-list(set1,set2)
pdf("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/all_cell_venn.pdf")
v1<-venn.diagram(all_set,category.names = c("Hela Nuc","Hela Cyto"),
                 fill = c("#E69F00", "#56B4E9"),
                 cat.pos = c(0, 0),filename = NULL,
                 cat.dist = c(0.05, 0.05))
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(v1)
upViewport()
upViewport()
dev.off()
### volcano
library(DESeq2)
'''
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
coldata<-data.frame("Cell"=rep(c("Hela","K562","MDA","MCF"),each=6),
                    "Frac"=rep(c("Total","Nuc","Cyto"),8),
                    "Sample"=colnames(dat)[5:28])
coldata$Type<-paste(coldata$Cell,coldata$Frac,sep="_")
dds<-DESeqDataSetFromMatrix(countData = dat[,5:28],colData = coldata,design=~Type)
dds<-DESeq(dds,parallel = T)
saveRDS(dds,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/all_genes.deseq")

Res<-c(list(data.frame(results(dds,contrast=c("Type","K562_Cyto","K562_Nuc")))),
       list(data.frame(results(dds,contrast=c("Type","Hela_Cyto","Hela_Nuc")))),
       list(data.frame(results(dds,contrast=c("Type","MDA_Cyto","MDA_Nuc")))),
       list(data.frame(results(dds,contrast=c("Type","MCF_Cyto","MCF_Nuc")))))
names(Res)<-c("K562","Hela","MDA","MCF")
saveRDS(Res,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/Res_CytovsNuc_allgenes")
'''
library(gplots)
plot_volcano<-function(dat,x_lim,y_lim,xlab_seq,
                       ylab_seq,title,leg_posx=-1,leg_posy=-1){
  plot(dat$log2FoldChange,dat$log10p,xlim=x_lim,ylim=y_lim,
       main=title, xlab=bquote(log[2](FC)),xaxt="n",bty="n", yaxt="n",
       ylab=bquote("-"*log[10]*"(adjusted p-value)"),col=pcol[dat$Col],
       pch=20, cex=1)
  axis(side = 1,at = xlab_seq,labels = xlab_seq)
  axis(side = 2,labels = ylab_seq,at = ylab_seq,las=1)
  abline(v=0, col="black", lty=3, lwd=1.0)
  if ((leg_posx != -1)&&(leg_posy!= -1)){
    par(xpd=F)
    legend(leg_posx,leg_posy,pch=21,bty="n",pt.bg=pcol,cex=1,
           pt.lwd = 0.5,
           legend =c("Non-significant","Nuclear locolized RNA","Cytoplasm locolized RNA"))
  }
}

dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
Res<-readRDS("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/Res_CytovsNuc_allgenes")
pcol<-c("gray80","skyblue","salmon")
pcol<-paste0(col2hex(pcol),"80")
for (i in 1:4){
  T_name<-names(Res)[i]
  tmp<-Res[[i]]
  tmp<-cbind(tmp,dat$Type)
  tmp<-tmp[complete.cases(tmp),]
  tmp$log10p<-log10(1/tmp$padj)
  tmp$Col<-1
  tmp$Col[tmp$`dat$Type`=="snRNA"]<-2
  tmp$Col[tmp$`dat$Type`=="snoRNA"]<-2
  tmp$Col[tmp$`dat$Type`=="scaRNA"]<-2
  tmp$Col[tmp$`dat$Type`=="tRNA"]<-3
  tmp$Col[tmp$`dat$Type`=="Mt_protein_coding"]<-3
  tmp$Col[tmp$`dat$Type`=="Mt_rRNA"]<-3
  tmp$Col[tmp$`dat$Type`=="Mt_tRNA"]<-3
  x_ran<-range(tmp$log2FoldChange)%/%2
  # left end +1
  x_ran[2]<-x_ran[2]+1
  x_ran<-x_ran*2
  y_ran<-max(tmp$log10p)%/%30+1
  y_ran<-c(0,y_ran*30)
  plot_volcano(tmp,x_ran,y_ran,seq(x_ran[1],x_ran[2],2),seq(y_ran[1],y_ran[2],3),
               "Protein coding",x_ran[1],y_ran[2])
}

### scatter plots
dat$K562_Cyto<-dat$K562_Cyto_1+dat$K562_Cyto_2
dat$K562_Nuc<-dat$K562_Nuc_1+dat$K562_Nuc_2
dat$HeLa_Cyto<-dat$HeLa_Cyto_1+dat$HeLa_Cyto_2
dat$HeLa_Nuc<-dat$HeLa_Nuc_1+dat$HeLa_Nuc_2
dat$MDA_Cyto<-dat$MDA_Cyto_1+dat$MDA_Cyto_2
dat$MDA_Nuc<-dat$MDA_Nuc_1+dat$MDA_Nuc_2
dat$MCF_Cyto<-dat$MCF7_Cyto_1+dat$MCF7_Cyto_2
dat$MCF_Nuc<-dat$MCF7_Nuc_1+dat$MCF7_Nuc_2
dat<-dat[,c(1:4,29:36)]
dat<-dat[dat$Type!="ERCC",]
dat<-dat[dat$Type!="28S_rRNA",]
dat<-dat[dat$Type!="18S_rRNA",]
dat<-dat[dat$Type!="5.8S_rRNA",]
dat<-dat[dat$Type!="5S_rRNA",]
dat<-dat[dat$Type!="Mt_rRNA",]
dat$Col<-1
dat$Col[dat$Type=="snRNA"]<-2
dat$Col[dat$Type=="snoRNA"]<-2
dat$Col[dat$Type=="scaRNA"]<-2
dat$Col[dat$Type=="tRNA"]<-3
dat$Col[dat$Type=="Mt_protein_coding"]<-3
dat$Col[dat$Type=="Mt_tRNA"]<-3
dat[,5:12]<-1e6*prop.table(as.matrix(dat[,5:12]),2)
dat[,5:12]<-log2(dat[,5:12]+1)
T_name<-c("K-562","HeLa S3","MDA-MB-231","MCF7")
pdf("Frac_scatter.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s")
for (i in c(5,7,9,11)){
  tmp<-dat[,c(i,i+1,13)]
  plot(tmp[tmp$Col==1,1:2],pch=19,col=pcol[1],xlim=c(0,20),ylim=c(0,20),
       main=T_name[(i-3)/2],xlab="Cytoplasm RNA","Nuclear RNA")
  points(tmp[tmp$Col==2,1:2],pch=19,col=pcol[2])
  points(tmp[tmp$Col==3,1:2],pch=19,col=pcol[3])
  abline(0,1)
  if(i==5){
    legend(bty="n","topleft",pch=19,col=pcol,
          legend = c("Other","Nuclear locolized RNA","Cytoplasm locolized RNA"))
  }
}
dev.off()



### cluster I FLEXIs
FLEXI_RBP<-read.table("/stor/work/Lambowitz/ref/TGSEQ/RBP/FLEXI_RBP.bed12")
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
dat<-dat[dat$Type=="FLEXI",]
dat<-dat[rowSums(dat[,5:28])>0,]
Clusters<-list(Cluster1=c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"),
               Cluster2=c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24","FXR2"),
               Cluster3=c("TIA1","TIAL1"),
               Cluster4=c("U2AF1","U2AF2","KHSRP"),
               Cluster5=c("AATF","DKC1","NOLC1","SMNDC1"),
               Cluster6=c("AGO","DICER"))

Clus_list<-c(list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]]])),
             list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]]])),
             list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[3]]])),
             list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[4]]])),
             list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[5]]])),
             list(unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[6]]])))
names(Clus_list)<-paste0("Cluster",1:6)
FLEXI<-read.delim("/stor/work/Lambowitz/ref/TGSEQ/full_length_intron/GRCh38.93.intron_deduped.tsv")
library(ComplexHeatmap)
pdf("Cluster_upset.pdf",height=4,width=11)
set1<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][1]])
set2<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][2]])
set3<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][3]])
set4<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][4]])
set5<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][5]])
set6<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][6]])
set7<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][7]])
set8<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[1]][8]])

set <- list (set1,set2,set3,set4,set5,set6,set7,set8)
names(set)<-Clusters[[1]]
m = make_comb_mat(set)
m<-m[comb_degree(m)>2]
ss<-set_size(m)
cs=comb_size(m)
print(sum(cs))
od<-1:55
UpSet(m,
      comb_col = c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2],
      column_title="FLEXI RNAs",
      comb_order=od,
      top_annotation = HeatmapAnnotation(
        "Counts" = anno_barplot(cs, 
                                ylim = c(0, max(cs)*1.1),
                                border = F,
                                gp = gpar(border =NA,lty=0,
                                          fill =c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2]), 
                                height = unit(5, "cm")), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {grid.text(cs[od], x = 1:55, y = unit(cs[od], "native") + unit(2, "pt"), 
                                         default.units = "native", just = c("left", "bottom"), 
                                         gp = gpar(fontsize = 8,
                                                   col=c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2]),
                                         rot = 45)})

### cluster 2
set1<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][1]])
set2<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][2]])
set3<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][3]])
set4<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][4]])
set5<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][5]])
set6<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][6]])
set7<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][7]])
set8<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][8]])
set9<-unique(FLEXI_RBP$ID[FLEXI_RBP$RBP%in%Clusters[[2]][9]])

set <- list (set1,set2,set3,set4,set5,set6,set7,set8,set9)
names(set)<-Clusters[[2]]
m = make_comb_mat(set)
m<-m[comb_degree(m)>2]
print(sum(comb_size(m)))
ss<-set_size(m)
cs=comb_size(m)
od<-1:126
UpSet(m,
      comb_col = c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2],
      column_title="FLEXI RNAs",
      comb_order=od,
      top_annotation = HeatmapAnnotation(
        "Counts" = anno_barplot(cs, 
                                ylim = c(0, max(cs)*1.1),
                                border = F,
                                gp = gpar(border =NA,lty=0,
                                          fill =c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2]), 
                                height = unit(5, "cm")), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {grid.text(cs[od], x = 1:126, y = unit(cs[od], "native") + unit(2, "pt"), 
                                         default.units = "native", just = c("left", "bottom"), 
                                         gp = gpar(fontsize = 8,
                                                   col=c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2]),
                                         rot = 45)})
dev.off()

### individual cell line upset
## Venn

### get some numbers
old_FLEXI<-read.delim("../Paper/FLEXI_git/all.FLEXI")
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
dat<-dat[dat$Type=="FLEXI",]
dat<-dat[rowSums(dat[,5:28])>0,]
all_FLEXI_ID<-union(old_FLEXI$ID,dat$Ensbl)
tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%all_FLEXI_ID,]


pdf("Cluster_upset_deected_FLEXI.pdf",height=4,width=11)
set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][1]])
set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][2]])
set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][3]])
set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][4]])
set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][5]])
set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][6]])
set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][7]])
set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][8]])

set <- list (set1,set2,set3,set4,set5,set6,set7,set8)
names(set)<-Clusters[[1]]
m = make_comb_mat(set)
print(sum((comb_size(m))))
m<-m[comb_degree(m)>1]
print(sum(comb_size(m)))
m<-m[comb_degree(m)>2]
print(sum(comb_size(m)))

ss<-set_size(m)
cs=comb_size(m)

od<-1:34
UpSet(m,
      comb_col = c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2],
      column_title="FLEXI RNAs",
      comb_order=od,
      top_annotation = HeatmapAnnotation(
        "Counts" = anno_barplot(cs, 
                                ylim = c(0, max(cs)*1.1),
                                border = F,
                                gp = gpar(border =NA,lty=0,
                                          fill =c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2]), 
                                height = unit(5, "cm")), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {grid.text(cs[od], x = od, y = unit(cs[od], "native") + unit(2, "pt"), 
                                         default.units = "native", just = c("left", "bottom"), 
                                         gp = gpar(fontsize = 8,
                                                   col=c("black","tomato","royalblue1","goldenrod","orchid")[comb_degree(m)-2]),
                                         rot = 45)})

### cluster 2
set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][1]])
set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][2]])
set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][3]])
set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][4]])
set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][5]])
set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][6]])
set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][7]])
set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][8]])
set9<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][9]])
set <- list (set1,set2,set3,set4,set5,set6,set7,set8,set9)
names(set)<-Clusters[[2]]
m = make_comb_mat(set)
print(sum((comb_size(m))))
m<-m[comb_degree(m)>1]
print(sum(comb_size(m)))
m<-m[comb_degree(m)>2]
print(sum(comb_size(m)))
ss<-set_size(m)
cs=comb_size(m)
od<-1:82
UpSet(m,
      comb_col = c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2],
      column_title="FLEXI RNAs",
      comb_order=od,
      top_annotation = HeatmapAnnotation(
        "Counts" = anno_barplot(cs, 
                                ylim = c(0, max(cs)*1.1),
                                border = F,
                                gp = gpar(border =NA,lty=0,
                                          fill =c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2]), 
                                height = unit(5, "cm")), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      right_annotation = rowAnnotation("Size"=anno_text(formatC(ss,big.mark = ","))))
decorate_annotation("Counts", {grid.text(cs[od], x = od, y = unit(cs[od], "native") + unit(2, "pt"), 
                                         default.units = "native", just = c("left", "bottom"), 
                                         gp = gpar(fontsize = 8,
                                                   col=c("black","tomato","royalblue1","goldenrod","orchid","skyblue")[comb_degree(m)-2]),
                                         rot = 45)})
dev.off()

### combination counts
### cluster1
set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][1]])
set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][2]])
set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][3]])
set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][4]])
set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][5]])
set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][6]])
set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][7]])
set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][8]])
set <- list (set1,set2,set3,set4,set5,set6,set7,set8)
names(set)<-Clusters[[1]]
for (i in 1:(length(set)-1)){
  for (j in (i+1):length(set)){
    name1<-names(set)[i]
    name2<-names(set)[j]
    unit_num<-length(intersect(set[[i]],set[[j]]))
    if (i==1){
      Clu1_comb<-data.frame(name1,name2,unit_num)
    } else (
      Clu1_comb<-rbind(Clu1_comb,data.frame(name1,name2,unit_num))
    )
  }
}

### cluster2
set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][1]])
set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][2]])
set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][3]])
set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][4]])
set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][5]])
set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][6]])
set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][7]])
set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][8]])
set9<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][9]])
set <- list (set1,set2,set3,set4,set5,set6,set7,set8,set9)
names(set)<-Clusters[[2]]
for (i in 1:(length(set)-1)){
  for (j in (i+1):length(set)){
    name1<-names(set)[i]
    name2<-names(set)[j]
    unit_num<-length(intersect(set[[i]],set[[j]]))
    if (i==1){
      Clu2_comb<-data.frame(name1,name2,unit_num)
    } else (
      Clu2_comb<-rbind(Clu2_comb,data.frame(name1,name2,unit_num))
    )
  }
}

### >=2, top popular combinations

'''
library(stringdist)
# calculate Levenshtein Distance
## make a master distance matrix
tmp<-FLEXI[FLEXI$ID%in%dat$Ensbl,c("ID","Seq","Len")]
rownames(tmp)<-tmp$ID

for (i in 1:dim(tmp)[1]){
  print (i)
  tmp$LV<-sapply(tmp$Seq,
                     function(x){stringdist(x, tmp$Seq[i], method = "lv")/max(nchar(x),nchar(tmp$Seq[i]))})
  colnames(tmp)[dim(tmp)[2]]<-rownames(tmp)[i]
}
tmp<-tmp[,c(-1:-3)]
saveRDS(tmp,"/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/FLEXI_seq-LVdis.matrix")
'''
FLEXI_dis<-readRDS("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/FLEXI_seq-LVdis.matrix")
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)
pdf("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/FLEXI_distance_by_cluster.pdf",height=8,width=8,onefile = T)
for (i in 1:6){
  pick_by_ID<-rownames(FLEXI_dis)%in%Clus_list[[i]]
  tmp<-FLEXI_dis[pick_by_ID,pick_by_ID]
  heatmap.2(as.matrix(tmp),labRow =NA,dendrogram = "none",
            scale="none",breaks=seq(0,1,0.01),margins = c(1, 1),
            main=paste0("Cluster ",i),
            density.info = "none",trace="none",symm=F,key=T,
            symbreaks=F,keysize=1,symkey=F,labCol =NA,
            col =heatPalette)
}
dev.off()

### RBP sites overlap
RBP_over<-read.table("/stor/work/Lambowitz/ref/TGSEQ/RBP/RBP_152.overlap10percent.info",
                     col.names=c("ID","Overlap"))
RBP_name<-sort(unique(RBP_over$ID))
### make tables
for (i in 1:152){
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
#RBP_matrix<-apply(RBP_matrix,2,function(x){x/max(x)})
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)
L_col<-c("black","red","orange","skyblue","orchid","blue1","goldenrod4")
Core_splicesome<-c("PRPF8","SF3B4","AQR","EFTUD2","BUD13","PPIG")
Clusters<-c(Clusters,list("Cluster7"=c("G3BP1", "IGF2BP2",  "YBX3")))
pdf("RBP_site_overlap.pdf",height=11,width=8)
for (i in 1:7){
  tmp<-data.frame(t(RBP_matrix[rownames(RBP_matrix)%in%Clusters[[i]],]))
  tmp<-tmp[!rownames(tmp)%in%Core_splicesome,]
  tmp$Col<-1
  tmp$Col[rownames(tmp)%in%Clusters[[1]]]<-2
  tmp$Col[rownames(tmp)%in%Clusters[[2]]]<-3
  tmp$Col[rownames(tmp)%in%Clusters[[3]]]<-4
  tmp$Col[rownames(tmp)%in%Clusters[[4]]]<-5
  tmp$Col[rownames(tmp)%in%Clusters[[5]]]<-6
  tmp$Col[rownames(tmp)%in%Clusters[[6]]]<-7
  if (i<7){
    heatmap.2(as.matrix(tmp[,1:(dim(tmp)[2]-1)]),labRow =rownames(tmp),
              dendrogram = "none",
              scale="none",breaks=seq(0,5,0.05),margins = c(10, 10),
              main=paste0("Cluster ",i),colRow=L_col[tmp[,dim(tmp)[2]]],
              density.info = "none",trace="none",symm=F,key=F,
              key.title = "log10(binding sites)",key.xlab = NA,key.ylab=NA,
              keysize = 0.5,colCol = rep(L_col[i+1],length(Clusters[[i]])),
              symbreaks=F,labCol =colnames(tmp),cexRow=0.5,
              col =heatPalette)
  } else {
    heatmap.2(as.matrix(tmp[,1:(dim(tmp)[2]-1)]),labRow =rownames(tmp),
              dendrogram = "none",
              scale="none",breaks=seq(0,5,0.05),margins = c(10, 10),
              main=paste0("Cluster ",i),colRow=L_col[tmp[,dim(tmp)[2]]],
              density.info = "none",trace="none",symm=F,key=T,
              key.title = "log10(binding sites)",key.xlab = NA,key.ylab=NA,
              keysize = 1,colCol = rep(L_col[i+1],length(Clusters[[i]])),
              symbreaks=F,labCol =colnames(tmp),cexRow=0.5,
              col =heatPalette)
  }

}
dev.off()

pdf("RBP_site_overlap_allRBPs.pdf",height=11,width=11)
tmp<-data.frame(t(RBP_matrix))
tmp<-tmp[!rownames(tmp)%in%Core_splicesome,]
tmp<-tmp[,!colnames(tmp)%in%Core_splicesome]
tmp$Col<-1
tmp$Col[rownames(tmp)%in%Clusters[[1]]]<-2
tmp$Col[rownames(tmp)%in%Clusters[[2]]]<-3
tmp$Col[rownames(tmp)%in%Clusters[[3]]]<-4
tmp$Col[rownames(tmp)%in%Clusters[[4]]]<-5
tmp$Col[rownames(tmp)%in%Clusters[[5]]]<-6
tmp$Col[rownames(tmp)%in%Clusters[[6]]]<-7
tmp<-rbind(tmp,rep(1,147))
rownames(tmp)[147]<-"Col"
tmp["Col",colnames(tmp)%in%Clusters[[1]]]<-2
tmp["Col",colnames(tmp)%in%Clusters[[2]]]<-3
tmp["Col",colnames(tmp)%in%Clusters[[3]]]<-4
tmp["Col",colnames(tmp)%in%Clusters[[4]]]<-5
tmp["Col",colnames(tmp)%in%Clusters[[5]]]<-6
tmp["Col",colnames(tmp)%in%Clusters[[6]]]<-7
heatmap.2(as.matrix(tmp[1:(dim(tmp)[1]-1),1:(dim(tmp)[2]-1)]),labRow =rownames(tmp),
          dendrogram = "none",
          scale="none",breaks=seq(0,5,0.05),margins = c(10, 10),
          colRow=L_col[tmp[1:146,dim(tmp)[2]]],
          density.info = "none",trace="none",symm=F,key=F,
          key.title = "log10(binding sites)",key.xlab = NA,key.ylab=NA,
          keysize = 0.5,colCol = L_col[unlist(tmp["Col",1:146])],
          symbreaks=F,labCol =colnames(tmp),cexRow=0.5,cexCol=0.5,
          col =heatPalette)

dev.off()

RBP4<-c("YBX3", "IGF2BP2", "G3BP1", "RBM15")
pdf("RBP_site_overlap_4RPBs.pdf",height=11,width=6)
tmp<-data.frame(t(RBP_matrix[rownames(RBP_matrix)%in%RBP4,]))
tmp<-tmp[!rownames(tmp)%in%Core_splicesome,]
tmp$Col<-1
tmp$Col[rownames(tmp)%in%Clusters[[1]]]<-2
tmp$Col[rownames(tmp)%in%Clusters[[2]]]<-3
tmp$Col[rownames(tmp)%in%Clusters[[3]]]<-4
tmp$Col[rownames(tmp)%in%Clusters[[4]]]<-5
tmp$Col[rownames(tmp)%in%Clusters[[5]]]<-6
tmp$Col[rownames(tmp)%in%Clusters[[6]]]<-7
heatmap.2(as.matrix(tmp[,1:(dim(tmp)[2]-1)]),labRow =rownames(tmp),
          dendrogram = "none",
          scale="none",breaks=seq(0,5,0.05),margins = c(10, 10),
          main="4 RBPs",colRow=L_col[tmp[,dim(tmp)[2]]],
          density.info = "none",trace="none",symm=F,key=T,
          key.title = "log10(binding sites)",key.xlab = NA,key.ylab=NA,
          keysize = 1,colCol = rep(L_col[i+1],length(RBP4)),
          symbreaks=F,labCol =colnames(tmp),cexRow=0.5,
          col =heatPalette)

dev.off()
### scatter plots , CPM of frac
'''
dat<-read.csv("/stor/work/Lambowitz/yaojun/Work/Documents/NGS/full_length_intron/Cellular_frac/raw_counts.csv",row.names=1)
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
library(RColorBrewer)
pcol<-c("gray80","red","skyblue","orange","brown","blue1")
pcol<-col2hex(pcol)
pcol<-paste0(pcol,"80")
### 
pdf("Cell_frac.scatter.pdf",height=8,width=8)
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
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
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
png("Cell_frac.scatter1.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter2.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter3.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter4.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter5.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter6.png",height=3,width=12,units = "in",res = 600)
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

png("Cell_frac.scatter7.png",height=3,width=12,units = "in",res = 600)
par(mfrow=c(1,4),pty="s",pch=16,mar=c(1,1,1,1))

## MT genes
tmp<-dat[dat$Type2=="MT",]
tmp$Type<-factor(tmp$Type)
plot(tmp$HeLa_Nuc~tmp$HeLa_Cyto,xlim=c(0,20),ylim=c(0,20),col=ncol[tmp$Type],
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
pdf("Cell_frac.scatter_legend.pdf",height=9,width=6)
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
legend("bottomright",legend = levels(tmp$Type),col=ncol,pch=16,bty="n")

dev.off()


### RBP overlap???
### all 4 cell cytoplasm FLEXIs
### Res names: "K562" "Hela" "MDA"  "MCF"  "All" 
#barplot of Cluster 1-6 FLEXIs in cyto and nuc
for (i in 1:4){
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  bar_tmp<-data.frame(c(length(unique(tmp$ID[tmp$RBP%in%Clusters[[1]]])),
                               length(unique(tmp$ID[tmp$RBP%in%Clusters[[2]]])),
                               length(unique(tmp$ID[tmp$RBP%in%Clusters[[3]]])),
                               length(unique(tmp$ID[tmp$RBP%in%Clusters[[4]]])),
                               length(unique(tmp$ID[tmp$RBP%in%Clusters[[5]]])),
                               length(unique(tmp$ID[tmp$RBP%in%Clusters[[6]]])),0))
  bar_tmp[7,1]<-length(unique(tmp$ID))-sum(bar_tmp[,1])
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange<0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  bar_tmp<-cbind(bar_tmp,data.frame(c(length(unique(tmp$ID[tmp$RBP%in%Clusters[[1]]])),
                        length(unique(tmp$ID[tmp$RBP%in%Clusters[[2]]])),
                        length(unique(tmp$ID[tmp$RBP%in%Clusters[[3]]])),
                        length(unique(tmp$ID[tmp$RBP%in%Clusters[[4]]])),
                        length(unique(tmp$ID[tmp$RBP%in%Clusters[[5]]])),
                        length(unique(tmp$ID[tmp$RBP%in%Clusters[[6]]])),0)))
  bar_tmp[7,2]<-length(unique(tmp$ID))-sum(bar_tmp[,2])
  bar_tmp<-bar_tmp[c(7,1:6),]
  barplot(prop.table(as.matrix(bar_tmp[-1,]),2),beside=F,col=col[-1])
}

pdf("RBP_overlaps_bar.pdf",width=11,height=8)
par(mfrow=c(2,3))
### overalpping in cluster 1
for (k in 1:4){
  tmp<-Res[[k]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  # Cluster 1
  tmp1<-tmp[tmp$RBP%in%Clusters[[1]],]
  F_id<-unique(tmp1$ID)
  for (i in 1:length(F_id)){
    if (i==1){bar_tmp<-rep(0,5)}
    id<-F_id[i]
    FLEXI_tmp<-tmp1[tmp1$ID%in%id,c(1:4,10,8:9)]
    FLEXI_tmp<-FLEXI_tmp[FLEXI_tmp$RBP%in%Clusters[[1]],]
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
barplot(prop.table(bar_4cell,2),beside=F,col=col)
### overalpping in cluster 2
for (k in 1:4){
  tmp<-Res[[k]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  # Cluster 1
  tmp1<-tmp[tmp$RBP%in%Clusters[[2]],]
  F_id<-unique(tmp1$ID)
  for (i in 1:length(F_id)){
    print (i)
    if (i==1){bar_tmp<-rep(0,5)}
    id<-F_id[i]
    FLEXI_tmp<-tmp1[tmp1$ID%in%id,c(1:4,10,8:9)]
    FLEXI_tmp<-FLEXI_tmp[FLEXI_tmp$RBP%in%Clusters[[2]],]
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
barplot(prop.table(bar_4cell,2),beside=F,col=col)
dev.off()

pdf("RBP_overlaps_bar.pdf",width=11,height=8)
par(mfrow=c(2,3))
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
  barplot(prop.table(bar_4cell,2),beside=F,col=col,main=paste0("Cluster ",m))
}
dev.off()

### upset of only cytoplasm FLEXI
### get some numbers
pdf("Cluster_combination_CytoFLEXI.pdf",height=11,width=10)
par(mfrow=c(4,2))
for (i in 1:4){
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][1]])
  set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][2]])
  set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][3]])
  set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][4]])
  set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][5]])
  set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][6]])
  set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][7]])
  set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][8]])
  set <- list (set1,set2,set3,set4,set5,set6,set7,set8)
  for (j in 1:(length(set)-1)){
    for (k in (j+1):length(set)){
      Over<-0
      RBP1<-Clusters[[1]][j]
      RBP2<-Clusters[[1]][k]
      sub_tmp<-tmp[tmp$RBP%in%c(RBP1,RBP2),]
      for (l in 1:dim(sub_tmp)[1]){
        sub_tmp$range[l]<-list(c(sub_tmp$RBP_st[l]:sub_tmp$RBP_ed[l]))
      }
      Unique_id<-unique(sub_tmp$ID)
      for (l in 1:length(Unique_id)){
        F_id<-Unique_id[l]
        F_tmp<-sub_tmp[sub_tmp$ID%in%F_id,]
        range1<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP1]))
        range2<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP2]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          Over<-Over+1
        }
      }
      
      
      if (j==1 &k==2){
        Clu1_up<-data.frame(Clusters[[1]][j],
                            Clusters[[1]][k],
                            Over,
                            length(intersect(set[[j]],set[[k]]))-Over)
      } else {
        Clu1_up<-rbind(Clu1_up,
                       data.frame(Clusters[[1]][j],
                                  Clusters[[1]][k],
                                  Over,
                                  length(intersect(set[[j]],set[[k]]))-Over))
      }
    }
  }
  Clu1_up$Total<-rowSums(Clu1_up[,3:4])
  Clu1_up<-Clu1_up[Clu1_up$Total>0,]
  Clu1_up[,1]<-factor(Clu1_up[,1],levels=names(sort(table(rep(Clu1_up[,1],Clu1_up[,5])),decreasing = T)))
  Clu1_up[,2]<-factor(Clu1_up[,2],levels=names(sort(table(rep(Clu1_up[,2],Clu1_up[,5])),decreasing = T)))
  Clu1_up<-Clu1_up[order(Clu1_up[,1],-Clu1_up[,5]),]
  
  bp<-barplot(t(as.matrix(Clu1_up[,3:4])),beside=F,xaxt='n',
              main=paste0(names(Res)[i],": Cluster 1"))
  axis(1,at=bp,labels = Clu1_up[,1],cex.axis=0.5,las=2)
  axis(1,at=bp,labels = Clu1_up[,2],line = 2,tick = F,cex.axis=0.5,las=2)
  
  ###cluster 2
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][1]])
  set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][2]])
  set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][3]])
  set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][4]])
  set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][5]])
  set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][6]])
  set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][7]])
  set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][8]])
  set9<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][9]])
  set <- list (set1,set2,set3,set4,set5,set6,set7,set8,set9)
  for (j in 1:(length(set)-1)){
    for (k in (j+1):length(set)){
      Over<-0
      RBP1<-Clusters[[2]][j]
      RBP2<-Clusters[[2]][k]
      sub_tmp<-tmp[tmp$RBP%in%c(RBP1,RBP2),]
      for (l in 1:dim(sub_tmp)[1]){
        sub_tmp$range[l]<-list(c(sub_tmp$RBP_st[l]:sub_tmp$RBP_ed[l]))
      }
      Unique_id<-unique(sub_tmp$ID)
      for (l in 1:length(Unique_id)){
        F_id<-Unique_id[l]
        F_tmp<-sub_tmp[sub_tmp$ID%in%F_id,]
        range1<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP1]))
        range2<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP2]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          Over<-Over+1
        }
      }
      
      
      if (j==1 &k==2){
        Clu1_up<-data.frame(Clusters[[2]][j],
                            Clusters[[2]][k],
                            Over,
                            length(intersect(set[[j]],set[[k]]))-Over)
      } else {
        Clu1_up<-rbind(Clu1_up,
                       data.frame(Clusters[[2]][j],
                                  Clusters[[2]][k],
                                  Over,
                                  length(intersect(set[[j]],set[[k]]))-Over))
      }
    }
  }
  Clu1_up$Total<-rowSums(Clu1_up[,3:4])
  Clu1_up<-Clu1_up[Clu1_up$Total>0,]
  Clu1_up[,1]<-factor(Clu1_up[,1],levels=names(sort(table(rep(Clu1_up[,1],Clu1_up[,5])),decreasing = T)))
  Clu1_up[,2]<-factor(Clu1_up[,2],levels=names(sort(table(rep(Clu1_up[,2],Clu1_up[,5])),decreasing = T)))
  Clu1_up<-Clu1_up[order(Clu1_up[,1],-Clu1_up[,5]),]
  
  bp<-barplot(t(as.matrix(Clu1_up[,3:4])),beside=F,xaxt='n',
              main=paste0(names(Res)[i],": Cluster 2"))
  axis(1,at=bp,labels = Clu1_up[,1],cex.axis=0.5,las=2)
  axis(1,at=bp,labels = Clu1_up[,2],line = 2,tick = F,cex.axis=0.5,las=2)
}
dev.off()

pdf("Cluster_combination_NucFLEXI.pdf",height=11,width=10)
par(mfrow=c(4,2))
for (i in 1:4){
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange<0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][1]])
  set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][2]])
  set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][3]])
  set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][4]])
  set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][5]])
  set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][6]])
  set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][7]])
  set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[1]][8]])
  set <- list (set1,set2,set3,set4,set5,set6,set7,set8)
  for (j in 1:(length(set)-1)){
    for (k in (j+1):length(set)){
      Over<-0
      RBP1<-Clusters[[1]][j]
      RBP2<-Clusters[[1]][k]
      sub_tmp<-tmp[tmp$RBP%in%c(RBP1,RBP2),]
      for (l in 1:dim(sub_tmp)[1]){
        sub_tmp$range[l]<-list(c(sub_tmp$RBP_st[l]:sub_tmp$RBP_ed[l]))
      }
      Unique_id<-unique(sub_tmp$ID)
      for (l in 1:length(Unique_id)){
        F_id<-Unique_id[l]
        F_tmp<-sub_tmp[sub_tmp$ID%in%F_id,]
        range1<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP1]))
        range2<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP2]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          Over<-Over+1
        }
      }
      
      
      if (j==1 &k==2){
        Clu1_up<-data.frame(Clusters[[1]][j],
                            Clusters[[1]][k],
                            Over,
                            length(intersect(set[[j]],set[[k]]))-Over)
      } else {
        Clu1_up<-rbind(Clu1_up,
                       data.frame(Clusters[[1]][j],
                                  Clusters[[1]][k],
                                  Over,
                                  length(intersect(set[[j]],set[[k]]))-Over))
      }
    }
  }
  Clu1_up$Total<-rowSums(Clu1_up[,3:4])
  Clu1_up<-Clu1_up[Clu1_up$Total>0,]
  Clu1_up[,1]<-factor(Clu1_up[,1],levels=names(sort(table(rep(Clu1_up[,1],Clu1_up[,5])),decreasing = T)))
  Clu1_up[,2]<-factor(Clu1_up[,2],levels=names(sort(table(rep(Clu1_up[,2],Clu1_up[,5])),decreasing = T)))
  Clu1_up<-Clu1_up[order(Clu1_up[,1],-Clu1_up[,5]),]
  
  bp<-barplot(t(as.matrix(Clu1_up[,3:4])),beside=F,xaxt='n',
              main=paste0(names(Res)[i],": Cluster 1"))
  axis(1,at=bp,labels = Clu1_up[,1],cex.axis=0.5,las=2)
  axis(1,at=bp,labels = Clu1_up[,2],line = 2,tick = F,cex.axis=0.5,las=2)
  
  ###cluster 2
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  set1<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][1]])
  set2<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][2]])
  set3<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][3]])
  set4<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][4]])
  set5<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][5]])
  set6<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][6]])
  set7<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][7]])
  set8<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][8]])
  set9<-unique(tmp$ID[tmp$RBP%in%Clusters[[2]][9]])
  set <- list (set1,set2,set3,set4,set5,set6,set7,set8,set9)
  for (j in 1:(length(set)-1)){
    for (k in (j+1):length(set)){
      Over<-0
      RBP1<-Clusters[[2]][j]
      RBP2<-Clusters[[2]][k]
      sub_tmp<-tmp[tmp$RBP%in%c(RBP1,RBP2),]
      for (l in 1:dim(sub_tmp)[1]){
        sub_tmp$range[l]<-list(c(sub_tmp$RBP_st[l]:sub_tmp$RBP_ed[l]))
      }
      Unique_id<-unique(sub_tmp$ID)
      for (l in 1:length(Unique_id)){
        F_id<-Unique_id[l]
        F_tmp<-sub_tmp[sub_tmp$ID%in%F_id,]
        range1<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP1]))
        range2<-unique(unlist(F_tmp$range[F_tmp$RBP==RBP2]))
        inter_tmp<-intersect(range1,range2)
        if (length(inter_tmp)>0){
          Over<-Over+1
        }
      }
      
      
      if (j==1 &k==2){
        Clu1_up<-data.frame(Clusters[[2]][j],
                            Clusters[[2]][k],
                            Over,
                            length(intersect(set[[j]],set[[k]]))-Over)
      } else {
        Clu1_up<-rbind(Clu1_up,
                       data.frame(Clusters[[2]][j],
                                  Clusters[[2]][k],
                                  Over,
                                  length(intersect(set[[j]],set[[k]]))-Over))
      }
    }
  }
  Clu1_up$Total<-rowSums(Clu1_up[,3:4])
  Clu1_up<-Clu1_up[Clu1_up$Total>0,]
  Clu1_up[,1]<-factor(Clu1_up[,1],levels=names(sort(table(rep(Clu1_up[,1],Clu1_up[,5])),decreasing = T)))
  Clu1_up[,2]<-factor(Clu1_up[,2],levels=names(sort(table(rep(Clu1_up[,2],Clu1_up[,5])),decreasing = T)))
  Clu1_up<-Clu1_up[order(Clu1_up[,1],-Clu1_up[,5]),]
  
  bp<-barplot(t(as.matrix(Clu1_up[,3:4])),beside=F,xaxt='n',
              main=paste0(names(Res)[i],": Cluster 2"))
  axis(1,at=bp,labels = Clu1_up[,1],cex.axis=0.5,las=2)
  axis(1,at=bp,labels = Clu1_up[,2],line = 2,tick = F,cex.axis=0.5,las=2)
}
dev.off()

### COMBination of Cluster 1
library(RColorBrewer)
R_col<-c("gray80",brewer.pal(8,"Set1"))
R_col<-col2hex(R_col)
R_col<-paste0(R_col,"80")
B_col<-c("gray80",brewer.pal(9,"Set3"))
B_col<-col2hex(B_col)
B_col<-paste0(B_col,"80")
pdf("Cluster1_2_sites.pdf",height=11,width=10)
par(mfrow=c(5,2))
for (i in 1:4){
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  tmp<-tmp[tmp$RBP%in%Clusters[[1]],]
  Unique_id<-unique(tmp$ID)
  ### retract FLEXIs by id
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%Unique_id,]
  tmp$Col<-1
  tmp$Col[tmp$RBP%in%Clusters[[1]][1]]<-2
  tmp$Col[tmp$RBP%in%Clusters[[1]][2]]<-3
  tmp$Col[tmp$RBP%in%Clusters[[1]][3]]<-4
  tmp$Col[tmp$RBP%in%Clusters[[1]][4]]<-5
  tmp$Col[tmp$RBP%in%Clusters[[1]][5]]<-6
  tmp$Col[tmp$RBP%in%Clusters[[1]][6]]<-7
  tmp$Col[tmp$RBP%in%Clusters[[1]][7]]<-8
  tmp$Col[tmp$RBP%in%Clusters[[1]][8]]<-9
  plot(NULL,xlim=c(-25,125),ylim=c(1,length(Unique_id)*3),
       bty="n",axes=F,xlab=NA,ylab="FLEXIs",main=paste0(names(Res)[i],": Cluster I"))
  axis(1,at=seq(-25,125,25),labels = seq(-25,125,25))
  print (length(Unique_id))
  for (j in 1:(length(Unique_id))){
    F_tmp<-tmp[tmp$ID==Unique_id[j],]
    F_tmp<-F_tmp[order(F_tmp$Col,F_tmp$RBP),]
    # ignore strandness
    F_tmp$Len<-(F_tmp$FLEXI_ed-F_tmp$FLEXI_st)/100
    F_tmp$Left<-(F_tmp$RBP_st-F_tmp$FLEXI_st)/F_tmp$Len
    F_tmp$Right<-(F_tmp$RBP_ed-F_tmp$FLEXI_st)/F_tmp$Len
    F_tmp$Left[F_tmp$Left< -25]<- -25
    F_tmp$Right[F_tmp$Right> 125]<- 125
    for (k in 1:dim(F_tmp)[1]){
      if (k==1){
        segments(0,j*3,100,j*3)
      }
      rect(F_tmp$Left[k],(j*3)-0.75,F_tmp$Right[k],(j*3)+0.75,col=R_col[F_tmp$Col[k]])
      
    }
  }
  
  
  ###cluster 2
  tmp<-Res[[i]]
  tmp<-sub("_FLEXI","",rownames(tmp)[tmp$log2FoldChange>0])
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%tmp,]
  tmp<-tmp[tmp$RBP%in%Clusters[[2]],]
  Unique_id<-unique(tmp$ID)
  print(length(Unique_id))
  ### retract FLEXIs by id
  tmp<-FLEXI_RBP[FLEXI_RBP$ID%in%Unique_id,]
  tmp$Col<-1
  tmp$Col[tmp$RBP%in%Clusters[[2]][1]]<-2
  tmp$Col[tmp$RBP%in%Clusters[[2]][2]]<-3
  tmp$Col[tmp$RBP%in%Clusters[[2]][3]]<-4
  tmp$Col[tmp$RBP%in%Clusters[[2]][4]]<-5
  tmp$Col[tmp$RBP%in%Clusters[[2]][5]]<-6
  tmp$Col[tmp$RBP%in%Clusters[[2]][6]]<-7
  tmp$Col[tmp$RBP%in%Clusters[[2]][7]]<-8
  tmp$Col[tmp$RBP%in%Clusters[[2]][8]]<-9
  tmp$Col[tmp$RBP%in%Clusters[[2]][9]]<-10
  plot(NULL,xlim=c(-25,125),ylim=c(1,length(Unique_id)*3),
       bty="n",axes=F,xlab=NA,ylab="FLEXIs",main=paste0(names(Res)[i],": Cluster II"))
  axis(1,at=seq(-25,125,25),labels = seq(-25,125,25))
  for (j in 1:(length(Unique_id))){
    F_tmp<-tmp[tmp$ID==Unique_id[j],]
    F_tmp<-F_tmp[order(F_tmp$Col,F_tmp$RBP),]
    # ignore strandness
    F_tmp$Len<-(F_tmp$FLEXI_ed-F_tmp$FLEXI_st)/100
    F_tmp$Left<-(F_tmp$RBP_st-F_tmp$FLEXI_st)/F_tmp$Len
    F_tmp$Right<-(F_tmp$RBP_ed-F_tmp$FLEXI_st)/F_tmp$Len
    F_tmp$Left[F_tmp$Left< -25]<- -25
    F_tmp$Right[F_tmp$Right> 125]<- 125
    for (k in 1:dim(F_tmp)[1]){
      if (k==1){
        segments(0,j*3,100,j*3)
      }
      rect(F_tmp$Left[k],(j*3)-0.75,F_tmp$Right[k],(j*3)+0.75,col=B_col[F_tmp$Col[k]])
      
    }
  }
}
plot.new()
legend("bottom",fill=R_col,legend = c("Other",Clusters[[1]]),bty="n")
plot.new()
legend("bottom",fill=B_col,legend = c("Other",Clusters[[2]]),bty="n")
dev.off()
