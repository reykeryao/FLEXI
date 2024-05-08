library(gplots)
library(RColorBrewer)
library(openxlsx)
Clusters<-list(Cluster1=c("LARP4","PABPC4","SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"),
               Cluster2=c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24","FXR2"),
               Cluster3=c("TIA1","TIAL1"),
               Cluster4=c("U2AF1","U2AF2","KHSRP"),
               Cluster5=c("AATF","DKC1","NOLC1","SMNDC1"),
               Cluster6=c("AGO","DICER"))
'''
## preprare ENSEMBL ID for DAVID
#### for DAVID, FLEXI/Other shirt intron/long intron ; host gene go
### read FLEXI, 4 cell lines only
dat<-read.delim("all.FLEXI")
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
FLEXI<-dat$ID[rowSums(dat[,73:76])>0]
FLEXI_list<-list("FLEXI_bg"=unique(dat$GID[rowSums(dat[,73:76])>0]))

FLEXI<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI,]
FLEXI<-unique(FLEXI)
#### read all intron intersect info
all_intron<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron<-unique(all_intron)
### remove all FLEXI introns
all_intron<-all_intron[!all_intron$FLEXI%in%unique(FLEXI$ID),]
### make long introns
long_intron<-all_intron[all_intron$Len>300,c(1,3)]
#make Other short introns
other_short_intron<-all_intron[all_intron$Len<=300,c(1,3)]

# all intron
tmp<-read.table(gzfile("all_intron_ID_len.txt.gz"),col.names=c("ID","Len"))
tmp<-separate(tmp,"ID",into =c("IID","GID","TID","Gtype","TType"),sep="___",remove=F)
tmp_long<-tmp[tmp$Len>300,]
tmp_ID<-dat$ID[rowSums(dat[,73:76])>0]
tmp_short<-tmp[tmp$Len<=300 & !tmp$ID%in%tmp_ID,]
OtherShortIntron_list<-list("OSI_bg"=unique(tmp_short$GID))
LongIntron_list<-list("LongIntron_bg"=unique(tmp_long$GID))
LongIntron_match<-list("LongIntron_bg"=sample(unique(tmp_long$GID),length(OtherShortIntron_list[[1]]),replace = F))
### make gene list and background

### convert intron ID to host genes
FLEXI<-separate(FLEXI,ID, into =c("IID","GID","TID","Gtype","TType"),sep="___")
FLEXI<-unique(FLEXI[,c("GID","RBP")])
other_short_intron<-separate(other_short_intron,FLEXI, into =c("IID","GID","TID","Gtype","TType"),sep="___")
other_short_intron<-unique(other_short_intron[,c("GID","RBP")])
long_intron<-separate(long_intron,FLEXI, into =c("IID","GID","TID","Gtype","TType"),sep="___")
long_intron<-unique(long_intron[,c("GID","RBP")])

for (i in 1:length(Clusters)){
  cluster_RBP<-Clusters[[i]]
  for (j in 1:length(cluster_RBP)){
    RBP_name<-cluster_RBP[j]
    FLEXI_list<-c(FLEXI_list,setNames(list(unique(FLEXI$GID[FLEXI$RBP==RBP_name])),RBP_name))
    OtherShortIntron_list<-c(OtherShortIntron_list,setNames(list(unique(other_short_intron$GID[other_short_intron$RBP==RBP_name])),RBP_name))
    pick_len<-length(unique(other_short_intron$GID[other_short_intron$RBP==RBP_name]))
    LongIntron_list<-c(LongIntron_list,setNames(list(unique(long_intron$GID[long_intron$RBP==RBP_name])),RBP_name))
  }
}

write.xlsx(FLEXI_list,"DAVID_hostgene/FLEXI_cluster1To6.hostgene.xlsx")
write.xlsx(OtherShortIntron_list,"DAVID_hostgene/OtherSHortIntron_cluster1To6.hostgene.xlsx")
write.xlsx(LongIntron_list,"DAVID_hostgene/LongIntron_cluster1To6.hostgene.xlsx")

## write the rest 25 from 53 RBPs

RBP53<-read.delim("53_RBP_info_fig4_7.txt")
RBP53<-RBP53$RBP.name
RBP53<-RBP53[!RBP53%in%unlist(Clusters)]
FLEXI_list<-list("FLEXI_bg"=unique(dat$GID[rowSums(dat[,73:76])>0]))
for (i in 1:length(RBP53)){
  RBP_name<-RBP53[[i]]
  FLEXI_list<-c(FLEXI_list,setNames(list(unique(FLEXI$GID[FLEXI$RBP==RBP_name])),RBP_name))
}
write.xlsx(FLEXI_list,"DAVID_hostgene/FLEXI_other25.hostgene.xlsx")

### DAVID results cleanup
RBP_name<-unlist(Clusters)
for (i in 1:length(RBP_name)){
  RBP<-RBP_name[i]
  F_res<-read.delim(paste0("DAVID_hostgene/FLEXI/",RBP,".txt"))
  F_res<-F_res[F_res$FDR<0.1 ,c(1,2,13)]
  S_res<-read.delim(paste0("DAVID_hostgene/OtherShortIntron/",RBP,".txt"))
  S_res<-S_res[S_res$FDR<0.1 ,c(1,2,13)]
  L_res<-read.delim(paste0("DAVID_hostgene/LongIntron/",RBP,".txt"))
  L_res<-L_res[L_res$FDR<0.1 ,c(1,2,13)]
  LM_res<-read.delim(paste0("DAVID_hostgene/LongIntron_matched/",RBP,".txt"))
  LM_res<-LM_res[LM_res$FDR<0.1 ,c(1,2,13)]
  colnames(F_res)[3]<-colnames(S_res)[3]<-colnames(L_res)[3]<-colnames(LM_res)[3]<-RBP
  if(i==1){
    res1<-F_res
    res2<-S_res
    res3<-L_res
    res4<-LM_res
  }else {
    res1<-merge(res1,F_res,by=1:2,all=T)
    res2<-merge(res2,S_res,by=1:2,all=T)
    res3<-merge(res3,L_res,by=1:2,all=T)
    res4<-merge(res4,LM_res,by=1:2,all=T)
  }
}
res1[is.na(res1)]<-0.1
res2[is.na(res2)]<-0.1
res3[is.na(res3)]<-0.1
res4[is.na(res4)]<-0.1
res1[,-1:-2]<- -log10(res1[,-1:-2])
res2[,-1:-2]<- -log10(res2[,-1:-2])
res3[,-1:-2]<- -log10(res3[,-1:-2])
res4[,-1:-2]<- -log10(res4[,-1:-2])
res<-list("FLEXI"=res1,"OSI"=res2,"Longintron"=res3,"LongIntronMatched"=res4)
saveRDS(res,"DAVID_hostgene/combined_cleaned.res")
'''
RBP53<-read.delim("53_RBP_info_fig4_7.txt")
RBP53<-RBP53$RBP.name
for (i in 1:length(RBP53)){
  RBP<-RBP53[i]
  A_res<-read.delim(paste0("DAVID_hostgene/FLEXI/",RBP,".txt"))
  A_res<-A_res[A_res$FDR<0.1 ,c(1,2,13)]
  colnames(A_res)[3]<-RBP
  if(i==1){
    res1<-A_res
  }else {
    res1<-merge(res1,A_res,by=1:2,all=T)
  }
}
res1[is.na(res1)]<-0.1
res1[,-1:-2]<- -log10(res1[,-1:-2])
res1<-separate(res1,Term,into=c("GO","Term"),sep="~")
### all 53
### top 50
res1$Sum<-rowSums(res1[,-1:-3])
res1<-res1[order(res1$Sum,decreasing = T),]
heatPalette = colorRampPalette(c("dodgerblue4",  "white", "orangered"))(100)
pdf("Figures/Fig4B.pdf",width=8,height=9)
ht<-heatmap.2(as.matrix(res1[1:50,c(-1:-3,-57)]),labRow = res1$Term,dendrogram = "none",
              scale="none",margins = c(10, 20),breaks=seq(1,5,0.04),
              density.info = "none",trace="none",symm=F,key=T,
              symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=F,cexRow=0.2,cexCol=0.2,
              col =heatPalette)
dev.off()

### 
res<-readRDS("DAVID_hostgene/combined_cleaned.res")
res1<-res[[1]]
res2<-res[[2]]
res3<-res[[3]]
res4<-res[[4]]
### find common terms
terms<-intersect(intersect(res1$Term,res2$Term),res3$Term)
res1<-res1[rowMeans(res1[,-1:-2])>1.1,]
terms<-intersect(res1$Term,terms)
res1<-res1[res1$Term%in%terms,]
res2<-res2[res2$Term%in%terms,]
res3<-res3[res3$Term%in%terms,]
res4<-res4[res4$Term%in%terms,]

res1<-separate(res1,Term,into=c("GO","Term"),sep="~")
res2<-separate(res2,Term,into=c("GO","Term"),sep="~")
res3<-separate(res3,Term,into=c("GO","Term"),sep="~")
res4<-separate(res4,Term,into=c("GO","Term"),sep="~")

### sort res1,2,3 by terms
res1<-res1[order(res1$Term),]
res2<-res2[order(res2$Term),]
res3<-res3[order(res3$Term),]
## heatmap
heatPalette = colorRampPalette(c("dodgerblue4",  "white", "orangered"))(100)
pdf("Figures/Fig7B.pdf",width=8,height=9)
ht<-heatmap.2(as.matrix(res1[,-1:-3]),labRow = res1$Term,dendrogram = "none",
              scale="none",margins = c(10, 20),breaks=seq(1,5,0.04),
              density.info = "none",trace="none",symm=F,key=T,
              symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=T,cexRow=0.7,cexCol=0.7,
              col =heatPalette)
dev.off()
### make OSI teh same order
res2<-res2[order(match(res2$Term,res1$Term[ht$rowInd])),]
res3<-res3[order(match(res3$Term,res1$Term[ht$rowInd])),]
pdf("Figures/Fig9_1.pdf",width=8,height=9)
heatmap.2(as.matrix(res2[,-1:-3]),labRow = res2$Term,dendrogram = "none",
          scale="none",margins = c(10, 20),breaks=seq(1,5,0.04),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=F,cexRow=0.7,cexCol=0.7,
          col =heatPalette)
dev.off()
pdf("Figures/Fig9_2.pdf",width=8,height=9)
heatmap.2(as.matrix(res3[,-1:-3]),labRow = res3$Term,dendrogram = "none",
          scale="none",margins = c(10, 20),breaks=seq(1,30,0.29),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=F,cexRow=0.7,cexCol=0.7,
          col =heatPalette)
dev.off()

pdf("DAVID_hostgene/LI_match_DAVID_heat.pdf",width=8,height=9)
heatmap.2(as.matrix(res4[,-1:-3]),labRow = res4$Term,dendrogram = "none",
          scale="none",margins = c(10, 20),breaks=seq(1,5,0.04),
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,keysize=1,symkey=F,Colv=F,Rowv=T,cexRow=0.7,cexCol=0.7,
          col =heatPalette)
dev.off()
