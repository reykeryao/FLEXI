dat<-read.delim("all.FLEXI")
mapped_reads<-c(207.491024,692.091831,666.341854,713.775291,715.241521,768.433748,71.116246)
names(mapped_reads)<-colnames(dat[,71:77])
gene_counts<-read.delim("combined_counts.tsv")
FLEXI_CPM<-dat[,c(1,8,73:76)]
Unfrag_total<-mapped_reads[3:6]
FLEXI_CPM[,3:6]<-t(t(FLEXI_CPM[,3:6])/Unfrag_total)

rownames(gene_counts)<-gene_counts$ID
gene_counts$ID<-as.character(gene_counts$ID)
#Subset only FLEXI host genes
GID_list<-unique(FLEXI_CPM$GID)
gene_counts<-gene_counts[gene_counts$ID%in%GID_list,c(1,6:9)]
gene_counts[,2:5]<-t(t(gene_counts[,2:5])/Unfrag_total)
colnames(gene_counts)[-1]<-paste0(colnames(gene_counts)[-1],"_gene")
colnames(gene_counts)[1]<-"GID"
tmp<-merge(FLEXI_CPM,gene_counts,by="GID",all=T)
tmp[is.na(tmp)]<-0
tmp<-tmp[,c(1,2,4,5,6,3,8,9,10,7)]


### process non-FLEXI
file_name<-c("K562","HEK","Hela","UHRR")
for (i in 1:4){
  per_dat<-read.delim(gzfile(paste0(file_name[i],".per.info.gz")))
  per_dat<-per_dat[per_dat$nonFLEXI==1,]
  SI<-aggregate(Per~ID,per_dat,sum)
  colnames(SI)<-c("ID",file_name[i])
  if (i==1){
    SI_tmp<-SI
  } else {
    SI_tmp<-merge(SI_tmp,SI,by=1,all=T)
  }
}
SI_tmp[is.na(SI_tmp)]<-0
SI_tmp[,-1]<-SI_tmp[,-1]/100
non_FLEXI<-SI_tmp[!SI_tmp$ID%in%FLEXI_CPM$ID,]
non_FLEXI<-separate(non_FLEXI,"ID",into=c("IID","GID","TID","GTYpe","TTYpe"),sep="___",remove = F)
non_FLEXI<-non_FLEXI[,c(1,3,7:10)]
non_FLEXI[,3:6]<-t(t(non_FLEXI[,3:6])/Unfrag_total)

gene_counts<-read.delim("combined_counts.tsv")
rownames(gene_counts)<-gene_counts$ID
gene_counts$ID<-as.character(gene_counts$ID)
#Subset only non-FLEXI host genes
GID_list<-unique(non_FLEXI$GID)
gene_counts<-gene_counts[gene_counts$ID%in%GID_list,c(1,6:9)]
gene_counts[,2:5]<-t(t(gene_counts[,2:5])/Unfrag_total)
colnames(gene_counts)[-1]<-paste0(colnames(gene_counts)[-1],"_gene")
colnames(gene_counts)[1]<-"GID"
tmp1<-merge(non_FLEXI,gene_counts,by="GID",all=T)
tmp1[is.na(tmp1)]<-0
tmp1<-tmp1[,c(1:6,8,9,10,7)]

### FLEXI and non FLEXI on teh same gene
GID_list<-unique(FLEXI_CPM$GID)
tmp2<-SI_tmp
tmp2<-separate(tmp2,"ID",into=c("IID","GID","TID","GTYpe","TTYpe"),sep="___",remove = F)
tmp2<-tmp2[,c(1,3,7:10)]
tmp2[,3:6]<-t(t(tmp2[,3:6])/Unfrag_total)
tmp2<-tmp2[tmp2$GID%in%GID_list,]
colnames(tmp2)[-1:-2]<-paste0(colnames(tmp2)[-1:-2],"_nonFLEXI")
tmp2<-merge(FLEXI_CPM,tmp2,by="GID",all=T)
tmp2<-tmp2[,c(1,2,7,4,5,6,3,8:11)]
tmp2[is.na(tmp2)]<-0

pdf("temp_fig/FigS7A.pdf",width = 12,height=9)
pck_name<-c("K-562","HEK-293T","HeLa S3","UHRR")
par(pch=16,mfcol=c(3,4),pty="s")
for (i in 1:4){
  tmp_for_plot<-tmp[,c(i+6,i+2)]
  tmp_for_plot<-tmp_for_plot[rowSums(tmp_for_plot)>0,]
  plot(log2(tmp_for_plot),xlim=c(-15,15),ylim=c(-15,15),main=pck_name[i],
       ylab="log2RPM (FLEXI)",xlab="log2RPM (FLEXI host gene)")
  abline(0,1,col="red")
  cor_p=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "pearson"),digits=3, format="f")
  cor_s=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "spearman"),digits=2, format="f")
  text(-7,13,bquote(atop(italic(r)== .(cor_p)~phantom())))
  text(-7,11,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
  
  ###non-FLEXI OSI
  tmp_for_plot<-tmp1[,c(i+6,i+2)]
  tmp_for_plot<-tmp_for_plot[rowSums(tmp_for_plot)>0,]
  plot(log2(tmp_for_plot),xlim=c(-15,15),ylim=c(-15,15),main=pck_name[i],
       ylab="log2RPM (Other short intron)",xlab="log2RPM (Other short intron host gene)")
  abline(0,1,col="red")
  cor_p=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "pearson"),digits=3, format="f")
  cor_s=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "spearman"),digits=2, format="f")
  text(-7,13,bquote(atop(italic(r)== .(cor_p)~phantom())))
  text(-7,11,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
  
  ###FLEXI vs non-FLEXI OSI on teh same gene
  tmp_for_plot<-tmp2[,c(i+7,i+3)]
  tmp_for_plot<-tmp_for_plot[rowSums(tmp_for_plot)>0,]
  plot(log2(tmp_for_plot),xlim=c(-15,5),ylim=c(-15,5),main=pck_name[i],
       xlab="log2RPM (Other short intron)",ylab="log2RPM (FLEXI)")
  abline(0,1,col="red")
  cor_p=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "pearson"),digits=3, format="f")
  cor_s=formatC(cor(tmp_for_plot[,1],tmp_for_plot[,2],method = "spearman"),digits=2, format="f")
  text(3,-11,bquote(atop(italic(r)== .(cor_p)~phantom())))
  text(3,-13,bquote(atop(italic(r[s])== .(cor_s)~phantom())))
}
dev.off()

#ccor list for test
cocor_dat<-list(F=data.frame("F1"=FLEXI_by_GID$K562,"F2"=FLEXI_by_GID$Hela))
cocor_dat<-c(cocor_dat,list(G=data.frame("G1"=FLEXI_by_GID$K562,"G2"=FLEXI_by_GID$Hela)))
cocor(~F1 + F2 | G1 + G2, cocor_dat)
#p<0.0001


### FLEXI heatmap
FLEXI<-read.delim("all.FLEXI")
Total<-c(109.704771,97.786253,70.563696,83.971287,93.530133,96.188889,84.195614,95.958499,90.154958,
         77.528755,87.438356,80.725874,84.193148,84.695000,79.516751,94.126936,80.171312,75.474477,
         87.388516,82.600414,84.957488,88.141802,95.927677,113.807174,87.425522,73.526698,96.566478,
         71.023706,90.569418,84.985720,107.210447,92.506910,88.853873,83.524969,91.268315,83.392799,
         82.750963,62.981188,60.057226,67.704164,86.607803,72.132500,96.955981,64.582809)
FLEXI[,27:70]<-t(t(FLEXI[,27:70])/Total)
FLEXI_ID<-FLEXI$ID[rowSums(FLEXI[,73:76])>0]

RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI_ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
FLEXI_ID<-unique(RBP_4cell_plasma$ID)
tmp<-FLEXI[FLEXI$ID%in%FLEXI_ID,27:70]
tmp<-tmp[rowSums(tmp)>0,]
tmp[tmp==0]<-2^-7
tmp<-log2(tmp)

library(gplots)
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)
pdf("temp_fig/FLEXI_heat.pdf",height=4,width=6)
ht<-heatmap.2(as.matrix(tmp),labRow = F,dendrogram = "none",labCol=colnames(tmp),
          scale="none",breaks=seq(-7,1,0.08),
          density.info = "none",trace="none",symm=F,key=F,
          symbreaks=F,Colv=F,Rowv=T,cexCol=0.5,
          col =heatPalette)
pick<-rev(ht$rowInd)[1:50]
tmp1<-tmp[pick,]
heatmap.2(as.matrix(tmp1),labRow = F,dendrogram = "none",labCol=colnames(tmp1),
          scale="none",breaks=seq(-7,1,0.08),
          density.info = "none",trace="none",symm=F,key=F,
          symbreaks=F,Colv=F,Rowv=T,cexCol=0.5,
          col =heatPalette)
dev.off()
### snoRNA vs snoRNA FLEXI
snoFLEX<-dat[dat$Has_snoRNA!=".",]
mapped_reads<-c(207.491024,692.091831,666.341854,713.775291,715.241521,768.433748,71.116246)
names(mapped_reads)<-colnames(dat[,71:77])
gene_counts<-read.delim("combined_counts.tsv")
gene_counts[6:9]<-t(t(gene_counts[6:9])/Unfrag_total)
gene_counts<-gene_counts[,c(1,6:9)]
snoFLEX<-snoFLEX[,c(1:8,15,73:76)]
Unfrag_total<-mapped_reads[3:6]
snoFLEX[,10:13]<-t(t(snoFLEX[,10:13])/Unfrag_total)

snoBed<-read.table("/stor/home/yaojun/Ref/TGSEQ/combined_bed/intron_embedded_sno_scaRNA.bed")
colnames(snoBed)<-c("Chr","St","Ed","Has_snoRNA","Score","Strand","TYpe","GID")
snoBed<-snoBed[,c(8,1:3,6,4)]
snoBed<-merge(snoBed,gene_counts,by=1)
snoBed<-unique(snoBed)

tmp<-unique(merge(snoFLEX,snoBed,by=c("Chr","Strand","Has_snoRNA"),all=T))

tmp2<-tmp[is.na(tmp$ID),]
tmp2[is.na(tmp2)]<-0
tmp1<-tmp[!is.na(tmp$ID),]

tmp1<-tmp1[tmp1$St.x<tmp1$St.y,]
tmp<-rbind(tmp1,tmp2)

tmp<-tmp[,c(3,10:13,17:20)]
tmp[,-1]<-log10(tmp[,-1])

pdf("temp_fig/sno_snoFLEXI.pdf",height=8,width=8)
par(mfrow=c(2,2),psy="s")
for (i in c(2:5)){
  tmp1<-tmp[,c(i,i+4)]
  tmp2<-tmp1[is.finite(tmp1[,1]),]
  tmp1$col="red"
  r<-cor(tmp2)[1,2]
  rs<-cor(tmp2,method = "spearman")[1,2]
  tmp1$col[!is.finite(tmp1[,1])]<-"black"
  tmp1[,1][!is.finite(tmp1[,1])]<- -4
  
  plot(x=tmp1[,1],y=tmp1[,2],xlim=c(-4,0),ylim=c(-2,6),pch=19,main=colnames(tmp1)[1],col=tmp1$col,
      xlab="FLEXI (log10RPM)",ylab="snoRNA (log10RPM)")
  text(-3,2,paste0("r = ",round(r,2)))
  text(-3,1,paste0("rs = ",round(rs,2)))
}
dev.off()

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
### mark reproducibility
dat$MCF_repo<-apply(dat[,2:9],1,FUN=function(x){sum(x>0)})
dat$Hela_repo<-apply(dat[,10:19],1,FUN=function(x){sum(x>0)})
dat$MDA1_repo<-apply(dat[,20:21],1,FUN=function(x){sum(x>0)})
dat$UHRR1_repo<-apply(dat[,22:29],1,FUN=function(x){sum(x>0)})
dat$K5621_repo<-apply(dat[,30:37],1,FUN=function(x){sum(x>0)})
dat$HEK1_repo<-apply(dat[,38:45],1,FUN=function(x){sum(x>0)})
dat$MDA2_repo<-apply(dat[,46:50],1,FUN=function(x){sum(x>0)})
dat$UHRR2_repo<-apply(dat[,51:60],1,FUN=function(x){sum(x>0)})
dat$K5622_repo<-apply(dat[,61:67],1,FUN=function(x){sum(x>0)})
dat$HEK2_repo<-apply(dat[,68:69],1,FUN=function(x){sum(x>0)})

dat1<-data.frame(label=c(rep("MCF",8),rep("Hela",10),
              rep("MDA_Bio1",2),rep("UHRR_Bio1",8),rep("K562_Bio1",8),rep("HEK_Bio1",8),
              rep("MDA_Bio2",5),rep("UHRR_Bio2",10),rep("K562_Bio2",7),rep("HEK_Bio2",2)))
#### 


## Hela marker
dat2<-dat[dat$Hela_repo>9,]
pdf("temp_fig/Hela_rep_box.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=3,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

## K562
dat2<-dat[dat$K5621_repo>7 & dat$K5622_repo>6,]
pdf("temp_fig/K562_rep_box.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=3,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

## HEK
dat2<-dat[dat$HEK1_repo>7 &dat$HEK2_repo>1,]
pdf("temp_fig/HEK_rep_box.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=3,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

## MCF
dat2<-dat[dat$MCF_repo>7,]
pdf("temp_fig/MCF_rep_box.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=3,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

## MDA
dat2<-dat[dat$MDA1_repo>1 &dat$MDA2_repo>4,]
pdf("temp_fig/MDA_rep_box.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=3,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()
### pick examples.

## MCF:23I_DOCK6___ENSG00000130158___ENST00000587656___protein_coding___protein_coding
## Hela: 11I_SGK1___ENSG00000118515___ENST00000237305___protein_coding___protein_coding;1I_FTH1___ENSG00000167996___ENST00000532601___protein_coding___protein_coding
## MDA: none
## UHRR: 11I_THBS1___ENSG00000137801___ENST00000260356___protein_coding___protein_coding
## K562: 13I_CCHCR1___ENSG00000204536___ENST00000396268___protein_coding___protein_coding
## HEK: 2I_ACTG1___ENSG00000184009___ENST00000575842___protein_coding___protein_coding
pick<-c("23I_DOCK6___ENSG00000130158___ENST00000587656___protein_coding___protein_coding",
        "11I_SGK1___ENSG00000118515___ENST00000237305___protein_coding___protein_coding",
        "1I_FTH1___ENSG00000167996___ENST00000532601___protein_coding___protein_coding",
        "11I_THBS1___ENSG00000137801___ENST00000260356___protein_coding___protein_coding",
        "1I_RPS2___ENSG00000140988___ENST00000526586___protein_coding___protein_coding",
        "2I_ACTG1___ENSG00000184009___ENST00000575842___protein_coding___protein_coding")
dat2<-dat[dat$ID%in%pick,]
pdf("temp_fig/FLEXI_box.pdf",height=8,width=11)
par(mfrow=c(2,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[grep("UHRR",tmp$label,invert = T),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,cex=2.5,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()

### pick from 3B
pick<-c("2I_RPL22L1___ENSG00000163584___ENST00000475836___protein_coding___lncRNA",
        "7I_SGK1___ENSG00000118515___ENST00000367857___protein_coding___protein_coding",
        "3I_FOS___ENSG00000170345___ENST00000303562___protein_coding___protein_coding" ,
        "1I_TIAL1___ENSG00000151923___ENST00000369086___protein_coding___protein_coding",
        "10I_ITGA5___ENSG00000161638___ENST00000293379___protein_coding___protein_coding",
        "15I_CDC45___ENSG00000093009___ENST00000407835___protein_coding___protein_coding",
        "18I_ANK1___ENSG00000029534___ENST00000520299___protein_coding___protein_coding" )
dat2<-dat[dat$ID%in%pick,]
pdf("temp_fig/FLEXI_box_2.pdf",height=8,width=11)
par(mfrow=c(3,3),cex=0.5)
for(i in 1:dim(dat2)[1]){
  tmp<-cbind(t(dat2[i,2:69]),dat1)
  tmp<-tmp[tmp$label%in%c("Hela","K562_Bio1","HEK_Bio1"),]
  colnames(tmp)[1]<-"FLEXI"
  boxplot(FLEXI~label,tmp,las=2,xlab=NA,ylab="RPM",main=dat2$ID[i],cex.main=0.5)
  stripchart(FLEXI~label,tmp,pch=19,method = "jitter",vertical = TRUE,add = TRUE)
}
dev.off()


### hetamap
library(gplots)
library(stringr)
heatPalette = colorRampPalette(c("dodgerblue4", "skyblue", "white",
                                 "goldenrod", "orangered"))(100)

dat1<-dat
rownames(dat1)<-dat1$ID
dat1<-dat1[,-1]
dat1<-dat1[,order(colnames(dat1))]
dat1<-dat1[rowSums(dat1)>0,]
dat1[dat1==0]<-2^-7
dat1<-log2(dat1)

ht<-heatmap.2(as.matrix(dat1))
pick<-rev(ht$rowInd)[1:100]
tmp1<-dat1[pick,]
rowLab<-str_split_i(rownames(tmp1),pattern = "___",i=1)
pdf("temp_fig/FLEXI_heat2.pdf",height=11,width=8.5)
heatmap.2(as.matrix(tmp1),labRow = rowLab,
          dendrogram = "none",labCol=colnames(tmp1),
          scale="none",breaks=seq(-7,1,0.08),cexRow = 0.1,
          density.info = "none",trace="none",symm=F,key=T,
          symbreaks=F,Colv=F,Rowv=T,cexCol=0.25,
          col =heatPalette)
dev.off()
