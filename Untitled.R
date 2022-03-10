
pdf("temp_fig/Fig1A_cutoff.pdf",height=4,width=8,onefile = T)
for (i in c(1,2,3,6)){
  cut_off=i
  set_1 <- as.character(dat$ID[dat$UHRR>=cut_off])
  set_2 <- as.character(dat$ID[dat$K562>=cut_off])
  set_3 <- as.character(dat$ID[dat$HEK>=cut_off])
  set_4 <- as.character(dat$ID[dat$Hela>=cut_off])
  set_5 <- as.character(dat$ID[dat$Plasma>=cut_off])
  set <- list ("UHRR"=set_1,
               "K-562"=set_2,"HEK 293T"=set_3,
               "Hela S3"=set_4,"Plasma"=set_5)
  m = make_comb_mat(set)
  print(UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
        comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)]))
}
dev.off()

i=5
pdf(paste0("temp_fig//Fig1D_cutoff",i,".pdf"),width=12,height=12)
col=c("red","orchid","black")
dat<-temp[rowMaxs(as.matrix(temp[,88:91]))>i,]
#length
pdf(NULL)
dev.control(displaylist="enable")
dat1<-dat[rowSums(dat[,88:91])>0,]
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

ggarrange(ggarrange(length.line,GC.line,MFE.line,ncol=2,nrow = 2))
dev.off()


plot(RBP_fre[RBP_fre$col==0,c(4,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.01 & (RBP_fre[,4]>=4 | RBP_fre[,3]>=4)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.01 & (RBP_fre[,4]>=4 | RBP_fre[,3]>=4)])




RBP_clus<-read.delim(paste0(cell_list[i],"_RBP_clus_froGower.txt"))
marker_text<-RBP_clus
marker_text[marker_text=="U"]<-"X"
marker_text[marker_text=="N"]<-""
marker_text[marker_text=="D"]<-"X"
for (j in 1:dim(marker_text)[2]){
  colname<-rownames(marker_text)[j]
  colname<-paste0(colname,"BF")
  marker_text[j,colname]<-""
}

RBP_clus<-data.frame(t(RBP_clus))
for (j in 1:dim(RBP_clus)[2]){
  RBP_clus[,j]<-factor(RBP_clus[,j],levels = c("D","N","U"))
}
RBP_clus_log10<-read.delim(paste0(cell_list[i],"_RBP_clus_froHeatmap.txt"))

RBP_clus<-RBP_clus[RBP_47_name]
RBP_clus<-RBP_clus[RBP_47,]

marker_text.47<-marker_text[RBP_47]
marker_text.47<-marker_text.47[RBP_47_name,]

RBP_clus_log10<-RBP_clus_log10[RBP_47]
RBP_clus_log10<-RBP_clus_log10[RBP_47_name,]

gower.dist <- daisy(t(RBP_clus_log10), metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}


#p 005 fig4B
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>0,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
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
colnames(RBP_fre)<-c("RBP","Cell","Other SHort introns","Long introns")
RBP_fre<-merge(RBP_fre,RBP[,c(1,4,5,11)],by=1)
RBP_fre$col<-(RBP_fre$Splicing.regulation+RBP_fre$Spliceosome)/3+RBP_fre$microRNA.processing
RBP_fre<-RBP_fre[,c(1:4,8)]
RBP_fre$Bpvlue<-1
RBP_fre$Cpvlue<-1
R_sum<-colSums(RBP_fre[,2:4])
for (i in 1:152){
  RBP_fre[i,6]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,3)],R_sum[c(1,2)])))$p.value
  RBP_fre[i,7]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,4)],R_sum[c(1,3)])))$p.value
}
RBP_fre[,2:4]<-data.frame(prop.table(as.matrix(RBP_fre[,2:4]),margin = 2)*100)
RBP_fre$Bpvlue<-p.adjust(RBP_fre$Bpvlue,method="fdr")
RBP_fre$Cpvlue<-p.adjust(RBP_fre$Cpvlue,method="fdr")

RBP_fre[RBP_fre$col>1,5]<-4
RBP_fre[RBP_fre$col==1,5]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,5]<-2
RBP_fre[RBP_fre$col==0,5]<-1

percent_cutoff<-1
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")

pdf("temp_fig/Fig4B.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)

plot(RBP_fre[RBP_fre$col==1,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)],
     col=col[RBP_fre$col[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4| RBP_fre[,3]>=4)]])
#subpanel
plot(RBP_fre[RBP_fre$col==1,c(3,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,3]>=1)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,3]>=1)],
     col=col[RBP_fre$col[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=1 | RBP_fre[,3]>=1)]])


plot(RBP_fre[RBP_fre$col==1,c(4,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)],
     col=col[RBP_fre$col[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)]])

plot(RBP_fre[RBP_fre$col==1,c(4,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=1 | RBP_fre[,2]>=1)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=1 | RBP_fre[,2]>=1)],
     col=col[RBP_fre$col[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=1 | RBP_fre[,2]>=1)]])

dev.off()

pdf("temp_fig/Fig4B_p005.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)

plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)])
#subpanel
plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2)])


plot(RBP_fre[RBP_fre$col==0,c(4,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)])

plot(RBP_fre[RBP_fre$col==0,c(4,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2)])

dev.off()


#FIg4B without core 6 RBPS
dat<-read.delim("all.FLEXI")
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>0,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
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
RBP_fre$Bpvlue<-1
RBP_fre$Cpvlue<-1
core_list<-c("PRPF8","SF3B4","AQR","EFTUD2","BUD13")
RBP_fre<-RBP_fre[!RBP_fre$RBP%in%core_list,]
R_sum<-colSums(RBP_fre[,2:4])
for (i in 1:147){
  RBP_fre[i,6]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,3)],R_sum[c(1,2)])))$p.value
  RBP_fre[i,7]<-fisher.test(as.matrix(rbind(RBP_fre[i,c(2,4)],R_sum[c(1,3)])))$p.value
}
RBP_fre[,2:4]<-data.frame(prop.table(as.matrix(RBP_fre[,2:4]),margin = 2)*100)
RBP_fre$Bpvlue<-p.adjust(RBP_fre$Bpvlue,method="fdr")
RBP_fre$Cpvlue<-p.adjust(RBP_fre$Cpvlue,method="fdr")

RBP_fre[RBP_fre$col>1,5]<-4
RBP_fre[RBP_fre$col==1,5]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,5]<-2
RBP_fre[RBP_fre$col==0,5]<-1

percent_cutoff<-1
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")

pdf("temp_fig/Fig4B_no6Core.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)

plot(RBP_fre[RBP_fre$col==1,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4)],
     col=col[RBP_fre$col[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4| RBP_fre[,3]>=4)]])
#subpanel
plot(RBP_fre[RBP_fre$col==1,c(3,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2)),c(3,2)],
     labels = RBP_fre$RBP[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2)],
     col=col[RBP_fre$col[RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2)]])


plot(RBP_fre[RBP_fre$col==1,c(4,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)],
     col=col[RBP_fre$col[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4)]])

plot(RBP_fre[RBP_fre$col==1,c(4,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
     ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
points(RBP_fre[RBP_fre$col==2 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==3,c(4,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col==4,c(4,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2)),c(4,2)],
     labels = RBP_fre$RBP[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2)],
     col=col[RBP_fre$col[RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2)]])

dev.off()



# FLEXIs with RPM cutoff
RPM_cufoff<-c(0.2,0.3,0.4,0.5)
all_intron_RBP<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron_RBP<-unique(all_intron_RBP)
all_intron_RBP_short<-all_intron_RBP[all_intron_RBP$Len<=300,c(1,3)]
all_intron_RBP<-all_intron_RBP[all_intron_RBP$Len>300,c(1,3)]
all_intron_RBP_short<-all_intron_RBP_short[!all_intron_RBP_short$FLEXI%in%FLEXI$ID,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)

pdf("temp_fig/RPM_RBP1.pdf",height=12,width=12)
par(mfrow=c(4,4),pty="s",pch=16)
for (i in 1:4){
  temp<-dat_CPM$ID[rowMaxs(as.matrix(dat_CPM[,33:36]))>RPM_cufoff[i]]
  temp<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%temp,]
  RBP_fre<-data.frame(table(temp$RBP))
  long_fre<-data.frame(table(all_intron_RBP$RBP))
  OtherShort_fre<-data.frame(table(all_intron_RBP_short$RBP))
  RBP_fre<-merge(RBP_fre,OtherShort_fre,by=1,all=T)
  RBP_fre<-merge(RBP_fre,long_fre,by=1,all=T)
  RBP_fre[is.na(RBP_fre)]<-0
  colnames(RBP_fre)<-c("RBP","Cell","Other SHort introns","Long introns")
  RBP_fre<-merge(RBP_fre,RBP[,c(1,4,5,11)],by=1)
  RBP_fre$col<-(RBP_fre$Splicing.regulation+RBP_fre$Spliceosome)/3+RBP_fre$microRNA.processing
  RBP_fre<-RBP_fre[,c(1:4,8)]
  RBP_fre$Bpvlue<-1
  RBP_fre$Cpvlue<-1
  R_sum<-colSums(RBP_fre[,2:4])
  for (j in 1:152){
    RBP_fre[j,6]<-fisher.test(as.matrix(rbind(RBP_fre[j,c(2,3)],R_sum[c(1,2)])))$p.value
    RBP_fre[j,7]<-fisher.test(as.matrix(rbind(RBP_fre[j,c(2,4)],R_sum[c(1,3)])))$p.value
  }
  RBP_fre[,2:4]<-data.frame(prop.table(as.matrix(RBP_fre[,2:4]),margin = 2)*100)
  RBP_fre$Bpvlue<-p.adjust(RBP_fre$Bpvlue,method="fdr")
  RBP_fre$Cpvlue<-p.adjust(RBP_fre$Cpvlue,method="fdr")
  
  
  plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",main=paste0(RPM_cufoff[i]," RPM"),
       ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
  points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
  points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
  points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=4 | RBP_fre[,3]>=4))
  if (sum(sig)>0){
    text(RBP_fre[sig,c(3,2)],
         labels = RBP_fre$RBP[sig])
  }
  #subpanel
  plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
       ylab="FLEXIs (% RBP sites)",xlab="Other short introns (% RBP sites)")
  points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
  points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
  points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-(RBP_fre$Bpvlue<=0.05 & (RBP_fre[,2]>=2 | RBP_fre[,3]>=2))
  if (sum(sig)>0){
    text(RBP_fre[sig,c(3,2)],
         labels = RBP_fre$RBP[sig])
  }
  

  plot(RBP_fre[RBP_fre$col==0,c(4,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
       ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
  points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
  points(RBP_fre[RBP_fre$col==1,c(4,2)],col="orange",cex=1.5)
  points(RBP_fre[RBP_fre$col>1,c(4,2)],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=4 | RBP_fre[,2]>=4))
  if(sum(sig)>0){
    text(RBP_fre[sig,c(4,2)],labels = RBP_fre$RBP[sig])
  }
  
  plot(RBP_fre[RBP_fre$col==0,c(4,2)],xlim=c(0,4),ylim=c(0,4),cex=1.5,bty="n",
       ylab="FLEXIs (% RBP sites)",xlab="Long introns (% RBP sites)")
  points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(4,2)],col="red",cex=1.5)
  points(RBP_fre[RBP_fre$col==1,c(4,2)],col="orange",cex=1.5)
  points(RBP_fre[RBP_fre$col>1,c(4,2)],col="skyblue",cex=1.5)
  abline(0,1,col="red")
  sig<-(RBP_fre$Cpvlue<=0.05 & (RBP_fre[,4]>=2 | RBP_fre[,2]>=2))
  if(sum(sig)>0){
    text(RBP_fre[sig,c(4,2)],labels = RBP_fre$RBP[sig])
  }

}
dev.off()

dPCR<-read.table("dPCR.txt",sep="\t",
                 col.names=c("ID","HEK_TG_RPM","HEK_TG_lm","HEK_dPCR_lm","HEK_dPCR_U12","HEK_dPCR_U14",
                                        "S3_TG_RPM","S3_TG_lm","S3_dPCR_lm","S3_dPCR_U12","S3_dPCR_U14",
                                        "UHRR_TG_RPM","UHRR_TG_lm","UHRR_dPCR_lm","UHRR_dPCR_U12","UHRR_dPCR_U14"))
dPCR[is.na(dPCR)]<-0
cell_name<-c("HEK-293T","HeLa S3","UHRR")
pdf("temp_fig/dPCR_cor.pdf",width=8,height=8)
par(mfrow=c(3,3),pty="s",pch=16)
for (i in 1:3){
  plot(dPCR[,c((i-1)*5+2,(i-1)*5+4)],main=cell_name[i],xlab="TGIRT-seq (RPM)",ylab=c("dPCR (lm)"))
  plot(dPCR[,c((i-1)*5+2,(i-1)*5+5)],main=cell_name[i],xlab="TGIRT-seq (RPM)",ylab=c("dPCR (U12)"))
  plot(dPCR[,c((i-1)*5+2,(i-1)*5+6)],main=cell_name[i],xlab="TGIRT-seq (RPM)",ylab=c("dPCR (SNORD14)"))
}
dev.off()


ECLIP<-read.delim("~/Desktop/tempFLEXI/eCLIP/eCLIP.metadata")
CONTR<-read.delim("~/Desktop/tempFLEXI/eCLIP/Control.metadata")
SIreads<-read.delim("~/Desktop/tempFLEXI/eCLIP/SIreads_by_type.txt")
eclip_sum<-c(sum(ECLIP$Total[ECLIP$Cell=="K562"]),
             sum(ECLIP$Total[ECLIP$Cell=="HepG2"]),
             sum(CONTR$Total[CONTR$Cell=="K562"]),
             sum(CONTR$Total[CONTR$Cell=="HepG2"]))/1e6
K562_FELXI<-FourCell$ID[FourCell$K562>0]
Cell_FLEXI<-FourCell$ID
library(tidyr)
#k562 scatter plot from eCLIPvs Control
pdf("temp_fig/eCLIP-scatter.pdf",width=8,height=8)
par(mfrow=c(2,2),pty="s")
tmp<-SIreads[rowSums(SIreads[,c(2,6)])>0,c(1,2,6)]
tmp$eCLIP_K562_SI<-tmp$eCLIP_K562_SI/eclip_sum[1]
tmp$Ctl_K562_SI<-tmp$Ctl_K562_SI/eclip_sum[3]
tmp[,2:3][tmp[,2:3]==0]<-2^-10
tmp[,2:3]<-log2(tmp[,2:3])
plot(tmp[!tmp$ID%in%K562_FELXI,2:3],xlim=c(-10,10),ylim=c(-10,10),pch=16,
     xlab="eCLIP (log2RPM)",ylab="Control IgG (log2RPM)",
     main="K-562 (intronic reads)")
abline(0,1,col="red")
points(tmp[tmp$ID%in%K562_FELXI,2:3],pch=16,col="red")
legend("topleft",legend = c("K-562 FLEXI","Other short introns"),
       col=c("red","black"),pch=16,bty="n")

tmp<-SIreads[rowSums(SIreads[,c(3,7)])>0,c(1,3,7)]
tmp$eCLIP_K562_FLEXI<-tmp$eCLIP_K562_FLEXI/eclip_sum[1]
tmp$Ctl_K562_FLEXI<-tmp$Ctl_K562_FLEXI/eclip_sum[3]
tmp[,2:3][tmp[,2:3]==0]<-2^-10
tmp[,2:3]<-log2(tmp[,2:3])
plot(tmp[!tmp$ID%in%K562_FELXI,2:3],xlim=c(-10,0),ylim=c(-10,0),pch=16,
     xlab="eCLIP (log2RPM)",ylab="Control IgG (log2RPM)",
     main="K-562 (FLEXI reads)")
abline(0,1,col="red")
points(tmp[tmp$ID%in%K562_FELXI,2:3],pch=16,col="red")

tmp<-SIreads[rowSums(SIreads[,c(4,8)])>0,c(1,4,8)]
tmp$eCLIP_HepG2_SI<-tmp$eCLIP_HepG2_SI/eclip_sum[2]
tmp$Ctl_HepG2_SI<-tmp$Ctl_HepG2_SI/eclip_sum[4]
tmp[,2:3][tmp[,2:3]==0]<-2^-10
tmp[,2:3]<-log2(tmp[,2:3])
plot(tmp[!tmp$ID%in%K562_FELXI,2:3],xlim=c(-10,10),ylim=c(-10,10),pch=16,
     xlab="eCLIP (log2RPM)",ylab="Control IgG (log2RPM)",
     main="HepG2 (intronic reads)")
abline(0,1,col="red")
points(tmp[tmp$ID%in%K562_FELXI,2:3],pch=16,col="red")
legend("topleft",legend = c("K-562 FLEXI","Other short introns"),
       col=c("red","black"),pch=16,bty="n")

tmp<-SIreads[rowSums(SIreads[,c(5,9)])>0,c(1,5,9)]
tmp$eCLIP_HepG2_FLEXI<-tmp$eCLIP_HepG2_FLEXI/eclip_sum[2]
tmp$Ctl_HepG2_FLEXI<-tmp$Ctl_HepG2_FLEXI/eclip_sum[4]
tmp[,2:3][tmp[,2:3]==0]<-2^-10
tmp[,2:3]<-log2(tmp[,2:3])
plot(tmp[!tmp$ID%in%K562_FELXI,2:3],xlim=c(-10,0),ylim=c(-10,0),pch=16,
     xlab="eCLIP (log2RPM)",ylab="Control IgG (log2RPM)",
     main="K-562 (FLEXI reads)")
abline(0,1,col="red")
points(tmp[tmp$ID%in%Cell_FLEXI,2:3],pch=16,col="red")
dev.off()



library(UpSetR)
library(ComplexHeatmap)
set_1 <- as.character(SIreads$ID[SIreads$eCLIP_K562_FLEXI>0])
set_2 <- as.character(SIreads$ID[SIreads$Ctl_K562_FLEXI>0])
set_3 <- as.character(SIreads$ID[SIreads$eCLIP_HepG2_FLEXI>0])
set_4 <- as.character(SIreads$ID[SIreads$Ctl_HepG2_FLEXI>0])
set_5 <- as.character(FourCell$ID[FourCell$K562>0])
set <- list ("eCLIP-K562"=set_1,
             "Control-K562"=set_2,"eCLIP-HepG2"=set_3,
             "Control-HepG2"=set_4,"K562"=set_5)
m = make_comb_mat(set)
pdf("temp_fig/upset_eclip_FLEXI.pdf")
UpSet(m,set_order=c("K562","eCLIP-K562","Control-K562","eCLIP-HepG2","Control-HepG2"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)],
      top_annotation = upset_top_annotation(m, add_numbers = T),
      right_annotation = upset_right_annotation(m, add_numbers = T))
dev.off()
library(VennDiagram)
library(ggpubr)
pdf("temp_fig/eCLIP_K562_venn.pdf",height=4,width=8)
col=c("tomato","royalblue1","greenyellow","goldenrod","orchid","black")
vennplot1 <- venn.diagram (set, filename=NULL,category.names=names(set),
                           cat.col = col[1:5],
                           fill = col[1:5],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(0,0,180,180,0),
                           main.cex=1, cat.cex = 1) 
vennplot2 <- venn.diagram (set[c(1,2,5)], filename=NULL,
                           category.names=names(set)[c(1,2,5)],
                           cat.col = col[c(1,2,5)],
                           fill = col[c(1,2,5)],
                           height = 300, width = 300, units = "px",
                           cex = 1,cat.pos=c(-30,30,180),
                           main.cex=1, cat.cex = 1) 
ggarrange(vennplot1,vennplot2,nrow = 1,ncol=2)
dev.off()


#new figS9



x<-log10(gene_for_correlation[,i])
lm.out <- lm(y ~ x)
newx = log10(c(1,0.01))
conf_interval <- data.frame(predict(lm.out, newdata=data.frame(x=newx)))
t(t(round(10^conf_interval,0)))


FLEXI_CPM<-dat[,c(1,8,89:91)]
Unfrag_total<-mapped_reads[8:10]
FLEXI_CPM[,3:5]<-t(t(FLEXI_CPM[,3:5])/Unfrag_total)

FLEXI_CPM$K562_repo<-apply(dat[,56:63],1,FUN=function(x){sum(x>0)})
FLEXI_CPM$HEK_repo<-apply(dat[,64:71],1,FUN=function(x){sum(x>0)})
FLEXI_CPM$Hela_repo<-apply(dat[,72:81],1,FUN=function(x){sum(x>0)})
FLEXI_CPM<-FLEXI_CPM[rowSums(FLEXI_CPM[,3:5])>0,]
temp<-FLEXI_CPM[,3:5]
temp[temp==0]<-2^-10
FLEXI_CPM[,3:5]<-temp
FLEXI_CPM<-FLEXI_CPM[,2:8]
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

FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562> 2^-10 &
                      FLEXI_by_GID$Hela==2^-10 &
                      FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela>2^-10 &
                      FLEXI_by_GID$K562==2^-10 &
                      FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(2,4)]),col="red")
points(log2(sig.x[,c(2,4)]),col="blue")




Repo<-dat[,c(1,32:34,29:31,26:28,35:37,38:47)]
Repo<-Repo[rowSums(Repo[,2:23])>0,]
Repo$BCH3_repro<-apply(Repo[,2:4],1,FUN=function(x){sum(x>0)})
Repo$BCH4_repro<-apply(Repo[,5:7],1,FUN=function(x){sum(x>0)})
Repo$BC3_repro<-apply(Repo[,8:10],1,FUN=function(x){sum(x>0)})
Repo$BC4_repro<-apply(Repo[,11:13],1,FUN=function(x){sum(x>0)})
Repo$MDA_repro<-apply(Repo[,14:15],1,FUN=function(x){sum(x>0)})
Repo$MCF_repro<-apply(Repo[,16:23],1,FUN=function(x){sum(x>0)})
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
Repo$Combined_H_repo<-Repo$BCH3_repro + Repo$BCH4_repro
Repo<-Repo[,c(1,15,2:7,14,8:13)]
Repo$ID<-as.character(Repo$ID)
Repo<-separate(Repo,ID,into=c("IID","GID"),sep = "___",remove = F,extra="drop")
Repo<-separate(Repo,IID,into=c("Intron","GName"),sep = "_",remove = T,extra="drop")
colnames(Repo)[6:9]<-c("PatientA_H_repro","PatientB_H_repro","PatientA_C_repro","PatientB_C_repro")


sig.x<-Repo[Repo$PatientA_C>=0.01 &
              Repo$Combined_H==2^-10 &
              Repo$PatientA_C_repro>=2,]
sig.y<-Repo[Repo$PatientB_C>=0.01 &
              Repo$Combined_H==2^-10 &
              Repo$PatientB_C_repro>=2,]



FLEXI_by_GID<-Repo[!(Repo[,13]==2^-10 & Repo[,15]==2^-10),]
plot(log2(FLEXI_by_GID[,c(13,15)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientA_H==2^-10 & FLEXI_by_GID$PatientA_C>2^-10  & 
                           FLEXI_by_GID$PatientA_C_repro>=3,c(13,15)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientA_H==2^-10 & 
                                    FLEXI_by_GID$PatientA_C>=0.01 & 
                                    FLEXI_by_GID$PatientA_C_repro>=3])
print(unique(FLEXI_by_GID$ID[FLEXI_by_GID$PatientA_H==2^-10 & 
                               FLEXI_by_GID$PatientA_C>2^-10  & 
                               FLEXI_by_GID$PatientA_C_repro>=3]))



x<-log10(gene_for_correlation[,2*i-1])
lm.out <- lm(y ~ x)
newx = log10(1)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                         level = 0.95)
10^conf_interval[1]

FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562>2^-10  &
                      FLEXI_by_GID$HEK==2^-10 &
                      FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$HEK>2^-10  &
                      FLEXI_by_GID$K562==2^-10 &
                      FLEXI_by_GID$HEK_repo>=8,]
points(log2(sig.y[,c(2,3)]),col="red")
points(log2(sig.x[,c(2,3)]),col="yellow")





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
                            data.frame("U11"=colSums(temp[grep("U11\\b|U11[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U12"=colSums(temp[grep("U12\\b|U12[A-Z]",temp$Name),2:9])))
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
temp<-gene_counts[gene_counts$Type=="snoRNA",c(4,9:16)]
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD3"=colSums(temp[grep("SNORD3\\b|SNORD3[A-Z]|U3",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD118"=colSums(temp[grep("SNORD118\\b|SNORD118[A-Z]|U8",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD13"=colSums(temp[grep("SNORD13\\b|SNORD13[A-Z]|U13",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD14"=colSums(temp[grep("SNORD14\\b|SNORD14[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("SNORD22"=colSums(temp[grep("SNORD22\\b|SNORD22[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RMRP"=colSums(gene_counts[gene_counts$Name=="RMRP",9:16])))
#mappedreads and corrected mapped reads
correct_mapped_reads<-c(715.24152,277.8381,768.43375,368.4548,666.34185,362.6804,713.77529,230.0604)
names(correct_mapped_reads)<-rownames(gene_for_correlation)
gene_for_correlation<-gene_for_correlation/correct_mapped_reads
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",2),"7SK","7SL","RNase P","MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)


#####
# all replcate FLEXIs in one, all replicates i
pdf("Figures/Fig1B.pdf",width = 10,height=6)
par(pch=16,mfcol=c(2,3),pty="s")
#FLEXIs scatter
#FELXIs between K562 and HeLa
col<-col2hex(c("red","blue"))
col<-paste0(col,"80")

FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
sig.x<-rownames(FLEXI_by_GID)[FLEXI_by_GID$K562>2^-10 & FLEXI_by_GID$K562_repo>=8]
sig.y<-rownames(FLEXI_by_GID)[FLEXI_by_GID$Hela>2^-10 & FLEXI_by_GID$Hela_repo>=10]

plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray80")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[sig.y,c(2,4)]),col=col[1])
points(log2(FLEXI_by_GID[sig.x,c(2,4)]),col=col[2])
#Unfrag scatter
#FELXIs between K562 and HeLa
sig.x<-unique(FLEXI_by_GID[sig.x,1])
sig.y<-unique(FLEXI_by_GID[sig.y,1])
sig.all<-intersect(sig.x,sig.y)
sig.x<-setdiff(sig.x,sig.all)
sig.y<-setdiff(sig.y,sig.all)
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,4]==2^-10),]
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")

plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray80")
abline(0,1,col="red")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,4)]),col=col[1])
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,4)]),col=col[2])

#FLEXIs scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,3]==2^-10 & FLEXI_CPM[,4]==2^-10),]
sig.x<-rownames(FLEXI_by_GID)[FLEXI_by_GID$HEK>2^-10 & FLEXI_by_GID$HEK_repo>=8]
sig.y<-rownames(FLEXI_by_GID)[FLEXI_by_GID$Hela>2^-10 & FLEXI_by_GID$Hela_repo>=10]

plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray80")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[sig.y,c(3,4)]),col=col[1])
points(log2(FLEXI_by_GID[sig.x,c(3,4)]),col=col[2])

#Unfrag scatter
#FELXIs between HEK and HeLa
sig.x<-unique(FLEXI_by_GID[sig.x,1])
sig.y<-unique(FLEXI_by_GID[sig.y,1])
sig.all<-intersect(sig.x,sig.y)
#sig.x<-setdiff(sig.x,sig.all)
#sig.y<-setdiff(sig.y,sig.all)
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,3]==2^-10 & Cell_counts[,4]==2^-10),]
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")

plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray80")
abline(0,1,col="red")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(3,4)]),col=col[1])
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(3,4)]),col=col[2])

#FELXIs between K562 and HEK
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
sig.x<-rownames(FLEXI_by_GID)[FLEXI_by_GID$K562>2^-10 & FLEXI_by_GID$K562_repo>=8]
sig.y<-rownames(FLEXI_by_GID)[FLEXI_by_GID$HEK>2^-10 & FLEXI_by_GID$HEK_repo>=8]

plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5),col="gray80")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[sig.y,c(2,3)]),col=col[1])
points(log2(FLEXI_by_GID[sig.x,c(2,3)]),col=col[2])

#Unfrag scatter
#FELXIs between K562 and HEK
sig.x<-unique(FLEXI_by_GID[sig.x,1])
sig.y<-unique(FLEXI_by_GID[sig.y,1])
sig.all<-intersect(sig.x,sig.y)
#sig.x<-setdiff(sig.x,sig.all)
#sig.y<-setdiff(sig.y,sig.all)
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,3]==2^-10),]
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")

plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15),col="gray80")
abline(0,1,col="red")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,3)]),col=col[1])
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,3)]),col=col[2])

dev.off()


# get copy number of agotron, mirtron snoRNA_FLEIX etc from 4 cells
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
#mappedreads and corrected mapped reads
correct_mapped_reads<-c(715.24152,277.8381,768.43375,368.4548,666.34185,362.6804,713.77529,230.0604)
names(correct_mapped_reads)<-rownames(gene_for_correlation)
gene_for_correlation<-gene_for_correlation/correct_mapped_reads
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",2),"7SK","7SL","RNase P","MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)

y<-log10(gene_for_correlation[,9])
FLEXI_for_test<-data.frame((dat[dat$Is_agotron!=".",c(90,91,88,89)])/correct_mapped_reads[c(1,3,5,7)]
FLEXI_index<-c(90,91,88,89)
#HEk,Hela,UHRR,K562
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i-1])
  lm.out <- lm(y ~ x)
  tmp<-dat[dat[,FLEXI_index[i]]>0,]
  newx<-c(range(tmp[tmp$Is_agotron!=".",FLEXI_index[i]])/correct_mapped_reads[2*i-1],
          range(tmp[tmp$Is_mirtron!=".",FLEXI_index[i]])/correct_mapped_reads[2*i-1],
          range(tmp[tmp$Has_snoRNA!=".",FLEXI_index[i]])/correct_mapped_reads[2*i-1])
  newx = log10(newx)
  conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence",
                           level = 0.95)
  print(round(10^conf_interval[,1]),0)
}  
  
# get # of ribosomal protrin

Clus<-list(group1=c("LARP4","PABPC4", "SUB1","DDX3X","RPS3","NCBP2","DDX55","METAP2"))
Clus$group2<-c("BCLAF1","UCHL5","ZNF622","TRA2A","ZNF800","GRWD1","PUM1","DDX24")
Clus$group3<-c("TIA1","TIAL1")
Clus$group4<-c("U2AF1","U2AF2","KHSRP")
Clus$group5<-c("AATF","NOLC1","DKC1","SMNDC1")
Clus$group6<-c("AGO","DICER")
Clus$group5M<-c("AATF","NOLC1","DKC1")
Ribo<-data.frame()
for (i in 1:7){
    RBP_list_sig<-Clus[[i]]
    FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_list_sig])
    FLEXI_dat<-FourCell[FourCell$ID%in%FLEXI_list,]
    temp2<-FourCell[!FourCell$ID%in%FLEXI_list,]
    cont<-matrix((c(length(FLEXI_dat$GName),length(FLEXI_dat$GName[grep("^RPS|^RPL",FLEXI_dat$GName)]),
          length(temp2$GName),length(temp2$GName[grep("^RPS|^RPL",temp2$GName)]))),2,2)
    print (cont)
    print(fisher.test(cont)$p.value)
}

#make upset in clusters
size_max=c(200,400,100,400,150,400,50)
for (i in 1:7){
  RBP_num<-length(Clus[[i]])
  RBP_name<-sort(Clus[[i]])
  RBP_set<-list()
  for (j in 1:RBP_num){
    RBP_set[[j]]<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_name[j]])
  }
  names(RBP_set)<-RBP_name
  pdf(paste0("temp_fig/Cluster",i,".pdf"))
  print(upset(fromList(RBP_set),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = size_max[i]))
  dev.off()
}


RBP_set<-list()
for (i in 1:6){
  RBP_name<-sort(Clus[[i]])
  RBP_set[[i]]<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_name])
}
names(RBP_set)<-paste0("Cluster ",c("I","II","III","IV","V","VI"))
pdf("temp_fig/Cluster.pdf")
upset(fromList(RBP_set),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 1250,
      sets =names(RBP_set),keep.order=T)
dev.off()


#check snoRNAin RBOP
AATF_ID<-RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP=="AATF"]
table(dat$Has_snoRNA[dat$ID%in%AATF_ID])
DKC1_ID<-RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP=="DKC1"]
table(dat$Has_snoRNA[dat$ID%in%DKC1_ID])
NOLC1_ID<-RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP=="NOLC1"]
table(dat$Has_snoRNA[dat$ID%in%NOLC1_ID])
SMNDC1_ID<-RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP=="SMNDC1"]
table(dat$Has_snoRNA[dat$ID%in%SMNDC1_ID])

Clu5<-list("AATF"=AATF_ID,"DKC1"=DKC1_ID,"NOLC1"=NOLC1_ID,"SMNDC1"=SMNDC1_ID)
pdf(paste0("temp_fig/Cluster",5,".pdf"))
upset(fromList(Clu5),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 150)
dev.off()




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
  
  pvalue<-wilcox.test(FLEXI_dat$Len,temp2$Len,exact = F)$p.value
  print (paste("length:",pvalue,(Len_t/rep_times)))
  #GC
  pvalue<-wilcox.test(FLEXI_dat$GC,temp2$GC,exact = F)$p.value
  print (paste("GC:",pvalue,(GC_t/rep_times)))
  #MEF
  pvalue<-wilcox.test(FLEXI_dat$MFE,temp2$MFE,exact = F)$p.value
  print (paste("MFE:",pvalue, (MFE_t/rep_times)))
}

#snoRNA upset
snoA<-FourCellPlasma[grep("^SNORA",FourCellPlasma$Has_snoRNA),]
snoD<-FourCellPlasma[grep("^SNORD",FourCellPlasma$Has_snoRNA),]
#snoRA upset
set1<-snoA$ID[grep("DKC1",snoA$RBP)]
set2<-snoA$ID[grep("NOLC1",snoA$RBP)]
set3<-snoA$ID[grep("AATF",snoA$RBP)]
snoA<-list("DKC1"=set1,"NOLC1"=set2)
pdf("temp_fig/snoRNAA_upset.pdf")
upset(fromList(snoA),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 15,
      sets =names(snoA),keep.order=T)
dev.off()
pdf("temp_fig/snoRNAD_upset.pdf")
set1<-snoD$ID[grep("DKC1",snoD$RBP)]
set2<-snoD$ID[grep("NOLC1",snoD$RBP)]
set3<-snoD$ID[grep("AATF",snoD$RBP)]
snoD<-list("AATF"=set3,"DKC1"=set1,"NOLC1"=set2)
upset(fromList(snoD),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 40,
      sets =names(snoD),keep.order=T)
dev.off()

#set2
snoA<-FourCellPlasma[grep("^SNORA",FourCellPlasma$Has_snoRNA),]
snoD<-FourCellPlasma[grep("^SNORD",FourCellPlasma$Has_snoRNA),]
#snoRA upset
set1<-snoA$ID[grep("DKC1",snoA$RBP)]
set2<-snoA$ID[grep("NOLC1",snoA$RBP)]
set3<-snoA$ID[grep("AATF",snoA$RBP)]
set4<-snoA$ID[grep("SMNDC1",snoA$RBP)]
snoA<-list("DKC1"=set1,"NOLC1"=set2,"SMNDC1"=set4)
pdf("temp_fig/snoRNAA_upset2.pdf")
upset(fromList(snoA),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 15,
      sets =names(snoA),keep.order=T)
dev.off()
pdf("temp_fig/snoRNAD_upset2.pdf")
set1<-snoD$ID[grep("DKC1",snoD$RBP)]
set2<-snoD$ID[grep("NOLC1",snoD$RBP)]
set3<-snoD$ID[grep("AATF",snoD$RBP)]
set4<-snoD$ID[grep("SMNDC1",snoD$RBP)]
snoD<-list("AATF"=set3,"DKC1"=set1,"NOLC1"=set2,"SMNDC1"=set4)
upset(fromList(snoD),set_size.show = TRUE,mb.ratio = c(0.7, 0.3),set_size.scale_max = 40,
      sets =names(snoD),keep.order=T)
dev.off()



pdf("Figures/Fig3A.pdf",width = 9,height=6)
par(pch=16,mfcol=c(2,3),pty="s")
#FLEXIs scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(2,4)]),col="#FF000050")
points(log2(sig.x[,c(2,4)]),col="#0000FF50")
#Unfrag scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,4)]),col="#0000FF80")


#FLEXIs scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,3]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(3,4)]),col="#FF000080")
points(log2(sig.x[,c(3,4)]),col="#0000FF80")
#Unfrag scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,3]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(3,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(3,4)]),col="#0000FF80")

#FELXIs between K562 and HEK
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
points(log2(sig.y[,c(2,3)]),col="#FF000080")
points(log2(sig.x[,c(2,3)]),col="#0000FF80")

#Unfrag scatter
#FELXIs between K562 and HEK
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,3)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,3)]),col="#0000FF80")
dev.off()
#remove the temporay 