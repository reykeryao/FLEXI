rm(list=ls())
setwd("~/Desktop/tempFLEXI/RBP_KD/")
ID<-read.delim("FigS11.csv",sep=",")
NameID<-ID$ID
pcol<-c("gray","black","red")
for (i in 1:30){
  Acc<-ID[i,2]
  if (i==1 || i==30 || i==7){
    dat1<-read.delim(paste0(Acc,".tsv"))
    dat1<-dat1[,c(2,7,9)]
    dat2<-read.delim(paste0(Acc,"BoundHostExp.csv"),sep=",")
    dat2<-dat2[,c(1,7,9)]
    dat2$Type<-"3"
    dat<-merge(dat1,dat2,by=1:3,all=T)
    dat$Type[is.na(dat$Type)]<-"1"
    dat<-dat[complete.cases(dat),]
    dat$log10padj<-log10(1/dat$padj)
    dat$Type[dat$Type!="3" & dat$padj<=0.05]<-"2"
    dat$Type<-as.integer(dat$Type)
  }  else {
    dat1<-read.delim(paste0(Acc,".tsv"))
    dat1<-dat1[,c(2,6)]
    dat2<-read.delim(paste0(Acc,"BoundHostExp.csv"),sep=",",row.names=1)
    dat2<-dat2[,c(2,6)]
    dat2$Type<-"3"
    dat<-merge(dat1,dat2,by=0:2,all=T)
    dat$Type[is.na(dat$Type)]<-"1"
    dat<-dat[complete.cases(dat),]
    dat$log10padj<-log10(1/dat$padj)
    dat$Type[dat$Type!="3" & dat$padj<=0.05]<-"2"
    dat$Type<-as.integer(dat$Type)
  }
  png(paste0("tiff/",ID[i,1],".png"),width = 900,height=900,res=300)
  par(pch=20,pty="s",cex=0.5,mai=c(0.2,0.2,0.2,0.2))
  plot(dat[dat$Type==2,c(2,5)],xlab=NA,ylab=NA,labels = NA)
  points(dat[dat$Type==1,c(2,5)],col=pcol[1])
  points(dat[dat$Type==3,c(2,5)],col=pcol[3])
  dev.off()
}
