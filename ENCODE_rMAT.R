rm(list=ls())
set.seed(740714)
library(tidyverse)
### download rMAT normalized results from ENCODE, rMAT results are hg19 based
MATS_dir<-"MATS/RBP51/"
### read RBP and cell line info
RBP_meta<-read.table(paste0(MATS_dir,"51RBP.folder"),col.names = c("SubDir"))

####
'''
### create bed file with all introns of SE and RI event
### make sure AGO and DICER is not included
### creat SE and RI bed file with all introns
for (i in 1:length(RBP_meta$SubDir)){
  subDir<-RBP_meta$SubDir[i]
  RBP<-strsplit(subDir,split = "-")[[1]][1]
  print (paste0("Processing ",subDir))
  SE<-read.delim(gzfile(paste0(MATS_dir,subDir,"/MATS_Norm_output/SE.MATS.JunctionCountOnly.txt.gz")))
  ## in SE the introns around (up/down) the skipped exon is the one should check: if bound FLEXI|unbound FLEXI|short intron
  ## ie. intron start/end (hg19) should be upstreamEE to exon_Start_0base, or exonEnd to downstreamES
  colNames<-c("Chr","St","Ed","ID","Symbol","Strand","PValue","FDR","IncLevelDifference")
  bedFile<-rbind(setNames(SE[,c(4,9,6,2,3,5,19:23)],colNames),
                 setNames(SE[,c(4,7,10,2,3,5,19:23)],colNames))
  ### sort bed
  bedFile<-bedFile[order(bedFile$Chr,bedFile$St),]
  write.table(bedFile,paste0(MATS_dir,subDir,"/MATS_Norm_output/SE_intron.bed"),quote=F,sep="\t",row.names=F,col.names=F)
  
  RI<-read.delim(gzfile(paste0(MATS_dir,subDir,"/MATS_Norm_output/RI.MATS.JunctionCountOnly.txt.gz")))
  ## in RI the retained intron is the one should check: if bound FLEXI|unbound FLEXI|short intron
  ## ie. intron start/end (hg19) should be upstreamEE to downstreamES
  RI<-RI[,c(4,9,10,2,3,5,19:23)]
  colNames<-c("Chr","St","Ed","ID","Symbol","Strand","PValue","FDR","IncLevelDifference")
  ### sort bed
  colnames(RI)<-colNames
  RI<-RI[order(RI$Chr,RI$St),]
  write.table(RI,paste0(MATS_dir,subDir,"/MATS_Norm_output/RI_intron.bed"),quote=F,sep="\t",row.names=F,col.names=F)
}

### gzip bed file after the loop
'''
### gennerate RBP intersection with hg19 version of the IDR_RBP.bed
### for i in *HepG2 *K562;do 
### bedtools map -s -c 4 -o distinct -a $i/MATS_Norm_output/SE_intron.bed -b $REF/TGSEQ/RBP/hg19/RBP_IDR.bed.gz | \
### awk '{if ($12=="") $12=".";print}' FS=\\t OFS=\\t | \
### bedtools intersect -s -f 0.9 -r -wao -a - -b $i/MATS_Norm_output/FLEXI.bed | cut -f 1-12,16,17 | \
### sort -u > $i/MATS_Norm_output/SE.intersect;
### bedtools map -s -c 4 -o distinct -a $i/MATS_Norm_output/RI_intron.bed -b $REF/TGSEQ/RBP/hg19/RBP_IDR.bed.gz | \
### awk '{if ($12=="") $12=".";print}' FS=\\t OFS=\\t | \
### bedtools intersect -s -f 0.9 -r -wao -a - -b $i/MATS_Norm_output/FLEXI.bed | cut -f 1-12,16,17 | \
### sort -u > $i/MATS_Norm_output/RI.intersect;
### done
'''
### pre process RI/SE.intersect, add Type (FLEXI/OSI.LI) and Bound (T/F), 
RBP_dir<-"MATS/RBP51/"
RBP_folder<-read.table(paste0(RBP_dir,"51RBP.folder"))
colNames<-c("Chr","St","Ed","ID","Symbol","Strand","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference","RBP","FLEXI","Bound")
for (i in 1:dim(RBP_folder)[1]){
  RBP<-strsplit(RBP_folder$V1[i],split = "-")[[1]][1]
  SE<-read.table(gzfile(paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/SE.intersect.gz")),col.names = colNames)
  ### re-typing column 17 "FLEXI", into $Type: FLEXI, OSI and LI
  SE$Type<-"FLEXI"
  SE$Type[(SE$Ed-SE$St)<301 & SE$FLEXI=="." ]<-"OSI"
  SE$Type[(SE$Ed-SE$St)>300 & SE$FLEXI=="."]<-"LI"
  ### recalculate bound
  SE$Bound="F"
  SE$Bound[grepl(RBP,SE$RBP)]<-"T"
  write.table(SE,paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/SE.info"),quote=F,sep="\t",row.names=F)
  
  RI<-read.table(gzfile(paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/RI.intersect.gz")),col.names = colNames)
  ### re-typing column 17 "FLEXI", into $Type: FLEXI, OSI and LI
  RI$Type<-"FLEXI"
  RI$Type[(RI$Ed-RI$St)<301 & RI$FLEXI=="." ]<-"OSI"
  RI$Type[(RI$Ed-RI$St)>300 & RI$FLEXI=="."]<-"LI"
  ### recalculate bound
  RI$Bound="F"
  RI$Bound[grepl(RBP,RI$RBP)]<-"T"
  write.table(RI,paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/RI.info"),quote=F,sep="\t",row.names=F)
}

## gzip info file after the loop
'''
RBP_dir<-"MATS/RBP51/"
RBP_folder<-read.table(paste0(RBP_dir,"51RBP.folder"))
Types<-c("FLEXI","OSI","LI")
###plot SE
for (j in 1:3){
  pdf_name<-paste0("Figures/FigS15/",Types[j],"_SE_ecdf.pdf")
  pdf(pdf_name,width=8,height=12)
  par(mfrow=c(6,4),pty="s",mar=c(2,2,3,2))
  SE_list<-""
  for (i in 1:dim(RBP_folder)[1]){
    RBP<-strsplit(RBP_folder$V1[i],split = "-")[[1]][1]
    SE<-read.delim(gzfile(paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/SE.info.gz")))
    ### make ecdf of sub types
    SE<-SE[SE$Type==Types[j],]
    pV<-ks.test(SE$IncLevelDifference[SE$Bound],SE$IncLevelDifference[!SE$Bound])$p.value
    if (pV<=0.05){
      SE_list<-c(SE_list,strsplit(RBP_folder$V1[i],"-")[[1]][1])
      dens = split(SE, SE$Bound) %>% 
        map_df(function(d) {
          dens = density(d$IncLevelDifference, adjust=1, from=min(SE$IncLevelDifference) - 0.05*diff(range(SE$IncLevelDifference)), 
                         to=max(SE$IncLevelDifference) + 0.05*diff(range(SE$IncLevelDifference)))
          data.frame(x=dens$x, y=dens$y, cd=cumsum(dens$y)/sum(dens$y), group=d$Bound[1])
        })
      plot(dens$cd[dens$group]~dens$x[dens$group],type="l",col="red",axes=F,main=RBP_folder$V1[i],xlab=NA,ylab=NA,xlim=c(-0.5,0.5),ylim=c(0,1))
      axis(1,at=seq(-0.5,0.5,0.25),labels = NA)
      axis(2,at=seq(0,1,0.2),labels = NA)
      lines(dens$cd[!dens$group]~dens$x[!dens$group],col="gray50")
    }
  }
  par(mar=c(5,5,2,2))
  plot(x=NA,y=NA,xlim=c(-0.5,0.5),ylim=c(0,1),axes=F,xlab="IncLevelDifference",ylab="ECDF")
  axis(1,at=seq(-0.5,0.5,0.25))
  axis(2,at=seq(0,1,0.2),las=2)
  legend("center",legend = c("Bound","Unbound"),col=c("red","gray50"),lty=1,bty="n",title = "Skipped exon")
  dev.off()
  if(j==1){
    SE_name<-list("FLEXI"=unique(SE_list[SE_list!=""]))
  } else {
    SE_name<-c(SE_name,setNames(list("FLEXI"=unique(SE_list[SE_list!=""])),Types[j]))
  }
}
## plot RI
for (j in 1:3){
  pdf_name<-paste0("Figures/FigS15/",Types[j],"_RI_ecdf.pdf")
  pdf(pdf_name,width=8,height=12)
  par(mfrow=c(6,4),pty="s",mar=c(2,2,3,2))
  RI_list<-""
  for (i in 1:dim(RBP_folder)[1]){
    RBP<-strsplit(RBP_folder$V1[i],split = "-")[[1]][1]
    RI<-read.delim(gzfile(paste0(RBP_dir,RBP_folder$V1[i],"/MATS_Norm_output/RI.info.gz")))
    ### make ecdf of sub types
    RI<-RI[RI$Type==Types[j],]
    pV<-ks.test(RI$IncLevelDifference[RI$Bound],RI$IncLevelDifference[!RI$Bound])$p.value
    if (pV<=0.05){
      RI_list<-c(RI_list,strsplit(RBP_folder$V1[i],"-")[[1]][1])
      dens = split(RI, RI$Bound) %>% 
        map_df(function(d) {
          dens = density(d$IncLevelDifference, adjust=1, from=min(RI$IncLevelDifference) - 0.05*diff(range(RI$IncLevelDifference)), 
                         to=max(RI$IncLevelDifference) + 0.05*diff(range(RI$IncLevelDifference)))
          data.frame(x=dens$x, y=dens$y, cd=cumsum(dens$y)/sum(dens$y), group=d$Bound[1])
        })
      plot(dens$cd[dens$group]~dens$x[dens$group],type="l",col="red",axes=F,main=RBP_folder$V1[i],xlab=NA,ylab=NA,xlim=c(-0.5,0.5),ylim=c(0,1))
      axis(1,at=seq(-0.5,0.5,0.25),labels = NA)
      axis(2,at=seq(0,1,0.2),labels = NA)
      lines(dens$cd[!dens$group]~dens$x[!dens$group],col="gray50")
    }
  }
  par(mar=c(5,5,2,2))
  plot(x=NA,y=NA,xlim=c(-0.5,0.5),ylim=c(0,1),axes=F,xlab="IncLevelDifference",ylab="ECDF")
  axis(1,at=seq(-0.5,0.5,0.25))
  axis(2,at=seq(0,1,0.2),las=2)
  legend("center",legend = c("Bound","Unbound"),col=c("red","gray50"),lty=1,bty="n",title = "Retained intron")
  dev.off()
  if(j==1){
    RI_name<-list("FLEXI"=unique(RI_list[RI_list!=""]))
  } else {
    RI_name<-c(RI_name,setNames(list("FLEXI"=unique(RI_list[RI_list!=""])),Types[j]))
  }
}
### previous results
print(SE_name)
#$FLEXI
#[1] "AATF"   "AQR"    "BUD13"  "DKC1"   "EFTUD2" "KHSRP"  "LARP4"  "LIN28B" "NOLC1"  "PABPC4" "PCBP1"  "PCBP2"  "PRPF8"  "RBM22"  "RPS3"   "SF3A3"  "SF3B4"  "SMNDC1"
#[19] "SND1"   "SUB1"   "TIA1"   "TIAL1"  "U2AF1"  "XRN2"   "ZNF800"
#$OSI
#[1] "AATF"    "AQR"     "BCLAF1"  "BUD13"   "DDX3X"   "EFTUD2"  "G3BP1"   "IGF2BP1" "KHSRP"   "LARP4"   "LIN28B"  "LSM11"   "NCBP2"   "PCBP2"   "PPIG"    "PRPF8"  
#[17] "PUM1"    "RBM15"   "RBM22"   "RPS3"    "SF3A3"   "SF3B4"   "SND1"    "SUB1"    "U2AF1"   "U2AF2"   "XRN2"   
#$LI
#[1] "AATF"    "AQR"     "BCLAF1"  "BUD13"   "DDX24"   "DDX3X"   "DDX55"   "DKC1"    "EFTUD2"  "FXR2"    "G3BP1"   "GEMIN5"  "GPKOW"   "GRWD1"   "IGF2BP1" "KHSRP"  
#[17] "LARP4"   "LIN28B"  "METAP2"  "NCBP2"   "NOLC1"   "PABPC4"  "PCBP1"   "PCBP2"   "PPIG"    "PRPF4"   "PRPF8"   "PUM1"    "RBFOX2"  "RBM15"   "RBM22"   "RPS3"   
#[33] "SF3A3"   "SF3B4"   "SMNDC1"  "SND1"    "SUB1"    "TIA1"    "TIAL1"   "TRA2A"   "U2AF1"   "U2AF2"   "UCHL5"   "XRN2"    "YBX3"    "ZNF622"  "ZNF800" 
print (RI_name)
#$FLEXI
#[1] "AQR"    "DKC1"   "EFTUD2" "LIN28B" "LSM11"  "METAP2" "RBFOX2" "RBM15"  "RPS3"   "SF3B4"  "SND1"   "TIA1"   "U2AF1"  "XRN2"   "ZNF622" "ZNF800"
#$OSI
#[1] "AQR"     "DDX24"   "IGF2BP1" "KHSRP"   "LIN28B"  "PABPC4"  "RBM15"   "RPS3"    "SF3A3"   "SF3B4"   "YBX3"    "ZNF622" 
#$LI
#[1] "AQR"     "BUD13"   "DDX24"   "DDX3X"   "DKC1"    "GRWD1"   "IGF2BP1" "KHSRP"   "LARP4"   "LSM11"   "METAP2"  "PABPC4"  "PCBP1"   "PPIG"    "PRPF8"   "RBM15"  
#[17] "RBM22"   "RPS3"    "SF3B4"   "SMNDC1"  "SND1"    "SUB1"    "TIA1"    "U2AF1"   "UCHL5"   "YBX3"   

### update SE RI results in 53 RBP table
RBPinfo<-read.delim("53_RBP_info_fig4_7.txt")
RBPinfo$AltSPl_SE<-0
RBPinfo$AltSPl_SE[RBPinfo$RBP.name%in%SE_name[[1]]]<-1
RBPinfo$AltSPl_RI<-0
RBPinfo$AltSPl_RI[RBPinfo$RBP.name%in%RI_name[[1]]]<-1

RBPinfo$AltSPl_SE_OSI<-0
RBPinfo$AltSPl_SE_OSI[RBPinfo$RBP.name%in%SE_name[[2]]]<-1
RBPinfo$AltSPl_RI_OSI<-0
RBPinfo$AltSPl_RI_OSI[RBPinfo$RBP.name%in%RI_name[[2]]]<-1

RBPinfo$AltSPl_SE_LI<-0
RBPinfo$AltSPl_SE_LI[RBPinfo$RBP.name%in%SE_name[[3]]]<-1
RBPinfo$AltSPl_RI_LI<-0
RBPinfo$AltSPl_RI_LI[RBPinfo$RBP.name%in%RI_name[[3]]]<-1

write.table(RBPinfo,"53_RBP_info_fig4_7.txt",quote=F,sep="\t",row.names=F)
### curate position of RBP site (mid-point) in excel
### add mRNA (DE) changes in excel too
