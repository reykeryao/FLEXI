rm(list=ls())
set.seed(740714)
library(tidyverse)
### download rMAT normalized results from ENCODE, rMAT results are hg19 based
DE_dir<-"DEseq2_KD/RBP51/"
### read RBP and cell line info
DE_meta<-read.table(paste0(DE_dir,"RBP51.info"),col.names = c("F_name","RBP","Cell"))
### read FLEXIs (only from 4 cell types), keep in mind this is hg38 based
FLEXI<-read.delim("all.FLEXI")
FLEXI<-FLEXI[rowSums(FLEXI[,73:76])>0,]
### get hg19 version of ENSG ID
hg19_bed<-read.table("4cell_plasma_hg19.bed",col.names=c("Chr_hg19","St_hg19","Ed_hg19","ID","Score","Strand_hg19"))
hg19_bed<-hg19_bed[,c(4,1,2,3,6)]
FLEXI<-merge(FLEXI,hg19_bed,by="ID")
FLEXI<-FLEXI[,c(1:5,7,8,78:81,23)]

FLEXI$GID_hg19<-FLEXI$GID
ENSG_dic<-read.delim(paste0(DE_dir,"liftedhg38ENSG_ID_table_to_hg19ENSG_ID"))
for (i in 1:dim(FLEXI)[1]){
  hg38ID<-FLEXI$GID[i]
  if (hg38ID%in%ENSG_dic$hg38_ENSG){
    FLEXI$GID_hg19[i]<-ENSG_dic$hg19_ENSG[ENSG_dic$hg38_ENSG==hg38ID]
  }
}
### 60 FLEXIs has different GID.
### make the same table for all intron/other short intron
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
tmp_short<-tmp[tmp$Len<=300 & !tmp$ID%in%FLEXI$ID,]

### merge long/OSI intronID to RBP
long_intron<-merge(long_intron,tmp_long[,c(1,3)],by=1,all=T)
long_intron[is.na(long_intron)]<-"."
other_short_intron<-merge(other_short_intron,tmp_short[,c(1,3)],by=1,all=T)
other_short_intron[is.na(other_short_intron)]<-"."

intersect_list<-list("FLEXI"=FLEXI[,c(1,12,7)],"OSI"=other_short_intron,"LI"=long_intron)
### loop through DE files and do fisher exact test on sig genes, and sig Up/Down genes
f_list<-list("FLEXI"="","OSI"="","LI"="")
p_col=c("gray80","black","red")
for (i in 1:dim(DE_meta)[1]){
  file_name<-DE_meta$F_name[i]
  RBP<-DE_meta$RBP[i]
  Cell<-DE_meta$Cell[i]
  ### read DE file
  DE_file<-read.delim(paste0(DE_dir,file_name))
  DE_file<-DE_file[complete.cases(DE_file),]
  ### remove log2foldchange is inf or -inf
  DE_file<-DE_file[DE_file$log2FoldChange!=Inf,]
  DE_file<-DE_file[DE_file$log2FoldChange!= -Inf,]
  DE_file<-separate(DE_file,id,into="GID",sep="\\.",extra="drop")
  ### calculate Down/up genes in total, FLEXI host genes, OSI host genes, LI host genes, the last three need to be bound by RBP
  Types<-c("FLEXI","OSI","LI")
  up<-max(0,sum(DE_file$log2FoldChange>=1 & DE_file$padj<=0.05),NA,na.rm = T)
  down<-max(0,sum(DE_file$log2FoldChange<= -1 & DE_file$padj<=0.05),NA,na.rm = T)
  sig<-up+down
  non_sig<-sum(abs(DE_file$log2FoldChange)<0 | DE_file$padj>0.05)
  Sig_test<-matrix(c(up,down,sig,non_sig),nrow = 1)
  colnames(Sig_test)<-c("Up","Down","Sig","Non_sig")
  for ( j in 1:3){
    tmp<-intersect_list[[j]]
    ### get OSI host gene contain RBP binding site
    bound_GID<-unique(tmp$GID[grep(RBP,tmp$RBP)])
    bound_DE<-DE_file[DE_file$GID%in%bound_GID,]
    up<-max(0,sum(bound_DE$log2FoldChange>=1 & bound_DE$padj<=0.05),NA,na.rm = T)
    down<-max(0,sum(bound_DE$log2FoldChange<= -1 & bound_DE$padj<=0.05),NA,na.rm = T)
    sig<-up+down
    non_sig<-sum(abs(bound_DE$log2FoldChange)<0 | bound_DE$padj>0.05)
    Sig_tmp<-rbind(Sig_test,matrix(c(up,down,sig,non_sig),nrow = 1))
    Sig_tmp[1,]<-Sig_tmp[1,]-Sig_tmp[2,]
    ### test the sig/non-sig in one of the group to total, if non-sig, skip
    ### significantly diff in number of DE genes in subgroup,all three possibilities 
    test_pV<-min(c(fisher.test(Sig_tmp[,3:4])$p.value,
                 fisher.test(Sig_tmp[,3:4],alternative = "less")$p.value,
                 fisher.test(Sig_tmp[,3:4],alternative = "greater")$p.value),na.rm = T)
    ### only care if the DE genes are significantly diff from the background
    if (test_pV>0.05){
      next
    } else {
      ### alt less is more Up DE genes, greater is more down DE genes
      test_dir<-c(fisher.test(Sig_tmp[,1:2],alternative = "less")$p.value<=0.05,
                     fisher.test(Sig_tmp[,1:2],alternative = "greater")$p.value<=0.05)
      ### add name to the list with up/down/bidirection info (U/D/B)
      if (test_dir[1] & !test_dir[2]){
        flag="U"
      } else if (!test_dir[1] & test_dir[2]) {
        flag="D"
      } else {
        flag="B"
      }
      f_name<-paste0(RBP,"_",Cell,"_",flag)
      f_list[[j]]<-c(f_list[[j]],f_name)
      ### plot FLEXI volcano plots, png format
      dat<-DE_file[,c(1,6,8)]
      dat$padj<-log10(1/dat$padj)
      ### set up type
      ### 1 : non-sig, gray
      ### 2 : sig non-bound, black
      ### 3: sig bound FLEXI, red
      dat$Type<-1
      dat$Type[dat$padj>log10(1/0.05) & abs(dat$log2FoldChange)>=1]<-2
      dat$Type[dat$Type==2 & dat$GID%in%bound_GID]<-3
      png(paste0("Figures/FigS14/",Types[j],"/",f_name,".png"),width = 900,height=900,res=300)
      par(pch=20,pty="s",cex=0.5,mai=c(0.2,0.2,0.2,0.2))
      plot(dat[dat$Type==2,c(2,3)],xlab=NA,ylab=NA,labels=F,col=p_col[2])
      points(dat[dat$Type==1,c(2,3)],col=p_col[1])
      points(dat[dat$Type==3,c(2,3)],col=p_col[3])
      abline(h=log10(1/0.05),lty=3)
      abline(v=0,lty=1)
      abline(v=1,lty=2)
      abline(v= -1,lty=2)
      dev.off()
    }
  }
}

lapply(f_list,function(x){x[x!=""]})

#$FLEXI
#[1] "RBM15_K562_B"  "RBM22_HepG2_U" "SF3A3_HepG2_B" "SF3B4_HepG2_D" "TIAL1_HepG2_U" "U2AF2_HepG2_U" "U2AF2_K562_B" 

#$OSI
#[1] "BCLAF1_HepG2_B" "BCLAF1_K562_B"  "BUD13_HepG2_B"  "BUD13_K562_B"   "DDX24_HepG2_B"  "DDX3X_HepG2_B"  "DDX3X_K562_B"   "LARP4_K562_B"   "LIN28B_HepG2_B"
#[10] "PABPN1_K562_D"  "PCBP2_K562_B"   "PPIG_HepG2_D"   "PRPF8_HepG2_U"  "RBFOX2_HepG2_B" "RBFOX2_K562_B"  "RBM15_HepG2_B"  "SF3A3_HepG2_U"  "SF3B4_HepG2_D" 
#[19] "SMNDC1_HepG2_D" "SND1_HepG2_B"   "SRSF1_HepG2_U"  "SRSF1_K562_U"   "TIA1_K562_B"    "TIAL1_HepG2_U"  "U2AF2_HepG2_U"  "UCHL5_K562_D"   "XRN2_HepG2_B"  

#$LI
#[1] "BCLAF1_HepG2_U"  "BCLAF1_K562_U"   "BUD13_HepG2_B"   "BUD13_K562_B"    "DDX24_HepG2_U"   "DDX3X_HepG2_U"   "DDX3X_K562_D"    "DKC1_HepG2_U"    "EFTUD2_HepG2_B" 
#[10] "GEMIN5_HepG2_B"  "GPKOW_HepG2_U"   "GPKOW_K562_B"    "GRWD1_HepG2_B"   "IGF2BP1_HepG2_B" "KHSRP_HepG2_U"   "KHSRP_K562_B"    "LARP4_HepG2_U"   "LARP4_K562_B"   
#[19] "LIN28B_HepG2_B"  "LSM11_HepG2_B"   "PCBP1_HepG2_B"   "PCBP2_HepG2_U"   "PCBP2_K562_U"    "PPIG_HepG2_B"    "PPIG_K562_B"     "PRPF8_HepG2_U"   "RBFOX2_HepG2_B" 
#[28] "RBM15_HepG2_U"   "SF3A3_HepG2_U"   "SF3B4_HepG2_B"   "SF3B4_K562_B"    "SMNDC1_HepG2_D"  "SND1_HepG2_B"    "SRSF1_HepG2_U"   "SRSF1_K562_U"    "TIA1_HepG2_B"   
#[37] "TIAL1_HepG2_U"   "TIAL1_K562_B"    "TRA2A_HepG2_B"   "U2AF1_K562_B"    "U2AF2_HepG2_U"   "UCHL5_HepG2_B"   "UCHL5_K562_B"    "XRN2_HepG2_B"    "YBX3_HepG2_B"   
#[46] "YBX3_K562_B"     "ZNF622_HepG2_U" 

lapply(f_list,function(x){unique(unlist(lapply(strsplit(x[x!=""],split = "_"),`[`, 1)))})

## assign value, if U/D in one cell line and B in other, use U/D, if U in one ,D in another, use B
#$FLEXI
#[1] "RBM15" "RBM22" "SF3A3" "SF3B4" "TIAL1" "U2AF2"
#.     B       U        B       D       U      U
#$OSI
#[1] "BCLAF1" "BUD13"  "DDX24"  "DDX3X"  "LARP4"  "LIN28B" "PABPN1" "PCBP2"  "PPIG"   "PRPF8"  "RBFOX2" "RBM15"  "SF3A3"  "SF3B4"  "SMNDC1" "SND1"   "SRSF1"  "TIA1"  
#.      B       B        B         B        B         B       D       B        D        U         B       B        U        D.        D.      B.        U.      B
#[19] "TIAL1"  "U2AF2"  "UCHL5"  "XRN2"  
#.       U.      U.       D.       B
#$LI
#[1] "BCLAF1"  "BUD13"   "DDX24"   "DDX3X"   "DKC1"    "EFTUD2"  "GEMIN5"  "GPKOW"   "GRWD1"   "IGF2BP1" "KHSRP"   "LARP4"   "LIN28B"  "LSM11"   "PCBP1"   "PCBP2"  
#       U.         B        U        B          U         B            B      U        B           B        U       U          B        B         B          U
#[17] "PPIG"    "PRPF8"   "RBFOX2"  "RBM15"   "SF3A3"   "SF3B4"   "SMNDC1"  "SND1"    "SRSF1"   "TIA1"    "TIAL1"   "TRA2A"   "U2AF1"   "U2AF2"   "UCHL5"   "XRN2"   
#       B          U         B        U           U       B          D         B          U         B       U         B         B           U       B        B
#[33] "YBX3"    "ZNF622" 
#       B         U

## update mRNA level change
## 0 non-sig, 1"blue",2,B:bidirectional (gray), 3 (sig bias towards LFC<0,down),red, 4 (sig bias towards LFC>0,up), skyblue
RBPinfo<-read.delim("53_RBP_info_fig4_7.txt")
RBPinfo$mRNA<-0
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("RBM15","SF3A3")]<-2
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("SF3B4")]<-3
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("RBM22","TIAL1","U2AF2")]<-4

RBPinfo$mRNA_OSI<-0
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("BCLAF1","BUD13","DDX24","DDX3X","LARP4","LIN28B","PCBP2","RBFOX2","RBM15","SND1","TIA1","XRN2")]<-2
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("PABPN1","PPIG","SF3B4","SMNDC1","UCHL5")]<-3
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("PRPF8","SF3A3","SRSF1","TIAL1","U2AF2")]<-4

RBPinfo$mRNA_LI<-0
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("BUD13","DDX3X","EFTUD2","GEMIN5","GRWD1","IGF2BP1","LIN28B","LSM11","PCBP1","PPIG","RBFOX2","SF3B4","SND1","TIA1",
                                      "TRA2A","U2AF1","UCHL5","XRN2","YBX3")]<-2
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("SMNDC1")]<-3
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("BCLAF1","DDX24","DKC1","GPKOW","KHSRP","LARP4","PCBP2","PRPF8","RBM15","SF3A3","SRSF1","TIAL1","U2AF2","ZNF622")]<-4

write.table(RBPinfo,"53_RBP_info_fig4_7.txt",quote=F,sep="\t",row.names=F)
