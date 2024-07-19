rm(list=ls())
set.seed(740714)
library(tidyverse)
### download DESeq2 results from ENCODE, hg19 based
DE_dir<-"DEseq2_KD/RBP51/"
### read RBP and cell line info
DE_meta<-read.table(paste0(DE_dir,"RBP51.info"),col.names = c("F_name","RBP","Cell"))

'''
### read FLEXIs (only from 4 cell types), keep in mind this is hg38 based
FLEXI<-read.delim("all.FLEXI")
FLEXI<-FLEXI[rowSums(FLEXI[,73:76])>0,]

### map hg19 GID to hg38 GID in all introns
ENSG_dic<-read.delim(paste0(DE_dir,"liftedhg38ENSG_ID_table_to_hg19ENSG_ID"))
#### read all intron intersect info
all_intron<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron<-unique(all_intron)
all_intron<-separate(all_intron,"FLEXI",into =c("IID","GID","TID","Gtype","TType"),sep="___",remove=F)
all_intron$GID_hg19<-all_intron$GID
for (i in 1:dim(all_intron)[1]){
  hg38ID<-all_intron$GID[i]
  if (hg38ID%in%ENSG_dic$hg38_ENSG){
    all_intron$GID_hg19[i]<-ENSG_dic$hg19_ENSG[ENSG_dic$hg38_ENSG==hg38ID]
  }
}
### all FLEXI
RBP_4cell<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell<-RBP_4cell[RBP_4cell$ID%in%FLEXI$ID,]
RBP_4cell<-unique(RBP_4cell)
RBP_4cell<-separate(RBP_4cell,"ID",into =c("IID","GID","TID","Gtype","TType"),sep="___",remove=F)
RBP_4cell<-RBP_4cell[,c("GID","RBP")]
RBP_4cell$GID_hg19<-RBP_4cell$GID
for (i in 1:dim(RBP_4cell)[1]){
  hg38ID<-RBP_4cell$GID[i]
  if (hg38ID%in%ENSG_dic$hg38_ENSG){
    RBP_4cell$GID_hg19[i]<-ENSG_dic$hg19_ENSG[ENSG_dic$hg38_ENSG==hg38ID]
  }
}
RBP_4cell<-unique(RBP_4cell)
### remove all FLEXI introns
all_intron<-all_intron[!all_intron$FLEXI%in%unique(FLEXI$ID),]
### make long introns
long_intron<-unique(all_intron[all_intron$Len>300,c("GID","RBP","GID_hg19")])
#make Other short introns
other_short_intron<-unique(all_intron[all_intron$Len<=300,c("GID","RBP","GID_hg19")])

### get GID and RBP of genes with any RBP sites in any part of the gene
G_by_RBP<-read.table("genes_with_152RBPsite.info",col.names = c("GID","RBP"))
G_by_RBP$GID_hg19<-G_by_RBP$GID
for (i in 1:dim(G_by_RBP)[1]){
  hg38ID<-G_by_RBP$GID[i]
  if (hg38ID%in%ENSG_dic$hg38_ENSG){
    G_by_RBP$GID_hg19[i]<-ENSG_dic$hg19_ENSG[ENSG_dic$hg38_ENSG==hg38ID]
  }
}


intersect_list<-list("FLEXI"=RBP_4cell,"OSI"=other_short_intron,"LI"=long_intron,"Gene"=G_by_RBP)
saveRDS(intersect_list,"FigS14.obj")
'''



intersect_list<-readRDS("FigS14.obj")
### loop through DE files and do fisher exact test on sig genes, and sig Up/Down genes
f_list<-list("FLEXI"="","OSI"="","LI"="")
Gene<-intersect_list$Gene
p_col<-c("gray","black","red")
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
  ### calculate Down/up genes in total (no RBP sites), FLEXI host genes, OSI host genes, LI host genes,
  ### the last three need to be bound by RBP
  Types<-c("FLEXI","OSI","LI")
  ### up/down/sig/non-sign in genes without the RBP binding sites in any part 
  tmp_GID<-unique(Gene$GID_hg19[Gene$RBP==RBP])
  tmp<-DE_file[!DE_file$GID%in%tmp_GID,]
  up<-max(0,sum(tmp$log2FoldChange>=1 & tmp$padj<=0.05),NA,na.rm = T)
  down<-max(0,sum(tmp$log2FoldChange<= -1 & tmp$padj<=0.05),NA,na.rm = T)
  sig<-up+down
  non_sig<-length(tmp$GID)-sig
  Sig_test<-matrix(c(up,down,sig,non_sig),nrow = 1)
  colnames(Sig_test)<-c("Up","Down","Sig","Non_sig")
  for ( j in 1:3){
    tmp<-intersect_list[[j]]
    ### get intron host gene contain RBP binding site
    bound_GID<-unique(tmp$GID_hg19[tmp$RBP==RBP])
    bound_DE<-DE_file[DE_file$GID%in%bound_GID,]
    up<-max(0,sum(bound_DE$log2FoldChange>=1 & bound_DE$padj<=0.05),NA,na.rm = T)
    down<-max(0,sum(bound_DE$log2FoldChange<= -1 & bound_DE$padj<=0.05),NA,na.rm = T)
    sig<-up+down
    non_sig<-length(bound_DE$GID)-sig
    Sig_tmp<-rbind(Sig_test,matrix(c(up,down,sig,non_sig),nrow = 1))
    #Sig_tmp[1,]<-Sig_tmp[1,]-Sig_tmp[2,]
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
      ### set up type
      ### 1 : non-sig, gray
      ### 2 : sig in genes w/o binding sites, black
      ### 3: sig bound FLEXI, red
      dat$Type<-1
      dat$Type[!dat$GID%in%tmp_GID & abs(dat$log2FoldChange)>=1 & dat$padj<=0.05]<-2
      dat$Type[dat$GID%in%bound_GID & abs(dat$log2FoldChange)>=1 & dat$padj<=0.05]<-3
      
      dat$padj<-log10(1/dat$padj)
      png(paste0("Figures/FigS14/",Types[j],"/",f_name,".png"),width = 900,height=900,res=300)
      par(pch=20,pty="s",cex=0.5,mai=c(0.2,0.2,0.2,0.2))
      plot(dat[dat$Type==1,c(2,3)],xlab=NA,ylab=NA,labels=F,col=p_col[1])
      points(dat[dat$Type==2,c(2,3)],col=p_col[2])
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
#[1] "AQR_K562_B"     "BCLAF1_HepG2_B" "BUD13_HepG2_B"  "EFTUD2_K562_B"  "LARP4_HepG2_B"  "LIN28B_HepG2_B"
#[7] "PRPF8_HepG2_U"  "RBM15_K562_B"   "SF3A3_K562_B"   "SF3B4_K562_D"   "TIAL1_HepG2_U"  "UCHL5_HepG2_B" 
#[13] "UCHL5_K562_B"  

#$OSI
#[1] "BCLAF1_HepG2_B" "BCLAF1_K562_B"  "BUD13_HepG2_B"  "BUD13_K562_B"   "EFTUD2_HepG2_U" "EFTUD2_K562_B" 
#[7] "G3BP1_HepG2_B"  "GRWD1_HepG2_B"  "KHSRP_HepG2_B"  "KHSRP_K562_B"   "LIN28B_HepG2_B" "NCBP2_HepG2_B" 
#[13] "PABPN1_K562_D"  "PPIG_HepG2_D"   "PPIG_K562_B"    "PRPF8_HepG2_U"  "RBFOX2_HepG2_B" "RBFOX2_K562_B" 
#[19] "RBM15_HepG2_B"  "RBM22_HepG2_U"  "SF3A3_HepG2_U"  "SF3B4_K562_D"   "SND1_HepG2_B"   "SRSF1_K562_U"  
#[25] "TIA1_K562_B"    "TIAL1_HepG2_U"  "U2AF2_K562_B"   "UCHL5_HepG2_B"  "UCHL5_K562_D"   "ZNF622_HepG2_U"

#$LI
#[1] "AATF_HepG2_B"    "AQR_K562_D"      "BCLAF1_HepG2_U"  "BCLAF1_K562_B"   "BUD13_HepG2_B"   "BUD13_K562_B"   
#[7] "DDX3X_HepG2_U"   "DDX55_HepG2_D"   "DKC1_HepG2_U"    "EFTUD2_HepG2_B"  "EFTUD2_K562_B"   "G3BP1_HepG2_B"  
#[13] "GEMIN5_HepG2_B"  "GEMIN5_K562_B"   "GPKOW_HepG2_U"   "GPKOW_K562_B"    "GRWD1_HepG2_B"   "IGF2BP1_HepG2_B"
#[19] "KHSRP_HepG2_U"   "KHSRP_K562_B"    "LARP4_K562_B"    "LIN28B_HepG2_B"  "LSM11_HepG2_B"   "METAP2_HepG2_B" 
#[25] "NCBP2_HepG2_B"   "PABPC4_HepG2_U"  "PCBP2_HepG2_U"   "PPIG_HepG2_B"    "PPIG_K562_B"     "PRPF4_HepG2_U"  
#[31] "PRPF8_HepG2_U"   "PUM1_HepG2_B"    "RBFOX2_HepG2_B"  "RBM15_HepG2_U"   "RBM22_HepG2_U"   "RBM22_K562_B"   
#[37] "SF3A3_HepG2_U"   "SF3B4_HepG2_B"   "SMNDC1_HepG2_D"  "SRSF1_K562_U"    "TIA1_HepG2_B"    "TIAL1_HepG2_U"  
#[43] "TIAL1_K562_B"    "TRA2A_HepG2_B"   "U2AF1_HepG2_B"   "U2AF2_HepG2_U"   "U2AF2_K562_B"    "UCHL5_HepG2_B"  
#[49] "UCHL5_K562_B"    "XRN2_HepG2_B"    "YBX3_HepG2_B"    "YBX3_K562_B"     "ZNF622_HepG2_U"  

lapply(f_list,function(x){unique(unlist(lapply(strsplit(x[x!=""],split = "_"),`[`, 1)))})

## assign value, if U/D in one cell line and B in other, B, if U in one ,D in another, use B
#$FLEXI
#[1] "AQR"    "BCLAF1" "BUD13"  "EFTUD2" "LARP4"  "LIN28B" "PRPF8"  "RBM15"  "SF3A3"  "SF3B4"  "TIAL1"  "UCHL5" 
#.     B         B        B         B       B        B       U         B        B       D        U         B     

#$OSI
# [1] "BCLAF1" "BUD13"  "EFTUD2" "G3BP1"  "GRWD1"  "KHSRP"  "LIN28B" "NCBP2"  "PABPN1" "PPIG"   "PRPF8"  "RBFOX2" "RBM15" 
#.      B         B        B         B        B       B         B       B        D        B        U       B        B
#[14] "RBM22"  "SF3A3"  "SF3B4"  "SND1"   "SRSF1"  "TIA1"   "TIAL1"  "U2AF2"  "UCHL5"  "ZNF622"
#.       U.       U.       D       B        U         B        U         B         B      U 

#$LI
#[1] "AATF"    "AQR"     "BCLAF1"  "BUD13"   "DDX3X"   "DDX55"   "DKC1"    "EFTUD2"  "G3BP1"   "GEMIN5"  "GPKOW"  
#       B        D          B         B         U         D        U           B        B        B           B
#[12] "GRWD1"   "IGF2BP1" "KHSRP"   "LARP4"   "LIN28B"  "LSM11"   "METAP2"  "NCBP2"   "PABPC4"  "PCBP2"   "PPIG" 
#       B          B         B        B          B        B          B         B          U        U         B
#[23] "PRPF4"   "PRPF8"   "PUM1"    "RBFOX2"  "RBM15"   "RBM22"   "SF3A3"   "SF3B4"   "SMNDC1"  "SRSF1"   "TIA1"  
#       U         U          B        B          U        B          U         B          D        U         B
#[34] "TIAL1"   "TRA2A"   "U2AF1"   "U2AF2"   "UCHL5"   "XRN2"    "YBX3"    "ZNF622" 
#        B         B          B        B         B        B          B         U   

## update mRNA level change
## 0 non-sig, 1"blue",2,B:bidirectional (gray), 3 (sig bias towards LFC<0,down),red, 4 (sig bias towards LFC>0,up), skyblue
RBPinfo<-read.delim("53_RBP_info_fig4_7.txt")
RBPinfo$mRNA<-0
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("AQR","BCLAF1","BUD13",
                                   "EFTUD2","LARP4","LIN28B",
                                   "RBM15","SF3A3","UCHL5")]<-2
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("SF3B4")]<-3
RBPinfo$mRNA[RBPinfo$RBP.name%in%c("PRPF8","TIAL1")]<-4

RBPinfo$mRNA_OSI<-0
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("BCLAF1","BUD13","EFTUD2",
                                       "G3BP1" ,"GRWD1" , "KHSRP",
                                       "LIN28B" ,"NCBP2","PPIG",
                                       "RBFOX2", "RBM15", "SND1",
                                       "TIA1" ,"U2AF2",  "UCHL5"
                                       )]<-2
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("PABPN1","SF3B4")]<-3
RBPinfo$mRNA_OSI[RBPinfo$RBP.name%in%c("PRPF8","RBM22","SF3A3","SRSF1","TIAL1","ZNF622")]<-4

RBPinfo$mRNA_LI<-0
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("AATF","BCLAF1" , "BUD13" ,"EFTUD2",
                                      "G3BP1" , "GEMIN5","GPKOW","GRWD1" ,
                                      "IGF2BP1", "KHSRP", "LARP4","LIN28B",
                                      "LSM11", "METAP2","NCBP2","PPIG",
                                      "PUM1" , "RBFOX2", "RBM22","SF3B4",
                                      "TIA1" ,
                                      "TIAL1", "TRA2A","U2AF1","U2AF2",
                                      "UCHL5", "XRN2" , "YBX3")]<-2
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("AQR","DDX55", "SMNDC1" )]<-3
RBPinfo$mRNA_LI[RBPinfo$RBP.name%in%c("DDX3X","DKC1","PABPC4","PCBP2","PRPF4",
                                      "PRPF8","RBM15","SF3A3","SRSF1","ZNF622" )]<-4

write.table(RBPinfo,"53_RBP_info_fig4_7_new.txt",quote=F,sep="\t",row.names=F)
