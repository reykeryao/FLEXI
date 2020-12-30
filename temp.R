library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(stringi)
library(matrixStats)
library(ComplexHeatmap)
library(tidyr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggfortify)
library(UpSetR)
pdf("RBP_by_category.pdf",width=8,height=12)
par(mfrow=c(3,1),mar = c(3,15,2.5,2))
pick<-c(5,96,151,102,12,9,32,94,113,111,117,68,131,134,119,130,132,135,150,11,136,47,
        53,101,145,93,105,29,74,87,30,81,124,49,69,85,76,2,24,20,146,48,71,19,91,86,63,103,118,51,100,72,88)
image(1:53,1:20,as.matrix(RBP[pick,22:3]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(2,las=2,at = 1:20,labels = colnames(RBP)[22:3],tick = FALSE)

image(1:53,1:14,as.matrix(RBP[pick,36:23]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(2,las=2,at = 1:14,labels = colnames(RBP)[36:23],tick = FALSE)
image(1:53,1:7,as.matrix(RBP[pick,43:37]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(2,las=2,at = 1:7,labels = colnames(RBP)[43:37],tick = FALSE)
axis(1,las=2,at = 1:53,labels = RBP$RBP.name[pick],tick = FALSE)
dev.off()


sp_name<-RBP$RBP.name[RBP$RBP.name%in%RBP_list & (RBP$Spliceosome==1 | RBP$Splicing.regulation==1)]
Fun[Fun$Name%in%sp_name,]
mi_name<-RBP$RBP.name[RBP$RBP.name%in%RBP_list & (RBP$microRNA.processing==1)]
Fun[Fun$Name%in%mi_name,]

test<-RBP$RBP.name[RBP$RBP.name%in%RBP_list & RBP$Splicing.regulation==0 & RBP$Spliceosome==0 &RBP$microRNA.processing==0]
Fun[Fun$Name%in%test,]

sig_RBP<-RBP[RBP$RBP.name%in%RBP_list,]
apply(sig_RBP[,4:22],2,sum)

apply(RBP[RBP$RBP.name%in%test,4:48],2,sum)
RBP[RBP$RBP.name%in%test,c(1,40)]
RBP[RBP$RBP.name==test[7],1:22]

#module for upset
RBP_wish_lis<-c("KHSRP", "TIAL1", "TRA2A", "U2AF2","BCLAF1")

postscript(paste0("Figures/Upset_lowGC.eps"),height=4,width=8)
set_1 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[1],FourCellPlasma$RBP)])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set_4 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[4],FourCellPlasma$RBP)])
set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
#set_6 <- as.character(dat$ID[dat$MCF7>=cut_off])
set <- list ("KHSRP"=set_1,"TIAL1"=set_2,"TRA2A"=set_3,"U2AF2"=set_4,"BCLAF1"=set_5)
upset(fromList(set))
dev.off()

#upset for AQR, BUD3, EFTUD, PPIG, PRPF8, SF3B4

RBP_wish_lis<-c("AQR", "BUD3", "EFTUD", "PPIG", "PRPF8", "SF3B4")
postscript("Figures/Upset_6RBPs.eps",height=4,width=8)
set_1 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[1],FourCellPlasma$RBP)])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set_4 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[4],FourCellPlasma$RBP)])
set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
set_6 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[6],FourCellPlasma$RBP)])
set <- list ("AQR"=set_1,"BUD3"=set_2,"EFTUD"=set_3,
             "PPIG"=set_4,"PRPF8"=set_5,"SF3B4"=set_6)
upset(fromList(set))
dev.off()


RBP_wish_lis<-c("AATF", "DKC1", "NOLC1")
postscript("Figures/Upset_3RBPs.eps",height=4,width=8)
set_1 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[1],FourCellPlasma$RBP)])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set <- list ("AATF"=set_1,"DKC1"=set_2,"NOLC1"=set_3)
upset(fromList(set))
dev.off()

RBP_wish_lis<-c("DICER","NCBP2","AATF","SUB1","DDX55","RBM15","XRN2","DDX3X","NOLC1",
                "METAP2","PCBP1","PRPF4","DKC1","PABPC4","LARP4","GEMIN5","YBX3","RPS3")
postscript("Figures/Upset_18RBPs.eps",height=8,width=16)
set_1 <- as.character(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set_4 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[4],FourCellPlasma$RBP)])
set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
set_6 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[6],FourCellPlasma$RBP)])
set_7 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[7],FourCellPlasma$RBP)])
set_8 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[8],FourCellPlasma$RBP)])
set_9 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[9],FourCellPlasma$RBP)])
set_10 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[10],FourCellPlasma$RBP)])
set_11<- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[11],FourCellPlasma$RBP)])
set_12<- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[12],FourCellPlasma$RBP)])
set_13 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[13],FourCellPlasma$RBP)])
set_14<- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[14],FourCellPlasma$RBP)])
set_15 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[15],FourCellPlasma$RBP)])
set_16 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[16],FourCellPlasma$RBP)])
set_17 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[17],FourCellPlasma$RBP)])
set_18 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[18],FourCellPlasma$RBP)])
set <- list (set_1,set_2,set_3,set_4,set_5,set_6,set_7,set_8,set_9,
             set_10,set_11,set_12,set_13,set_14,set_15,set_16,set_17,set_18)
names(set)<-RBP_wish_lis
upset(fromList(set),nsets = 18,nintersects = 100)
dev.off()


RBP_wish_lis<-c("GRWD1","BCLAF1","U2AF2","DDX24","SRSF1","TRA2A","U2AF1","UCHL5",
                "ZNF622","PPIG","TIA1","FXR2")
postscript("Figures/Upset_12RBPs.eps",height=6,width=12)
set_1 <- as.character(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set_4 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[4],FourCellPlasma$RBP)])
set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
set_6 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[6],FourCellPlasma$RBP)])
set_7 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[7],FourCellPlasma$RBP)])
set_8 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[8],FourCellPlasma$RBP)])
set_9 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[9],FourCellPlasma$RBP)])
set_10 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[10],FourCellPlasma$RBP)])
set_11<- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[11],FourCellPlasma$RBP)])
set_12<- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[12],FourCellPlasma$RBP)])
set <- list (set_1,set_2,set_3,set_4,set_5,set_6,set_7,set_8,set_9,
             set_10,set_11,set_12)
names(set)<-RBP_wish_lis
upset(fromList(set),nsets = 12,nintersects = 100)
dev.off()

RBP_wish_lis<-c("G3BP1","LIN28B","GPKOW","RBFOX2","RBM5","PCBP2","LSM11","SND1","IGF2BP1","PABPN1")
postscript("Figures/Upset_10RBPs.eps",height=6,width=12)
set_1 <- as.character(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
set_2 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[2],FourCellPlasma$RBP)])
set_3 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[3],FourCellPlasma$RBP)])
set_4 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[4],FourCellPlasma$RBP)])
set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
set_6 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[6],FourCellPlasma$RBP)])
set_7 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[7],FourCellPlasma$RBP)])
set_8 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[8],FourCellPlasma$RBP)])
set_9 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[9],FourCellPlasma$RBP)])
set_10 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[10],FourCellPlasma$RBP)])
set <- list (set_1,set_2,set_3,set_4,set_5,set_6,set_7,set_8,set_9,
             set_10)
names(set)<-RBP_wish_lis
upset(fromList(set),nsets = 10,nintersects = 100)
dev.off()

Black_lis<-c("ZNF800",
                "SUB1","LARP4","PABPC4","DKC1","METAP2","NOLC1","RPS3","DDX3X",
                "XRN2","DDX55","AATF","YBX3",
                "ZNF622","UCHL5","FXR2","DDX24","BCLAF1","GRWD1",
                "PABPN1","IGF2BP1","LSM11","G3BP1")
Splicesome<-c("PRPF8","AQR","EFTUD2","BUD13","SF3B4")
MIR<-c("LIN28B","PUM1","SND1","RBFOX2","KHSRP")
RBP_info_4cell<-RBP_info[RBP_info$ID%in%FourCellPlasma$ID,]
set_1 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome])
set_2 <- as.character(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
set_3 <- as.character(FourCellPlasma$ID[FourCellPlasma$AGO_CCR!="."])
set_4 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%MIR])
union_set1_4<-union(set_1,union(set_2,union(set3,set_4)))

plot_l <- list()
set_5<-unique(FourCellPlasma$ID[FourCellPlasma$AGO_CCR!="."])
pec<-round(100*length(setdiff(set_5,set_1))/length(set_5),0)
pec_dat<-data.frame(Name="AGO1-4",Pec=pec)
set_5<-unique(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."])
pec<-round(100*length(setdiff(set_5,set_1))/length(set_5),0)
pec_dat<-rbind(pec_dat,data.frame(Name="DICER",Pec=pec))
for (i in 1:23){
  set_5 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Black_lis[i]])
  set <- list (set_1,set_2,set_3,set_4,set_5)
  names(set)<-c("Splicesome","DICER","AGO1-4","miRNA related",Black_lis[i])
  vp <-upset(fromList(set),keep.order=T,number.angles = 45,set_size.show = T,
        mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",sets=names(set))
  pec<-round(100*length(setdiff(set_5,set_1))/length(set_5),0)
  pec_dat<-rbind(pec_dat,data.frame(Name=Black_lis[i],Pec=pec))
# grid.text(paste0(Black_lis[i]," (", pec, "%)"),x = 0.65, y=0.95, gp=gpar(fontsize=12))
  vp <- plot_grid(NULL, vp$Main_bar, vp$Sizes, vp$Matrix,
                             nrow=2, align='hv', rel_heights = c(3,1),
                             rel_widths = c(2,3))
  plot_l[[Black_lis[i]]] <-vp
}

pdf("Figures/23RBPs.pdf",onefile = T,width=16,height=24)
pec_dat<-pec_dat[order(pec_dat$Pec,decreasing = T),]
pec_dat$Title<-paste0(pec_dat$Name," (",pec_dat$Pec,"%)")
for (i in 0:3){
  if (i<3){
    grid.arrange(arrangeGrob(plot_l[[pec_dat$Name[6*i+1]]],top=pec_dat$Title[6*i+1]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+2]]],top=pec_dat$Title[6*i+2]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+3]]],top=pec_dat$Title[6*i+3]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+4]]],top=pec_dat$Title[6*i+4]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+5]]],top=pec_dat$Title[6*i+5]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+6]]],top=pec_dat$Title[6*i+6]), 
                 ncol = 2,nrow=3)
  } else {
    grid.arrange(arrangeGrob(plot_l[[pec_dat$Name[6*i+1]]],top=pec_dat$Title[6*i+1]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+2]]],top=pec_dat$Title[6*i+2]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+3]]],top=pec_dat$Title[6*i+3]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+4]]],top=pec_dat$Title[6*i+4]),
                 arrangeGrob(plot_l[[pec_dat$Name[6*i+5]]],top=pec_dat$Title[6*i+5]),
                 ncol = 2,nrow=3)
  }
}
dev.off()

Splicesome<-c("SF3B4","PRPF8","EFTUD2","BUD13","AQR")

set_1 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[1]])
set_2 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[2]])
set_3 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[3]])
set_4 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[4]])
set_5 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%Splicesome[5]])
RBP_30<-c("DICER","AGO1-4","BCLAF1","FXR2","PABPN1","IGF2BP1","ZNF622","G3BP1")
plot_l <- list()
for (i in 1:8){
  if (i==1){
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
  } else if (i==2){
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
pdf("Figures/RBPs_30.pdf",width=16,height=22)
grid.arrange(arrangeGrob(plot_l[[RBP_30[1]]],top=RBP_30[1]),
             arrangeGrob(plot_l[[RBP_30[2]]],top=RBP_30[2]),
             arrangeGrob(plot_l[[RBP_30[3]]],top=RBP_30[3]),
             arrangeGrob(plot_l[[RBP_30[4]]],top=RBP_30[4]),
             arrangeGrob(plot_l[[RBP_30[5]]],top=RBP_30[5]),
             arrangeGrob(plot_l[[RBP_30[6]]],top=RBP_30[6]),
             arrangeGrob(plot_l[[RBP_30[7]]],top=RBP_30[7]),
             arrangeGrob(plot_l[[RBP_30[8]]],top=RBP_30[8]),
             ncol = 2,nrow=4)
dev.off()
#other 5
RBP_5<-c("DDX3X","LARP4","RPS3","YBX3","DDX55")
plot_l <- list()
for (i in 1:5){
    set_6 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_5[i]])
    set <- list (set_1,set_2,set_3,set_4,set_5,set_6)
    names(set)<-c(Splicesome,RBP_5[i])
    vp <-upset(fromList(set),keep.order=T,set_size.show = T,
               mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",mainbar.y.max=600,
               sets=names(set),nintersects = 1000,nsets = 6)
    vp <- plot_grid(NULL, vp$Main_bar, vp$Sizes, vp$Matrix,
                    nrow=2, align='hv', rel_heights = c(3,1),
                    rel_widths = c(1,3))
    plot_l[[RBP_5[i]]] <-vp
}
pdf("Figures/RBPs_5.pdf",width=16,height=22)
grid.arrange(arrangeGrob(plot_l[[RBP_5[1]]],top=RBP_5[1]),
             arrangeGrob(plot_l[[RBP_5[2]]],top=RBP_5[2]),
             arrangeGrob(plot_l[[RBP_5[3]]],top=RBP_5[3]),
             arrangeGrob(plot_l[[RBP_5[4]]],top=RBP_5[4]),
             arrangeGrob(plot_l[[RBP_5[5]]],top=RBP_5[5]),
             ncol = 2,nrow=4)
dev.off()
#3 snoRNA RBPs
pdf("Figures/RBPs_3.pdf",width=8,height=5)
RBP_3<-c("AATF","DKC1","NOLC1")
set_6 <- union(set_1,union(set_2,union(set_3,union(set_4,set_5))))
set_7 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[1]])
set_8 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[2]])
set_9 <- unique(RBP_info_4cell$ID[RBP_info_4cell$RBP%in%RBP_3[3]])
set <- list (set_6,set_7,set_8,set_9)
names(set)<-c("Splicesome RBP",RBP_3)
upset(fromList(set),keep.order=T,set_size.show = T,
               mainbar.y.label = "FLEXI RNAs",sets.x.label = "FLEXI RNAs",
               sets=names(set),nintersects = 1000,nsets = 4)
dev.off()



#
FourCellPlasma_snoRNA<-FourCellPlasma[FourCellPlasma$Has_snoRNA!=".",]
postscript(paste0("Figures/Upset_snoRNA.eps"),height=4,width=8)
set_1 <- as.character(FourCellPlasma_snoRNA$ID[grep(RBP_wish_lis[1],FourCellPlasma_snoRNA$RBP)])
set_2 <- as.character(FourCellPlasma_snoRNA$ID[grep(RBP_wish_lis[2],FourCellPlasma_snoRNA$RBP)])
set_3 <- as.character(FourCellPlasma_snoRNA$ID[grep(RBP_wish_lis[3],FourCellPlasma_snoRNA$RBP)])
set_4 <- as.character(FourCellPlasma_snoRNA$ID[grep(RBP_wish_lis[4],FourCellPlasma_snoRNA$RBP)])
#set_5 <- as.character(FourCellPlasma$ID[grep(RBP_wish_lis[5],FourCellPlasma$RBP)])
#set_6 <- as.character(dat$ID[dat$MCF7>=cut_off])
set <- list ("DKC1"=set_1,"NOLC1"=set_2,"AATF"=set_3,"AQR"=set_4)
upset(fromList(set))
dev.off()

FourCellPlasma_snoRNA<-FourCellPlasma[FourCellPlasma$Has_snoRNA!=".",]
sno=c()
for (i in 1:53){
  if (length(grep(RBP_list[i],FourCellPlasma_snoRNA$RBP))) {
    sno=c(sno,RBP_list[i])} 
}
FourCellPlasma_nonsnoRNA<-FourCellPlasma[FourCellPlasma$Has_snoRNA==".",]
nonSno=c()
for (i in 1:53){
  if (length(grep(RBP_list[i],FourCellPlasma_nonsnoRNA$RBP))) {
    nonSno=c(nonSno,RBP_list[i])} 
  }

set <- list ("snoRNA FLEXI"=sno,"Others"=nonSno)
vennplot1<-venn.diagram (set, category.names=names(set),
                           cat.col = col[1:2], 
                           fill = col[1:2],
                           height = 300, width = 300, units = "px",
                           cex = 1,filename=NULL,
                           main.cex=1, cat.cex = 1) 
unlink("*.log")
ggarrange(vennplot1)
dev.off()

pdf("Figures/Density_by_phastCon.pdf",width=12,height=12)
par(mfrow=c(2,2))
#length
FourCellPlasma<-read.delim("4cell_plasma_FLEXI.tsv")
High_p<-FourCellPlasma[FourCellPlasma$PhastCon30>=0.8,]
Low_p<-FourCellPlasma[FourCellPlasma$PhastCon30<0.8,]
plot(density(High_p$Len),bty="n",xlim=c(0,350),lwd=1.5,
     ylim=c(0,0.015),main=NA,xlab="Intron length (bp)",col="red")
lines(density(Low_p$Len),xlim=c(0,350),lwd=1.5,col="black")
legend(150,0.012,lty=c(1,1),lwd=1.5,col=c("black","red"),
       legend = c("phastCons < 0.8","phastCons >= 0.8"),bty="n")
#GC
plot(density(High_p$GC),bty="n",xlim=c(0,100),lwd=1.5,
     ylim=c(0,0.05),main=NA,xlab="GC (%)",col="red")
lines(density(Low_p$GC),xlim=c(0,100),lwd=1.5,col="black")

#MEF
plot(density(High_p$MFE),bty="n",xlim=c(-150,0),lwd=1.5,
     ylim=c(0,0.03),main=NA,xlab="Minimal free energy (MFE; kcal/mol)",col="red")
lines(density(Low_p$MFE),xlim=c(-150,0),lwd=1.5,col="black")
dev.off()


table(FourCellPlasma$Len[FourCellPlasma$PhastCon30>=0.80 & FourCellPlasma$GC<50]%%3)
table(FourCellPlasma$Len[FourCellPlasma$PhastCon30>=0.80]%%3)

RBP_wish_lis<-c("KHSRP", "TIAL1", "TRA2A", "U2AF2","BCLAF1")
AU_rich<-FourCellPlasma[grep(paste(RBP_wish_lis,collapse="|"),FourCellPlasma$RBP),]
table(AU_rich$Len%%3)
BCL_seq<-FourCellPlasma$Seq[grep("BCLAF1",FourCellPlasma$RBP)]
BCL_pattern<-"GGAGGCUGGGGC"
BCL_seq[agrep(BCL_pattern,BCL_seq)]
length(FourCellPlasma$Has_snoRNA[grep("AATF",FourCellPlasma$RBP)])

snoRNA_list<-unique(FourCellPlasma$Has_snoRNA)[2:37]

for (i in c(1:36)){
  print (1e6*sum(rowSums(gene_counts[gene_counts$Name==snoRNA_list[i],26:59]))/2743033649)
}

pdf("Figures/RBP_4cell_log-all.pdf",height=6,width=12)
Fun<-read.table("4cell_plasma_RBP_by_FLEXI.counts",col.names=c("Name","RBP_by_FLEXI"))
Fun<-Fun[order(Fun$RBP_by_FLEXI,decreasing = T),]
barplot(Fun$RBP_by_FLEXI,ylim=c(1,10000),log="y",
        names.arg =Fun$Name,las=2)
dev.off()

pdf("Figures/RBP_4cell_log_over30.pdf",height=6,width=12)
barplot(Fun$RBP_by_FLEXI[Fun$RBP_by_FLEXI>=30],ylim=c(1,10000),log="y",
        names.arg =Fun$Name[Fun$RBP_by_FLEXI>=30],las=2)
dev.off()

AATF_fleix<-FourCellPlasma[grep("AATF",FourCellPlasma$RBP) ,c(1:7,93)]
unique(AATF_fleix$GName[AATF_fleix$Has_snoRNA=="."])


temp<-c("G3BP1","LIN28B","GPKOW","RBFOX2","RBM5","PCBP2","LSM11","SND1","IGF2BP1","PABPN1",
        "GRWD1","BCLAF1","U2AF2","DDX24","PPIG","SRSF1","TRA2A","TIA1","U2AF1","FXR2",
        "UCHL5","ZNF622","DICER","YBX3","NCBP2","AATF","DDX55","RBM15","XRN2","DDX3X",
        "RPS3","NOLC1","METAP2","PCBP1","PRPF4","DKC1","PABPC4","LARP4","GEMIN5","SUB1",
        "RBM22","TIAL1","PUM1","ZNF800","KHSRP","SF3A3","SMNDC1","AGO","PRPF8","AQR",
        "EFTUD2","BUD13","SF3B4")

pick<-c(48,71,51,100,103,88,72,118,63,86,
        53,11,135,19,91,119,132,130,134,47,
        136,150,29,146,76,2,24,101,145,20,
        105,81,74,87,93,30,85,69,49,124,102,131,96,151,68,111,117,5,94,9,32,12,113)

dat<-read.delim("all.FLEXI")
dat_CPM<-dat
for (i in 82:92){
  dat_CPM[,i]<-dat_CPM[,i]/mapped_reads[i-81]
}
apply(dat_CPM[,82:92],MARGIN = 2,FUN=function(x){min(x[x>0])})
dat_CPM<-dat_CPM[,c(1:25,93,82:92)]
dat_CPM[dat_CPM==0]<-2^-10
KRAS_up<-c("ANGPTL4","ITGA2","SPRY2","HBEGF","RBP4","HSD11B1","ETV4","GLRX","DUSP6","SCG5",
           "ETV5","ITGB2","AKT2","PPBP","G0S2","GABRA3","IRF8","BIRC3","FGF9","DCBLD2","INHBA",
           "TFPI","TSPAN1","ADAM8","SLPI","PRKG2","MMP11","MMP10","TMEM158","TNFAIP3","PRDM1",
           "GALNT3","ETS1","MMP9","WNT7A","IGFBP3","SPP1","ETV1","CLEC4A","CCND2","TSPAN7",
           "ITGBL1","EMP1","CDADC1","KIF5C","TRIB2","SDCCAG8","PCP4","CFHR2","ALDH1A2","NR0B2",
           "ALDH1A3","AMMECR1","SATB1","GUCY1A1","CSF2","APOD","TOR1AIP2","CMKLR1","TMEM176B",
           "ADGRA2","LAPTM5","CD37","CAB39L","CIDEA","ZNF639","IL1B","GYPC","LY96","FLT4","SPON1",
           "BMP2","PLEK2","IGF2","NR1H4","SNAP25","ACE","PRRX1","C3AR1","TRAF1","TLR8","ID2",
           "TMEM100","PLAUR","GADD45G","CBX8","SCN1B","PTBP2","NAP1L2","AKAP12","PLAT","SCG3",
           "ANO1","IL1RL2","CXCL10","ATG10","YRDC","HDAC9","PEG3","SEMA3B","TNNT2","LIF","CFB",
           "BTC","PPP1R15A","PTPRR","CCL20","ARG1","RETN","KLF4","MMD","PDCD1LG2","H2BC3","HOXD11",
           "TRIB1","F2RL1","ANXA10","TSPAN13","MTMR10","CFH","LAT2","ERO1A","RELN","KCNN4",
           "TMEM176A","MAP4K1","PTGS2","IL33","MAFB","LCP1","NGF","CA2","SERPINA3","RGS16",
           "CTSS","USP12","CPE","SPARCL1","ABCB1","USH1C","CSF2RA","BTBD3","IL2RG","DNMBP",
           "IL10RA","EREG","PRELID3B","EPHB2","FBXO4","CROT","MPZL2","ANKH","CBR4","DOCK2",
           "GPRC5B","RABGAP1L","MALL","STRN","ST6GAL1","PIGR","VWA5A","PSMB8","F13A1","NRP1",
           "SOX9","JUP","ADGRL4","ZNF277","EPB41L3","PCSK1N","FUCA1","PLVAP","ADAM17","AVL9",
           "ADAMDEC1","HKDC1","MAP7","IL7R","RBM4","BPGM","ENG","GFPT2","PLAU","GNG11","PTCD2",
           "MAP3K1","CBL","CXCR4","NIN","IKZF1","WDR33","MYCN","FCER1G","PECAM1","CCSER2",
           "SNAP91","EVI5","TNFRSF1B","GPNMB","TPH1")
KRAS_dn<-c("CDH16","SPTBN2","FGFR3","NOS1","PDE6B","SIDT1","NTF3","SCN10A","TAS2R4","DTNB",
           "HTR1D","MAST3","SLC6A3","CLDN8","TGM1","PCDHB1","THRB","YPEL1","RYR2","GDNF",
           "CAMK1D","AMBN","LFNG","SERPINA10","ALOX12B","CALCB","FGGY","SPRR3","ATP6V1B1",
           "EDN1","UPK3B","GAMT","PRODH","RYR1","GPR19","SLC29A3","EDN2","GPRC5C","C5",
           "ZC2HC1C","EDAR","KCNMB1","TENT5C","STAG3","KRT4","SLC25A23","INSL5","GP1BA",
           "KMT2D","ABCG4","MYH7","BRDT","PTGFR","KCNN1","SELENOP","IFI44L","PROP1",
           "LGALS7","BMPR1B","PDK2","DCC","SNCB","COL2A1","UGT2B17","NPY4R","DLK2","WNT16",
           "SLC6A14","ITGB1BP2","SLC12A3","YBX2","CKM","CPB1","MAGIX","ARHGDIG","CALML5",
           "KRT1","MTHFR","ABCB11","MYOT","CELSR2","KLK7","THNSL2","CHRNG","NR4A2","CD40LG",
           "KRT13","GP2","CD207","CCR8","ZBTB16","PRKN","CPA2","MEFV","CCNA1","SLC38A3",
           "KLK8","CYP39A1","KCNQ2","SLC30A3","CHST2","ITIH3","EGF","MSH5","TEX15","CLPS",
           "ENTPD7","CLSTN3","CYP11B2","CLDN16","HSD11B2","COQ8A","HNF1A","FGF22","TFCP2L1",
           "OXT","KCND1","MACROH2A2","NRIP2","RGS11","CPEB3","TG","PTPRJ","NUDT11","FSHB",
           "IDUA","TNNI3","GTF3C5","TFF2","SSTR4","COPZ2","PDCD1","SLC16A7","KCNE2","MFSD6",
           "KLHDC8A","CNTFR","IL5","NPHS1","SCGB1A1","CCDC106","PAX3","VPREB1","TSHB",
           "CACNA1F","AKR1B10","P2RX6","KRT15","PAX4","BTG2","GPR3","GRID2","TCL1A",
           "PLAG1","PKP1","IRS4","IL12B","EPHA5","SOX10","SNN","CACNG1","SPHK2","IGFBP2",
           "CAPN9","TCF7L1","TGFB2","SMPX","LYPD3","PNMT","SYNPO","MX1","IFNG","NR6A1",
           "ACTC1","RSAD2","ADRA2C","BARD1","HTR1B","FGF16","TENM2","CDKAL1","SHOX2",
           "SGK1","RIBC2","SKIL","NGB","ASB7","MYO15A","SLC5A5","KRT5","ZNF112","TFAP2B",
           "VPS50","CD80","ATP4A","ARPP21","SERPINB2","TLX1","EFHD1","P2RY4")
KRASUP<-dat_CPM[dat_CPM$GName%in%KRAS_up,]

KRASDN<-dat_CPM[dat_CPM$GName%in%KRAS_dn,]

pdf("Figures/scatter_KRAS.pdf",height=12,width=6)
par(mfrow=c(2,1),pty="s")
p_col<-alpha("black",alpha = .50)
plot(log2(KRASUP$MDA)~log2(KRASUP$MCF7),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="MCF7 log2(RPM)",main="KRAS(UP, 10/13/37/200)",col=p_col,pch=16,bty="n")
abline(0,1,col="red")
plot(log2(KRASDN$MDA)~log2(KRASDN$MCF7),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="MCF7 log2(RPM)",main="KRAS(DOWN, 8/14/88/200)",col=p_col,pch=16,bty="n")
abline(0,1,col="red")
dev.off()

pdf("Figures/scatter_KRAS_vs_Healthy.pdf",height=12,width=6)
par(mfrow=c(2,1),pty="s")
p_col<-alpha("black",alpha = .50)
plot_set<-KRASUP
plot_set$Healthy<-(plot_set$BCH3+plot_set$BCH4)/2
plot_set<-plot_set[,c(7,38,31)]
plot_set<-plot_set[(plot_set$Healthy+plot_set$MDA)>2^-9,]
plot(log2(plot_set$MDA)~log2(plot_set$Healthy),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="Combined Healthy log2(RPM)",main="KRAS(UP, 10/12/37/200)",col=p_col,pch=16,bty="n")
abline(0,1,col="red")
plot_set<-KRASDN
plot_set$Healthy<-(plot_set$BCH3+plot_set$BCH4)/2
plot_set<-plot_set[,c(7,38,31)]
plot_set<-plot_set[(plot_set$Healthy+plot_set$MDA)>2^-9,]
plot(log2(plot_set$MDA)~log2(plot_set$Healthy),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="Combined Healthy log2(RPM)",main="KRAS(DOWN, 12/18/88/200)",col=p_col,pch=16,bty="n")
abline(0,1,col="red")
dev.off()

MMC<-dat_CPM[,c(7,32,33)]
MDA_vs_MCF_up<-c("RHGEF3","ATP7A","CCNA2","CD74","CDH11","F2R","RAC2","ARHGEF18","CUL5","LOX",
                 "MSN","RCN1","RUNX2","S100A8","TIMP2","VIM","BCL2L1","CAV2","CCNB1","F2RL3",
                 "AXL","MMP14","FLNB","CCNE1","COL6A1")
MDA_vs_MCF_dn<-c("CDH1","RARA","FHL1","AGR2","OGFOD3","KRT18","ESRP1","CDH3","TPI1","JUP",
                 "ANXA9","SPINT1","AURKB")

pdf("Figures/scatter_MDA_MCF.pdf",height=12,width=6)
par(mfrow=c(2,1),pty="s")
p_col<-alpha("black",alpha = .50)
plot_set<-dat_CPM[dat_CPM$GName%in%MDA_vs_MCF_up,c(7,32,31)]
plot(log2(plot_set$MDA)~log2(plot_set$MCF7),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="Combined Healthy log2(RPM)",main="Up in MDA-MB-231",col=p_col,pch=16,bty="n")
abline(0,1,col="red")

plot_set<-dat_CPM[dat_CPM$GName%in%MDA_vs_MCF_dn,c(7,32,31)]
plot_set<-KRASDN
plot_set$Healthy<-(plot_set$BCH3+plot_set$BCH4)/2
plot_set<-plot_set[,c(7,38,31)]
plot_set<-plot_set[(plot_set$Healthy+plot_set$MDA)>2^-9,]
plot(log2(plot_set$MDA)~log2(plot_set$Healthy),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 log2(RPM)",
     xlab="Combined Healthy log2(RPM)",main="KRAS(DOWN, 12/18/88/200)",col=p_col,pch=16,bty="n")
abline(0,1,col="red")
dev.off()

pdf("Figures/RBP_scatter.pdf",height=8,width=8)
par(mfrow=c(2,2),pty="s",pch=16)
RBP_fre<-RBP[,c(46,48,49)]
RBP_fre$All_Intron<-RBP$All_Intron-RBP$All_FLEXI
RBP_fre<-data.frame(prop.table(as.matrix(RBP_fre),margin = 2)*100)
RBP_fre$Name<-RBP$RBP.name
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-RBP_fre[,c(2,1,3,4,5)]
plot(RBP_fre[RBP_fre$col==0,1:2],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (%)",xlab="Long introns (%)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,1:2],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,1:2],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,1:2],col="skyblue",cex=1.5)
text(RBP_fre[RBP_fre$Cells>4,1:2],labels = RBP_fre$Name[RBP_fre$Cells>4])
abline(0,1,col="red")

plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="GRCh38 (%)",ylab="FLEXIs (%)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(RBP_fre[RBP_fre$Cells>4,c(3,2)],labels = RBP_fre$Name[RBP_fre$Cells>4])

RBP_fre<-RBP_fre[RBP_fre$Name%in% RBP_list,]
plot(RBP_fre[RBP_fre$col==0,1:2],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     ylab="FLEXIs (%)",xlab="Long introns (%)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,1:2],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,1:2],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,1:2],col="skyblue",cex=1.5)
#text(RBP_fre[RBP_fre$Cells>4,1:2],labels = RBP_fre$Name[RBP_fre$Cells>4])
abline(0,1,col="red")

plot(RBP_fre[RBP_fre$col==0,c(3,2)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="GRCh38 (%)",ylab="FLEXIs (%)")
points(RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,c(3,2)],col="red",cex=1.5)
points(RBP_fre[RBP_fre$col==1,c(3,2)],col="orange",cex=1.5)
points(RBP_fre[RBP_fre$col>1,c(3,2)],col="skyblue",cex=1.5)
abline(0,1,col="red")

dev.off()

Phast99<-read.table("44phastCon99.RBP.counts",col.names=c("ID","P99"))
Phast80<-read.table("120phastCon80.RBP.counts",col.names=c("ID","P80"))
Phast<-merge(Cell,Phast80,by=1,all=T)
Phast<-merge(Phast,Phast99,by=1,all=T)
Phast[is.na(Phast)]<-0
rownames(Phast)<-Phast$ID
Phast<-data.frame(prop.table(as.matrix(Phast[,2:4]),margin = 2)*100)
Phast$Name<-rownames(Phast)
Phast<-Phast[,c(4,1:3)]
Phast<-merge(Phast,RBP[,c(1,4,5,11)],by=1)
Phast$col<-(Phast$Splicing.regulation+Phast$Spliceosome)/3+Phast$microRNA.processing
Phast<-Phast[,c(1:4,8)]
pdf("Figures/RBP_byPhast_scatter.pdf",height=11,width=4)
par(mfrow=c(3,1),pty="s",pch=16)
plot(Phast[Phast$col==0,2:3],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="Cellular RNA (%)",ylab="phastCons >= 0.80 (%)")
points(Phast[Phast$col<1 & Phast$col>0,2:3],col="red",cex=1.5)
points(Phast[Phast$col==1,2:3],col="orange",cex=1.5)
points(Phast[Phast$col>1,2:3],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Phast[Phast$P80/Phast$Cells>2 & Phast$P80>4,c(2,3)],pos = 3,
     labels = Phast$Name[Phast$P80/Phast$Cells>2 & Phast$P80>4])

plot(Phast[Phast$col==0,c(2,4)],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="Cellular RNA (%)",ylab="phastCons >= 0.99 (%)")
points(Phast[Phast$col<1 & Phast$col>0,c(2,4)],col="red",cex=1.5)
points(Phast[Phast$col==1,c(2,4)],col="orange",cex=1.5)
points(Phast[Phast$col>1,c(2,4)],col="skyblue",cex=1.5)
abline(0,1,col="red")
text(Phast[Phast$P99/Phast$Cells>2 & Phast$P99>3,c(2,4)],pos = 3,
     labels = Phast$Name[Phast$P99/Phast$Cells>2 & Phast$P99>3])

plot(Phast[Phast$col==0,3:4],xlim=c(0,20),ylim=c(0,20),cex=1.5,bty="n",
     xlab="phastCons >= 0.80 (%)",ylab="phastCons >= 0.99 (%)")
points(Phast[Phast$col<1 & Phast$col>0,3:4],col="red",cex=1.5)
points(Phast[Phast$col==1,3:4],col="orange",cex=1.5)
points(Phast[Phast$col>1,3:4],col="skyblue",cex=1.5)
abline(0,1,col="red")
dev.off()



postscript("Figures/H_Vs_C_RPM005.eps",width = 15,height=10,horizontal = F)

par(pch=16,mfrow=c(2,3),pty="s")
A_unique<-unique(dat[dat$PatientA_H==2^-10 & dat$PatientA_C>=0.05 & dat$PatientA_C_repro,4])

plot(log2(CPM_FLEXI[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,2],CPM_FLEXI[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,2],CPM_FLEXI[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(CPM_FLEXI[CPM_FLEXI$ID%in%A_unique,c(2,4)]),col="red")

plot(log2(CPM_FLEXI[,c(8,10)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,8],CPM_FLEXI[,10],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,8],CPM_FLEXI[,10],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(CPM_FLEXI[CPM_FLEXI$ID%in%A_unique,c(8,10)]),col="red")

plot(log2(dat[,c(13,15)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(dat[,11],dat[,13],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,11],dat[,13],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$PatientA_H==2^-10 & dat$PatientA_C>=0.05 & dat$PatientA_C_repro,c(11,13)]),col="red")
legend(2,-7,col = c("red"),pch=16,legend = c("Cancer only"),bty="n")

B_unique<-unique(dat[dat$PatientB_H==2^-10 & dat$PatientB_C>=0.05 & dat$PatientB_C_repo,4])

plot(log2(CPM_FLEXI[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,3],CPM_FLEXI[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,3],CPM_FLEXI[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(CPM_FLEXI[CPM_FLEXI$ID%in%B_unique,c(3,5)]),col="red")

plot(log2(CPM_FLEXI[,c(9,11)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,9],CPM_FLEXI[,11],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,9],CPM_FLEXI[,11],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(CPM_FLEXI[CPM_FLEXI$ID%in%B_unique,c(9,11)]),col="red")

plot(log2(dat[,c(12,14)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(dat[,12],dat[,14],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,12],dat[,14],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$PatientB_H==2^-10 & dat$PatientB_C>=0.05 & dat$PatientB_C_repo,c(12,14)]),col="red")

dev.off()



postscript("Figures/BC_cellline_RPM005.eps",width = 10,height=15,horizontal = F)
BC_cell<-CPM_FLEXI[,c(1,6,7)]
BC_cell$Healthy<-rowMeans(CPM_FLEXI[,2:3])
BC_cell<-BC_cell[,c(1,4,2,3)]
par(pch=16,mfrow=c(3,2),pty="s")
A_unique<-unique(dat[dat$Combined_H==2^-10 & dat$MCF>=0.01 ,4])

plot(log2(BC_cell[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(BC_cell[,2],BC_cell[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(BC_cell[,2],BC_cell[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(BC_cell[BC_cell$ID%in%A_unique,c(2,4)]),col="red")

plot(log2(dat[,c(10,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(dat[,10],dat[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,10],dat[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$Combined_H==2^-10 & dat$MCF>=0.05 & dat$MCF_repro,c(10,16)]),col="red")
#legend(2,-7,col = c("red"),pch=16,legend = c("Cancer only"),bty="n")

B_unique<-unique(dat[dat$Combined_H==2^-10 & dat$MDA>=0.05 & dat$MDA_repo,4])

plot(log2(BC_cell[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(BC_cell[,2],BC_cell[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(BC_cell[,2],BC_cell[,3],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(BC_cell[BC_cell$ID%in%B_unique,c(2,3)]),col="red")

plot(log2(dat[,c(10,15)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(dat[,10],dat[,15],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,10],dat[,15],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$Combined_H==2^-10 & dat$MDA>=0.05 & dat$MDA_repo,c(10,15)]),col="red")

C_unique<-unique(dat[dat$MDA==2^-10,4])
D_unique<-unique(dat[dat$MCF==2^-10,4])

plot(log2(BC_cell[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(BC_cell[,3],BC_cell[,3],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(BC_cell[,3],BC_cell[,3],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))

plot(log2(dat[,c(15,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(dat[,15],dat[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,15],dat[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))

dev.off()

dat<-read.delim("all.FLEXI")
dat<-dat[,c(1:25,93,82:92)]
dat[,27:37]<-t(t(dat[,27:37])/mapped_reads)
dat$Healthy<-(dat$BCH3+dat$BCH4)/2
temp<-dat[,27:38]
temp[temp==0]<-2^-10
dat[,27:38]<-log2(temp)
rm(temp)

postscript("Figures/23RBP_scatter_BC.eps",width = 10,height=10,horizontal = F)
par(mfrow=c(2,2),pty="s")
RBP23<-c("AATF","BCLAF1","DDX24","DDX3X","DDX55","DKC1","FXR2","G3BP1","GRWD1","IGF2BP1",
         "LARP4","LSM11","METAP2","NOLC1","PABPC4","PABPN1","RPS3","SUB1","UCHL5","XRN2",
         "YBX3","ZNF622","ZNF800")
dat1<-dat[grep(paste(RBP23,collapse = "|"),dat$RBP),]
dat2<-dat[grep(paste(RBP23,collapse = "|"),dat$RBP,invert = T),]

plot(dat2[rowSums(dat2[,c(27,29)])!=-20,c(27,29)],pch=16,xlim=c(-10,5),ylim=c(-10,5),
     xlab="PatientA (Healthy)",ylab="PatientA (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(27,29)])!=-20,c(27,29)],col="red",pch=16)

plot(dat2[rowSums(dat2[,c(28,30)])!=-20,c(28,30)],pch=16,xlim=c(-10,5),ylim=c(-10,5),
     xlab="PatientB (Healthy)",ylab="PatientB (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(28,30)])!=-20,c(28,30)],col="red",pch=16)

plot(dat2[rowSums(dat2[,c(31,38)])!=-20,c(31,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),
     xlab="MDA-MB-231)",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(31,38)])!=-20,c(31,38)],col="red",pch=16)

plot(dat2[rowSums(dat2[,c(32,38)])!=-20,c(32,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),
     xlab="MCF7",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(32,38)])!=-20,c(32,38)],col="red",pch=16)
dev.off()

#postscript("Figures/MDA_MCF_scatter_BC.eps",width = 10,height=10,horizontal = F)
pdf("Figures/MDA_MCF_scatter_BC.pdf",width = 10,height=10)
par(mfrow=c(3,2),pty="s")
MDA_vs_MCF_up<-c("RHGEF3","ATP7A","CCNA2","CD74","CDH11","F2R","RAC2","ARHGEF18","CUL5","LOX",
                 "MSN","RCN1","RUNX2","S100A8","TIMP2","VIM","BCL2L1","CAV2","CCNB1","F2RL3",
                 "AXL","MMP14","FLNB","CCNE1","COL6A1")
MDA_vs_MCF_dn<-c("CDH1","RARA","FHL1","AGR2","OGFOD3","KRT18","ESRP1","CDH3","TPI1","JUP",
                 "ANXA9","SPINT1","AURKB")

dat1<-dat[grep(paste(MDA_vs_MCF_up,collapse = "|"),dat$GName),]
dat2<-dat[grep(paste(MDA_vs_MCF_up,collapse = "|"),dat$GName,invert = T),]
dat3<-dat2[grep(paste(MDA_vs_MCF_dn,collapse = "|"),dat2$GName),]
dat2<-dat2[grep(paste(MDA_vs_MCF_dn,collapse = "|"),dat2$GName,invert = T),]

plot(dat2[rowSums(dat2[,c(27,29)])!=-20,c(27,29)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientA (Healthy)",ylab="PatientA (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(27,29)])!=-20,c(27,29)],col="red",pch=16)
points(dat3[rowSums(dat3[,c(27,29)])!=-20,c(27,29)],col="blue",pch=16)
text(dat1[rowSums(dat1[,c(27,29)])!=-20,c(27,29)],labels =dat1$GName[rowSums(dat1[,c(27,29)])!=-20],
     pos = 3,cex=0.7,col="red")
text(dat3[rowSums(dat3[,c(27,29)])!=-20,c(27,29)],labels =dat3$GName[rowSums(dat3[,c(27,29)])!=-20],
     pos = 3,cex=0.7,col="blue")

plot(dat2[rowSums(dat2[,c(28,30)])!=-20,c(28,30)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientB (Healthy)",ylab="PatientB (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(28,30)])!=-20,c(28,30)],col="red",pch=16)
points(dat3[rowSums(dat3[,c(28,30)])!=-20,c(28,30)],col="blue",pch=16)
text(dat1[rowSums(dat1[,c(28,30)])!=-20,c(28,30)],labels =dat1$GName[rowSums(dat1[,c(28,30)])!=-20],
     pos = 3,cex=0.7,col="red")
text(dat3[rowSums(dat3[,c(28,30)])!=-20,c(28,30)],labels =dat3$GName[rowSums(dat3[,c(28,30)])!=-20],
     pos = 3,cex=0.7,col="blue")

plot(dat2[rowSums(dat2[,c(31,38)])!=-20,c(31,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MDA-MB-231)",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(31,38)])!=-20,c(31,38)],col="red",pch=16)
points(dat3[rowSums(dat3[,c(31,38)])!=-20,c(31,38)],col="blue",pch=16)
text(dat1[rowSums(dat1[,c(31,38)])!=-20,c(31,38)],labels =dat1$GName[rowSums(dat1[,c(31,38)])!=-20],
     pos = 3,cex=0.7,col="red")
text(dat3[rowSums(dat3[,c(31,38)])!=-20,c(31,38)],labels =dat3$GName[rowSums(dat3[,c(31,38)])!=-20],
     pos = 3,cex=0.7,col="blue")

plot(dat2[rowSums(dat2[,c(32,38)])!=-20,c(32,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MCF7",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(32,38)])!=-20,c(32,38)],col="red",pch=16)
points(dat3[rowSums(dat3[,c(32,38)])!=-20,c(32,38)],col="blue",pch=16)
text(dat1[rowSums(dat1[,c(32,38)])!=-20,c(32,38)],labels =dat1$GName[rowSums(dat1[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="red")
text(dat3[rowSums(dat3[,c(32,38)])!=-20,c(32,38)],labels =dat3$GName[rowSums(dat3[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="blue")

plot(dat2[rowSums(dat2[,c(32,33)])!=-20,c(32,33)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MCF7",ylab="MDA-MB-231")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(32,33)])!=-20,c(32,33)],col="red",pch=16)
points(dat3[rowSums(dat3[,c(32,33)])!=-20,c(32,33)],col="blue",pch=16)
text(dat1[rowSums(dat1[,c(32,33)])!=-20,c(32,33)],labels =dat1$GName[rowSums(dat1[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="red")
text(dat3[rowSums(dat3[,c(32,33)])!=-20,c(32,33)],labels =dat3$GName[rowSums(dat3[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="blue")

dev.off()

#BC marker
pdf("Figures/BCmarker_scatter_BC.pdf",width = 10,height=10)
par(mfrow=c(3,2),pty="s")
BC_marker<-c("ESR1", "PGR","BCL2","SCUBE2","ERBB2","GRB7","MKI67","AURKA","BIRC5",
             "CCNB1","MYBL2","MMP11","CTSV","GSTM1","CD68","BAG1","BRCA1","BRCA2")

dat1<-dat[grep(paste(BC_marker,collapse = "|"),dat$GName),]
dat2<-dat[grep(paste(BC_marker,collapse = "|"),dat$GName,invert = T),]

plot(dat2[rowSums(dat2[,c(27,29)])!=-20,c(27,29)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientA (Healthy)",ylab="PatientA (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(27,29)])!=-20,c(27,29)],col="red",pch=16)
text(dat1[rowSums(dat1[,c(27,29)])!=-20,c(27,29)],labels =dat1$GName[rowSums(dat1[,c(27,29)])!=-20],
     pos = 3,cex=0.7,col="red")

plot(dat2[rowSums(dat2[,c(28,30)])!=-20,c(28,30)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientB (Healthy)",ylab="PatientB (Cancer)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(28,30)])!=-20,c(28,30)],col="red",pch=16)
text(dat1[rowSums(dat1[,c(28,30)])!=-20,c(28,30)],labels =dat1$GName[rowSums(dat1[,c(28,30)])!=-20],
     pos = 3,cex=0.7,col="red")

plot(dat2[rowSums(dat2[,c(31,38)])!=-20,c(31,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MDA-MB-231)",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(31,38)])!=-20,c(31,38)],col="red",pch=16)
text(dat1[rowSums(dat1[,c(31,38)])!=-20,c(31,38)],labels =dat1$GName[rowSums(dat1[,c(31,38)])!=-20],
     pos = 3,cex=0.7,col="red")

plot(dat2[rowSums(dat2[,c(32,38)])!=-20,c(32,38)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MCF7",ylab="Healthy (A+B)")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(32,38)])!=-20,c(32,38)],col="red",pch=16)
text(dat1[rowSums(dat1[,c(32,38)])!=-20,c(32,38)],labels =dat1$GName[rowSums(dat1[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="red")

plot(dat2[rowSums(dat2[,c(32,33)])!=-20,c(32,33)],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="MCF7",ylab="MDA-MB-231")
abline(0,1,col="red")
points(dat1[rowSums(dat1[,c(32,33)])!=-20,c(32,33)],col="red",pch=16)
text(dat1[rowSums(dat1[,c(32,33)])!=-20,c(32,33)],labels =dat1$GName[rowSums(dat1[,c(32,38)])!=-20],
     pos = 3,cex=0.7,col="red")

dev.off()

#BC marker remaster
BC_FLEXI<-dat
BC_marker<-c("ESR1", "PGR","BCL2","SCUBE2","ERBB2","GRB7","MKI67","AURKA","BIRC5",
             "CCNB1","MYBL2","MMP11","CTSV","GSTM1","CD68","BAG1","BRCA1","BRCA2")
Oncotype<-c("BAG1","BCL2","BIRC5","CCNB1","CD68","CTSV","ESR1","GRB7",
            "GSTM1","HER2","MKI67","MMP11","MYBL2","PGR","SCUBE2","AURKA")
Onco<-read.delim("OncoGenes.table")
Onco$Onco<-1
Onco<-Onco[,c(2,8)]
colnames(Onco)[1]<-"GName"
TSG<-read.delim("TSGs.table")
TSG$TSG<-1
TSG<-TSG[,c(2,9)]
colnames(TSG)[1]<-"GName"
BC_FLEXI<-merge(BC_FLEXI,Onco,by="GName",all.x=T)
BC_FLEXI<-merge(BC_FLEXI,TSG,by="GName",all.x=T)
BC_FLEXI[is.na(BC_FLEXI)]<-0

pdf("Figures/BC_scatter_005RPM.pdf",width = 10,height=10)
par(mfrow=c(3,2),pty="s")
#patient A Cancer vs Healthy
pick1<-27
pick2<-29
pick=c(pick1,pick2)
cutOff<-log2(0.05)

#BC >= 0.05 RPM, n = 49, onco=1, TSG=2
BC_Sig<- BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
BC_Onco<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
BC_TSG<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(BC_Sig+BC_Onco+BC_TSG)
#BCH >= 0.05 RPM, n = 19, onco=5, TSG=3
H_Sig<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
H_Onco<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
H_TSG<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(H_Sig+H_Onco+H_TSG)
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)-BC_Sig-BC_Onco-BC_TSG-H_Sig-H_Onco-H_TSG

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientA (Healthy)",ylab="PatientA (Cancer)")
abline(0,1,col="red")
points(BC_FLEXI[BC_Sig,pick],col="red",pch=16)
points(BC_FLEXI[BC_Onco,pick],col="black",pch=10)
text(BC_FLEXI[BC_Onco,pick],labels =BC_FLEXI$GName[BC_Onco],
     pos = 4,cex=0.5,col="red")
points(BC_FLEXI[BC_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[BC_TSG,pick],labels =BC_FLEXI$GName[BC_TSG],
     pos = 4,cex=0.5,col="blue")

points(BC_FLEXI[H_Sig,pick],col="red",pch=16)
points(BC_FLEXI[H_Onco,pick],col="black",pch=10)
text(BC_FLEXI[H_Onco,pick],labels =BC_FLEXI$GName[H_Onco],srt=90,
     pos = 3,cex=0.5,col="red")
points(BC_FLEXI[H_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[H_TSG,pick],labels =BC_FLEXI$GName[H_TSG],srt=90,
     pos = 3,cex=0.5,col="blue")

#patient B Cancer vs Healthy
pick1<-28
pick2<-30
pick=c(pick1,pick2)

#BC >= 0.05 RPM, n = 19, onco=1, TSG=3
BC_Sig<- BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
BC_Onco<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
BC_TSG<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(BC_Sig+BC_Onco+BC_TSG)
#BCH >= 0.05 RPM, n = 7, onco=2, TSG=1
H_Sig<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
H_Onco<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
H_TSG<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(H_Sig+H_Onco+H_TSG)
#non-sig, n=2557
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)-BC_Sig-BC_Onco-BC_TSG-H_Sig-H_Onco-H_TSG

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientB (Healthy)",ylab="PatientB (Cancer)")
abline(0,1,col="red")
points(BC_FLEXI[BC_Sig,pick],col="red",pch=16)
points(BC_FLEXI[BC_Onco,pick],col="black",pch=10)
text(BC_FLEXI[BC_Onco,pick],labels =BC_FLEXI$GName[BC_Onco],
     pos = 4,cex=0.5,col="red")
points(BC_FLEXI[BC_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[BC_TSG,pick],labels =BC_FLEXI$GName[BC_TSG],
     pos = 4,cex=0.5,col="blue")

points(BC_FLEXI[H_Sig,pick],col="red",pch=16)
points(BC_FLEXI[H_Onco,pick],col="black",pch=10)
text(BC_FLEXI[H_Onco,pick],labels =BC_FLEXI$GName[H_Onco],srt=90,
     pos = 3,cex=0.5,col="red")
points(BC_FLEXI[H_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[H_TSG,pick],labels =BC_FLEXI$GName[H_TSG],srt=90,
     pos = 3,cex=0.5,col="blue")

#MDA Cancer vs Healthy (Combined)
pick1<-38
pick2<-31
pick=c(pick1,pick2)

#BC >= 0.05 RPM, n = 32
BC_Sig<- BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
BC_Onco<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
BC_TSG<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(BC_Sig+BC_Onco+BC_TSG)
#BCH >= 0.05 RPM, n = 6
H_Sig<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
H_Onco<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
H_TSG<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(H_Sig+H_Onco+H_TSG)
#non-sig, n=2557
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)-BC_Sig-BC_Onco-BC_TSG-H_Sig-H_Onco-H_TSG

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="Combined healthy",ylab="MDA-MB-231")
abline(0,1,col="red")
points(BC_FLEXI[BC_Sig,pick],col="red",pch=16)
points(BC_FLEXI[BC_Onco,pick],col="black",pch=10)
text(BC_FLEXI[BC_Onco,pick],labels =BC_FLEXI$GName[BC_Onco],
     pos = 4,cex=0.5,col="red")
points(BC_FLEXI[BC_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[BC_TSG,pick],labels =BC_FLEXI$GName[BC_TSG],
     pos = 4,cex=0.5,col="blue")

points(BC_FLEXI[H_Sig,pick],col="red",pch=16)
points(BC_FLEXI[H_Onco,pick],col="black",pch=10)
text(BC_FLEXI[H_Onco,pick],labels =BC_FLEXI$GName[H_Onco],srt=90,
     pos = 3,cex=0.5,col="red")
points(BC_FLEXI[H_TSG,pick],col="blue",pch=9)
text(BC_FLEXI[H_TSG,pick],labels =BC_FLEXI$GName[H_TSG],srt=90,
     pos = 3,cex=0.5,col="blue")

#MCF7 Cancer vs Healthy (Combined)
pick1<-38
pick2<-32
pick=c(pick1,pick2)

#BC >= 0.05 RPM, n = 2
BC_Sig<- BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
BC_Onco<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
BC_TSG<-BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(BC_Sig+BC_Onco+BC_TSG)
#BCH >= 0.05 RPM, n = 3
H_Sig<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==0
H_Onco<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==1 & BC_FLEXI$TSG==0
H_TSG<-BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]>=cutOff & BC_FLEXI$Onco==0 & BC_FLEXI$TSG==1
sum(H_Sig+H_Onco+H_TSG)
#non-sig, n=2557
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)-BC_Sig-BC_Onco-BC_TSG-H_Sig-H_Onco-H_TSG

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="Combined healthy",ylab="MDA-MB-231")
abline(0,1,col="red")
points(BC_FLEXI[BC_Sig,pick],col="red",pch=16)
points(BC_FLEXI[BC_Onco,pick],col="black",pch=10)
#text(BC_FLEXI[BC_Onco,pick],labels =BC_FLEXI$GName[BC_Onco],
#     pos = 4,cex=0.5,col="red")
points(BC_FLEXI[BC_TSG,pick],col="blue",pch=9)
#text(BC_FLEXI[BC_TSG,pick],labels =BC_FLEXI$GName[BC_TSG],
#     pos = 4,cex=0.5,col="blue")

points(BC_FLEXI[H_Sig,pick],col="red",pch=16)
points(BC_FLEXI[H_Onco,pick],col="black",pch=10)
text(BC_FLEXI[H_Onco,pick],labels =BC_FLEXI$GName[H_Onco],srt=90,
     pos = 3,cex=0.5,col="red")
points(BC_FLEXI[H_TSG,pick],col="blue",pch=9)
#text(BC_FLEXI[H_TSG,pick],labels =BC_FLEXI$GName[H_TSG],srt=90,
#     pos = 3,cex=0.5,col="blue")

dev.off()

#Top30 pick in BC
pdf("Figures/BC_scatter_top30.pdf",width=10,height=10)
BC_marker<-c("ESR1", "PGR","BCL2","SCUBE2","ERBB2","GRB7","MKI67","AURKA","BIRC5",
             "CCNB1","MYBL2","MMP11","CTSV","GSTM1","CD68","BAG1","BRCA1","BRCA2")
MDA_vs_MCF_up<-c("RHGEF3","ATP7A","CCNA2","CD74","CDH11","F2R","RAC2","ARHGEF18","CUL5","LOX",
                 "MSN","RCN1","RUNX2","S100A8","TIMP2","VIM","BCL2L1","CAV2","CCNB1","F2RL3",
                 "AXL","MMP14","FLNB","CCNE1","COL6A1")
MDA_vs_MCF_dn<-c("CDH1","RARA","FHL1","AGR2","OGFOD3","KRT18","ESRP1","CDH3","TPI1","JUP",
                 "ANXA9","SPINT1","AURKB")
par(mfrow=c(3,2),pty="s",mar=c(4.1,4.1,1.1,10.1))
#patient A Cancer vs Healthy
pick1<-27
pick2<-29
pick=c(pick1,pick2)
Top_pick<-30
P_col<-c("black","red","blue")
P_col[2:3]<-col2rgb(P_col[2:3])

Sig<- BC_FLEXI[BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]> -10,]
Sig<- Sig[order(Sig[,pick2],decreasing = T),]
Sig<-Sig[1:Top_pick,]
BC_Sig<- Sig$ID[Sig$Onco==0 & Sig$TSG==0]
BC_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
BC_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
BC_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
BC_name<-aggregate(Counts~GName+Color,data=Sig,sum)
BC_name$GName<-paste0(BC_name$GName," (",BC_name$Counts,")")
BC_name$Color<-factor(BC_name$Color)

Sig<- BC_FLEXI[BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]> -10,]
Sig<- Sig[order(Sig[,pick1],decreasing = T),]
Sig<-Sig[1:Top_pick,]
H_Sig<-Sig$ID[Sig$Onco==0 & Sig$TSG==0]
H_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
H_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
H_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
H_name<-aggregate(Counts~GName+Color,data=Sig,sum)
H_name$GName<-paste0(H_name$GName," (",H_name$Counts,")")
H_name$Color<-factor(H_name$Color)

BC_name<-BC_name[order(BC_name$GName),]
H_name<-H_name[order(H_name$GName),]
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientA (Healthy)",ylab="PatientA (Cancer)")
abline(0,1,col="red")
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_TSG,pick],col="blue",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_TSG,pick],col="blue",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$Onco==1,pick],col="red",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$TSG==1,pick],col="blue",pch=16)
legend(-10,4,legend = c("FLEXI in oncogene","FLEXI in TSG","Other"),
       text.col=c(P_col[2:3],"gray75"),col=c(P_col[2:3],"gray75"),pch=16,bty="n",cex=0.7)
par(xpd=T)
legend("topright",inset=c(-0.25,0),legend = BC_name$GName,text.col=P_col[BC_name$Color],bty="n",cex=0.7)
legend("topright",inset=c(-0.6,0),legend = H_name$GName,text.col=P_col[H_name$Color],bty="n",cex=0.7)
par(xpd=F)
#patient B Cancer vs Healthy
pick1<-28
pick2<-30
pick=c(pick1,pick2)
Sig<- BC_FLEXI[BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]> -10,]
Sig<- Sig[order(Sig[,pick2],decreasing = T),]
Sig<-Sig[1:Top_pick,]
BC_Sig<- Sig$ID[Sig$Onco==0 & Sig$TSG==0]
BC_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
BC_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
BC_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
BC_name<-aggregate(Counts~GName+Color,data=Sig,sum)
BC_name$GName<-paste0(BC_name$GName," (",BC_name$Counts,")")
BC_name$Color<-factor(BC_name$Color)

Sig<- BC_FLEXI[BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]> -10,]
Sig<- Sig[order(Sig[,pick1],decreasing = T),]
Sig<-Sig[1:Top_pick,]
H_Sig<-Sig$ID[Sig$Onco==0 & Sig$TSG==0]
H_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
H_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
H_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
H_name<-aggregate(Counts~GName+Color,data=Sig,sum)
H_name$GName<-paste0(H_name$GName," (",H_name$Counts,")")
H_name$Color<-factor(H_name$Color)

BC_name<-BC_name[order(BC_name$GName),]
H_name<-H_name[order(H_name$GName),]
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="PatientB (Healthy)",ylab="PatientB (Cancer)")
abline(0,1,col="red")
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_TSG,pick],col="blue",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$Onco==1,pick],col="red",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$TSG==1,pick],col="blue",pch=16)

#points(BC_FLEXI[BC_FLEXI$ID%in%H_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_TSG,pick],col="blue",pch=16)
par(xpd=T)
legend("topright",inset=c(-0.25,0),legend = BC_name$GName,text.col=P_col[BC_name$Color],bty="n",cex=0.7)
legend("topright",inset=c(-0.6,0),legend = H_name$GName,text.col=P_col[H_name$Color],bty="n",cex=0.7)
par(xpd=F)

#MDA Cancer vs Healthy (Combined)
pick1<-38 
pick2<-31
pick=c(pick1,pick2)
Sig<- BC_FLEXI[BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]> -10,]
Sig<- Sig[order(Sig[,pick2],decreasing = T),]
Sig<-Sig[1:Top_pick,]
BC_Sig<- Sig$ID[Sig$Onco==0 & Sig$TSG==0]
BC_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
BC_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
BC_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
BC_name<-aggregate(Counts~GName+Color,data=Sig,sum)
BC_name$GName<-paste0(BC_name$GName," (",BC_name$Counts,")")
BC_name$Color<-factor(BC_name$Color)

Sig<- BC_FLEXI[BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]> -10,]
Sig<- Sig[order(Sig[,pick1],decreasing = T),]
Sig<-Sig[1:Top_pick,]
H_Sig<-Sig$ID[Sig$Onco==0 & Sig$TSG==0]
H_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
H_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
H_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
H_name<-aggregate(Counts~GName+Color,data=Sig,sum)
H_name$GName<-paste0(H_name$GName," (",H_name$Counts,")")
H_name$Color<-factor(H_name$Color)

BC_name<-BC_name[order(BC_name$GName),]
H_name<-H_name[order(H_name$GName),]
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="Combined healthy",ylab="MDA-MB-231")
abline(0,1,col="red")
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_TSG,pick],col="blue",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$Onco==1,pick],col="red",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$TSG==1,pick],col="blue",pch=16)

#points(BC_FLEXI[BC_FLEXI$ID%in%H_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_TSG,pick],col="blue",pch=16)
par(xpd=T)
legend("topright",inset=c(-0.25,0),legend = BC_name$GName,text.col=P_col[BC_name$Color],bty="n",cex=0.7)
legend("topright",inset=c(-0.6,0),legend = H_name$GName,text.col=P_col[H_name$Color],bty="n",cex=0.7)
par(xpd=F)
#MCF7 Cancer vs Healthy (Combined)
pick1<-38
pick2<-32
pick=c(pick1,pick2)
Sig<- BC_FLEXI[BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]> -10,]
Sig<- Sig[order(Sig[,pick2],decreasing = T),]
Sig<-Sig[1:Top_pick,]
BC_Sig<- Sig$ID[Sig$Onco==0 & Sig$TSG==0]
BC_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
BC_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
BC_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
BC_name<-aggregate(Counts~GName+Color,data=Sig,sum)
BC_name$GName<-paste0(BC_name$GName," (",BC_name$Counts,")")
BC_name$Color<-factor(BC_name$Color)

Sig<- BC_FLEXI[BC_FLEXI[,pick2]==-10 & BC_FLEXI[,pick1]> -10,]
Sig<- Sig[order(Sig[,pick1],decreasing = T),]
Sig<-Sig[1:Top_pick,]
H_Sig<-Sig$ID[Sig$Onco==0 & Sig$TSG==0]
H_Onco<-Sig$ID[Sig$Onco==1 & Sig$TSG==0]
H_TSG<-Sig$ID[Sig$Onco==0 & Sig$TSG==1]
H_M<-Sig$ID[Sig$ID%in%BC_marker]
Sig$Color<-Sig$Onco+2*Sig$TSG
Sig$Counts<-1
H_name<-aggregate(Counts~GName+Color,data=Sig,sum)
H_name$GName<-paste0(H_name$GName," (",H_name$Counts,")")
H_name$Color<-factor(H_name$Color)

BC_name<-BC_name[order(BC_name$GName),]
H_name<-H_name[order(H_name$GName),]
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     xlab="Combined healthy",ylab="MDA-MB-231")
abline(0,1,col="red")
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%BC_TSG,pick],col="blue",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$Onco==1,pick],col="red",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$TSG==1,pick],col="blue",pch=16)

#points(BC_FLEXI[BC_FLEXI$ID%in%H_Sig,pick],col="black",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_Onco,pick],col="red",pch=16)
#points(BC_FLEXI[BC_FLEXI$ID%in%H_TSG,pick],col="blue",pch=16)
par(xpd=T)
legend("topright",inset=c(-0.25,0),legend = BC_name$GName,text.col=P_col[BC_name$Color],bty="n",cex=0.7)
legend("topright",inset=c(-0.6,0),legend = H_name$GName,text.col=P_col[H_name$Color],bty="n",cex=0.7)
par(xpd=F)

#MDA vs MCF
pick1<-31
pick2<-32
pick=c(pick1,pick2)
#non-sig
non_Sig<-!(BC_FLEXI[,pick1]==-10 & BC_FLEXI[,pick2]==-10)

plot(BC_FLEXI[non_Sig,pick],pch=16,xlim=c(-10,5),ylim=c(-10,5),col="gray75",
     ylab="MCF7",xlab="MDA-MB-231")
abline(0,1,col="red")
points(BC_FLEXI[non_Sig & BC_FLEXI$GName%in% MDA_vs_MCF_up,pick],col="red",pch=16)
points(BC_FLEXI[non_Sig & BC_FLEXI$GName%in% MDA_vs_MCF_dn,pick],col="blue",pch=16)
legend(-10,4,legend = c("Up in MDA-MB-231","Up in MCF7","Other"),
       text.col=c(P_col[2:3],"gray75"),col=c(P_col[2:3],"gray75"),pch=16,bty="n",cex=0.7)
dev.off()

FourCellPlasma_snoRNA<-FourCellPlasma[FourCellPlasma$Has_snoRNA!=".",]
sno=c()
for (i in 1:53){
  if (length(grep(RBP_list[i],FourCellPlasma_snoRNA$RBP))) {
    sno=c(sno,RBP_list[i])} 
}
FourCellPlasma_nonsnoRNA<-FourCellPlasma[FourCellPlasma$Has_snoRNA==".",]
nonSno=c()
for (i in 1:53){
  if (length(grep(RBP_list[i],FourCellPlasma_nonsnoRNA$RBP))) {
    nonSno=c(nonSno,RBP_list[i])} 
}

Vcol=c("tomato","royalblue1","greenyellow","goldenrod","orchid","black")

pdf("Figures/Venn_up_BC.pdf",height=8, width=8)
set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set5<-BC_FLEXI[(BC_FLEXI$BCH3-BC_FLEXI$BC3>=1) & BC_FLEXI$Onco==1,2]
set6<-BC_FLEXI[(BC_FLEXI$BCH4-BC_FLEXI$BC4>=1) & BC_FLEXI$Onco==1,2]
set7<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MDA>=1) & BC_FLEXI$Onco==1,2]
set8<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MCF7>=1) & BC_FLEXI$Onco==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
vennplot1<-venn.diagram (set, category.names=names(set),
                         cat.col = Vcol[1:4], main = "Oncogene (Up)",
                         fill = Vcol[1:4],
                         height = 300, width = 300, units = "px",
                         cex = 1,filename=NULL,
                         main.cex=1, cat.cex = 1) 
set <- list ("Patient A"=set5,"Patient B"=set6,
             "MDA-MB-231"=set7,"MCF7"=set8)
vennplot2<-venn.diagram (set, category.names=names(set),
                         cat.col = Vcol[1:4], main = "Oncogene (Down)",
                         fill = Vcol[1:4],
                         height = 300, width = 300, units = "px",
                         cex = 1,filename=NULL,
                         main.cex=1, cat.cex = 1) 

set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$TSG==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$TSG==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$TSG==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$TSG==1,2]
set5<-BC_FLEXI[(BC_FLEXI$BCH3-BC_FLEXI$BC3>=1) & BC_FLEXI$TSG==1,2]
set6<-BC_FLEXI[(BC_FLEXI$BCH4-BC_FLEXI$BC4>=1) & BC_FLEXI$TSG==1,2]
set7<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MDA>=1) & BC_FLEXI$TSG==1,2]
set8<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MCF7>=1) & BC_FLEXI$TSG==1,2]
set <- list ("Patient A"=set1,"Patient B"=set2,
             "MDA-MB-231"=set3,"MCF7"=set4)
vennplot3<-venn.diagram (set, category.names=names(set),
                         cat.col = Vcol[1:4], main = "TSG (Up)",
                         fill = Vcol[1:4],
                         height = 300, width = 300, units = "px",
                         cex = 1,filename=NULL,
                         main.cex=1, cat.cex = 1) 
set <- list ("Patient A"=set5,"Patient B"=set6,
             "MDA-MB-231"=set7,"MCF7"=set8)
vennplot4<-venn.diagram (set, category.names=names(set),
                         cat.col = Vcol[1:4], main = "TSG (Down)",
                         fill = Vcol[1:4],
                         height = 300, width = 300, units = "px",
                         cex = 1,filename=NULL,
                         main.cex=1, cat.cex = 1) 
ggarrange(vennplot1,vennplot2,vennplot3,vennplot4)
unlink("*.log")
dev.off()

#upset of Onco in BC
library(UpSetR)

pdf("Figures/BC_onco_upset.pdf",height=8,width=11)
set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]

A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))
Any_cancer<-union(set1,union(set2,union(set3,set4)))

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

A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))

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
A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))

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

A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))

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



snoRNA_FLEXI_RB<-read.table(pipe("pbpaste"),col.names="RBP")

#Top20 RBPs by occurrence in FLEXIs
pdf("Figures/RBP_4cells_log-all.pdf",height=12,width=12)
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
par(xpd=T)
legend("topright",legend = Fun$Name,text.col=B_col[Fun$Color],cex=0.3)
par(xpd=F)
dev.off()

#RBPs in Oncogene
RBP_info<-read.delim(pipe("cut -f 1,6 all_FLEXI_RBP_intersect.info"),col.names=c("ID","RBP"))
RBP_info<-RBP_info[RBP_info$RBP!=".",]
RBP_info<-separate(RBP_info,col="ID",into=c("IID","GID","TID","GType","TType"),
                   sep="___",remove=F)
RBP_info<-separate(RBP_info,col="IID",into=c("IID","GName"),sep="_")
#all_onco_RBP<-data.frame(table(RBP_info$RBP[RBP_info$GName%in%Onco$GName]))

set1<-BC_FLEXI[BC_FLEXI$MDA>BC_FLEXI$Healthy & BC_FLEXI$Onco==1 & BC_FLEXI$MDA>-10,2]
set2<-BC_FLEXI[BC_FLEXI$MCF>BC_FLEXI$Healthy & BC_FLEXI$Onco==1 & BC_FLEXI$MCF>-10,2]
set3<-BC_FLEXI[BC_FLEXI$BC3>BC_FLEXI$BCH3 & BC_FLEXI$Onco==1 & BC_FLEXI$BC3>-10,2]
set4<-BC_FLEXI[BC_FLEXI$BC4>BC_FLEXI$BCH4 & BC_FLEXI$Onco==1 & BC_FLEXI$BC4>-10,2]
set<-intersect(set1,set2)
set<-intersect(set,set3)
set<-intersect(set,set4)

Four_cancer_RBP<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%set]))
Four_cancer_RBP<-merge(Four_cancer_RBP,Cell,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,FLEXI,by=1,all=T)
Four_cancer_RBP[is.na(Four_cancer_RBP)]<-0
Four_cancer_RBP$Var1<-as.character(Four_cancer_RBP$Var1)
Four_cancer_RBP[Four_cancer_RBP$Var1=="AGO",][2]<-1
Four_cancer_RBP[Four_cancer_RBP$Var1=="DICER",][2]<-4

#77, 19008, 63568
for (i in 1:146){
  if (rowSums(Four_cancer_RBP[i,c(2,4)]>0)){
    set_num<-matrix(unlist(c(Four_cancer_RBP[i,c(2,4)],77,63568)),2,2)
    if (fisher.test(set_num)$p.value<0.05) {
      print (c(Four_cancer_RBP[i,1],fisher.test(set_num)$p.value))
    }
  }
}

Four_cancer_RBP[,2:4]<-prop.table(as.matrix(Four_cancer_RBP[,2:4]),margin = 2)

pdf("Figures/RBP_by_cancer.pdf",height=8,width=4)
par(mfrow=c(2,1),pty="s",pch=20)
plot(100*Four_cancer_RBP[,c(3,2)],xlim=c(0,25),ylim=c(0,25),cex=1.5,bty="n",
     xlab="FELXIs (% RBP sites)",ylab="Oncogene FLEXIs (% RBP sites)")
abline(0,1,col="red")
sig<-c("BUD13","FTO","GPKOW","NCBP2","RBM15","RBM22","XPO5","DICER")
points(100*Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,c(3,2)],col="red")
text(100*Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,c(3,2)],
     labels = Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,1],pos = 3,cex=0.5)


plot(100*Four_cancer_RBP[,c(4,2)],xlim=c(0,25),ylim=c(0,25),cex=1.5,bty="n",
     xlab="All short introns (% RBP sites)",ylab="Oncogene FLEXIs (% RBP sites)")
abline(0,1,col="red")
sig<-c("BUD13","GPKOW","NCBP2","RBM15","RBM22","XPO5","DICER")
points(100*Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,c(4,2)],col="red")
text(100*Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,c(4,2)],cex=0.5,
     labels = Four_cancer_RBP[Four_cancer_RBP$Var1%in%sig,1],pos = 3)

dev.off()

#RBP of BC(onco)
set1<-BC_FLEXI[(BC_FLEXI$BC3-BC_FLEXI$BCH3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BC4-BC_FLEXI$BCH4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$MDA-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$MCF7-BC_FLEXI$Healthy>=1) & BC_FLEXI$Onco==1,2]

A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
Any_cancer<-union(set1,union(set2,union(set3,set4)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))
AGO_sites<-c(dim(dat2[dat2$ID%in%A_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%B_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MDA_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MCF_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%Any_cancer & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%All_cancer & dat2$AGO_CCR!=".",])[1])
DICER_sites<-c(dim(dat2[dat2$ID%in%A_unique & dat2$DICER_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%B_unique & dat2$DICER_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MDA_unique & dat2$DICER_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MCF_unique & dat2$DICER_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%Any_cancer & dat2$DICER_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%All_cancer & dat2$DICER_CCR!=".",])[1])

A_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%A_unique]))
colnames(A_unique)<-c("RBP","A_unique")
B_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%B_unique]))
colnames(B_unique)<-c("RBP","B_unique")
MDA_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MDA_unique]))
colnames(MDA_unique)<-c("RBP","MDA_unique")
MCF_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MCF_unique]))
colnames(MCF_unique)<-c("RBP","MCF_unique")
Any_cancer<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%Any_cancer]))
colnames(Any_cancer)<-c("RBP","Any_cancer")
All_cancer<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%All_cancer]))
colnames(All_cancer)<-c("RBP","All_cancer")

Four_cancer_RBP<-merge(A_unique,B_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,MDA_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,MCF_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Any_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,All_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Cell,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,FLEXI,by=1,all=T)

Four_cancer_RBP[is.na(Four_cancer_RBP)]<-0
Four_cancer_RBP$RBP<-as.character(Four_cancer_RBP$RBP)

Four_cancer_RBP[Four_cancer_RBP$RBP=="AGO",][2:7]<-AGO_sites
Four_cancer_RBP[Four_cancer_RBP$RBP=="DICER",][2:7]<-DICER_sites

pdf("Figures/RBP_by_cancer.pdf",height=6,width=8)
par(mfrow=c(3,4),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
for (i in 2:7){
  Name=colnames(Four_cancer_RBP)[i]
  dat_set<-Four_cancer_RBP[,c(1,i,8,9)]
  RBP_sum<-colSums(dat_set[,2:4])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
    FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)]),2,2))$p.value})
  dat_set$F_sig<-apply(dat_set[,c(2,4)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,3)]),2,2))$p.value})
  dat_set[,2:4]<-100*prop.table(as.matrix(dat_set[,2:4]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
  text(dat_set[dat_set$C_sig<=0.05,c(3,2)],
       labels = dat_set[dat_set$C_sig<=0.05,1],pos = 3,cex=0.7)
  if (i<6) {
    axis(1,at=seq(0,25,5),labels = NA)
  } else {
    axis(1,at=seq(0,25,5),labels = seq(0,25,5))
    title(xlab="FLEXIs (% RBP sites)")
  }
  if (i%%2==0){
    axis(2,at=seq(0,25,5),labels = seq(0,25,5))
    title(ylab="Oncogene FLEXIs (% RBP sites)")
  } else {
    axis(2,at=seq(0,25,5),labels = NA)
  }
  plot(dat_set[,c(4,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$F_sig<=0.05,c(4,2)],col="red")
  text(dat_set[dat_set$F_sig<=0.05,c(4,2)],
       labels = dat_set[dat_set$F_sig<=0.05,1],pos = 3,cex=0.7)
  if (i<6) {
    axis(1,at=seq(0,25,5),labels = NA)
  } else {
    axis(1,at=seq(0,25,5),labels = seq(0,25,5))
    title(xlab="All short introns (% RBP sites)")
  }
  axis(2,at=seq(0,25,5),labels = NA)
}
dev.off()

pdf("Figures/RBP_by_cancer_0_6.pdf",height=7.5,width=10)
par(mfrow=c(3,4),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
for (i in 2:7){
  Name=colnames(Four_cancer_RBP)[i]
  dat_set<-Four_cancer_RBP[,c(1,i,8,9)]
  RBP_sum<-colSums(dat_set[,2:4])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)]),2,2))$p.value})
  dat_set$F_sig<-apply(dat_set[,c(2,4)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,3)]),2,2))$p.value})
  dat_set[,2:4]<-100*prop.table(as.matrix(dat_set[,2:4]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,6),ylim=c(0,6),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
  text(dat_set[dat_set$C_sig<=0.05,c(3,2)],
       labels = dat_set[dat_set$C_sig<=0.05,1],pos = 3,cex=0.7)
  if (i<6) {
    axis(1,at=seq(0,6,2),labels = NA)
  } else {
    axis(1,at=seq(0,6,2),labels = seq(0,6,2))
    title(xlab="FLEXIs (% RBP sites)")
  }
  if (i%%2==0){
    axis(2,at=seq(0,6,2),labels = seq(0,6,2))
    title(ylab="Oncogene FLEXIs (% RBP sites)")
  } else {
    axis(2,at=seq(0,6,2),labels = NA)
  }
  plot(dat_set[,c(4,2)],xlim=c(0,6),ylim=c(0,6),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  points(dat_set[dat_set$F_sig<=0.05,c(4,2)],col="red")
  text(dat_set[dat_set$F_sig<=0.05,c(4,2)],
       labels = dat_set[dat_set$F_sig<=0.05,1],pos = 3,cex=0.7)
  if (i<6) {
    axis(1,at=seq(0,6,2),labels = NA)
  } else {
    axis(1,at=seq(0,6,2),labels = seq(0,6,2))
    title(xlab="All short introns (% RBP sites)")
  }
  axis(2,at=seq(0,6,2),labels = NA)
}
dev.off()

set1<-BC_FLEXI[(BC_FLEXI$BCH3-BC_FLEXI$BC3>=1) & BC_FLEXI$Onco==1,2]
set2<-BC_FLEXI[(BC_FLEXI$BCH4-BC_FLEXI$BC4>=1) & BC_FLEXI$Onco==1,2]
set3<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MDA>=1) & BC_FLEXI$Onco==1,2]
set4<-BC_FLEXI[(BC_FLEXI$Healthy-BC_FLEXI$MCF7>=1) & BC_FLEXI$Onco==1,2]

A_unique<-setdiff(set1,union(set2,union(set3,set4)))
B_unique<-setdiff(set2,union(set1,union(set3,set4)))
MDA_unique<-setdiff(set3,union(set1,union(set2,set4)))
MCF_unique<-setdiff(set4,union(set1,union(set2,set3)))
Any_cancer<-union(set1,union(set2,union(set3,set4)))
All_cancer<-intersect(set1,intersect(set2,intersect(set3,set4)))
AGO_sites<-c(dim(dat2[dat2$ID%in%A_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%B_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MDA_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%MCF_unique & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%Any_cancer & dat2$AGO_CCR!=".",])[1],
             dim(dat2[dat2$ID%in%All_cancer & dat2$AGO_CCR!=".",])[1])
DICER_sites<-c(dim(dat2[dat2$ID%in%A_unique & dat2$DICER_CCR!=".",])[1],
               dim(dat2[dat2$ID%in%B_unique & dat2$DICER_CCR!=".",])[1],
               dim(dat2[dat2$ID%in%MDA_unique & dat2$DICER_CCR!=".",])[1],
               dim(dat2[dat2$ID%in%MCF_unique & dat2$DICER_CCR!=".",])[1],
               dim(dat2[dat2$ID%in%Any_cancer & dat2$DICER_CCR!=".",])[1],
               dim(dat2[dat2$ID%in%All_cancer & dat2$DICER_CCR!=".",])[1])

A_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%A_unique]))
colnames(A_unique)<-c("RBP","A_unique")
B_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%B_unique]))
colnames(B_unique)<-c("RBP","B_unique")
MDA_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MDA_unique]))
colnames(MDA_unique)<-c("RBP","MDA_unique")
MCF_unique<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%MCF_unique]))
colnames(MCF_unique)<-c("RBP","MCF_unique")
Any_cancer<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%Any_cancer]))
colnames(Any_cancer)<-c("RBP","Any_cancer")
All_cancer<-data.frame(table(RBP_info$RBP[RBP_info$ID%in%All_cancer]))
colnames(All_cancer)<-c("RBP","All_cancer")

Four_cancer_RBP<-merge(A_unique,B_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,MDA_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,MCF_unique,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Any_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,All_cancer,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,Cell,by=1,all=T)
Four_cancer_RBP<-merge(Four_cancer_RBP,FLEXI,by=1,all=T)

Four_cancer_RBP[is.na(Four_cancer_RBP)]<-0
Four_cancer_RBP$RBP<-as.character(Four_cancer_RBP$RBP)

Four_cancer_RBP[Four_cancer_RBP$RBP=="AGO",][2:7]<-AGO_sites
Four_cancer_RBP[Four_cancer_RBP$RBP=="DICER",][2:7]<-DICER_sites

pdf("Figures/RBP_by_cancer_down.pdf",height=7.5,width=10)
par(mfrow=c(3,4),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
for (i in 2:7){
  Name=colnames(Four_cancer_RBP)[i]
  dat_set<-Four_cancer_RBP[,c(1,i,8,9)]
  RBP_sum<-colSums(dat_set[,2:4])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)]),2,2))$p.value})
  dat_set$F_sig<-apply(dat_set[,c(2,4)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,3)]),2,2))$p.value})
  dat_set[,2:4]<-100*prop.table(as.matrix(dat_set[,2:4]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  if (sum(dat_set$C_sig<=0.05)>0){
  points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
  text(dat_set[dat_set$C_sig<=0.05,c(3,2)],
         labels = dat_set[dat_set$C_sig<=0.05,1],pos = 3,cex=0.7)
  }
  if (i<6) {
    axis(1,at=seq(0,25,5),labels = NA)
  } else {
    axis(1,at=seq(0,25,5),labels = seq(0,25,5))
    title(xlab="FLEXIs (% RBP sites)")
  }
  if (i%%2==0){
    axis(2,at=seq(0,25,5),labels = seq(0,25,5))
    title(ylab="Oncogene FLEXIs (% RBP sites)")
  } else {
    axis(2,at=seq(0,25,5),labels = NA)
  }
  plot(dat_set[,c(4,2)],xlim=c(0,25),ylim=c(0,25),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  if (sum(dat_set$F_sig<=0.05)>0){
    points(dat_set[dat_set$F_sig<=0.05,c(4,2)],col="red")
    text(dat_set[dat_set$F_sig<=0.05,c(4,2)],
         labels = dat_set[dat_set$F_sig<=0.05,1],pos = 3,cex=0.7)
  }
  if (i<6) {
    axis(1,at=seq(0,25,5),labels = NA)
  } else {
    axis(1,at=seq(0,25,5),labels = seq(0,25,5))
    title(xlab="All short introns (% RBP sites)")
  }
  axis(2,at=seq(0,25,5),labels = NA)
}
dev.off()

pdf("Figures/RBP_by_cancer_down_0_6.pdf",height=7.5,width=10)
par(mfrow=c(3,4),pty="s",pch=16,cex=0.7,mai=c(0.3,0.3,0.3,0.3))
for (i in 2:7){
  Name=colnames(Four_cancer_RBP)[i]
  dat_set<-Four_cancer_RBP[,c(1,i,8,9)]
  RBP_sum<-colSums(dat_set[,2:4])
  dat_set$C_sig<-apply(dat_set[,c(2,3)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,2)]),2,2))$p.value})
  dat_set$F_sig<-apply(dat_set[,c(2,4)],MARGIN = 1,
                       FUN=function(x){fisher.test(as.matrix(rbind(x,RBP_sum[c(1,3)]),2,2))$p.value})
  dat_set[,2:4]<-100*prop.table(as.matrix(dat_set[,2:4]),margin = 2)
  plot(dat_set[,c(3,2)],xlim=c(0,6),ylim=c(0,6),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  if (sum(dat_set$C_sig<=0.05)>0){
    points(dat_set[dat_set$C_sig<=0.05,c(3,2)],col="red")
    text(dat_set[dat_set$C_sig<=0.05,c(3,2)],
         labels = dat_set[dat_set$C_sig<=0.05,1],pos = 3,cex=0.7)
  }
  if (i<6) {
    axis(1,at=seq(0,6,2),labels = NA)
  } else {
    axis(1,at=seq(0,6,2),labels = seq(0,6,2))
    title(xlab="FLEXIs (% RBP sites)")
  }
  if (i%%2==0){
    axis(2,at=seq(0,6,2),labels = seq(0,6,2))
    title(ylab="Oncogene FLEXIs (% RBP sites)")
  } else {
    axis(2,at=seq(0,6,2),labels = NA)
  }
  plot(dat_set[,c(4,2)],xlim=c(0,6),ylim=c(0,6),bty="n",main=Name,axes=F,frame.plot=TRUE,
       xlab="",ylab="")
  abline(0,1,col="red")
  if (sum(dat_set$F_sig<=0.05)>0){
    points(dat_set[dat_set$F_sig<=0.05,c(4,2)],col="red")
    text(dat_set[dat_set$F_sig<=0.05,c(4,2)],
         labels = dat_set[dat_set$F_sig<=0.05,1],pos = 3,cex=0.7)
  }
  if (i<6) {
    axis(1,at=seq(0,6,2),labels = NA)
  } else {
    axis(1,at=seq(0,6,2),labels = seq(0,6,2))
    title(xlab="All short introns (% RBP sites)")
  }
  axis(2,at=seq(0,6,2),labels = NA)
}
dev.off()


#Fig 6B gene list
library(matrixStats)
dat<-read.delim("all.FLEXI")
dat<-dat[,c(1,32:34,29:31,26:28,35:37,38:47)]
dat<-dat[rowSums(dat[,2:23])>0,]
dat$BCH3_repro<-rowMedians(as.matrix(dat[,2:4]))>0
dat$BCH4_repro<-rowMedians(as.matrix(dat[,5:7]))>0
dat$BC3_repro<-rowMedians(as.matrix(dat[,8:10]))>0
dat$BC4_repro<-rowMedians(as.matrix(dat[,11:13]))>0
dat$MDA_repro<-rowMins(as.matrix(dat[,14:15]))>0
dat$MCF_repro<-apply(dat[,16:23],MARGIN = 1,FUN=function(x){sort(x)[5]>0})

dat$PatientA_H<-rowSums(dat[,c(2:4)])
dat$PatientB_H<-rowSums(dat[,c(5:7)])
dat$PatientA_C<-rowSums(dat[,c(8:10)])
dat$PatientB_C<-rowSums(dat[,c(11:13)])
dat$MDA<-rowSums(dat[,c(14:15)])
dat$MCF<-rowSums(dat[,c(16:23)])
Unfrag_total<-c(305.069837,251.558067,268.210336,477.543790,207.491024,692.091831)
dat[,30:35]<-t(t(dat[,30:35])/Unfrag_total)
dat1<-dat[,30:35]
dat1[dat1==0]<-2^-10
dat[,30:35]<-dat1
dat$Combined_H<-rowMeans(dat[,30:31])
dat$Combined_H_repo<-dat$BCH3_repro | dat$BCH4_repro
dat<-dat[,c(1,37,26:29,36,30:35)]
dat$ID<-as.character(dat$ID)
dat<-separate(dat,ID,into=c("IID","GID"),sep = "___",remove = F,extra="drop")
dat<-separate(dat,IID,into=c("Intron","GName"),sep = "_",remove = T,extra="drop")
colnames(dat)[6:7]<-c("PatientA_C_repro","PatientB_C_repo")

postscript("Figures/H_Vs_C_RPM005.eps",width = 15,height=10,horizontal = F)
par(pch=16,mfrow=c(2,3),pty="s")
plot(log2(CPM_FLEXI[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,2],CPM_FLEXI[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,2],CPM_FLEXI[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))

plot(log2(CPM_FLEXI[,c(6,8)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,6],CPM_FLEXI[,8],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,6],CPM_FLEXI[,8],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
text(log2(CPM_FLEXI[1942,c(6,8)])+1,"CYP2B7P")

plot(log2(dat[,c(11,13)]),xlim=c(-10,5),ylim=c(-10,5),col="gray50")
abline(0,1,col="red")
cor_s=formatC(cor(dat[,11],dat[,13],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,11],dat[,13],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$PatientA_H==2^-10 & dat$PatientA_C>=0.05 & dat$PatientA_C_repro,c(11,13)]),col="red")
legend(-2,-7,col = c("red"),pch=16,legend = c("Cancer only"),bty="n")

#points(log2(dat[dat$GName%in%Oncotype,c(11,13)]),col="red")


plot(log2(CPM_FLEXI[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,3],CPM_FLEXI[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,3],CPM_FLEXI[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))

plot(log2(CPM_FLEXI[,c(7,9)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(CPM_FLEXI[,7],CPM_FLEXI[,9],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(CPM_FLEXI[,7],CPM_FLEXI[,9],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))

plot(log2(dat[,c(12,14)]),xlim=c(-10,5),ylim=c(-10,5),col="gray50")
abline(0,1,col="red")
cor_s=formatC(cor(dat[,12],dat[,14],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(dat[,12],dat[,14],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(dat[dat$PatientB_H==2^-10 & dat$PatientB_C>=0.05 & dat$PatientB_C_repo,c(12,14)]),col="red")
dev.off()

RBP_list<-unique(RBP_info_4cell$RBP)
RBP_4cell_FLEXI_matrix<-data.frame(table(FourCellPlasma$ID[FourCellPlasma$DICER_CCR!="."]))
colnames(RBP_4cell_FLEXI_matrix)<-c("ID","DICER")
temp<-data.frame(table(FourCellPlasma$ID[FourCellPlasma$AGO_CCR!="."]))
colnames(temp)<-c("ID","AGO1-4")
RBP_4cell_FLEXI_matrix<-merge(RBP_4cell_FLEXI_matrix,temp,by=1,all=T)
for (i in 1:124){
  name_temp<-unique(RBP_info_4cell$ID[RBP_info_4cell$RBP==RBP_list[i]])
  temp<-data.frame(ID=name_temp,Counts=rep(1,length(name_temp)))
  colnames(temp)[2]<-RBP_list[i]
  RBP_4cell_FLEXI_matrix<-merge(RBP_4cell_FLEXI_matrix,temp,by=1,all=T)
}
RBP_4cell_FLEXI_matrix[is.na(RBP_4cell_FLEXI_matrix)]<-0
rownames(RBP_4cell_FLEXI_matrix)<-RBP_4cell_FLEXI_matrix$ID
RBP_4cell_FLEXI_matrix<-RBP_4cell_FLEXI_matrix[,2:127]
RBP_4cell_FLEXI_matrix<-t(RBP_4cell_FLEXI_matrix)
pca_res <- prcomp(RBP_4cell_FLEXI_matrix)
autoplot(pca_res)

#reads cov density in intron
AGO<-read.table(pipe("cut -f 4 AGO_FLEXI.bed"),col.names="ID")
DICER<-read.table(pipe("cut -f 4 DICER_FLEXI.bed"),col.names="ID")
SPLICE<-read.table(pipe("cut -f 4 Splice5core_FLEXI.bed"),col.names="ID")
RBP_other<-read.table(pipe("cut -f 4 Other_RBPs_FLEXI.bed"),col.names="ID")
nonRBP<-read.table(pipe("cut -f 4 non_RBP_snoRNA.bed"),col.names="ID")

name<-c("UHRR","K562","HEK293T","Hela")
file_name<-paste0(name,".filtered.wao.gz")
out_name<-paste0(name,".percentile.info")
for (i in 1:4){
  con<-gzfile(file_name[i],open="r")
  write(c("ID","Per","AGO","DICER","Core5","OtherRBP","nonRBP"),out_name[i],sep="\t",ncolumns=7)
  while (TRUE){
    line<-readLines(con,n=1)
    if (length(line) == 0) {
      break
    }
    line<-strsplit(line,split="\t")
    L1<-line[[1]][10]
    L2<-100*as.numeric(line[[1]][13])/(as.numeric(line[[1]][9])-as.numeric(line[[1]][8]))
    L3<-sum(line[[1]][10]%in%AGO$ID)
    L4<-sum(line[[1]][10]%in%DICER$ID)
    L5<-sum(line[[1]][10]%in%SPLICE$ID)
    L6<-sum(line[[1]][10]%in%RBP_other$ID)
    L7<-sum(line[[1]][10]%in%nonRBP$ID)
    write(c(L1,L2,L3,L4,L5,L6,L7),out_name[i],append=T,sep="\t",ncolumns=7)
  }
  close(con)
}
