library(matrixStats)
library(tidyr)
library(ggpubr)
library("Rtsne")
library(zinbwave)
args<-commandArgs(trailingOnly=TRUE)
cutoff<-args[1]
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
mapped_reads<-c(sum(sub_mapped_reads[19:20]),sum(sub_mapped_reads[21:28]),
                sum(sub_mapped_reads[29:36]),sum(sub_mapped_reads[37:44]),
                sum(sub_mapped_reads[45:49]),sum(sub_mapped_reads[50:59]),
                sum(sub_mapped_reads[60:66]),sum(sub_mapped_reads[67:68]))
names(mapped_reads)<-c("MDA1","UHRR1","K562-1","HEK1","MDA2","UHRR2","K562-2","HEK2")

dat1<-dat
dat1$MDA1<-rowSums(dat1[,20:21])
dat1$UHRR1<-rowSums(dat1[,22:29])
dat1$K562_1<-rowSums(dat1[,30:37])
dat1$HEK1<-rowSums(dat1[,38:45])
dat1$MDA2<-rowSums(dat1[,46:50])
dat1$UHRR2<-rowSums(dat1[,51:60])
dat1$K562_2<-rowSums(dat1[,61:67])
dat1$HEK2<-rowSums(dat1[,68:69])
dat1$MDA1_repo<-rowSums(dat1[,20:21]>0)
dat1$UHRR1_repo<-rowSums(dat1[,22:29]>0)
dat1$K5621_repo<-rowSums(dat1[,30:37]>0)
dat1$HEK1_repo<-rowSums(dat1[,38:45]>0)
dat1<-dat1[,c(1,70:81)]
dat1[,2:9]<-t(t(dat1[,2:9])/mapped_reads)

temp<-dat1[,2:9]
temp[temp==0]<-2^-10
dat1[,2:9]<-temp
dat1[,2:9]<-log2(dat1[,2:9])

dat_clus<-dat
# and save the object
dat<-dat_clus
obj_name<-paste0("44_zinbwave",cutoff,"RPM.pbj")
##apply rowmax,rowmin give back very low number of flexi
dat<-dat[rowMaxs(t(t(as.matrix(dat[,2:45]))/sub_mapped_reads[1:44]))>=cutoff,2:45]
#apply rowMeans
#obj_name<-paste0("44_zinbwave_rowMean",cutoff,"RPM.pbj")
#dat<-dat[rowMeans(t(t(as.matrix(dat[,2:45]))/sub_mapped_reads[1:44]))>=cutoff,2:45]
dat<-as.matrix(dat[rowSums(dat)>0,])
Fc<-SummarizedExperiment(assays = list(counts = dat),
                         colData = data.frame(label=c(rep(1,8),rep(2,10),rep(3,2),rep(4:6,8))))
zinb<-zinbwave(Fc,normalizedValues = TRUE)
saveRDS(zinb,obj_name)

# biological replicate ZINB batch correction
dat<-dat_clus
dat1<-data.frame(t(dat[,c(20:69)]))
dat1$label<-c(rep("MDA",2),rep("UHRR",8),rep("K562",8),rep("HEK",8),
              rep("MDA",5),rep("UHRR",10),rep("K562",7),rep("HEK",2))
dat1$batch<-c(rep("1",26),rep("2",24))
dat1$label<-factor(dat1$label)
dat1$batch<-factor(dat1$batch)

dat<-dat[rowMaxs(t(t(as.matrix(dat[,20:69]))/sub_mapped_reads[19:68]))>=cutoff,20:69]
#dat<-dat[rowMeans(t(t(as.matrix(dat[,30:54]))/sub_mapped_reads[29:53]))>=cutoff,30:54]
dat<-as.matrix(dat[rowSums(dat)>0,])
Fc<-SummarizedExperiment(assays = list(counts = dat),
                         colData = data.frame(label=dat1$label,batch=dat1$batch))
obj_name<-paste0("Bio_rep_zinbwave",cutoff,"RPM.pbj")
#obj_name<-paste0("Bio_rep_zinbwave_rowmean_",cutoff,"RPM.pbj")
zinb<-zinbwave(Fc,normalizedValues = TRUE)
saveRDS(zinb,obj_name)
obj_name<-paste0("Bio_rep_zinbwave_Batch_cor_",cutoff,"RPM.pbj")
#obj_name<-paste0("Bio_rep_zinbwave_Batch_cor_rowmean_",cutoff,"RPM.pbj")
zinb<-zinbwave(Fc,X="~batch",normalizedValues = TRUE)
saveRDS(zinb,obj_name)