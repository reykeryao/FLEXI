#Fig3A using stringtie read counts
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

Cell_counts<-read.delim("combined_gene_bystringtie.counts")
rownames(Cell_counts)<-Cell_counts$ID
Cell_counts$ID<-as.character(Cell_counts$ID)
#Subset only FLEXI host genes
GID_list<-unique(FLEXI_CPM$GID)
Cell_counts<-Cell_counts[Cell_counts$ID%in%GID_list,c(1,9,7,8)]
Cell_counts[,2:4]<-t(t(Cell_counts[,2:4])/Unfrag_total)
#Assign log2CPM of -10 for those with 0 count
Cell_counts[Cell_counts==0]<-2^-10

pdf("Figures/Fig3A_stringtie.pdf",width = 9,height=6)
par(pch=16,mfcol=c(2,3),pty="s")
#FLEXIs scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(2,4)]),col="#FF000050")
points(log2(sig.x[,c(2,4)]),col="#0000FF50")
#Unfrag scatter
#FELXIs between K562 and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,4)]),col="#0000FF80")


#FLEXIs scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,3]==2^-10 & FLEXI_CPM[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$Hela_repo>=10,]
points(log2(sig.y[,c(3,4)]),col="#FF000080")
points(log2(sig.x[,c(3,4)]),col="#0000FF80")
#Unfrag scatter
#FELXIs between HEK and HeLa
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,3]==2^-10 & Cell_counts[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,4)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(3,4)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(3,4)]),col="#0000FF80")

#FELXIs between K562 and HEK
FLEXI_by_GID<-FLEXI_CPM[!(FLEXI_CPM[,2]==2^-10 & FLEXI_CPM[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,5),ylim=c(-10,5),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(1,-7,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-FLEXI_by_GID[FLEXI_by_GID$K562_repo>=8,]
sig.y<-FLEXI_by_GID[FLEXI_by_GID$HEK_repo>=8,]
points(log2(sig.y[,c(2,3)]),col="#FF000080")
points(log2(sig.x[,c(2,3)]),col="#0000FF80")

#Unfrag scatter
#FELXIs between K562 and HEK
FLEXI_by_GID<-Cell_counts[!(Cell_counts[,2]==2^-10 & Cell_counts[,3]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,3)]),xlim=c(-10,15),ylim=c(-10,15),col="gray")
abline(0,1,col="red")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,3],method = "pearson"),digits=2, format="f")
text(10,-5,bquote(atop(italic(r)== .(cor_p)~phantom())))
sig.x<-unique(sig.x$GID)
sig.y<-unique(sig.y$GID)
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.y,c(2,3)]),col="#FF000080")
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%sig.x,c(2,3)]),col="#0000FF80")
dev.off()



#Fig 8B scatter plot
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

#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("IBC_frag_gene_bystringtie.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_total<-c(50.48402,33.084889,56.983289,54.093371)
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,]
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")


pdf("Figures/Fig8B_frag_stringtie.pdf",width = 10,height=10)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
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
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>2^-10  & 
                           FLEXI_by_GID$PatientB_C_repro>=3,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & 
                                    FLEXI_by_GID$PatientB_C>=0.01  &
                                    FLEXI_by_GID$PatientB_C_repro>=3])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()

#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("combined_run.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_by_FLEXI$BCH3<-rowSums(Frag_by_FLEXI[,10:12])
Frag_by_FLEXI$BCH4<-rowSums(Frag_by_FLEXI[,4:6])
Frag_by_FLEXI$BC3<-rowSums(Frag_by_FLEXI[,7:9])
Frag_by_FLEXI$BC4<-rowSums(Frag_by_FLEXI[,13:15])
Frag_by_FLEXI<-Frag_by_FLEXI[,c(1,60:63)]
Frag_total<-mapped_reads[1:4]
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,]
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")


pdf("Figures/Fig8B_unfrag.pdf",width = 10,height=10)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
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
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>2^-10  & 
                           FLEXI_by_GID$PatientB_C_repro>=3,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & 
                                    FLEXI_by_GID$PatientB_C>=0.01  &
                                    FLEXI_by_GID$PatientB_C_repro>=3])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()


#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("combined_gene_bystringtie.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_by_FLEXI<-Frag_by_FLEXI[,c(1,5,6,3,4)]
Frag_total<-mapped_reads[1:4]
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,]
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")


pdf("Figures/Fig8B_unfrag_stringtie.pdf",width = 10,height=10)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
FLEXI_by_GID<-Repo[!(Repo[,13]==2^-10 & Repo[,15]==2^-10),]
plot(log2(FLEXI_by_GID[,c(13,15)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,13],FLEXI_by_GID[,15],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.01 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientA_H==2^-10 & FLEXI_by_GID$PatientA_C>2^-10  & 
                           FLEXI_by_GID$PatientA_C_repro>=3,c(13,15)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientA_H==2^-10 & 
                                    FLEXI_by_GID$PatientA_C>=0.01 & 
                                    FLEXI_by_GID$PatientA_C_repro>=3])
print(unique(FLEXI_by_GID$ID[FLEXI_by_GID$PatientA_H==2^-10 & 
                               FLEXI_by_GID$PatientA_C>2^-10  & 
                               FLEXI_by_GID$PatientA_C_repro>=3]))
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>2^-10  & 
                           FLEXI_by_GID$PatientB_C_repro>=3,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & 
                                    FLEXI_by_GID$PatientB_C>=0.01  &
                                    FLEXI_by_GID$PatientB_C_repro>=3])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,15),ylim=c(-10,15))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-5,10,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()

#check correlation coefficiency diff
library(cocor)
cocor_dat<-Repo[,c(4,13:16)]
cocor_dat<-merge(cocor_dat,Frag_by_FLEXI,by=1,all=T)
cocor_dat[is.na(cocor_dat)]<-0
cocor_dat<-cocor_dat[rowMaxs(as.matrix(cocor_dat[,2:9]))>2^-10,] 
#test patientA
A<-cocor_dat[,c(2,4,6,8)]
A<-A[rowMaxs(as.matrix(A))>2^-10,]
cor(A)
cocor(~PatientA_H + PatientA_C | PatientA_Healthy + PatientA_Cancer, A)
B<-cocor_dat[,c(3,5,7,9)]
B<-B[rowMaxs(as.matrix(B))>2^-10,]
cor(B)
cocor(~PatientB_H + PatientB_C | PatientB_Healthy + PatientB_Cancer, B)


#downsampl
#FigS17A downsampled scatter as Fig8B
#new Fig8B alt
#Fig 8B scatter plot
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
Repo_total<-colSums(Repo[,30:33])
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
'''
#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("combined_run.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_by_FLEXI$BCH3<-rowSums(Frag_by_FLEXI[,10:12])
Frag_by_FLEXI$BCH4<-rowSums(Frag_by_FLEXI[,4:6])
Frag_by_FLEXI$BC3<-rowSums(Frag_by_FLEXI[,7:9])
Frag_by_FLEXI$BC4<-rowSums(Frag_by_FLEXI[,13:15])
Frag_by_FLEXI<-Frag_by_FLEXI[,c(1,60:63)]
Frag_total<-mapped_reads[1:4]
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
FLEXI_host_total<-colSums(Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(2:5)])
Frag_ratio<-colSums(Frag_by_FLEXI[,2:5])/FLEXI_host_total
DownSampling<-ceiling(Frag_ratio*Repo_total)

#down sampling
for (i in 1:4){
  tmp<-rep(Frag_by_FLEXI$ID,Frag_by_FLEXI[,1+i])
  tmp<-sample(tmp)
  tmp<-sample(tmp,DownSampling[i],replace = T)
  tmp<-data.frame(table(tmp))
  Frag_by_FLEXI<-merge(Frag_by_FLEXI,tmp,by=1,all=T)
}
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(1,6:9)]
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")
write.table(Frag_by_FLEXI,"Fig8B_unfrag_downsampled.counts",quote=F,sep="\t",row.names=F)
'''
Frag_by_FLEXI<-read.delim("Fig8B_unfrag_downsampled.counts")
pdf("Figures/FigS17A_unfrag.pdf",width = 10,height=10)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
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
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>2^-10  & 
                           FLEXI_by_GID$PatientB_C_repro>=3,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & 
                                    FLEXI_by_GID$PatientB_C>=0.01  &
                                    FLEXI_by_GID$PatientB_C_repro>=3])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()
'''
#downsample mRNA in cancer
#Calculate CPM of fragmented patientA/B cellular RNA
Frag_by_FLEXI<-read.delim("IBC_frag_transcriptome_mapping.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_total<-c(50.48402,33.084889,56.983289,54.093371)

#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
FLEXI_host_total<-colSums(Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(2:5)])
Frag_ratio<-Frag_total*1e6/FLEXI_host_total
DownSampling<-ceiling(Frag_ratio*Repo_total)

#down sampling
for (i in 1:4){
  tmp<-rep(Frag_by_FLEXI$ID,Frag_by_FLEXI[,1+i])
  tmp<-sample(tmp)
  tmp<-sample(tmp,DownSampling[i],replace = T)
  tmp<-data.frame(table(tmp))
  Frag_by_FLEXI<-merge(Frag_by_FLEXI,tmp,by=1,all=T)
}
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(1,6:9)]
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")
write.table(Frag_by_FLEXI,"Fig8B_frag_transcriptome_downsampled.counts",quote=F,sep="\t",row.names=F)
'''
Frag_by_FLEXI<-read.delim("Fig8B_frag_transcriptome_downsampled.counts")
pdf("Figures/FigS17A_frag_transcriptome.pdf",width = 10,height=10)
par(pch=16,mfrow=c(2,2),pty="s")
#patientA FLEXI scatter
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
#patientA fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,2]==2^-10 & Frag_by_FLEXI[,4]==2^-10),]
plot(log2(FLEXI_by_GID[,c(2,4)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,2],FLEXI_by_GID[,4],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(2,4)]),col="red")
#patientB FELXI scatter
FLEXI_by_GID<-Repo[!(Repo[,14]==2^-10 & Repo[,16]==2^-10),]
plot(log2(FLEXI_by_GID[,c(14,16)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,14],FLEXI_by_GID[,16],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
#cancer ≥ 0.05 RPM and reproducible FLEXIs (detected in ≥ 50% dataset)
points(log2(FLEXI_by_GID[FLEXI_by_GID$PatientB_H==2^-10 & FLEXI_by_GID$PatientB_C>2^-10  & 
                           FLEXI_by_GID$PatientB_C_repro>=3,c(14,16)]),col="red")
GID_list<-unique(FLEXI_by_GID$GID[FLEXI_by_GID$PatientB_H==2^-10 & 
                                    FLEXI_by_GID$PatientB_C>=0.01  &
                                    FLEXI_by_GID$PatientB_C_repro>=3])
#patientB fragmented cellular RNA scatter plot
FLEXI_by_GID<-Frag_by_FLEXI[!(Frag_by_FLEXI[,3]==2^-10 & Frag_by_FLEXI[,5]==2^-10),]
plot(log2(FLEXI_by_GID[,c(3,5)]),xlim=c(-10,5),ylim=c(-10,5))
abline(0,1,col="red")
cor_s=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "spearman"),digits=2, format="f")
cor_p=formatC(cor(FLEXI_by_GID[,3],FLEXI_by_GID[,5],method = "pearson"),digits=2, format="f")
text(-8,4,bquote(atop(italic(r)== .(cor_p)~phantom(),italic(r[s])== .(cor_s))))
points(log2(FLEXI_by_GID[FLEXI_by_GID$ID%in%GID_list,c(3,5)]),col="red")
dev.off()



Frag_by_FLEXI<-read.delim("IBC_frag_stringtie.counts")
Frag_by_FLEXI$ID<-as.character(Frag_by_FLEXI$ID)
Frag_total<-c(50.48402,33.084889,56.983289,54.093371)
Transcript_total<-colSums(Frag_by_FLEXI[,3:6])
#Subset only FLEXI host genes
GID_list<-unique(Repo$GID)
FLEXI_host_total<-colSums(Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(3:6)])
Frag_ratio<-Transcript_total/FLEXI_host_total
DownSampling<-ceiling(Frag_ratio*Repo_total)

#down sampling
for (i in 1:4){
  tmp<-rep(Frag_by_FLEXI$ID,Frag_by_FLEXI[,2+i])
  tmp<-sample(tmp)
  tmp<-sample(tmp,DownSampling[i],replace = T)
  tmp<-data.frame(table(tmp))
  Frag_by_FLEXI<-merge(Frag_by_FLEXI,tmp,by=1,all=T)
}
Frag_by_FLEXI[is.na(Frag_by_FLEXI)]<-0
Frag_by_FLEXI<-Frag_by_FLEXI[Frag_by_FLEXI$ID%in%GID_list,c(1,7:10)]
Frag_by_FLEXI[,2:5]<-t(t(Frag_by_FLEXI[,2:5])/Frag_total)
#Assign log2CPM of -10 for those with 0 count
Frag_by_FLEXI[Frag_by_FLEXI==0]<-2^-10
colnames(Frag_by_FLEXI)<-c("ID","PatientA_Healthy","PatientB_Healthy",
                           "PatientA_Cancer","PatientB_Cancer")


tmp<-dat
tmp$MDA<-tmp$MDA/mapped_reads[5]
tmp$MCF7<-tmp$MCF7/mapped_reads[6]
ks.test(tmp$MDA[tmp$MDA>0],tmp$MCF7[tmp$MCF7>0],alternative = "less")
wilcox.test(tmp$MDA[tmp$MDA>0],tmp$MCF7[tmp$MCF7>0],alternative = "greater")
tmp1<-tmp[tmp$Is_agotron!=".",]
ks.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "less")
wilcox.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "greater")
tmp1<-tmp[tmp$Is_mirtron!=".",]
ks.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "less")
wilcox.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "greater")
tmp1<-tmp[tmp$Has_snoRNA!=".",]
ks.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "less")
wilcox.test(tmp1$MDA[tmp1$MDA>0],tmp1$MCF7[tmp1$MCF7>0],alternative = "greater")

pdf("~/Desktop/scatter.pdf",height=4,width=8)
par(mfrow=c(1,2),pty="s",pch=16)
tmp1<-tmp[tmp$Is_agotron!=".",]
tmp1<-tmp1[tmp1$MCF7+tmp1$MDA>0,]
tmp1[tmp1==0]<-2^-10
plot(log2(tmp1$MCF7),log2(tmp1$MDA),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 (log2RPM)",xlab="MCF7 (log2RPM)",main="Agotron FLEXI")
abline(0,1,col="red")
points(log2(tmp1[tmp1$MDA>=0.01& tmp1$MCF7==2^-10,c(87,86)]),col="red")
points(log2(tmp1[tmp1$MCF7>=0.01& tmp1$MDA==2^-10,c(87,86)]),col="red")
sapply(strsplit(tmp1$ID[tmp1$MCF7>=0.01& tmp1$MDA==2^-10],split ="___"),"[[",1)
tmp1<-tmp[tmp$Is_mirtron!=".",]
tmp1<-tmp1[tmp1$MCF7+tmp1$MDA>0,]
tmp1[tmp1==0]<-2^-10
plot(log2(tmp1$MCF7),log2(tmp1$MDA),xlim=c(-10,0),ylim=c(-10,0),ylab="MDA-MB-231 (log2RPM)",xlab="MCF7 (log2RPM)",
     main="Mirtron FLEXI")
abline(0,1,col="red")
points(log2(tmp1[tmp1$MDA>=0.01& tmp1$MCF7==2^-10,c(87,86)]),col="red")
sapply(strsplit(tmp1$ID[tmp1$MDA>=0.01& tmp1$MCF7==2^-10],split ="___"),"[[",1)
points(log2(tmp1[tmp1$MCF7>=0.01& tmp1$MDA==2^-10,c(87,86)]),col="red")
sapply(strsplit(tmp1$ID[tmp1$MCF7>=0.01& tmp1$MDA==2^-10],split ="___"),"[[",1)
dev.off()