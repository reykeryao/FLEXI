RBP_fre<-RBP[,c(1,46)]
RBP_fre$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-RBP_fre[RBP_fre$Cells>0,]
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_list<-sort(unique(RBP_4cell_plasma$RBP))
for (i in 1:126) {
  RBP_name<-RBP_list[i]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name])
  FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Padj_withSign<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:126){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=3,7]<-"D"
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=3,7]<-"U"
  RBP_plot[RBP_plot$Cells<=RBP_plot$Freq,6]<- -log10(RBP_plot[RBP_plot$Cells<=RBP_plot$Freq,5])
  RBP_plot[RBP_plot$Cells>RBP_plot$Freq,6]<- log10(RBP_plot[RBP_plot$Cells>RBP_plot$Freq,5])
  if (i==1){
    RBP_clus<-RBP_plot[,c(1,2,7)]
    colnames(RBP_clus)[3]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-RBP_plot[,c(1,2,6)]
    colnames(RBP_clus_log10)[3]<-paste0(RBP_name,"BF")
    
  } else {
    RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,2,7)],by=1:2)
    colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,2,6)],by=1:2)
    colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
  }
}
rownames(RBP_clus)<-RBP_clus$RBP.name
RBP_clus<-RBP_clus[,c(3:128)]
write.table(RBP_clus,"RBP_clus_froGower_cutoff3.txt",quote=F,sep="\t")
rownames(RBP_clus_log10)<-RBP_clus_log10$RBP.name
RBP_clus_log10<-RBP_clus_log10[,3:128]
write.table(RBP_clus_log10,"RBP_clus_froHeatmap_cutoff3.txt",quote=F,sep="\t")



RBP_clus<-read.delim("RBP_clus_froGower_cutoff3.txt")
Sig_label<-RBP_clus
RBP_clus<-data.frame(t(RBP_clus))
for (i in 1:126){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}

RBP_clus_log10<-read.delim("RBP_clus_froHeatmap_cutoff3.txt")

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}
RBP_col<-RBP[,c(1,46)]
RBP_col$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_col[RBP_col$col>1,3]<-4
RBP_col[RBP_col$col==1,3]<-3
RBP_col[RBP_col$col<1 & RBP_col$col>0,3]<-2
RBP_col[RBP_col$col==0,3]<-1
RBP_col<-RBP_col[,c(1,3)]
rownames(RBP_col)<-RBP_col$RBP.name

col<-c("black","red","orange","skyblue")
pdf("temp_fig/126RBP_heatmap_cluster3cutoff.pdf",height=8,width=8)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 4,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        show_column_names = F,
        cluster_rows = F, cluster_columns =aggl.clust.r)
dev.off()



RBP_47<-paste0(Fun$Name[Fun$RBP_by_FLEXI>=30][7:53],"BF")
RBP_47_name<-Fun$Name[Fun$RBP_by_FLEXI>=30][7:53]

RBP_clus_log10.47<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
RBP_clus_log10.47<-RBP_clus_log10.47[rownames(RBP_clus_log10.47)%in%RBP_47_name,]

RBP_clus.47<-RBP_clus[rownames(RBP_clus)%in%RBP_47,]
RBP_clus.47<-RBP_clus.47[,colnames(RBP_clus.47)%in%RBP_47_name]

gower.dist <- daisy(RBP_clus.47, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10.47
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue>=5]<-5
RBP_clus_log10_convertedvalue[RBP_clus_log10_convertedvalue<=-5]<- -5
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

pdf("temp_fig/47RBP_heatmap_cluster3cutoff2.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
grid_text<-Sig_label[,colnames(Sig_label)%in%RBP_47]
grid_text<-grid_text[labels.r,]
grid_text[grid_text=="N"]<-""
grid_text[grid_text=="U"]<-"X"
grid_text[grid_text=="D"]<-"X"
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        show_column_names = F,
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        #cell_fun = function(j, i, x, y, w, h, col) {
        #  grid.text(grid_text[i, j], x, y,gp=gpar(fontsize=7))},
        cluster_rows = F, cluster_columns =aggl.clust.r)
dev.off()

RBP_47_info<-RBP[match(labels.r,RBP$RBP.name),1:43]

RBP_47_info$StressGranule<-0
SG<-read.table("SG_MSGP_RGD.list",col.names="ID")
RBP_47_info[RBP_47_info$RBP.name%in%SG$ID,44]<-1
RBP_47_info<-RBP_47_info[,c(1:31,44,32:43)]
RBP_47_info[is.na(RBP_47_info)]<-0
RBP_47_info<-RBP_47_info[,c(TRUE,TRUE,colSums(RBP_47_info[,3:44])>0)]
#AGO is in P body, stress granule, nuclei and cytoplasm
RBP_47_info[RBP_47_info$RBP.name=="AGO",c(19,23,26,27)]<-1
#DICER is in P body, ER, nuclei and cytoplasm
RBP_47_info[RBP_47_info$RBP.name=="DICER",c(19,23,26,28)]<-1

pdf("temp_fig/Fig5B.pdf",width=12,height=8)
par(mfrow=c(3,1),mar = c(5,2,2,20))
image(1:47,1:16,as.matrix(RBP_47_info[,18:3]),col=c("white","blue"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:16,labels = colnames(RBP_47_info)[18:3],tick = FALSE)
image(1:47,1:14,as.matrix(RBP_47_info[,32:19]),col=c("white","red"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:14,labels = colnames(RBP_47_info)[32:19],tick = FALSE)
image(1:47,1:6,as.matrix(RBP_47_info[38:33]),col=c("white","green"),bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:6,labels = colnames(RBP_47_info)[38:33],tick = FALSE)
axis(1,las=2,at = 1:47,labels = RBP_47_info$RBP.name,tick = FALSE,cex=0.5)
dev.off()


group_list<-list(group1=labels.r[c(1:7,18:20)])
group_list$group2<-labels.r[c(10:15)]
group_list$group3<-labels.r[c(24:25)]
group_list$group4<-labels.r[c(26:27)]
group_list$group5<-labels.r[c(28:31,45:47)]
group_list$group6<-labels.r[c(40:44)]

pdf("temp_fig/RBP_boundFLEXI_RBP_by_group.pdf",width=6,height=9)
col<-c("black","red","orange","skyblue")
par(mfrow=c(3,2))
par(pch=16,pty="s")
for (i in 1:6) {
  RBP_name<-group_list[[i]]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%RBP_name])
  FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:126){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=1,6]<-"D"
  RBP_plot[RBP_plot$Padj<=0.01 & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=1,6]<-"U"
  #RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="bonferroni")
  axis_max<-floor(max(RBP_plot[,3:4])/5)*5+5
  plot(RBP_plot[,c(3,4)],xlim=c(0,axis_max),ylim=c(0,axis_max),cex=1.5,
       bty="n",col=col[RBP_plot$col],
       ylab=paste0("Group ",i,"RBP bound FLEXIs (% RBP sites)"),xlab="ALl FLEXIs (% RBP sites)")
  #sig_cutoff<-(RBP_plot$Padj<=0.05 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  sig_cutoff<-(RBP_plot$Padj<=0.01 & (RBP_plot[,3]>=1 | RBP_plot[,4]>=1))
  sig_cutoffL<-(RBP_plot$Padj<=0.01 & (RBP_plot[,3]>=3 | RBP_plot[,4]>=3))
  text(RBP_plot[sig_cutoffL,3:4],cex=0.5,pos = 4,
       col=col[RBP_plot$col[sig_cutoffL]],
       labels = RBP_plot$RBP.name[sig_cutoffL])
  abline(0,1,col="red")
  legend(axis_max/2,axis_max,cex=0.5,
         text.col=col[RBP_plot$col[sig_cutoff& (RBP_plot$Cells<RBP_plot[,4])]],
         legend = RBP_plot$RBP.name[sig_cutoff & (RBP_plot$Cells<RBP_plot[,4])])
  legend(axis_max-2,axis_max-5,cex=0.5,
         text.col=col[RBP_plot$col[sig_cutoff& (RBP_plot$Cells>RBP_plot[,4])]],
         legend = RBP_plot$RBP.name[sig_cutoff & (RBP_plot$Cells>RBP_plot[,4])])
}
dev.off()
#Fig1G
pdf("temp_fig/Fig1G.pdf",height=8,width=6)

group_CPM<-dat_CPM[,c(1,33:36)]
group_FLEXI<-list(group1=unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[1]]]))
group_FLEXI$group2<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[2]]])
group_FLEXI$group3<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[3]]])
group_FLEXI$group4<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[4]]])
group_FLEXI$group5<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[5]]])
group_FLEXI$group6<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP%in%group_list[[6]]])
par(mfrow=c(3,2))
for (i in 1:6){
  plot(density(log10(rowMeans(group_CPM[group_CPM$ID%in%group_FLEXI[[i]],2:5]))),col="red",
       bty="n",xlim=c(-4,1),ylim=c(0,2),xlab="RPM",main=paste0("Complex ",i))
  lines(density(log10(rowMeans(group_CPM[!group_CPM$ID%in%group_FLEXI[[i]],2:5]))))
}
dev.off()