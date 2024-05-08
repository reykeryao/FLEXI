dat<-read.delim("4cell_plasma_FLEXI.tsv")
#this can be changed
cutoff<-0
percent_cutoff<-2
pvalue_cutoff<-0.05
col<-c("black","red","orange","skyblue")
#simple cutoff
FLEXI<-dat[rowMaxs(as.matrix(dat[,88:91]))>cutoff,]
RBP_4cell_plasma<-read.table("4cell_plasma_combined_RBP.info",col.names=c("ID","RBP"))
RBP_4cell_plasma<-RBP_4cell_plasma[RBP_4cell_plasma$ID%in%FLEXI$ID,]
RBP_4cell_plasma<-unique(RBP_4cell_plasma)
#make table for each subtype of cell
cell_list<-colnames(dat)[86:91]
RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP))
colnames(RBP_fre)<-c("RBP.name","Cells")
RBP$col<-(RBP$Splicing.regulation+RBP$Spliceosome)/3+RBP$microRNA.processing
RBP_fre<-merge(RBP_fre,RBP[,c(1,50)],by=1)
RBP_fre[RBP_fre$col>1,3]<-4
RBP_fre[RBP_fre$col==1,3]<-3
RBP_fre[RBP_fre$col<1 & RBP_fre$col>0,3]<-2
RBP_fre[RBP_fre$col==0,3]<-1
RBP_fre<-RBP_fre[,c(1,3,2)]
RBP_col<-RBP_fre[,1:2]
RBP_list<-sort(unique(RBP_4cell_plasma$RBP))
RBP_47<-data.frame(table(RBP_4cell_plasma$RBP))
RBP_47<-RBP_47[order(RBP_47$Freq,decreasing = T),]
RBP_47_name<-as.character(RBP_47$Var1[RBP_47$Freq>=30])
RBP_47_name<-RBP_47_name[7:length(RBP_47_name)]
RBP_47<-paste0(RBP_47_name,"BF")


all_intron_RBP<-read.table(gzfile("all_intron_RBP_inter.info.gz"),col.names=c("FLEXI","Len","RBP"))
all_intron_RBP<-unique(all_intron_RBP)
all_intron_RBP_short<-all_intron_RBP[all_intron_RBP$Len<=300,c(1,3)]
all_intron_RBP<-all_intron_RBP[all_intron_RBP$Len>300,c(1,3)]
#make table for other short intron ≤ 300 nt
All_FLEXI<-unique(dat$ID[rowMaxs(as.matrix(dat[,88:91]))>cutoff])
all_intron_RBP_short<-all_intron_RBP_short[!all_intron_RBP_short$FLEXI%in%All_FLEXI,]
RBP_list<-unique(all_intron_RBP_short$RBP)
RBP_fre<-data.frame(table(all_intron_RBP_short$RBP))
colnames(RBP_fre)<-c("RBP","AllShortIntron")
for (i in 1:length(RBP_list)) {
  RBP_name<-RBP_list[i]
  FLEXI_list<-unique(all_intron_RBP_short$FLEXI[all_intron_RBP_short$RBP==RBP_name])
  FLEXI_RBP<-data.frame(table(all_intron_RBP_short$RBP[all_intron_RBP_short$FLEXI%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Padj_withSign<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,2:3])
  # for larger number, fisher's exact won't work
  for (j in 1:length(RBP_list)){
    RBP_plot[j,4]<-fisher.test(as.matrix(rbind(RBP_plot[j,2:3],R_sum)))$p.value
  }
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$AllShortIntron>RBP_plot$Freq & RBP_plot$Cells>=percent_cutoff,6]<-"D"
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$AllShortIntron<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
  RBP_plot[RBP_plot$AllShortIntron<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$AllShortIntron<=RBP_plot$Freq,4])
  RBP_plot[RBP_plot$AllShortIntron>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$AllShortIntron>RBP_plot$Freq,4])
  if (i==1){
    RBP_clus<-RBP_plot[,c(1,6)]
    colnames(RBP_clus)[2]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-RBP_plot[,c(1,5)]
    colnames(RBP_clus_log10)[2]<-paste0(RBP_name,"BF")
    
  } else {
    RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,6)],by=1)
    colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,5)],by=1)
    colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
  }
}
rownames(RBP_clus)<-RBP_clus$RBP
RBP_clus<-RBP_clus[,c(2:dim(RBP_clus)[2])]
write.table(RBP_clus,"AllShortIntronRBP_clus_froGower.txt",quote=F,sep="\t")
rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
write.table(RBP_clus_log10,"AllShortIntronRBP_clus_froHeatmap.txt",quote=F,sep="\t")



RBP_fre<-data.frame(table(all_intron_RBP$RBP))
colnames(RBP_fre)<-c("RBP","AllLongIntron")
RBP_list<-unique(all_intron_RBP$RBP)
for (i in 1:length(RBP_list)) {
  RBP_name<-RBP_list[i]
  FLEXI_list<-unique(all_intron_RBP$FLEXI[all_intron_RBP$RBP==RBP_name])
  FLEXI_RBP<-data.frame(table(all_intron_RBP$RBP[all_intron_RBP$FLEXI%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Padj_withSign<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,2:3])
  # for larger number, fisher's exact won't work
  for (j in 1:152){
    RBP_plot[j,4]<-fisher.test(as.matrix(rbind(RBP_plot[j,2:3],R_sum)))$p.value
  }
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$AllLongIntron>RBP_plot$Freq & RBP_plot$Cells>=percent_cutoff,6]<-"D"
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$AllLongIntron<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
  RBP_plot[RBP_plot$AllLongIntron<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$AllLongIntron<=RBP_plot$Freq,4])
  RBP_plot[RBP_plot$AllLongIntron>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$AllLongIntron>RBP_plot$Freq,4])
  if (i==1){
    RBP_clus<-RBP_plot[,c(1,6)]
    colnames(RBP_clus)[2]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-RBP_plot[,c(1,5)]
    colnames(RBP_clus_log10)[2]<-paste0(RBP_name,"BF")
    
  } else {
    RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,6)],by=1)
    colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
    RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,5)],by=1)
    colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
  }
}
rownames(RBP_clus)<-RBP_clus$RBP
RBP_clus<-RBP_clus[,2:dim(RBP_clus)[2]]
write.table(RBP_clus,"AllLongIntronRBP_clus_froGower.txt",quote=F,sep="\t")
rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
write.table(RBP_clus_log10,"AllLongIntronRBP_clus_froHeatmap.txt",quote=F,sep="\t")

#make plots heatmap
RBP_clus<-read.delim("AllShortIntronRBP_clus_froGower.txt")
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
for (i in 1:dim(RBP_clus)[2]){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}
RBP_clus_log10<-read.delim("AllShortIntronRBP_clus_froHeatmap.txt")

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
col<-c("black","red","orange","skyblue")

for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

RBP_clus<-RBP_clus[RBP_47_name]
RBP_clus<-RBP_clus[RBP_47,]
marker_text.47<-marker_text[RBP_47]
marker_text.47<-marker_text.47[RBP_47_name,]
RBP_clus_log10<-RBP_clus_log10[RBP_47]
RBP_clus_log10<-RBP_clus_log10[RBP_47_name,]

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

pdf("temp_fig/AllShortIntron47RBP_heatmap.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
marker_text.47<-marker_text.47[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        column_names_gp = gpar(fontsize = 7,col=col[col.r]),
        cluster_rows = F, cluster_columns =aggl.clust.r,
        cell_fun = function(i, j, x, y,w, h, col) {
          grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))})
dev.off()
# make long

RBP_clus<-read.delim("AllLongIntronRBP_clus_froGower.txt")
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
for (i in 1:dim(RBP_clus)[2]){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}
RBP_clus_log10<-read.delim("AllLongIntronRBP_clus_froHeatmap.txt")

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

RBP_clus<-RBP_clus[RBP_47_name]
RBP_clus<-RBP_clus[RBP_47,]
marker_text.47<-marker_text[RBP_47]
marker_text.47<-marker_text.47[RBP_47_name,]

RBP_clus_log10<-RBP_clus_log10[RBP_47]
RBP_clus_log10<-RBP_clus_log10[RBP_47_name,]

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

pdf("temp_fig/AllLongIntron47RBP_heatmap.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
marker_text.47<-marker_text.47[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        column_names_gp = gpar(fontsize = 7,col=col[col.r]),
        cluster_rows = F
        ,cluster_columns =aggl.clust.r,
        cell_fun = function(i, j, x, y,w, h, col) {
        grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))})
dev.off()




#make plots heatmap
RBP_clus<-read.delim("AllShortIntronRBP_clus_froGower.txt")
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
for (i in 1:dim(RBP_clus)[2]){
  RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
}
RBP_clus_log10<-read.delim("AllShortIntronRBP_clus_froHeatmap.txt")

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
col<-c("black","red","orange","skyblue")

for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

RBP_clus<-RBP_clus[RBP_47_name]
RBP_clus<-RBP_clus[RBP_47,]
marker_text.47<-marker_text[RBP_47]
marker_text.47<-marker_text.47[RBP_47_name,]
RBP_clus_log10<-RBP_clus_log10[RBP_47]
RBP_clus_log10<-RBP_clus_log10[RBP_47_name,]

gower.dist <- daisy(RBP_clus, metric = c("gower"))
aggl.clust.r <- hclust(gower.dist, method = "complete")
RBP_clus_log10_convertedvalue<-RBP_clus_log10
for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
  colname<-rownames(RBP_clus_log10_convertedvalue)[j]
  colname<-paste0(colname,"BF")
  RBP_clus_log10_convertedvalue[j,colname]<-NA
}

pdf("temp_fig/AllShortIntron47RBP_heatmap.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
marker_text.47<-marker_text.47[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        column_names_gp = gpar(fontsize = 7,col=col[col.r]),
        cluster_rows = F, cluster_columns =aggl.clust.r,
        cell_fun = function(i, j, x, y,w, h, col) {
          grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))})
dev.off()


# #make random k=10 pick n=10k other short intron heatmap
# random samplying 3000 introns and repeat 10 times
# get all intron (short intron) list
#ShortIntronList<-read.delim("GRCh38.93.intron_deduped.tsv")
#ShortIntronList<-setdiff(ShortIntronList$ID,All_FLEXI)
#ShortIntronList<-unique(ShortIntronList)

# other strategy, pick 2k other short introns with RBP site
ShortIntronList<-unique(all_intron_RBP_short$FLEXI)
random_pick<-2000
for (k in 1:10){
  sample8000<-sample(ShortIntronList,random_pick)
  sample8000<-all_intron_RBP_short[all_intron_RBP_short$FLEXI%in%sample8000,]
  RBP_fre<-data.frame(table(sample8000$RBP))
  colnames(RBP_fre)<-c("RBP","Sample")
  RBP_list<-as.character(RBP_fre$RBP)
  for (i in 1:length(RBP_list)) {
    RBP_name<-RBP_list[i]
    FLEXI_list<-unique(sample8000$FLEXI[sample8000$RBP==RBP_name])
    FLEXI_RBP<-data.frame(table(sample8000$RBP[sample8000$FLEXI%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    RBP_plot$Padj_withSign<-1
    RBP_plot$Type<-"N"
    R_sum<-colSums(RBP_plot[,2:3])
    # for larger number, fisher's exact won't work
    for (j in 1:length(RBP_list)){
      RBP_plot[j,4]<-fisher.test(as.matrix(rbind(RBP_plot[j,2:3],R_sum)))$p.value
    }
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample>RBP_plot$Freq & RBP_plot$Sample>=percent_cutoff,6]<-"D"
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
    RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,4])
    RBP_plot[RBP_plot$Sample>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$Sample>RBP_plot$Freq,4])
    if (i==1){
      RBP_clus<-RBP_plot[,c(1,6)]
      colnames(RBP_clus)[2]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-RBP_plot[,c(1,5)]
      colnames(RBP_clus_log10)[2]<-paste0(RBP_name,"BF")
      
    } else {
      RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,6)],by=1)
      colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,5)],by=1)
      colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
    }
  }
  rownames(RBP_clus)<-RBP_clus$RBP
  RBP_clus<-RBP_clus[,c(2:dim(RBP_clus)[2])]
  rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
  RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
  if (k==1){
    clustObj<-list(RBP_clus)
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }  else {
    clustObj<-c(clustObj,list(RBP_clus))
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }
}
saveRDS(clustObj,paste0("OtherShortIntronRandompick",random_pick,".obj"))

clustObj<-readRDS(paste0("OtherShortIntronRandompick",random_pick,".obj"))
pdf(paste0("temp_fig/Randompick",random_pick,"Short_heatmap.pdf"),height=6,width=6,onefile = T)
for (k in 1:10){
  RBP_clus<-clustObj[[2*k-1]]
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
  for (i in 1:dim(RBP_clus)[2]){
    RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
  }
  RBP_clus_log10<-clustObj[[2*k]]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  
  RBP_clus<-RBP_clus[,colnames(RBP_clus)%in%RBP_47_name]
  RBP_clus<-RBP_clus[rownames(RBP_clus)%in%RBP_47,]
  marker_text.47<-marker_text[,colnames(marker_text)%in%RBP_47]
  marker_text.47<-marker_text.47[rownames(marker_text)%in%RBP_47_name,]
  
  RBP_clus_log10<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
  RBP_clus_log10<-RBP_clus_log10[rownames(RBP_clus_log10)%in%RBP_47_name,]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
  col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
  RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
  marker_text.47<-marker_text.47[labels.r,]
  print(Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
          row_names_gp = gpar(fontsize = 7,col=col[col.r]),
          col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
          column_names_gp = gpar(fontsize = 7,col=col[col.r]),
          cluster_rows = F
          ,cluster_columns =aggl.clust.r,
          cell_fun = function(i, j, x, y,w, h, col) {
            grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))}))
}
dev.off()

#LongIntronList<-read.delim(gzfile("allIntronOver300.list.gz"),header=F)
#LongIntronList<-LongIntronList$V1
LongIntronList<-unique(all_intron_RBP$FLEXI)
for (k in 1:10){
  sample8000<-sample(LongIntronList,random_pick)
  sample8000<-all_intron_RBP[all_intron_RBP$FLEXI%in%sample8000,]
  RBP_fre<-data.frame(table(sample8000$RBP))
  colnames(RBP_fre)<-c("RBP","Sample")
  RBP_list<-as.character(RBP_fre$RBP)
  for (i in 1:length(RBP_list)) {
    RBP_name<-RBP_list[i]
    FLEXI_list<-unique(sample8000$FLEXI[sample8000$RBP==RBP_name])
    FLEXI_RBP<-data.frame(table(sample8000$RBP[sample8000$FLEXI%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    RBP_plot$Padj_withSign<-1
    RBP_plot$Type<-"N"
    R_sum<-colSums(RBP_plot[,2:3])
    # for larger number, fisher's exact won't work
    for (j in 1:length(RBP_list)){
      RBP_plot[j,4]<-fisher.test(as.matrix(rbind(RBP_plot[j,2:3],R_sum)))$p.value
    }
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample>RBP_plot$Freq & RBP_plot$Sample>=percent_cutoff,6]<-"D"
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
    RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,4])
    RBP_plot[RBP_plot$Sample>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$Sample>RBP_plot$Freq,4])
    if (i==1){
      RBP_clus<-RBP_plot[,c(1,6)]
      colnames(RBP_clus)[2]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-RBP_plot[,c(1,5)]
      colnames(RBP_clus_log10)[2]<-paste0(RBP_name,"BF")
      
    } else {
      RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,6)],by=1)
      colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,5)],by=1)
      colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
    }
  }
  rownames(RBP_clus)<-RBP_clus$RBP
  RBP_clus<-RBP_clus[,c(2:dim(RBP_clus)[2])]
  rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
  RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
  if (k==1){
    clustObj<-list(RBP_clus)
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }  else {
    clustObj<-c(clustObj,list(RBP_clus))
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }
}
saveRDS(clustObj,paste0("LongIntronRandompick",random_pick,".obj"))

clustObj<-readRDS(paste0("LongIntronRandompick",random_pick,".obj"))
pdf(paste0("temp_fig/Randompick",random_pick,"Long_heatmap.pdf"),height=6,width=6,onefile = T)
for (k in 1:10){
  RBP_clus<-clustObj[[2*k-1]]
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
  for (i in 1:dim(RBP_clus)[2]){
    RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
  }
  RBP_clus_log10<-clustObj[[2*k]]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  
  RBP_clus<-RBP_clus[,colnames(RBP_clus)%in%RBP_47_name]
  RBP_clus<-RBP_clus[rownames(RBP_clus)%in%RBP_47,]
  marker_text.47<-marker_text[,colnames(marker_text)%in%RBP_47]
  marker_text.47<-marker_text.47[rownames(marker_text)%in%RBP_47_name,]
  
  RBP_clus_log10<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
  RBP_clus_log10<-RBP_clus_log10[rownames(RBP_clus_log10)%in%RBP_47_name,]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
  col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
  RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
  marker_text.47<-marker_text.47[labels.r,]
  print(Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
                row_names_gp = gpar(fontsize = 7,col=col[col.r]),
                col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
                column_names_gp = gpar(fontsize = 7,col=col[col.r]),
                cluster_rows = F
                ,cluster_columns =aggl.clust.r,
                cell_fun = function(i, j, x, y,w, h, col) {
                  grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))}))
}
dev.off()

##exon random pick
Exon<-read.delim(gzfile("All_exon_inter_RBP152_collapsed.info.gz"),header=F)
colnames(Exon)<-c("ID","RBP")
ExonList<-unique(Exon$ID)
for (k in 1:10){
  sample8000<-sample(ExonList,random_pick)
  sample8000<-Exon[Exon$ID%in%sample8000,]
  RBP_fre<-data.frame(table(sample8000$RBP))
  colnames(RBP_fre)<-c("RBP","Sample")
  RBP_list<-as.character(RBP_fre$RBP)
  for (i in 1:length(RBP_list)) {
    RBP_name<-RBP_list[i]
    FLEXI_list<-unique(sample8000$ID[sample8000$RBP==RBP_name])
    FLEXI_RBP<-data.frame(table(sample8000$RBP[sample8000$ID%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    RBP_plot$Padj_withSign<-1
    RBP_plot$Type<-"N"
    R_sum<-colSums(RBP_plot[,2:3])
    # for larger number, fisher's exact won't work
    for (j in 1:length(RBP_list)){
      RBP_plot[j,4]<-fisher.test(as.matrix(rbind(RBP_plot[j,2:3],R_sum)))$p.value
    }
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample>RBP_plot$Freq & RBP_plot$Sample>=percent_cutoff,6]<-"D"
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Sample<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
    RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$Sample<=RBP_plot$Freq,4])
    RBP_plot[RBP_plot$Sample>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$Sample>RBP_plot$Freq,4])
    if (i==1){
      RBP_clus<-RBP_plot[,c(1,6)]
      colnames(RBP_clus)[2]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-RBP_plot[,c(1,5)]
      colnames(RBP_clus_log10)[2]<-paste0(RBP_name,"BF")
      
    } else {
      RBP_clus<-merge(RBP_clus,RBP_plot[,c(1,6)],by=1)
      colnames(RBP_clus)[dim(RBP_clus)[2]]<-paste0(RBP_name,"BF")
      RBP_clus_log10<-merge(RBP_clus_log10,RBP_plot[,c(1,5)],by=1)
      colnames(RBP_clus_log10)[dim(RBP_clus_log10)[2]]<-paste0(RBP_name,"BF")
    }
  }
  rownames(RBP_clus)<-RBP_clus$RBP
  RBP_clus<-RBP_clus[,c(2:dim(RBP_clus)[2])]
  rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
  RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
  if (k==1){
    clustObj<-list(RBP_clus)
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }  else {
    clustObj<-c(clustObj,list(RBP_clus))
    clustObj<-c(clustObj,list(RBP_clus_log10))
  }
}
saveRDS(clustObj,paste0("ExonRandompick",random_pick,".obj"))

clustObj<-readRDS(paste0("ExonRandompick",random_pick,".obj"))
pdf(paste0("temp_fig/Randompick",random_pick,"Exon_heatmap.pdf"),height=6,width=6,onefile = T)
for (k in 1:10){
  RBP_clus<-clustObj[[2*k-1]]
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
  for (i in 1:dim(RBP_clus)[2]){
    RBP_clus[,i]<-factor(RBP_clus[,i],levels = c("D","N","U"))
  }
  RBP_clus_log10<-clustObj[[2*k]]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  
  RBP_clus<-RBP_clus[,colnames(RBP_clus)%in%RBP_47_name]
  RBP_clus<-RBP_clus[rownames(RBP_clus)%in%RBP_47,]
  marker_text.47<-marker_text[,colnames(marker_text)%in%RBP_47]
  marker_text.47<-marker_text.47[rownames(marker_text)%in%RBP_47_name,]
  
  RBP_clus_log10<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
  RBP_clus_log10<-RBP_clus_log10[rownames(RBP_clus_log10)%in%RBP_47_name,]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
  col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
  RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
  marker_text.47<-marker_text.47[labels.r,]
  print(Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
                row_names_gp = gpar(fontsize = 7,col=col[col.r]),
                col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
                column_names_gp = gpar(fontsize = 7,col=col[col.r]),
                cluster_rows = F
                ,cluster_columns =aggl.clust.r,
                cell_fun = function(i, j, x, y,w, h, col) {
                  grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))}))
}
dev.off()



#K562 image
library(readxl)
icol<-c("white",brewer.pal(8,"Set1"),brewer.pal(12,"Paired"))
pdf("temp_fig/K562_Over1_001RPM_image.pdf",width=12,height=8)
par(mfrow=c(2,1),mar = c(5,2,2,20))

RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 3))
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:36]
RBP_clus_image[is.na(RBP_clus_image)]<-0
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
image(1:35,1:21,as.matrix(RBP_clus_image[,21:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:21,labels = colnames(RBP_clus_image)[21:1],tick = FALSE,cex=0.5)
axis(1,las=2,at = 1:35,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)

RBP_clus_image<-data.frame(read_xlsx("RBP_clus_image.xlsx",sheet= 4))
rownames(RBP_clus_image)<-RBP_clus_image$Cell
RBP_clus_image<-RBP_clus_image[,2:26]
RBP_clus_image[is.na(RBP_clus_image)]<-0
RBP_clus_image<-RBP_clus_image+1
RBP_clus_image<-data.frame(t(RBP_clus_image))
RBP_clus_image$RBP.name<-rownames(RBP_clus_image)
image(1:25,1:21,as.matrix(RBP_clus_image[,21:1]),col=icol,bty="n",axes=F,xlab=NA,ylab=NA)
axis(4,las=2,at = 1:21,labels = colnames(RBP_clus_image)[21:1],tick = FALSE,cex=0.5)
axis(1,las=2,at = 1:25,labels = rownames(RBP_clus_image),tick = FALSE,cex=0.5)

dev.off()