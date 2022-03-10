library(tidyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(cluster)
library(purrr)
#read FLEIX--RBP table, set cutoff, rowMax ≤ cutoff will be discarded, remove duplicated RBP from teh same FLEXI
dat<-read.delim("4cell_plasma_FLEXI.tsv")
#this can be chanhged
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
RBP_fre<-merge(RBP_fre,RBP[,c(1,46)],by=1)
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

#make table for all FELXIs (>1 read in ≥ 1 dataset)
for (i in 1:length(RBP_list)) {
  RBP_name<-RBP_list[i]
  FLEXI_list<-unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name])
  FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
  RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
  RBP_plot[is.na(RBP_plot)]<-0
  RBP_plot$Padj<-1
  RBP_plot$Padj_withSign<-1
  RBP_plot$Type<-"N"
  R_sum<-colSums(RBP_plot[,3:4])
  for (j in 1:length(RBP_list)){
    RBP_plot[j,5]<-fisher.test(as.matrix(rbind(RBP_plot[j,3:4],R_sum)))$p.value
  }
  RBP_plot[,3:4]<-data.frame(prop.table(as.matrix(RBP_plot[,3:4]),margin = 2)*100)
  RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells>RBP_plot$Freq & RBP_plot$Cells>=percent_cutoff,7]<-"D"
  RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cells<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,7]<-"U"
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
RBP_clus<-RBP_clus[,3:dim(RBP_clus)[2]]
write.table(RBP_clus,"RBP_clus_froGower.txt",quote=F,sep="\t")
rownames(RBP_clus_log10)<-RBP_clus_log10$RBP.name
RBP_clus_log10<-RBP_clus_log10[,3:dim(RBP_clus_log10)[2]]
write.table(RBP_clus_log10,"RBP_clus_froHeatmap.txt",quote=F,sep="\t")

RBP_clus<-read.delim("RBP_clus_froGower.txt")
RBP_clus_log10<-read.delim("RBP_clus_froHeatmap.txt")
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


RBP_clus_log10.47<-RBP_clus_log10[,colnames(RBP_clus_log10)%in%RBP_47]
RBP_clus_log10.47<-RBP_clus_log10.47[rownames(RBP_clus_log10.47)%in%RBP_47_name,]
marker_text.47<-marker_text[,colnames(marker_text)%in%RBP_47]
marker_text.47<-marker_text.47[rownames(marker_text.47)%in%RBP_47_name,]
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

pdf("temp_fig/47RBP_heatmap.pdf",height=6,width=6)
labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
marker_text.47<-marker_text.47[labels.r,]
Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
        show_column_names = F,
        row_names_gp = gpar(fontsize = 7,col=col[col.r]),
        col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
        cluster_rows = F, cluster_columns =aggl.clust.r,
        cell_fun = function(i, j, x, y,w, h, col) {
          grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))})
dev.off()

#individual cell lines
for (i in 1:6){
  cell_specificFLEXI<-unique(dat$ID[dat[cell_list[i]]>cutoff])
  RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%cell_specificFLEXI]))
  colnames(RBP_fre)<-c("RBP","Cell")
  RBP_list<-as.character(RBP_fre$RBP)
  for (j in 1:length(RBP_list)) {
    RBP_name<-RBP_list[j]
    FLEXI_list<-intersect(unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name]),cell_specificFLEXI)
    FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    RBP_plot$Padj_withSign<-1
    RBP_plot$Type<-"N"
    R_sum<-colSums(RBP_plot[,2:3])
    # for larger number, fisher's exact won't work
    for (k in 1:length(RBP_list)){
      RBP_plot[k,4]<-fisher.test(as.matrix(rbind(RBP_plot[k,2:3],R_sum)))$p.value
    }
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cell>RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"D"
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cell<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
    RBP_plot[RBP_plot$Cell<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$Cell<=RBP_plot$Freq,4])
    RBP_plot[RBP_plot$Cell>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$Cell>RBP_plot$Freq,4])
    if (j==1){
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
  write.table(RBP_clus,paste0(cell_list[i],"_RBP_clus_froGower.txt"),quote=F,sep="\t")
  rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
  RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
  write.table(RBP_clus_log10,paste0(cell_list[i],"_RBP_clus_froHeatmap.txt"),quote=F,sep="\t")
}


#make plots heatmap
for (i in 1:6){
  
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
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  #gower.dist <- daisy(t(RBP_clus_log10), metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  pdf(paste0("temp_fig/",cell_list[i],"_47RBP_heatmap_cluster.pdf"),height=6,width=6)
  labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
  col.r<-RBP_col$col[match(labels.r,RBP_col$RBP)]
  RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
  marker_text.47<-marker_text.47[labels.r,]
  print(Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
          row_names_gp = gpar(fontsize = 7,col=col[col.r]),
          col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
          column_names_gp = gpar(fontsize = 7,col=col[col.r]),
          cluster_columns = aggl.clust.r,
          cluster_rows = F,
          cell_fun = function(i, j, x, y,w, h, col) {
            grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))}))
  dev.off()
}

###K562 different cutoff
cutoffOver1read_list<-list(unique(dat$ID[dat[cell_list[4]]>1]))
cutoffOver1read_list<-c(cutoffOver1read_list,list(unique(dat$ID[dat[cell_list[4]]/mapped_reads[8]>=0.01])))
name_list<-c("Over1","Over001RPM")
for (i in 1:2){
  cell_specificFLEXI<-cutoffOver1read_list[[i]]
  RBP_fre<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%cell_specificFLEXI]))
  colnames(RBP_fre)<-c("RBP","Cell")
  RBP_list<-as.character(RBP_fre$RBP)
  for (j in 1:length(RBP_list)) {
    RBP_name<-RBP_list[j]
    FLEXI_list<-intersect(unique(RBP_4cell_plasma$ID[RBP_4cell_plasma$RBP==RBP_name]),cell_specificFLEXI)
    FLEXI_RBP<-data.frame(table(RBP_4cell_plasma$RBP[RBP_4cell_plasma$ID%in%FLEXI_list]))
    RBP_plot<-merge(RBP_fre,FLEXI_RBP,by=1,all=T)
    RBP_plot[is.na(RBP_plot)]<-0
    RBP_plot$Padj<-1
    RBP_plot$Padj_withSign<-1
    RBP_plot$Type<-"N"
    R_sum<-colSums(RBP_plot[,2:3])
    # for larger number, fisher's exact won't work
    for (k in 1:length(RBP_list)){
      RBP_plot[k,4]<-fisher.test(as.matrix(rbind(RBP_plot[k,2:3],R_sum)))$p.value
    }
    RBP_plot$Padj<-p.adjust(RBP_plot$Padj,method="fdr")
    RBP_plot[,2:3]<-data.frame(prop.table(as.matrix(RBP_plot[,2:3]),margin = 2)*100)
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cell>RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"D"
    RBP_plot[RBP_plot$Padj<=pvalue_cutoff & RBP_plot$Cell<RBP_plot$Freq & RBP_plot$Freq>=percent_cutoff,6]<-"U"
    RBP_plot[RBP_plot$Cell<=RBP_plot$Freq,5]<- -log10(RBP_plot[RBP_plot$Cell<=RBP_plot$Freq,4])
    RBP_plot[RBP_plot$Cell>RBP_plot$Freq,5]<- log10(RBP_plot[RBP_plot$Cell>RBP_plot$Freq,4])
    if (j==1){
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
  write.table(RBP_clus,paste0(cell_list[4],name_list[i],"_RBP_clus_froGower.txt"),quote=F,sep="\t")
  rownames(RBP_clus_log10)<-RBP_clus_log10$RBP
  RBP_clus_log10<-RBP_clus_log10[,2:dim(RBP_clus_log10)[2]]
  write.table(RBP_clus_log10,paste0(cell_list[4],name_list[i],"_RBP_clus_froHeatmap.txt"),quote=F,sep="\t")
}


for (i in 1:2){
  RBP_clus<-read.delim(paste0(cell_list[4],name_list[i],"_RBP_clus_froGower.txt"))
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
  RBP_clus_log10<-read.delim(paste0(cell_list[4],name_list[i],"_RBP_clus_froHeatmap.txt"))
  
  RBP_clus<-RBP_clus[RBP_47_name]
  RBP_clus<-RBP_clus[RBP_47,]
  
  marker_text.47<-marker_text[RBP_47]
  marker_text.47<-marker_text.47[RBP_47_name,]
  
  RBP_clus_log10<-RBP_clus_log10[RBP_47]
  RBP_clus_log10<-RBP_clus_log10[RBP_47_name,]
  
  gower.dist <- daisy(RBP_clus, metric = c("gower"))
  #gower.dist <- daisy(t(RBP_clus_log10), metric = c("gower"))
  aggl.clust.r <- hclust(gower.dist, method = "complete")
  RBP_clus_log10_convertedvalue<-RBP_clus_log10
  for (j in 1:dim(RBP_clus_log10_convertedvalue)[2]){
    colname<-rownames(RBP_clus_log10_convertedvalue)[j]
    colname<-paste0(colname,"BF")
    RBP_clus_log10_convertedvalue[j,colname]<-NA
  }
  pdf(paste0("temp_fig/",cell_list[4],name_list[i],"_47RBP_heatmap_cluster.pdf"),height=6,width=6)
  labels.r<-substr(labels(aggl.clust.r),1,nchar(labels(aggl.clust.r))-2)
  col.r<-RBP_col$col[match(labels.r,RBP_col$RBP.name)]
  RBP_clus_log10_convertedvalue<-RBP_clus_log10_convertedvalue[labels.r,]
  marker_text.47<-marker_text.47[labels.r,]
  print(Heatmap(RBP_clus_log10_convertedvalue, na_col="black",
                row_names_gp = gpar(fontsize = 7,col=col[col.r]),
                col = colorRamp2(seq(-5,5,length = 3),c("blue", "#EEEEEE", "red")),
                column_names_gp = gpar(fontsize = 7,col=col[col.r]),
                cluster_columns = aggl.clust.r,
                cluster_rows = F,
                cell_fun = function(i, j, x, y,w, h, col) {
                  grid.text(marker_text.47[j, i], x, y,gp = gpar(fontsize = 7))}))
  dev.off()
}
