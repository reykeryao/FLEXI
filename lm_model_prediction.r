gene_counts<-read.delim("Combined.corrected.sncRNA.counts")
temp<-gene_counts[gene_counts$Type=="snRNA",c(4,9:16)]
gene_for_correlation<-data.frame("U1"=colSums(temp[grep("U1\\b|U1[A-Z]",temp$Name),2:9]))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U2"=colSums(temp[grep("U2\\b|U2[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4"=colSums(temp[grep("U4\\b|U4[B-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U5"=colSums(temp[grep("U5\\b|U5[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6"=colSums(temp[grep("U6\\b|U6[B-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U7"=colSums(temp[grep("U7\\b|U7[A-Z]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U4ATAC"=colSums(temp[grep("U4[ATAC|atac]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("U6ATAC"=colSums(temp[grep("U6[ATAC|atac]",temp$Name),2:9])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SK"=colSums(gene_counts[gene_counts$Type=="7SK",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RN7SL"=colSums(gene_counts[gene_counts$Type=="7SL",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RPPH1"=colSums(gene_counts[gene_counts$Name=="RPPH1",9:16])))
gene_for_correlation<-cbind(gene_for_correlation,
                            data.frame("RMRP"=colSums(gene_counts[gene_counts$Name=="RMRP",9:16])))

#mappedreads and corrected mapped reads
correct_mapped_reads<-c(715.24152,277.8381,768.43375,368.4548,666.34185,362.6804,713.77529,230.0604)
names(correct_mapped_reads)<-rownames(gene_for_correlation)
gene_for_correlation<-gene_for_correlation/correct_mapped_reads
gene_for_correlation<-rbind(gene_for_correlation,"Copy"=c(1e6,5e5,2e5,2e5,4e5,4e3,2e3,2e3,2e5,5e5,2e5,1e5))
gene_for_correlation<-data.frame(t(gene_for_correlation))
gene_for_correlation$Type=c(rep("GU-AG splicing",5),"U7",rep("AU-AC splicing",2),"7SK","7SL","RNase P","MRP")
gene_for_correlation$Type<-as.factor(gene_for_correlation$Type)

lm_list<-list()
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i-1])
  lm.out <- lm(y ~ x)
  lm_list<-c(lm_list,list(lm.out))
}
for (i in 1:4){
  x<-log10(gene_for_correlation[,2*i])
  lm.out <- lm(y ~ x)
  lm_list<-c(lm_list,list(lm.out))
}
names(lm_list)<-c(paste0(c("HEK","Hela","UHRR","K562"),"Before"),
                  paste0(c("HEK","Hela","UHRR","K562"),"After"))
for (i in c(1:4)){
  lm.out<-lm_list[[i]]
  x<-data.frame(x=log10(c(0.01,0.99,1)))
  predict_y<-predict(lm.out,newdata =x ,interval ="confidence" ,level = 0.95)
  print(10^predict_y)
}


