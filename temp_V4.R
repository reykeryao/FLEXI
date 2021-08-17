#fig1A various cutoff, 1,3,5,10
cut_off=1
set_1 <- as.character(dat$ID[dat$UHRR>=cut_off])
set_2 <- as.character(dat$ID[dat$K562>=cut_off])
set_3 <- as.character(dat$ID[dat$HEK>=cut_off])
set_4 <- as.character(dat$ID[dat$Hela>=cut_off])
set_5 <- as.character(dat$ID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK 293T"=set_3,
             "Hela S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
postscript("temp_fig//Fig1A_10.eps",height=4,width=8)
UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
dev.off()
#FLEXI host genes
set_1 <- unique(dat$GID[dat$UHRR>=cut_off])
set_2 <- unique(dat$GID[dat$K562>=cut_off])
set_3 <- unique(dat$GID[dat$HEK>=cut_off])
set_4 <- unique(dat$GID[dat$Hela>=cut_off])
set_5 <- unique(dat$GID[dat$Plasma>=cut_off])
set <- list ("UHRR"=set_1,
             "K-562"=set_2,"HEK 293T"=set_3,
             "Hela S3"=set_4,"Plasma"=set_5)
m = make_comb_mat(set)
postscript("temp_fig//Fig1AGene_10.eps",height=4,width=8)
UpSet(m,set_order=c("UHRR","K-562","HEK 293T","Hela S3","Plasma"),
      comb_col = c("tomato","royalblue1","goldenrod","orchid","black")[comb_degree(m)])
dev.off()
rm(list=c("set_1","set_2","set_3","set_4","set_5","set","m","cut_off"))


