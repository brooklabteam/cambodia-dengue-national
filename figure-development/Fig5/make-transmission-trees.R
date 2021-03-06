rm(list=ls())

library(plyr)
library(dplyr)
library(ape)
library(ggplot2)
library(treeio)
library(lubridate)
library(geosphere)
library(epitools)
library(ggtree)
library(rBt)
#get transmission trees

#and load the metadata
homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#get tree
tree1 <- read.annot.beast(file = paste0(homewd, "/BEAST-tree/denv1-out/DENV1avg.tree"))

#plot(tree)
#and load the metadata

#and load the metadata
dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
dat$date <- as.Date(dat$date)#, format="%m/%d/%y")
mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-07-13"

#plot timetree
p1 <- ggtree(tree1, mrsd= mrsd.denv1) + theme_tree2()
#extract the dates of the internal nodes
tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

# this matrix gives you the node that is the most recent common ancestor between
# each of the tips identified as the column and row names
mrca.matrix <- mrca(tree1)

#then, link time to mrca
#first, transform the matrix
combine.df <- cbind.data.frame(index = seq(1, (nrow(mrca.matrix)^2), 1), node = as.vector(mrca.matrix))
combine.df$beast_name1 <- rep(colnames(mrca.matrix), nrow(mrca.matrix))
combine.df$beast_name2 <- rep(rownames(mrca.matrix), each = nrow(mrca.matrix))
head(combine.df)
tail(combine.df)
combine.df$time1 <- as.Date(sapply(strsplit(combine.df$beast_name1, "_"), function(x) x[[length(x)]] ))
combine.df$time2 <- as.Date(sapply(strsplit(combine.df$beast_name2, "_"), function(x) x[[length(x)]] ))
combine.df$time1 <- year(combine.df$time1) + (yday(combine.df$time1)/365)
combine.df$time2 <- year(combine.df$time2) + (yday(combine.df$time2)/365)

#now remove any from before 2019
combine.df=subset(combine.df, time1>=2019)
combine.df=subset(combine.df, time2>=2019)
#now flag the repeats
combine.df$self <- 0
combine.df$self[which(combine.df$beast_name1==combine.df$beast_name2)] <-1
combine.df <- combine.df[combine.df$self<1,]
combine.df <- dplyr::select(combine.df, -(self), -(index))
#now add node date
combine.df <- merge(combine.df, node.sub, by="node", all.x=T, sort=F)
combine.df$merge_name <- paste(combine.df$beast_name1, combine.df$beast_name2, sep="_")
head(combine.df)

#now add in the gis data
tree.dat.merge = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1")
#tree.dat.merge <- dplyr::select(tree.dat.merge, -(X))
tree.dat.merge = tree.dat.merge[tree.dat.merge$date>"2018-12-31",]#48  sequences
head(tree.dat.merge)
#and distance matrix - here, we just choose Jess's sequences
xy.dat <- dplyr::select(tree.dat.merge, long, lat)
dist.mat <- distm(xy.dat)

head(tree.dat.merge)
rownames(dist.mat) <- colnames(dist.mat) <- tree.dat.merge$tip_name

#now move them to a data table
next.df <- cbind.data.frame(index = seq(1, (nrow(dist.mat)^2), 1), distance = as.vector(dist.mat))
next.df$beast_name1 <- rep(colnames(dist.mat), nrow(dist.mat))
next.df$beast_name2 <- rep(rownames(dist.mat), each = nrow(dist.mat))

head(next.df)
next.df$self <- 0
next.df$self[which(next.df$beast_name1==next.df$beast_name2)] <-1
next.df <- next.df[next.df$self<1,]
next.df <- dplyr::select(next.df, -(self), -(index))

next.df$merge_name <- paste(next.df$beast_name1, next.df$beast_name2, sep="_")

#now merge on merge name
next.df <- dplyr::select(next.df, -(beast_name2), -(beast_name1))

new.dat <- merge(next.df, combine.df, by="merge_name", all.x=T, sort=F)
head(new.dat)



##now that you have all pairs, find the earlier date for each and compute the time to mrca
split.pairs <- dlply(new.dat, .(rownames(new.dat)))
get.tMRCA <- function(df){
  if(df$time1>df$time2){
    df$early_time <- df$time2
    df$early_seq <- df$beast_name2
    df$late_time <- df$time1
    df$late_seq <- df$beast_name1
  }else{
    df$early_time <- df$time1
    df$early_seq <- df$beast_name1
    df$late_time <- df$time2
    df$late_seq <- df$beast_name2
  }
  
  df <- dplyr::select(df, -(merge_name), -(beast_name1), -(beast_name2), -(time1), -(time2))
  
  df$accession_early <- sapply(strsplit(df$early_seq, split="_"), "[",1)
  names(df) <-c("distance", "mrca_node", "mrcatime", "pairtime1", "pairseq1", "pairtime2", "pairseq2", "accession_early")
  df <- dplyr::select(df, accession_early, pairseq1, pairtime1, pairseq2, pairtime2,distance,  mrca_node, mrcatime)
  df$tMRCA <- df$pairtime1-df$mrcatime
  df$merge_name = paste(df$pairseq1, df$pairseq2,sep = "_")
  df$g1 <- df$pairtime1-df$mrcatime
  #df$g1 <- round(df$g1,0)
  df$g2 <- df$pairtime2-df$mrcatime
  #df$g2 <- round(df$g2,0)
  df$g1[df$g2<1]<- 0
  df$evol_time = ((df$g2-df$g1)/2) +df$g1
  
  #but if the cases are not in the same season, should throw them out
  if(df$pairtime2-df$pairtime1>.5){
    df$evol_time <- NA
  }
  return(df)
}

#and apply across all the pairs
pair.out <- lapply(split.pairs, get.tMRCA)

pair.df1 <- data.table::rbindlist(pair.out)
#and delete the duplicates
pair.df <- pair.df1[!duplicated(pair.df1$merge_name),]
pair.df <- dplyr::select(pair.df, -(merge_name))
##then link all the metadata for the early sequence and save to make transmission trees
head(pair.df)

#redo seq
length(unique(pair.df$pairseq1)) #48
##then link all the metadata for the early sequence and save to make transmission trees



head(dat)
merge.dat <- dplyr::select(dat, accession_num, age, sex,  DENV.serotype, DENV.subtype)
names(merge.dat)[names(merge.dat)=="accession_num"] <- "accession_early"
merge.dat <- merge.dat[!is.na(merge.dat$sex),]
head(merge.dat)
head(pair.df)

#and merge with pairs
pair.DENV1 <- merge(pair.df, merge.dat, by="accession_early", all.x=T, sort=F)
head(pair.DENV1)
tail(pair.DENV1)

#and the pair subtype
id.sub <- dplyr::select(dat, tip_name, DENV.subtype)
names(id.sub) <- c("pairseq2", "pair_subtype")

pair.DENV1 <- merge(pair.DENV1, id.sub, by="pairseq2", all.x=T, sort=F)

head(pair.DENV1)
#and save this for transmission trees
write.csv(pair.DENV1, file =paste0(homewd, "/data/DENV1transTreeDat.csv"), row.names = F)



#now do the same for DENV2

#get tree
rm(list=ls())
homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#get tree
tree2 <- read.annot.beast(file = paste0(homewd, "/BEAST-tree/denv2-out/DENV2avg.tree"))

#and load the metadata
dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
dat$date <- as.Date(dat$date)#, format="%m/%d/%y")
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"]) #"2020-09-23"

#plot timetree
p1 <- ggtree(tree2, mrsd= mrsd.denv2) + theme_tree2()
#extract the dates of the internal nodes
tree.dat <- p1$data
node.sub <- dplyr::select(tree.dat, node, x)
names(node.sub) <-  c("node", "nodetime")

# this matrix gives you the node that is the most recent common ancestor between
# each of the tips identified as the column and row names
mrca.matrix <- mrca(tree2)

#then, link time to mrca
#first, transform the matrix
combine.df <- cbind.data.frame(index = seq(1, (nrow(mrca.matrix)^2), 1), node = as.vector(mrca.matrix))
combine.df$beast_name1 <- rep(colnames(mrca.matrix), nrow(mrca.matrix))
combine.df$beast_name2 <- rep(rownames(mrca.matrix), each = nrow(mrca.matrix))
head(combine.df)
combine.df$time1 <- as.Date(sapply(strsplit(combine.df$beast_name1, "_"), function(x) x[[length(x)]]))
combine.df$time2 <- as.Date(sapply(strsplit(combine.df$beast_name2, "_"), function(x) x[[length(x)]]))
combine.df$time1 <- year(combine.df$time1) + (yday(combine.df$time1)/365)
combine.df$time2 <- year(combine.df$time2) + (yday(combine.df$time2)/365)

#now remove any from before 2019
combine.df=subset(combine.df, time1>=2019)
combine.df=subset(combine.df, time2>=2019)
#now flag the repeats
combine.df$self <- 0
combine.df$self[which(combine.df$beast_name1==combine.df$beast_name2)] <-1
combine.df <- combine.df[combine.df$self<1,]
combine.df <- dplyr::select(combine.df, -(self), -(index))
#now add node date
combine.df <- merge(combine.df, node.sub, by="node", all.x=T, sort=F)
combine.df$merge_name <- paste(combine.df$beast_name1, combine.df$beast_name2, sep="_")
head(combine.df)
length(combine.df$node[is.na(combine.df$nodetime)])#0

#now add in the gis data
tree.dat.merge = subset(dat, !is.na(lat) & DENV.serotype=="DENV-2")
tree.dat.merge = tree.dat.merge[tree.dat.merge$date>"2018-12-31",]#50  sequences. all but one have coord



#and distance matrix - here, we just choose Jess's sequences
xy.dat <- dplyr::select(tree.dat.merge, long, lat)
dist.mat <- distm(xy.dat)

head(tree.dat.merge)
#tree.dat.merge <- dplyr::select(tree.dat.merge, -(X))
rownames(dist.mat) <- colnames(dist.mat) <- tree.dat.merge$tip_name

#now move them to a data table
next.df <- cbind.data.frame(index = seq(1, (nrow(dist.mat)^2), 1), distance = as.vector(dist.mat))
next.df$beast_name1 <- rep(colnames(dist.mat), nrow(dist.mat))
next.df$beast_name2 <- rep(rownames(dist.mat), each = nrow(dist.mat))

head(next.df)
next.df$self <- 0
next.df$self[which(next.df$beast_name1==next.df$beast_name2)] <-1
next.df <- next.df[next.df$self<1,]
next.df <- dplyr::select(next.df, -(self), -(index))

next.df$merge_name <- paste(next.df$beast_name1, next.df$beast_name2, sep="_")

#now merge on merge name
next.df <- dplyr::select(next.df, -(beast_name2), -(beast_name1))

new.dat <- merge(next.df, combine.df, by="merge_name", all.x=T, sort=F)
head(new.dat)
#new.dat <- new.dat[!duplicated(new.dat),]

length(unique(new.dat$beast_name1))
length(unique(new.dat$beast_name2))

##now that you have all pairs, find the earlier date for each and compute the time to mrca
split.pairs <- dlply(new.dat, .(rownames(new.dat)))
get.tMRCA <- function(df){
  if(df$time1>df$time2){
    df$early_time <- df$time2
    df$early_seq <- df$beast_name2
    df$late_time <- df$time1
    df$late_seq <- df$beast_name1
  }else{
    df$early_time <- df$time1
    df$early_seq <- df$beast_name1
    df$late_time <- df$time2
    df$late_seq <- df$beast_name2
  }
  
  df <- dplyr::select(df, -(merge_name), -(beast_name1), -(beast_name2), -(time1), -(time2))
  
  df$accession_early <- sapply(strsplit(df$early_seq, split="_"), "[",1)
  names(df) <-c("distance", "mrca_node", "mrcatime", "pairtime1", "pairseq1", "pairtime2", "pairseq2", "accession_early")
  df <- dplyr::select(df, accession_early, pairseq1, pairtime1, pairseq2, pairtime2,distance,  mrca_node, mrcatime)
  df$tMRCA <- df$pairtime1-df$mrcatime
  df$merge_name = paste(df$pairseq1, df$pairseq2,sep = "_")
  df$g1 <- df$pairtime1-df$mrcatime
  #df$g1 <- round(df$g1,0)
  df$g2 <- df$pairtime2-df$mrcatime
  #df$g2 <- round(df$g2,0)
  df$g1[df$g2<1]<- 0
  df$evol_time = ((df$g2-df$g1)/2) +df$g1
  
  #but if the cases are not in the same season, should throw them out
  if(df$pairtime2-df$pairtime1>.5){
    df$evol_time <- NA
  }
  return(df)
}

#and apply across all the pairs
pair.out <- lapply(split.pairs, get.tMRCA)

pair.df1 <- data.table::rbindlist(pair.out)
#and delete the duplicates
pair.df <- pair.df1[!duplicated(pair.df1$merge_name),]
pair.df <- dplyr::select(pair.df, -(merge_name))
##then link all the metadata for the early sequence and save to make transmission trees
head(pair.df)
length(unique(pair.df$pairseq1)) #50

head(dat)
merge.dat <- dplyr::select(dat, accession_num, age, sex, DENV.serotype, DENV.subtype)
names(merge.dat)[names(merge.dat)=="accession_num"] <- "accession_early"
merge.dat <- merge.dat[!is.na(merge.dat$sex),]
head(merge.dat)

#and add in the subtype of the second strain too

#and merge with pairs
pair.DENV2 <- merge(pair.df, merge.dat, by="accession_early", all.x=T, sort=F)
head(pair.DENV2)
tail(pair.DENV2)

#and the pair subtype
id.sub <- dplyr::select(dat, tip_name, DENV.subtype)
names(id.sub) <- c("pairseq2", "pair_subtype")

pair.DENV2 <- merge(pair.DENV2, id.sub, by="pairseq2", all.x=T, sort=F)

write.csv(pair.DENV2, file =paste0(homewd, "/data/DENV2transTreeDat.csv"), row.names = F)

#and join together 
rm(list=ls())
homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

denv1 <- read.csv(file=paste0(homewd, "/data/DENV1transTreeDat.csv"), header = T, stringsAsFactors = F)
denv2 <- read.csv(file=paste0(homewd, "/data/DENV2transTreeDat.csv"), header = T, stringsAsFactors = F)


all.denv <- rbind(denv1, denv2)
unique(all.denv$DENV.serotype)
unique(all.denv$DENV.subtype)
all.denv$DENV.subtype[all.denv$DENV.subtype=="Cosmopolitan"]<- "DENV-2-Cosmopolitan"
all.denv$pair_subtype[all.denv$pair_subtype=="Cosmopolitan"]<- "DENV-2-Cosmopolitan"
all.denv$DENV.subtype[all.denv$DENV.subtype=="Asian-1"]<- "DENV-2-Asian-1"
all.denv$pair_subtype[all.denv$pair_subtype=="Asian-1"]<- "DENV-2-Asian-1"
all.denv$paired <- paste0(all.denv$DENV.subtype, "/", all.denv$pair_subtype)
unique(all.denv$paired)
#subset(all.denv, is.na(DENV.serotype))
write.csv(all.denv, file =paste0(homewd, "/data/AllDENVtransTreeDat.csv"), row.names = F)

