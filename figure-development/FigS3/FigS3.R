rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"

setwd(paste0(homewd, "/figure-development/FigS3"))

#load the DENV1 tree
treeD1 <-  read.tree(file = paste0(homewd, "figure-development/FigS3/denv1ML/T3.raxml.supportFBP"))

#and DENV2
treeD2 <-  read.tree(file = paste0(homewd, "figure-development/FigS3/denv2ML/T3.raxml.supportFBP"))
#root it


rooted.D1 <- root(treeD1, which(treeD1$tip.label == "NC_002640_DENV4"))
#take a quick look in base R
plot(rooted.D1)

#and DENV2
rooted.D2 <- root(treeD2, which(treeD2$tip.label == "NC_002640_DENV4"))
#take a quick look in base R
plot(rooted.D2)

#load tree data prepared from elsewhere
dat <- read.csv(file = paste0(homewd, "figure-development/FigS3/ML-Sequences.csv"))
head(dat)
#check subgroup names
unique(dat$Subclade[dat$Serotype=="DENV-1"])
unique(dat$Subclade[dat$Serotype=="DENV-2"])

#add denv4
dat.add <- c("NC_002640", NA, NA, "DENV-4", "outgroup")
dat <- rbind(dat, dat.add)

dat.denv1 <- subset(dat, Serotype!="DENV-2")

tree.denv1 <- cbind.data.frame(tip_label = rooted.D1$tip.label)
tree.denv1$Accession <- sapply(strsplit(tree.denv1$tip_label, split = "_"), function(x) x[[1]])

dat.all.denv1 <- merge(tree.denv1, dat.denv1, by ="Accession", all.x = T, sort=F)
dat.all.denv1$Accession[dat.all.denv1$tip_label=="NC_002640_DENV4"] <- "NC_002640"
dat.all.denv1$Serotype[dat.all.denv1$tip_label=="NC_002640_DENV4"] <- "DENV-4"

dat.all.denv1 <- dplyr::select(dat.all.denv1, tip_label, Accession, Locality, Serotype, Subclade)
unique(dat.all.denv1$Subclade)
subset(dat.all.denv1, is.na(Subclade))
dat.all.denv1$Subclade[dat.all.denv1$Serotype=="DENV-4"] <- "Outgroup"

colz1 = c('Genotype-I' = "tomato",
          'Genotype-II' = "royalblue", 
          'Genotype-III' ="darkgoldenrod1",
          'Old-Cambodia'= "magenta", 
          'New-Cambodia' = "darkorchid1", 
          'Outgroup' = "black")

#take a glance
p1 <- ggtree(rooted.D1) %<+% dat.all.denv1 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz1) + theme(legend.position = c(.2,.85), 
                                          legend.title = element_blank())

p1 #looks great - change to clade label and collapse the outgroup

#now shrink the long branch length for the outgroup
sort(rooted.D1$edge.length) #.03 is the next longest
rooted.D1$edge.length[rooted.D1$edge.length==max(rooted.D1$edge.length)] <- .15

#and add cladebars
#here for genotype II
gen2 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "DQ285561_Seychelles_2004" ),which(rooted.D1$tip.label == "AB204803_Micronesia_2004"))
gen3 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "EU081258_Singapore_2005" ),which(rooted.D1$tip.label == "EF122231_French_Guiana_2006"))
gen1 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "AY726550_Myanmar_2001" ),which(rooted.D1$tip.label == "OK159947_Cambodia_2019"))


dat.all.denv1$Subclade
#dat.all.denv1$Subclade[is.na(dat.all.denv1$Subclade)] <- "Outgroup"
dat.all.denv1$Subclade <- factor(dat.all.denv1$Subclade, 
                                 levels=c("Genotype-I", "Genotype-II", "Genotype-III",
                                          "Old-Cambodia", "New-Cambodia", "Outgroup")) 

#take a glance
pA <- ggtree(rooted.D1) %<+% dat.all.denv1 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz1) + theme(legend.position = c(.13,.83), 
                                          legend.background = element_rect(color="black"),
                                          legend.title = element_blank()) +
  geom_cladelabel(node=gen2, label="Genotype II", color="royalblue", offset = .0335) +
  geom_cladelabel(node=gen3, label="Genotype III", color="darkgoldenrod1",  offset = .082, extend = c(0,5)) +
  geom_cladelabel(node=gen1, label="Genotype I", color="tomato",  offset = .018) +
  xlim(c(0,.17)) +geom_treescale(fontsize=4, x=.005,y=110, linesize = .5) 

pA


#and do DENV2

dat.denv2 <- subset(dat, Serotype!="DENV-1")

tree.denv2 <- cbind.data.frame(tip_label = rooted.D2$tip.label)
tree.denv2$Accession <- sapply(strsplit(tree.denv2$tip_label, split = "_"), function(x) x[[1]])

dat.all.denv2 <- merge(tree.denv2, dat.denv2, by ="Accession", all.x = T, sort=F)
dat.all.denv2$Accession[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "NC_002640"
dat.all.denv2$Serotype[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "DENV-4"
dat.all.denv2$Subclade[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "Outgroup"

dat.all.denv2 <- dplyr::select(dat.all.denv2, tip_label, Accession, Locality, Serotype, Subclade)


#take a glance
p2 <- ggtree(rooted.D2) %<+% dat.all.denv2 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) + #scale_fill_manual(values=colz1) + 
  theme(legend.position = c(.2,.85), 
  legend.title = element_blank())

p2 #looks great - change to clade label and collapse the outgroup

unique(dat.all.denv2$Subclade)

#now shrink the long branch length for the outgroup
sort(rooted.D2$edge.length) #.07 is the next longest
rooted.D2$edge.length[rooted.D2$edge.length==max(rooted.D2$edge.length)] <- .14



dat.all.denv2$Subclade[dat.all.denv2$Subclade=="Cosmopolitan-I-C-1B" |dat.all.denv2$Subclade=="Cosmopolitan-I-C-1A"] <- "Cosmopolitan-I"

dat.all.denv2$Subclade <- factor(dat.all.denv2$Subclade, 
                                 levels=c("American", 
                                          "Asian-American", 
                                          "Asian-I", 
                                          "Asian-II", 
                                          "Cosmopolitan-I",
                                          "Cosmopolitan-II", 
                                          "Cosmopolitan-III",
                                          "Old-Cambodia",
                                          "New-Cambodia",
                                          "Outgroup")) 


colz2 = c('American'= "cornflowerblue", 
          'Asian-American' = "darkorange1", 
          'Asian-I' = "forestgreen", 
          'Asian-II' = "darkgoldenrod1", 
          'Cosmopolitan-I' = "tomato",
          'Cosmopolitan-II' = "mediumseagreen", 
          'Cosmopolitan-III' = "royalblue",
          'Old-Cambodia'="magenta",
          'New-Cambodia'="darkorchid1",
          'Outgroup' = "black")

cosmoIII <- MRCA(rooted.D2, which(rooted.D2$tip.label == "EU179857_Brunei_2005" ),which(rooted.D2$tip.label == "OL414753_Cambodia_2020"))
cosmoII <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ398259_Indonesia_1976" ),which(rooted.D2$tip.label == "GQ398258_Indonesia_1975"))
cosmoI <- MRCA(rooted.D2, which(rooted.D2$tip.label == "JX475906_India_2009" ),which(rooted.D2$tip.label == "KF041234_Pakistan_2011"))
American <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ868592_Columbia_1986" ),which(rooted.D2$tip.label == "HM582099_Fiji_1971"))
AsAm <- MRCA(rooted.D2, which(rooted.D2$tip.label == "HQ999999_Guatemala_2009" ),which(rooted.D2$tip.label == "FJ639700_Cambodia_2002"))
Asian1 <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ868591_Thailand_1964" ),which(rooted.D2$tip.label == "OL414721_Cambodia_2019"))
Asian2 <- MRCA(rooted.D2, which(rooted.D2$tip.label == "HQ891024_Taiwan_2008" ),which(rooted.D2$tip.label == "JF730050_Puerto_Rico_2007"))



#take a glance
pB <- ggtree(rooted.D2) %<+% dat.all.denv2 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz2) + theme(legend.position = c(.13,.76), 
                                          legend.background = element_rect(color="black"),
                                          legend.title = element_blank()) +
  geom_cladelabel(node=cosmoIII, label="Cosmopolitan III", color="royalblue", offset = .039, extend = c(0,8), fontsize = 3.5) +
  geom_cladelabel(node=cosmoII, label="Cosmopolitan II", color="mediumseagreen",  offset = .0795, extend = c(1,1), fontsize = 3.5) +
  geom_cladelabel(node=cosmoI, label="Cosmopolitan I", color="tomato",  offset = .0465, fontsize = 3.5) +
  geom_cladelabel(node=American, label="American", color="cornflowerblue",  offset = .013, fontsize = 3.5) +
  geom_cladelabel(node=AsAm , label="Asian-American", color="darkorange1",  offset = .057, fontsize = 3.5) +
  geom_cladelabel(node=Asian1 , label="Asian-I", color="forestgreen",  offset = .0565, fontsize = 3.5) +
  geom_cladelabel(node=Asian2 , label="Asian-II", color="darkgoldenrod1",  offset = .129, extend = c(1,1), fontsize = 3.5) +
  xlim(c(0,.16)) +geom_treescale( x=.004,y=74, linesize = .5, fontsize = 3.5) 

pB

#and together

FigS3 <- cowplot::plot_grid(pA, pB, ncol = 2, nrow = 1, labels = c("A", "B"),label_size = 18)


ggsave(file = paste0(homewd, "/final-figures/FigS3.png"),
       plot=FigS3,
       units="mm",  
       width=100, 
       height=60, 
       scale=2.8)#, 
