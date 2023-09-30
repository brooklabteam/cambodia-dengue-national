


rm(list=ls())


library(ggplot2)
library(ggtree)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(treeio)
library(ape)
library(sp)
library(dismo)
library(geosphere)
library(rgeos)
library(sf)




#make three timetrees and one map all together
#read in tree from beast
#first, make map

#and load the metadata
#homewd= "/home/rstudio"
homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)


#first, the colors
colorz = scales::hue_pal()(10)
#and the names
names(colorz) <- c("Singapore", "Brunei", "Vietnam","Malaysia", "Philippines", "Thailand", "Myanmar", "Cambodia", "Laos", "Indonesia")

#world map
world <- ne_countries(scale = "medium", returnclass = "sf")
#subset
SEA = subset(world, sovereignt=="Singapore"| sovereignt=="Brunei"| sovereignt=="Vietnam" | sovereignt==  "Cambodia"| sovereignt== "Philippines"| sovereignt== "Thailand"| sovereignt== "Myanmar"| sovereignt== "Malaysia"| sovereignt== "Laos"| sovereignt== "Indonesia")


#plot in color
# gene world map

pSEA <- ggplot(data = SEA) +
  geom_sf(aes(fill=sovereignt), show.legend = F, color="black", size=.2) +
  #labs( x = "Longitude", y = "Latitude") +
  scale_fill_manual(values = colorz) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines"),
        #plot.margin = unit(c(10,3,10,3),"lines"), 
        #panel.background  = element_rect(size=1, fill = NULL, color = "black"),
        panel.background  = element_rect(fill="white"),
        panel.border  = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank()) 
coord_sf(xlim=c(95,170), expand = T)
psub<-pSEA+
  theme(plot.margin = unit(c(0,0,0,-1.5), "cm"))

#now make into a grob, so it can get embedded as a subplot in the phylogeny



#now add the trees
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv1-out-final/DENV1avg.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv2-out-final/DENV2avg.tree"))

tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)

tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)


head(tree1dat)
head(tree2dat)

#and load the metadata


dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

#check the format
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")


mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-12-15"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"2020-09-23"

node.tree1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_2019-07-25" ),which(tree1@phylo$tip.label == "ON046271_2004-11-20"))


pA1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  geom_cladelab(node=node.tree1, label="Genotype I", textcolor="seagreen", barcolor="seagreen", fontsize=6,
                offset =-37, angle=270, offset.text = -12, vjust=2, hjust=.5)  +
  theme_tree2() + coord_cartesian(xlim=c(1930,2030), ylim=c(0,390)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  #scale_fill_continuous(low="yellow", high="red", limits=c(0,1))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  annotation_custom(ggplotGrob(psub), xmin = 1935, xmax = 1980,  ymin = 140, ymax = 385) #add map as embedded subplot



#node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-15" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
node.tree2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-23" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
#node.tree2.2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-15"),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))
node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-25" ),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))

pB2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + 
  coord_cartesian(xlim=c(1930,2030),  ylim=c(0,350))+
  geom_cladelab(node=node.tree2, label="Cosmopolitan I", textcolor="tomato",barcolor="tomato",
                offset =-37, angle=270, offset.text = -12, fontsize=6, vjust=2, hjust=.5)  +
  geom_cladelab(node=node.tree2.1, label="Asian I", textcolor="navy", barcolor="navy", fontsize=6,vjust=2, hjust=.5,
                offset =-37, angle=270, offset.text = -12)  +
  #geom_tiplab(size=1)+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red", limits=c(0,1)) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.15,.85), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))



tree1merge <- merge(x=tree1dat, y=dat, by="tip_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2dat, y=dat, by="tip_name", all.x = T, sort = F)
head(tree1merge)
tail(tree1merge)
head(tree2merge)



#and attach clusterID - for suppplot
#dat.clust.save <- dplyr::select(dat.clust.save, -(new_label))
#tree1merge <- merge(x=tree1merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)
#tree2merge <- merge(x=tree2merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)

tree1merge$new_label = sapply(strsplit(tree1merge$tip_name, "_"), function(x) x[[1]])
tree1merge$new_label <- paste0(tree1merge$new_label, " ", as.character(tree1merge$date))

tree1merge$new_seq = "no"
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$sex)] <- "yes"
tree1merge$new_seq <- as.factor(tree1merge$new_seq)

tree1merge$CambodiaSeq <- "no"
tree1merge$CambodiaSeq[tree1merge$country=="Cambodia"] <- "yes"



unique(tree2merge$country)
tree2merge$new_label = sapply(strsplit(tree2merge$tip_name, "_"), function(x) x[[1]])
tree2merge$new_label <- paste0(tree2merge$new_label, " ", as.character(tree2merge$date))

tree2merge$new_seq = "no"
tree2merge$new_seq[tree2merge$country=="Cambodia" & !is.na(tree2merge$sex)] <- "yes"
tree2merge$new_seq <- as.factor(tree2merge$new_seq)

tree2merge$CambodiaSeq <- "no"
tree2merge$CambodiaSeq[tree2merge$country=="Cambodia"] <- "yes"

shapez = c("yes"=21, "no"=24)


pA <- pA1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=country, shape=new_seq), 
                size = 2, show.legend = F, stroke=.1, color="black") +
  theme(legend.position = c(.08,.8), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_fill_manual(values=colorz) +
  scale_shape_manual(values=shapez) 




tree1.tiplabel<-tree1@phylo[["tip.label"]]




######################################################################
#################### add bar to subplot B tree #######################
######################################################################



pB <-pB2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  +
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18)) + 
  #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
  #legend.position = "none",
  #legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)


pAB <-cowplot::plot_grid(pA,pB,nrow=2,ncol=1,labels=c("A", "B"),label_size=22)+
  theme(plot.margin = unit(c(0,0,.5,0), "cm"))




######################################################################
#### now, combine with transmission trees plot (D) ###################
######################################################################

#load the transmission tree data
#here for both serotypes
all.denv <- read.csv(file=paste0(homewd, "/data/AllDENVtransTreeDat.csv"), header = T, stringsAsFactors = F)
all.denv$distance <- all.denv$distance/1000 #convert to km

#ggplot(dat) + geom_point(aes(x=evol_time,y=distance)) + facet_grid(~DENV.subtype) 

#are they in the same season?
all.denv$season <- "yes"
all.denv$season[all.denv$pairtime2-all.denv$pairtime1>.5] <- "no"

#now, only look at those within a season
all.denv <- subset(all.denv, season=="yes")
unique(all.denv$paired)


#dat=subset(all.denv,DENV.serotype=="DENV-1")
mrca_thresh=.5
#and write over to get the transmission trees
#geothresh <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
geothresh <- list(.2,.4,.6,.8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
#here's the threshold list of distances over which to compute

#here, decides whether each pair is in a transmission chain, based on the mrca_thresh (in years)
is.trans.chain <- function(dat, mrca_thresh){
  dat$trans_chain <- 0
  dat$trans_chain[dat$tMRCA <= mrca_thresh] <- 1
  return(dat)
  
}
get.prop.chain.discrete<- function(thresh, dat, thresh.bin){
  if(which(thresh.bin==thresh)==1){
    min.dist = 0
  }else{
    min.dist = thresh.bin[[(which(thresh.bin==thresh)-1)]]  
  }
  tot.chain <- sum(subset(dat, distance<=thresh & distance>min.dist)$trans_chain)
  N.pairs <- length(subset(dat, distance<=thresh& distance>min.dist)$trans_chain)
  if(N.pairs==0){
    prop.chain <-NA
    dat.new <- cbind.data.frame(distance = thresh, prop=NA, prop_lci = NA, prop_uci= NA, tot_pairs_in_chain=NA, tot_pairs=NA)
  }else{
    prop.chain <- tot.chain/N.pairs
    CI <- binom.test(x=tot.chain, n=N.pairs, alternative = "two.sided", conf.level = .95)
    dat.new <- cbind.data.frame(distance = thresh, prop=prop.chain, prop_lci = CI$conf.int[1], prop_uci= CI$conf.int[2], tot_pairs_in_chain=tot.chain, tot_pairs=N.pairs)
    
  }
  
  
  
  
  
  #and get CIs
  
  
  return(dat.new)
}
get.prop.chain.max <- function(thresh, dat){
  
  tot.chain <- sum(subset(dat, distance<=thresh)$trans_chain)
  N.pairs <- length(subset(dat, distance<=thresh)$trans_chain)
  prop.chain <- tot.chain/N.pairs
  
  CI <- binom.test(x=tot.chain, n=N.pairs, alternative = "two.sided", conf.level = .95)
  
  
  #and get CIs
  
  dat.new <- cbind.data.frame(distance = thresh, prop=prop.chain, prop_lci = CI$conf.int[1], prop_uci= CI$conf.int[2], tot_pairs_in_chain=tot.chain, tot_pairs=N.pairs)
  return(dat.new)
}
make.trans.chains <- function(dat, geothresh, mrca_thresh, character){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #and proportion of transmission chains
  dat.prop = lapply(geothresh, get.prop.chain.max, dat=dat1)
  #dat.prop = lapply(geothresh, get.prop.chain.discrete, dat=dat1, thresh.bin=geothresh)
  prop.df = data.table::rbindlist( dat.prop)
  prop.df$character = character
  prop.df$DENV.serotype=unique(dat$DENV.serotype)
  prop.df$transchain_threshold = mrca_thresh
  prop.df$character <- character
  return(prop.df)
  
}

combine.chain.prop<-function(mrca.thresh){
  out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Asian-1/DENV-2-Asian-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #head(out.prop1)
  out.prop1$DENV.subtype <- "DENV-1"
  #out.prop3$DENV.subtype <- "DENV-2-Cosmopolitan"
  out.prop7$DENV.subtype <- "DENV-2-All"
  #out.prop5$DENV.subtype <- "DENV-2-Asian-1"
  #out.all <-rbind(out.prop1, out.prop3, out.prop5, out.prop7)
  out.all <-rbind(out.prop1, out.prop7)
  return(out.all)
  
}
make.chain.diff <- function(mrca.thresh){
  out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Asian-1/DENV-2-Asian-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  head(out.prop1)
  
  out.prop <- out.prop3 #Cosmopolitan is dominant and we subtract from it
  out.prop$prop_diff <- out.prop3$prop - out.prop1$prop
  out.prop$prop_diff_lci <-  out.prop$prop_lci-out.prop1$prop
  out.prop$prop_diff_uci <-  out.prop$prop_uci-out.prop1$prop
  out.prop$comp <- "DENV-1"
  
  
  out.prop2 <- out.prop3 #Cosmopolitan is dominant and we subtract from it
  out.prop2$prop_diff <- out.prop2$prop - out.prop7$prop
  out.prop2$prop_diff_lci <-  out.prop2$prop_lci-out.prop7$prop
  out.prop2$prop_diff_uci <-  out.prop2$prop_uci-out.prop7$prop
  #out.prop2$comp <- "DENV-2-Asian-1"#
  out.prop2$comp <- "All-DENV-2"
  
  out.prop <- rbind(out.prop, out.prop2)
  
  out.prop <- dplyr::select(out.prop, distance, prop_diff, prop_diff_lci, prop_diff_uci, transchain_threshold, comp)
  return (out.prop)
}
mean.trans.chains <- function(dat, geothresh, mrca_thresh){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #get the probability that two cases are from the same chain
  prob <- sum(dat1$trans_chain)/length(dat1$trans_chain)
  prob.uci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[1]
  prob.lci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[2]
  
  eff.chains <- 1/prob #22.264 chains circulating in the region (Denv1)
  eff.chains.lci <- 1/prob.lci #17.019 chains circulating in the region (Denv1)
  eff.chains.uci <- 1/prob.uci #29.298 chains circulating in the region (Denv1)
  
  #and return
  chain.df <- cbind.data.frame(DENV.serotype=unique(dat1$DENV.serotype), DENV.subtype=unique(dat1$DENV.subtype), N_chains=eff.chains, N_chains_lci=eff.chains.lci, N_chains_uci = eff.chains.uci)
  
  return(chain.df)
  #pop size KP: 877,523... this is a little higher than the prediction for rural Thailand and lower
  #than the prediction for Bangkok
  #and get the mean number of transmission chains for this region that are in circulation
  
  
}


#out.pt5 = make.chain.diff(.5)
out.pt5 = combine.chain.prop(mrca.thresh=.5)
head(out.pt5)
tail(out.pt5)
#out.prop = rbind(out.pt3, out.pt5, out.1)#, out.3)
#head(out.prop)
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Cosmopolitan"))
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Asian-1", "DENV-2-Cosmopolitan"))
out.pt5$DENV.serotype<- factor(out.pt5$DENV.serotype, levels = c("DENV-1", "DENV-2"))


#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")
#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Asian-1"="cyan",  "DENV-2-Cosmopolitan"="dodgerblue")
colzB=c("DENV-1"="mediumseagreen", "DENV-2"="navy")


pC <- ggplot(data=out.pt5) + theme_bw()+
  #facet_grid(dummy_label~.) +
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=DENV.serotype, group=DENV.serotype), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=DENV.serotype, group=DENV.serotype)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        #axis.title.x = element_blank(), axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        legend.title = element_blank(),
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 18),
        axis.text = element_text(size=14), legend.text = element_text(size=12),
        legend.position = c(.73,.84)) + coord_cartesian(ylim=c(0,1), xlim=c(0,5), expand = F)+
  scale_color_manual(values=colzB) + 
  scale_fill_manual(values=colzB) + 
  #scale_color_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  #scale_fill_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  ylab("Proportion same transmission chain")  +
  xlab("Max distance between cases (km)")


######################################################################
################## and combine with data from Salje (C) ##############
######################################################################

#and plot with the salje data
salje.dat <- read.csv(file=paste0(homewd, "/data/salje_chains.csv"), header = T, stringsAsFactors = F)
head(salje.dat)


#and get mean chains
denv.1 = subset(all.denv, DENV.serotype=="DENV-1")
denv.1$DENV.subtype <- "DENV-1"
denv.1.mean = mean.trans.chains(dat=denv.1, geothresh, mrca_thresh)

denv.2 = subset(all.denv, DENV.serotype=="DENV-2")
denv.2$DENV.subtype <- "DENV-2"
denv.2.mean = mean.trans.chains(dat=denv.2, geothresh, mrca_thresh)


# denv.cosmo.2 = subset(all.denv, DENV.serotype=="DENV-2" & DENV.subtype=="DENV-2-Cosmopolitan")
# denv.cosmo.2$DENV.subtype <- "DENV-2-Cosmopolitan"
# denv.cosmo.2.mean = mean.trans.chains(dat=denv.cosmo.2, geothresh, mrca_thresh)

#all.denv.mean <- rbind(denv.1.mean, denv.2.mean, denv.cosmo.2.mean)
all.denv.mean <- rbind(denv.1.mean, denv.2.mean)
#fewer circulating chains for denv-2 vs. 1 and even fewer for co


#salje.dat$locale[salje.dat$locale=="Bangkok"] <- "Salje et al. Bangkok"
salje.dat$locale[salje.dat$locale=="Thai_Countryside"] <- "Rural Thailand"# "Salje et al.\nThai rural"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-1"] <- "Cambodia DENV-1"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2"] <- "Cambodia DENV-2"
#all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2-Cosmopolitan"] <- "Cambodia\nDENV-2-Cosmopolitan"

salje.dat$study <- "Salje et al. 2017"
all.denv.mean$study <- "Kampong Speu 2019-2020"

#colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia All-DENV-2' = "navy")
shapeznew <- c('Salje et al. 2017' = 21, "Kampong Speu 2019-2020" = 24)
colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia DENV-2' = "navy", 'Salje et al. 2017' = "black", "Kampong Speu 2019-2020" = "red")

#and plot
pD <- ggplot(data=salje.dat) + 
  geom_point(aes(x=pop_size, y=eff_chains, fill=locale, 
                 shape=study, color=locale), size=4, color="black", show.legend = F) +
  geom_errorbar(aes(x=pop_size, ymin=lci, ymax=uci, color=locale)) +
  scale_y_log10() + 
  scale_x_log10(breaks=c(1e+03, 1e+04, 1e+05, 1e+06, 1e+07), 
                labels=c(1,10,100,1000, 10000)) + theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        plot.margin = unit(c(.3,.2,.1,.5), "cm"),
        #legend.background = element_rect(color="black"),
        legend.position = c(.35,.88)) +
  ylab("Effective # Chains") +
  xlab("Population Size (x1000)") +
  geom_point(data=all.denv.mean,aes(x=877523., y=N_chains, fill=DENV.subtype, shape=study, color=study), size=5,  stroke=2) +
  scale_color_manual(values=colznew, guide=NULL) + scale_fill_manual(values=colznew, guide=NULL) +
  scale_shape_manual(values = shapeznew) 


#put C and D together


pCD <-cowplot::plot_grid(pC,pD,nrow=2,ncol=1,labels=c("C", "D"),label_size=22)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


#and all together

Fig4 <- cowplot::plot_grid(pAB, pCD, nrow=1, ncol = 2) + theme(plot.background = element_rect(fill ="white"))+ 
  theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                          axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())

ggsave(file = paste0(homewd, "/final-figures/fig4.png"),
       plot= Fig4,
       units="mm",  
       width=100, 
       height=85, 
       scale=3, 
       dpi=300)




