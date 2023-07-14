


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
pA<-pSEA+
  theme(plot.margin = unit(c(0,0,0,-1.5), "cm"))


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
dat$date <- as.Date(dat$date)#, format = "%m/%d/%y")
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")


mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-12-15"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"2020-09-23"

node.tree1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_2019-07-25" ),which(tree1@phylo$tip.label == "ON046271_2004-11-20"))


pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  geom_cladelab(node=node.tree1, label="Genotype I", textcolor="seagreen", barcolor="seagreen", fontsize=6,
              offset =-37, angle=270, offset.text = -12, vjust=2, hjust=.5)  +
  theme_tree2() + coord_cartesian(xlim=c(1930,2030), ylim=c(0,390)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red", limits=c(0,1))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))# + geom_tiplab() 



#node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-15" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
node.tree2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-23" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
#node.tree2.2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-15"),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))
node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-25" ),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))

pC2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + 
  coord_cartesian(xlim=c(1930,2030),  ylim=c(0,350))+
  geom_cladelab(node=node.tree2, label="Cosmopolitan I", textcolor="tomato",barcolor="tomato",
                offset =-37, angle=270, offset.text = -12, fontsize=6, vjust=2, hjust=.5)  +
  geom_cladelab(node=node.tree2.1, label="Asian I", textcolor="navy", barcolor="navy", fontsize=6,vjust=2, hjust=.5,
                  offset =-37, angle=270, offset.text = -12)  +
  #geom_tiplab(size=1)+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red", limits=c(0,1)) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))




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


pB <- pB1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=country, shape=new_seq), 
                size = 2, show.legend = F, stroke=.1, color="black") +
  theme(legend.position = c(.1,.75), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_fill_manual(values=colorz) +
  scale_shape_manual(values=shapez) #+
#scale_color_manual(values=colorz)#+



tree1.tiplabel<-tree1@phylo[["tip.label"]]




######################################################################
#################### add bar to subplot C tree #######################
######################################################################



pC <-pC2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  +
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)


pBC <-cowplot::plot_grid(pB,pC,nrow=2,ncol=1,labels=c("B", "C"),label_size=30)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

pALL <- cowplot::plot_grid(pA, pBC, nrow=2, ncol = 1, rel_widths = c(2,1), rel_heights = c(1,2.2),align = "w",labels=c("A"),label_size=30)

pALL_new <- pALL + theme(plot.background = element_rect(fill ="white"))+ theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                     axis.title.x=element_blank(),
                                                                                                     axis.title.y=element_blank())


 
ggsave(file = paste0(homewd, "/final-figures/fig4.png"),
       plot= pALL_new,
       units="mm",  
       width=55, 
       height=88, 
       scale=3.5, 
       dpi=300)




