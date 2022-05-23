

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
pA_new<-pSEA+
  theme(plot.margin = unit(c(0,4,0,-2), "cm"))


#now add the trees
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv1-small/DENV1avg.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv2-out/DENV2avg.tree"))

tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)

tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)


head(tree1dat)
head(tree2dat)

#and load the metadata


dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)

dat$date <- as.Date(dat$date)


mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-10-15"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#2020-09-15"


pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  theme_tree2() + coord_cartesian(xlim=c(1930,2021), ylim=c(0,350)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))# + geom_tiplab() 



pC2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
  coord_cartesian(xlim=c(1930,2021),  ylim=c(0,300))+
  geom_tiplab(size=1.5, offset = -2)+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))




tree1merge <- merge(x=tree1dat, y=dat, by="tip_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2dat, y=dat, by="tip_name", all.x = T, sort = F)
head(tree1merge)
tail(tree1merge)

#and attach clusterID
#dat.clust.save <- dplyr::select(dat.clust.save, -(new_label))
tree1merge <- merge(x=tree1merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)




head(tree1merge)

unique(tree1merge$country)
unique(tree2merge$country)

#add in cluster ID where relevant




#colorz = scales::hue_pal()(length(unique(c(tree1merge$country,tree2merge$country)))) #10 unique countries
#names(colorz) = unique(c(tree1merge$country,tree2merge$country))
# 
tree1merge$new_label = sapply(strsplit(tree1merge$beast_name, "_"), function(x) x[[1]])
tree1merge$new_label <- paste0(tree1merge$new_label, " ", as.character(tree1merge$date))

tree1merge$new_seq = "no"
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$sex)] <- "yes"
tree1merge$new_seq <- as.factor(tree1merge$new_seq)

tree1merge$CambodiaSeq <- "no"
tree1merge$CambodiaSeq[tree1merge$country=="Cambodia"] <- "yes"



unique(tree2merge$country)
tree2merge$new_label = sapply(strsplit(tree2merge$beast_name, "_"), function(x) x[[1]])
tree2merge$new_label <- paste0(tree2merge$new_label, " ", as.character(tree2merge$date))

# tree2merge$new_label = NA
# tree2merge$new_label[tree2merge$country=="Cambodia" & !is.na(tree2merge$sex)] <- tree2merge$year[tree2merge$country=="Cambodia" & !is.na(tree2merge$sex)]
# tree2merge$new_label <- as.factor(tree2merge$new_label)
tree2merge$new_seq = "no"
tree2merge$new_seq[tree2merge$country=="Cambodia" & !is.na(tree2merge$sex)] <- "yes"
tree2merge$new_seq <- as.factor(tree2merge$new_seq)

tree2merge$CambodiaSeq <- "no"
tree2merge$CambodiaSeq[tree2merge$country=="Cambodia"] <- "yes"




#colorz=c("no"="gray", "yes"="cornflowerblue")
shapez = c("yes"=21, "no"=24)
#shapez = c("yes"=19, "no"=15)


#colz=c('no'="gray", 'yes'="black")






pA <- pB1 %<+% tree1merge +
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
# tree1.tiplabel


######################################################################
#################### add bar to subplotB tree ########################
######################################################################


 
gen1_B_new <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_109-0011_2019-07-25" ),which(tree1@phylo$tip.label == "GQ357692_2008-07-31"))

pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen") +geom_cladelabel(node=gen1_B_new, label="Genotype I", color="seagreen",offset =-27, angle=270, offset.text = 4,vjust=10)  +#coord_cartesian(xlim=c(1930,2021), ylim=c(0,230))+ 
  theme_tree2() +   geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


pB_new<-pB1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  geom_tippoint(aes(fill=country, shape=new_seq), 
                size = 2, show.legend = F, stroke=.1, color="black") +
  theme(legend.position = c(.1,.75), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_fill_manual(values=colorz) +
  scale_shape_manual(values=shapez)


######################################################################
#################### add bar to subplot C tree #######################
######################################################################

gen1_C_new <- MRCA(tree4, which(tree4@phylo$tip.label== "OL414741_100-0277_2019-07-23" ),which(tree4@phylo$tip.label == "KU509277_2010-07-31"))
gen2_C_new <- MRCA(tree4, which(tree4@phylo$tip.label== "OL414721_100-0067_2019-07-25" ),which(tree4@phylo$tip.label == "FM210227_2002-07-31"))

pC2 <- ggtree(tree4, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
  # coord_cartesian(xlim=c(1930,2021),  ylim=c(0,230))+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=18)) # +geom_tiplab()
pC_new<-pC2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  +
  geom_cladelabel(node=gen1_C_new, label="Cosmopolitan I", color="tomato",offset =-35.5, angle=270, offset.text = 4,vjust=10, align=T)+
  geom_cladelabel(node=gen2_C_new, label="Asian I", color="navy",offset =-35.5, angle=270, offset.text = 4,vjust=10, align=T) +# ++
  
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)



 
pBC_new<-cowplot::plot_grid(pB_new,pC_new,nrow=2,ncol=1,labels=c("b", "c"),label_size=30)+
  theme(plot.margin = unit(c(0,4,0,0), "cm"))

pALL_new <- cowplot::plot_grid(pA_new, pBC_new, nrow=2, ncol = 1, rel_widths = c(2,1), rel_heights = c(1,2),align = "w",labels=c("a"),label_size=30)

pALL_new <- pALL_new + theme(plot.background = element_rect(fill ="white"))+ theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                     axis.title.x=element_blank(),
                                                                                                     axis.title.y=element_blank())


 
ggsave(file = paste0(homewd, "/fig4.png"),
       plot= pALL_new,
       units="mm",  
       width=55, 
       height=88, 
       scale=4, 
       dpi=300)


