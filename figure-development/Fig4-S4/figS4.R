
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
library(lubridate)




#make three timetrees and one map all together
#read in tree from beast
#first, make map

#and load the metadata
homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
#homewd= "/home/rstudio"
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
        axis.title = element_blank()) +
  geom_rect(xmin = 103.8, xmax = 104.8, ymin = 11, ymax = 12.2, 
            fill = NA, colour = "black", size = 1)  +
  coord_sf(xlim=c(95,170), expand = T)




#and get submap of kampong speu witht the points of cambodia sequences
cam = sf::st_read(paste0(homewd, "/data/shapefile/provinces.shp"))
sub = subset(cam, name=="Kampong Speu")


dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)

dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$year <- year(dat$date)
#dat = subset(dat, DENV.serotype=="DENV-2")
#names(dat)
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
#dat.plot$date <- as.Date(dat.plot$date)

head(dat.plot)
dat.plot$date <- as.character(as.Date(dat.plot$date))
#dat.plot$date <- as.character(as.Date(dat.plot$date, format = "%m/%d/%y"))
dat.plot$new_label = paste(dat.plot$NIH.ID, dat.plot$date, sep = " ")
dat.plot$new_label


#dat.plot  <- dplyr::select(dat.plot, new_label, unique.ID, NIH.ID, CZB.ID, beast_name, date, year, lat, long, accession_num, DENV.serotype)
head(dat.plot)

get.spatial.object <- function(dat1, dist.thresh, denv.serotype){
  
  #add in some local points and cluster them 
  xy <- SpatialPointsDataFrame(
    matrix(c(dat1$long, dat1$lat), ncol=2), data.frame(ID=seq(1:length(dat1$lat))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  # use the distm function to generate a geodesic distance matrix in meters
  mdist <- distm(xy)
  # cluster all points using a hierarchical clustering approach
  hc <- hclust(as.dist(mdist), method="complete")
  # define the distance threshold in meters. 
  #d=2000 #this is the typical
  d=dist.thresh #this is within 10km
  #d=500
  
  
  # define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
  xy$clust <- cutree(hc, h=d)
  
  # expand the extent of plotting frame
  xy@bbox[] <- as.matrix(extend(extent(xy),0.001))
  
  # get the centroid coords for each cluster
  cent <- matrix(ncol=2, nrow=max(xy$clust))
  for (i in 1:max(xy$clust)){
    # gCentroid from the rgeos package
    cent[i,] <- gCentroid(subset(xy, clust == i))@coords
  }
  
  
  # compute circles around the centroid coords using a radius of d meters
  # from the dismo package
  ci <- circles(cent, d=d, lonlat=T)
  
  dat1$cluster_ID <- NA
  dat1$cluster_ID[!is.na(dat1$lat)] <- xy$clust
  dat1$cluster_ID <- as.factor(dat1$cluster_ID)
  head(dat1)
  
  slim.gis <- dplyr::select(dat1, lat, long, cluster_ID, DENV.serotype)
  slim.gis  <- slim.gis[complete.cases(slim.gis),]
  #coord.dat <- st_as_sf(slim.gis, coords = c( "long", "lat"), crs = 4326)
  #circle.dat <- st_as_sf(ci@polygons, coords=c("X1", "X2"), crs=4326)
  
  #or the centroid dat
  if(denv.serotype=="all"){
    cent.gis <- ddply(slim.gis, .(cluster_ID), summarize, N=length(lat), DENV.serotype="all")
  }else if(denv.serotype!="all"){
    cent.gis <- ddply(slim.gis, .(cluster_ID), summarize, N=length(lat), DENV.serotype=unique( denv.serotype))
  }
  
  cent.gis$long <- cent[,1]
  cent.gis$lat <- cent[,2]
  circlpts.dat <- st_as_sf(cent.gis, coords = c( "long", "lat"), crs = 4326)
  
  return(list(circlpts.dat, dat1)) 
}

#now select those used for geospatial analysis
dat.plot = subset(dat.plot, !is.na(lat) & year >2018)#118 sequences
length(dat.plot$DENV.serotype[dat.plot$DENV.serotype=="DENV-1"]) #57
length(dat.plot$DENV.serotype[dat.plot$DENV.serotype=="DENV-2"]) #61

all.denv <- get.spatial.object(dat1=dat.plot, dist.thresh = 5000, denv.serotype = "all")
denv.map <- all.denv[[1]]
dat.plot.cluster <- all.denv[[2]]
dat.plot.cluster <- dplyr::select(dat.plot.cluster,new_label, cluster_ID)

#here, add cluster ID to the dataset
dat.plot <- merge(dat.plot, dat.plot.cluster, by="new_label", all.x=TRUE)

head(dat.plot)

# create your own color palette 
colorz1 = sample(rainbow(length(unique(dat.plot$cluster_ID))), size=length(unique(dat.plot$cluster_ID)), replace=F)
names(colorz1) <- sort(unique(dat.plot$cluster_ID))


# get map
pKPS <- ggplot(sub) + geom_sf(fill="#9590FF", color="black", size =.4) + 
  # theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  geom_sf(data=denv.map, aes(fill=cluster_ID), color="gray90",
          shape=21, size=2.5, stroke=1, show.legend = F) +
  scale_size_continuous(range=c(3,20)) +
  guides(fill="none") + scale_fill_manual(values=colorz1)




#now we add in the trees
dat.clust.save <- dplyr::select(dat.plot, tip_name, new_label, cluster_ID)

dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)

#dat = subset(dat, DENV.serotype=="DENV-2")
names(dat)
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
dat.plot$date <- as.Date(dat.plot$date)
dat.plot$accession_num <- c(unlist(sapply(strsplit(dat.plot$tip_name, "_"), function(x) x[[1]])))
head(dat.plot)
dat.plot$new_label = paste(dat.plot$accession_num, dat.plot$date, sep = " ")


#dat.plot  <- dplyr::select(dat.plot, new_label, unique.ID, NIH.ID, CZB.ID, beast_name, date, year, lat, long, accession_num, DENV.serotype)
head(dat.plot)

#now build the same trees as before, except select out the lineages of interest and color based on cluster instead
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


mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-07-13"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"2020-09-23"


node.tree1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_2019-07-25" ),which(tree1@phylo$tip.label == "ON046271_2004-11-20"))

#this is the large tree in panel A, Figure 4
pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  geom_cladelab(node=node.tree1, label="Genotype I", textcolor="seagreen", barcolor="seagreen", fontsize=6,
                offset =-37, angle=270, offset.text = -12, vjust=2, hjust=.5)  +
  theme_tree2() + coord_cartesian(xlim=c(1930,2030), ylim=c(0,350)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))# + geom_tiplab() 


#and here is panel B, Figure 4

#node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-15" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
node.tree2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-23" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
#node.tree2.2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-15"),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))
node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-25" ),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))

pC2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + 
  coord_cartesian(xlim=c(1930,2030),  ylim=c(0,300))+
  geom_cladelab(node=node.tree2, label="Cosmopolitan I", textcolor="tomato",barcolor="tomato",
                offset =-37, angle=270, offset.text = -12, fontsize=6, vjust=2, hjust=.5)  +
  geom_cladelab(node=node.tree2.1, label="Asian I", textcolor="navy", barcolor="navy", fontsize=6,vjust=2, hjust=.5,
                offset =-37, angle=270, offset.text = -12)  +
  #geom_tiplab(size=1)+
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
head(tree2merge)


tree1merge$new_label = sapply(strsplit(tree1merge$tip_name, "_"), function(x) x[[1]])
tree1merge$new_label <- paste0(tree1merge$new_label, " ", as.character(tree1merge$date))

tree1merge$new_seq = "no"
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$sex)] <- "yes"
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$lat)] <- "yes"
tree1merge$new_seq <- as.factor(tree1merge$new_seq)

tree1merge$CambodiaSeq <- "no"
tree1merge$CambodiaSeq[tree1merge$country=="Cambodia"] <- "yes"


unique(tree2merge$country)
tree2merge$new_label = sapply(strsplit(tree2merge$tip_name, "_"), function(x) x[[1]])
tree2merge$new_label <- paste0(tree2merge$new_label, " ", as.character(tree2merge$date))

tree2merge$new_seq = "no"
tree2merge$new_seq[tree2merge$country=="Cambodia" & !is.na(tree2merge$sex)] <- "yes"
tree2merge$new_seq[tree2merge$country=="Cambodia" & !is.na(tree2merge$lat)] <- "yes"
tree2merge$new_seq <- as.factor(tree2merge$new_seq)

tree2merge$CambodiaSeq <- "no"
tree2merge$CambodiaSeq[tree2merge$country=="Cambodia"] <- "yes"

shapez = c("yes"=21, "no"=24)

#Fig 4B
pB <- pB1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  geom_tiplab(size=1) +
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

#and mrca for cambodia with previously reported cambodian sequences
node1mrca <-MRCA(tree1, which(tree1@phylo$tip.label== "OL412140_2019-08-28" ),which(tree1@phylo$tip.label == "MW265679_2015-10-06"))

node.dat.denv1 <- pB$data
subset(node.dat.denv1, node==node1mrca) #node height is 7.86 years
mrsd.denv1 - (7.86*365) # "2013-02-06"
node.dat.denv1$height_0.95_HPD[[node1mrca]] #7.290718 8.450780

mrsd.denv1 - (7.290718*365) # "2013-09-01"
mrsd.denv1 - (8.450780*365) # "2012-07-05"


tree1.tiplabel<-tree1@phylo[["tip.label"]]


#Fig 4C
pC <-pC2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  +
  #geom_tiplab(size=1) +
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)


####################################################################
####################################################################
####################################################################

#now we subset the tree to just the recent Cambodia clades for analysis
#node2_sub <-MRCA(tree2, which(tree2@phylo$tip.label== "OL414746_2019-06-15" ),which(tree2@phylo$tip.label == "OL414724_2019-08-15"))
node2_sub <-MRCA(tree2, which(tree2@phylo$tip.label== "OL414739_2019-07-23" ),which(tree2@phylo$tip.label == "OL414736_2020-06-17"))
#node1_sub <-MRCA(tree1, which(tree1@phylo$tip.label== "OK159963_2019-07-15" ),which(tree1@phylo$tip.label == "OL412140_2019-08-15"))
node1_sub <-MRCA(tree1, which(tree1@phylo$tip.label== "OL412140_2019-08-28" ),which(tree1@phylo$tip.label == "OK159963_2019-07-21"))

#get
tree1sub <- tree_subset(tree=tree1, node =  node1_sub, levels_back=0)
tree2sub <- tree_subset(tree=tree2, node =  node2_sub, levels_back=0)


#and merge the data with the cluster ID
dat.clust.save <- dplyr::select(dat.clust.save, -(new_label))
tree1merge <- merge(tree1merge, dat.clust.save, by = "tip_name", all.x = T)
tree2merge <- merge(tree2merge, dat.clust.save, by = "tip_name", all.x = T)

#Fig S3 B.1
pS3B1 <- ggtree(tree1sub, mrsd=mrsd.denv1, color="forestgreen")  + 
  theme_tree2() + coord_cartesian(xlim=c(2016.5,2021.4), ylim=c(0,73), expand = T) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  #scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

#Fig S3C.1
pS3C1 <- ggtree(tree2sub, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
  coord_cartesian(xlim=c(2016.5,2021.4), ylim=c(0,55), expand = T) + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_fill_continuous(low="yellow", high="red") +
  #scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


pS3B <- pS3B1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  geom_label(aes(label=new_label), label.size = NA, size=1.5, hjust=-.1, fill=NA) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=cluster_ID, shape=new_seq), 
                size = 2, show.legend = F, stroke=1, color="#9590FF") +
  theme(legend.position = c(.1,.75), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        #plot.margin = unit(c(.2,5,.1,.1),"lines"), 
        # plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_x_continuous(breaks=c(2017, 2018, 2019, 2020)) +
  scale_fill_manual(values=colorz1) +
  scale_shape_manual(values=shapez) #+

#and the second half


pS3C <- pS3C1 %<+% tree2merge + 
  ggnewscale::new_scale_fill() +
  geom_label(aes(label=new_label), label.size = NA, size=1.5, hjust=-.1, fill=NA) +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=cluster_ID, shape=new_seq), size = 2, 
                show.legend = F, stroke=1, color="#9590FF") +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        #plot.margin = unit(c(.1,5,.1,.1),"lines"), 
        legend.title = element_blank()) + 
  scale_x_continuous(breaks=c( 2017, 2018, 2019, 2020)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz1)#



#and the cosmo 2 sequences

#and mrca for cambodia with previously reported cambodian sequences
node2cosmomrca <-MRCA(tree2, which(tree2@phylo$tip.label== "OL414736_2020-06-17" ),which(tree2@phylo$tip.label == "FJ639707_2004-07-31"))

node.dat.denv2 <- pC$data
subset(node.dat.denv2, node==node2cosmomrca) #node height is 88.5 years
mrsd.denv2 - (88.5 *365) # "1932-04-15"
node.dat.denv2$height_0.95_HPD[[node2cosmomrca]] #83.60570 93.28363

mrsd.denv2 - (83.6057*365) # "1937-03-06"
mrsd.denv2 - (93.28363*365) # "1927-07-05"



pBC <- cowplot::plot_grid(pS3B,  pS3C, labels=c("B", "C"), nrow=2, ncol=1, label_size = 22)



pS3_KPS<- cowplot::plot_grid(pKPS,  nrow=1, ncol=1, labels = c("A"), label_size = 22)


pALL <- cowplot::plot_grid(pS3_KPS, pBC, nrow=2, ncol = 1, rel_widths = c(1), rel_heights = c(1,2))



pALL <- pALL + theme(plot.background = element_rect(fill ="white"))

pALL<-pALL+ theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.title.y=element_blank())


ggsave(file = paste0(homewd, "/final-figures/figSuppPhyloTrees.png"),
       plot= pALL,
       units="mm",  
       width=55, 
       height=88, 
       scale=3, 
       dpi=300)




