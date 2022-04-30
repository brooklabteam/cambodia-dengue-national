



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
homewd= "/home/rstudio"
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
cam = sf::st_read(paste0(homewd, "/data/kampongspeu/provinces.shp"))
sub = subset(cam, name=="Kampong Speu")


dat <- read.csv(file = paste0(homewd, "/data/BEAST-seq-metadata-11-6.csv"), header = T, stringsAsFactors = F)

#dat = subset(dat, DENV.serotype=="DENV-2")
#names(dat)
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
#dat.plot$date <- as.Date(dat.plot$date)

head(dat.plot)
dat.plot$date <- as.character(as.Date(dat.plot$date, format = "%m/%d/%y"))
dat.plot$new_label = paste(dat.plot$seq.ID, dat.plot$date, sep = " ")
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
all.denv <- get.spatial.object(dat1=subset(dat.plot, !is.na(lat) & year >2018), dist.thresh = 12000, denv.serotype = "all")
denv.map <- all.denv[[1]]
dat.plot.cluster <- all.denv[[2]]
dat.plot.cluster <- dplyr::select(dat.plot.cluster,new_label, cluster_ID)

dat.plot <- merge(dat.plot, dat.plot.cluster, by="new_label", all.x=TRUE)


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





#and add in the points


# 
# pA <- pSEA + annotation_custom(grob = ggplotGrob(pKPS), 
#                                       xmin = 118, xmax = 145,
#                                       ymin = 12,
#                                       ymax = 32) +
#       geom_segment(aes(x=104.8, xend=121.5, y=12.2, yend=32), linetype=2)+
#       geom_segment(aes(x=104.8, xend=120, y=11, yend=12.2), linetype=2) +
#       theme(plot.background = element_rect(fill=NULL, color = "black", size=3))
#   
# 


pE <- pSEA + annotation_custom(grob = ggplotGrob(pKPS), 
                               xmin = 128, xmax = 175,
                               ymin = 5,
                               ymax = 43) +
  geom_segment(aes(x=104.8, xend=131, y=12.2, yend=42), linetype=2, size=1)+
  geom_segment(aes(x=104.8, xend=131, y=11, yend=6), linetype=2, size = 1) #+
#theme(plot.background = element_rect(fill=NULL, color = "black", size=3))


#and the other half

all.denv <- read.csv(file=paste0(homewd,"/data/AllDENVtransTreeDat.csv"), header = T, stringsAsFactors = F)
all.denv$distance <- all.denv$distance/1000 #convert to km


#are they in the same season?
all.denv$season <- "yes"
all.denv$season[all.denv$pairtime2-all.denv$pairtime1>.5] <- "no"


unique(all.denv$evol_time[all.denv$season=="no"]) #all NA
unique(all.denv$evol_time[all.denv$season=="yes"]) 
head(all.denv)
unique(all.denv$paired)

#all.denv = subset(all.denv, paired=="DENV-1/DENV-1" | paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"| paired=="DENV-2-Asian-1/DENV-2-Asian-1" | paired=="DENV-1-clade-2/DENV-1-clade-2"  )
all.denv = subset(all.denv, !is.na(evol_time))

colz= c("DENV-1" = "forestgreen", "DENV-2-Cosmopolitan" = "navy")
all.denv$evol_time <- round_any(all.denv$evol_time, .2)
#evo.sum <- ddply(subset(all.denv, distance<20), .(evol_time), summarise, spatial_med = median(distance), spatial_25 = quantile(distance)[2], spatial_75=quantile(distance)[4])

evo.sum <- ddply(all.denv, .(evol_time), summarise, spatial_med = median(distance), spatial_25 = quantile(distance)[2], spatial_75=quantile(distance)[4])

colz= c("DENV-1" = "forestgreen", "DENV-2" = "navy")
median.all.sub <- cbind.data.frame(med_dist=c(median(subset(all.denv, distance<20)$distance)))


pF <- ggplot(data=subset(evo.sum, spatial_med<20 & evol_time<10)) +
  geom_point(aes(x=evol_time,y=spatial_med), show.legend = F, size=5) +
  #geom_line(aes(x=evol_time,y=spatial_med), show.legend = F) +
  #geom_label(aes(x=2.3, y=48, label = DENV.serotype, size=14), label.size = NA, show.legend = F) +
  #scale_fill_manual(values = colz)+
  #scale_color_manual(values = colz)+
  geom_ribbon(aes(x=evol_time,ymin=spatial_25, ymax=spatial_75), alpha=.3, show.legend = F) + 
  geom_hline(data=median.all.sub, aes(yintercept = med_dist), linetype=2, size=1) +
  #coord_cartesian(ylim=c(0,15), xlim =c(0,10)) +
  #facet_grid(~DENV.serotype) + 
  ylab("Spatial distance (km)") +
  xlab("Evolutionary time (yrs)") + theme_bw() +
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(1,1,1,1),"lines"),
        strip.background = element_blank(),# element_rect(fill="white", color="white"),
        axis.title = element_text(size=14), 
        axis.text = element_text(size = 12),
        strip.text = element_blank())+
  scale_x_continuous(sec.axis=sec_axis(trans=~. *365/20, name="Transmission generations"))



pEF <- cowplot::plot_grid(pE, pF, ncol = 1, nrow = 2, rel_widths = c(1,1), labels = c("E", "F"), label_size = 22)

#pEF
#now we add in the big trees
dat.clust.save <- dplyr::select(dat.plot, beast_name, new_label, cluster_ID)

dat <- read.csv(file = paste0(homewd, "/data/BEAST-seq-metadata-11-6.csv"), header = T, stringsAsFactors = F)

#dat = subset(dat, DENV.serotype=="DENV-2")
names(dat)
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
dat.plot$date <- as.Date(dat.plot$date)
dat.plot$accession_num <- c(unlist(sapply(strsplit(dat.plot$beast_name, "_"), function(x) x[[1]])))
head(dat.plot)
dat.plot$new_label = paste(dat.plot$accession_num, dat.plot$date, sep = " ")


#dat.plot  <- dplyr::select(dat.plot, new_label, unique.ID, NIH.ID, CZB.ID, beast_name, date, year, lat, long, accession_num, DENV.serotype)
head(dat.plot)

#and the corresponding trees
#library(rBt)
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv1-strict/DENV1strictAVG.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv2-strict/DENV2strictAVG.tree"))

tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)
tree1dat$beast_name <-tree1dat$tip_name
tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)
tree2dat$beast_name <-tree2dat$tip_name

head(tree1dat)
head(tree2dat)
#dat$date <- as.Date(dat$date)

#dat <- read.csv(file = "/Users/caraebrook/Documents/R/R_repositories/cambodia-dengue/aug-2021/BEAST-seq-meta.csv", header = T, stringsAsFactors = F)


head(dat)
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")
mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2020-07-13"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"#"2020-09-23"



pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  theme_tree2() + coord_cartesian(xlim=c(1930,2021), ylim=c(0,230)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

# 
# #and collapse the clades with the sylvatic strains
# node1 <- MRCA(tree2, "FJ467493_2008-07-31", "OL414763_109-0039_2019-08-02")
# node2 <- MRCA(tree2, "KY923048_2015-07-31", "OL414763_109-0039_2019-08-02")
# node2 <- MRCA(tree2, "KY923048_2015-07-31", "FJ467493_2008-07-31")
# node2 <- nodeid(tree2, "KY923048_2015-07-31")
# node1<-  nodeid(tree2, "FJ467493_2008-07-31")

#p2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
# coord_cartesian( ylim=c(0,180)) #+geom_tiplab(size=2) #two very disparate lineages of denv2xlim=c(1900,2021),

tree3 <- tree_subset(tree=tree2, node = 177, levels_back=0)
tree4 <- tree_subset(tree=tree3, node = 176, levels_back=0)

pC2 <- ggtree(tree4, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
  coord_cartesian(xlim=c(1930,2021),  ylim=c(0,230))+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red") +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


tree1merge <- merge(x=tree1dat, y=dat, by="beast_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2dat, y=dat, by="beast_name", all.x = T, sort = F)


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


#
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




pB <- pC2 %<+% tree2merge + 
  ggnewscale::new_scale_fill() +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)# +
#scale_color_manual(values=colorz)#+


#and place above and below

pAB <- cowplot::plot_grid(pA, pB,  nrow=2, ncol=1, labels = c("a", "b"), label_size = 20)






#and the sublineages in the next panels
# #and collapse the clades with the sylvatic strains
node1_sub <- MRCA(tree1, "OK159963_100-0734_2019-07-21", "OL412678_109-0011_2019-07-25")
node2_sub <- MRCA(tree2, "OL414741_100-0277_2019-07-23", "OL414749_109-0293_2020-06-17")

#p2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
# coord_cartesian( ylim=c(0,180)) #+geom_tiplab(size=2) #two very disparate lineages of denv2xlim=c(1900,2021),

tree1sub <- tree_subset(tree=tree1, node =  node1_sub, levels_back=0)
tree2sub <- tree_subset(tree=tree2, node =  node2_sub, levels_back=0)



pD2 <- ggtree(tree1sub, mrsd=mrsd.denv1, color="forestgreen")  + 
  theme_tree2() + coord_cartesian(xlim=c(2016.3,2021.3), ylim=c(0,49), expand = T) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  #scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  scale_fill_continuous(low="yellow", high="red")+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))

pE2 <- ggtree(tree2sub, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + # nodelab() +
  coord_cartesian(xlim=c(2016.3,2021.3), ylim=c(0,49), expand = T) + 
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_fill_continuous(low="yellow", high="red") +
  #scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


pC <- pD2 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  geom_label(aes(label=new_label), label.size = NA, size=1.5, hjust=-.1, fill=NA) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=cluster_ID, shape=new_seq), 
                size = 2, show.legend = F, stroke=1, color="#9590FF") +
  theme(legend.position = c(.1,.75), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        plot.margin = unit(c(.2,5,.1,.1),"lines"), 
        # plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_x_continuous(breaks=c(2017, 2018, 2019, 2020)) +
  scale_fill_manual(values=colorz1) +
  scale_shape_manual(values=shapez) #+

#and the second half



pD <- pE2 %<+% tree2merge + 
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



pCD <- cowplot::plot_grid(pC,  pD, labels=c("b", "c"), nrow=2, ncol=1, label_size = 20)


#and stack

pALL <- cowplot::plot_grid(pAB, pCD, pEF, nrow=1, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1,1))

pALL <- pALL + theme(plot.background = element_rect(fill ="white"))




pS3_KPS<- cowplot::plot_grid(pKPS,  nrow=1, ncol=1, labels = c("a"), label_size = 20)


pALL <- cowplot::plot_grid(pAB, pCD, pEF, nrow=1, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1,1))

pALL <- pALL + theme(plot.background = element_rect(fill ="white"))


pALL <- cowplot::plot_grid(pS3_KPS, pCD, nrow=2, ncol = 1, rel_widths = c(1), rel_heights = c(1,2))





pALL <- pALL + theme(plot.background = element_rect(fill ="white"))

pALL<-pALL+ theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.title.y=element_blank())



ggsave(file = paste0(homewd, "/figS3.png"),
       plot= pALL,
       units="mm",  
       width=55, 
       height=88, 
       scale=3, 
       dpi=300)




