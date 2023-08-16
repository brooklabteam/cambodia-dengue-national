rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)
library(ggh4x)

#this one is just the plot of coherence with ONI (or lack thereof)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

oni.dat <- read.csv(file= paste0(homewd, "/data/oni_coherence.csv"), header = T, stringsAsFactors = F)

#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

setdiff(unique(oni.dat$provname), unique(centroid.prov$provname))
setdiff(unique(centroid.prov$provname), unique(oni.dat$provname))

oni.dat <- merge(oni.dat, centroid.prov, by="provname")

#and plot by latitude
oni.dat <- arrange(oni.dat, latitude, time)

oni.dat$provname <- factor(oni.dat$provname, levels=unique(oni.dat$provname))


#now make color bar
colz = scales::hue_pal()(length(unique((oni.dat$provname)))) #25
name.list <- sort(unique(as.character(oni.dat$provname))) #alphabetical
name.list <- c(name.list[name.list!="Ratanak Kiri" & name.list!="Mondul Kiri" & name.list!= "Tboung Khmum"], "Mondul Kiri", "Ratanak Kiri", "Tboung Khmum") #these three get moved to the end
names(colz) <- name.list

#now, reorder it in the same order that the provinces are plotted.
centroid.prov$color= NA
for(i in 1:length(colz)){
  centroid.prov$color[centroid.prov$provname==names(colz)[i]] <- colz[i]
}
centroid.prov <- arrange(centroid.prov, latitude)
colz = centroid.prov$color
names(colz) <- centroid.prov$provname

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

#and coherence
vert.df <- cbind.data.frame(xint = c(2007, 2008, 2012, 2013, 2019,2020))

#ggplot(subset(oni.dat, provname=="Battambang")) + geom_line(aes(x=time, y=avg_wave_coherence_oni))

#cross wavelet power
FigS13Bb <- ggplot(data=oni.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_oni, color=avg_cross_power_oni)) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt") +# limits=c(0.5,1)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt")+#, limits=c(0.5,1)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#coherence
FigS13Ab <- ggplot(data=oni.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_wave_coherence_oni, color=avg_wave_coherence_oni)) + 
  scale_fill_viridis_c( option="inferno", name="Average wave coherence\n (multiannual): dengue IR and ONI", trans="sqrt")+#, limits=c(0.5,1)) +
  scale_color_viridis_c( option="inferno", name="Average wave coherence\n (multiannual): dengue IR and ONI", trans="sqrt")+#, limits=c(0.5,1)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#and the distribution by month
dist.dat <- ddply(oni.dat, .(time, year, month), summarise, median_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["50%"],  min_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["25%"], max_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["75%"], median_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["50%"],  min_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["25%"], max_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["75%"])
head(dist.dat)

#and plot
FigS13Ba <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_oni_cross, ymax=max_oni_cross), alpha=.3) + ylab("distribution monthly\nmulti-annual dengue cycle cross\nwavelet power with ONI") +
  geom_line(aes(x=time, y=median_oni_cross), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS13Aa <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), ylim=c(.7,.9), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_oni_coherence, ymax=max_oni_coherence), alpha=.3) + ylab("distribution monthly\nmulti-annual dengue cycle\ncoherence with ONI") +
  geom_line(aes(x=time, y=median_oni_coherence), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)



FigS13A <- cowplot::plot_grid(FigS13Aa, FigS13Ab, ncol=1, nrow=2, rel_heights = c(.25,1))
FigS13B <- cowplot::plot_grid(FigS13Ba, FigS13Bb, ncol=1, nrow=2, rel_heights = c(.25,1))

FigS13 <- cowplot::plot_grid(FigS13A, FigS13B, ncol=2, nrow=1, labels = c("A", "B"), label_size = 22)


ggsave(file = paste0(homewd, "/final-figures/FigS13.png"),
       plot= FigS13,
       units="mm",  
       width=110, 
       height=70, 
       scale=3, 
       dpi=300)
