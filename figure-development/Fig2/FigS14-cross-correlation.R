rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(ggh4x)


homewd= "/Users/carabrook/Developer/cambodia-dengue-national"

#load data
pearsons.df <- read.csv(file = paste0(homewd, "/data/pearsons_correlations_provinces.csv"), header = T, stringsAsFactors = F)
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

#load and prep the climate data at the province level
temp.dat <- read.csv(file = paste0(homewd, "/data/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(homewd, "/data/biweek_ppt.csv"), header = T, stringsAsFactors = F )
head(temp.dat) # year and biweek
head(precip.dat) # year and biweek


temp.dat = subset(temp.dat, provname!="Administrative unit not available")
precip.dat = subset(precip.dat, provname!="Administrative unit not available")


#setdiff(unique(temp.dat$provname), unique(beta.df$provname)) #"Mondul Kiri"    "Oddar Meanchey" "Ratanak Kiri"   "Siemreap"       "Tboung Khmum"  
temp.dat$provname[temp.dat$provname=="Siemreap"] <- "Siem Reap"
precip.dat$provname[precip.dat$provname=="Siemreap"] <- "Siem Reap"
temp.dat$provname[temp.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
precip.dat$provname[precip.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"




pearsons.sum <- ddply(pearsons.df, .(provname, comp_prov), summarise, full_ts_corr=unique(full_ts_corr), dist_m = unique(dist_m), dist_km = unique(dist_km))
head(pearsons.sum)


colz = scales::hue_pal()(length(unique((tsir.dat$provname)))) #25
name.list <- sort(unique(as.character(tsir.dat$provname))) #alphabetical
name.list[name.list =="Otdar Meanchey"] <- "Oddar Meanchey"
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

p1 <- ggplot(pearsons.sum) + theme_bw() + ylab("full time series pearson's correlation between provinces") +
  xlab("distance between province centroids (km)") + coord_cartesian(ylim=c(0,1))+
  geom_point(aes(x=dist_km, y=full_ts_corr, color=comp_prov), size=3) +
  scale_color_manual(values=colz) +theme(panel.grid = element_blank(), legend.title = element_blank(),
                                         legend.position = "bottom") +facet_wrap(~provname)
p1

#then, get the average synchrony per province
pearsons.avg <- ddply(pearsons.sum, .(provname), summarise, mean_corr=mean(full_ts_corr ))

#then, look at avg synchrony per province per year
pearsons.avg.year <- ddply(pearsons.df, .(provname, year), summarise, mean_corr=mean(corr, na.rm=T))

pearsons.avg.year <- merge(pearsons.avg.year, centroid.prov, by="provname")

pearsons.avg.year <- arrange(pearsons.avg.year, latitude, year)
pearsons.avg.year$provname <- factor(pearsons.avg.year$provname, levels=unique(pearsons.avg.year$provname))


#and plot 
#vert.df <- cbind.data.frame(xint = c(2007, 2012,  2019))


vert.df <- cbind.data.frame(xint = c(2006.5, 2007.5, 2011.5, 2012.5,  2018.5, 2019.5))

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

# plot as heatmap - FigS14

FigS14b <- ggplot(data=pearsons.avg.year) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=year, y=provname, fill=mean_corr, color=mean_corr)) +
  scale_fill_viridis_c( option="inferno", name="Average annual\npairwise province Pearson's\ncorrelation coefficients", na.value = "black") +
  scale_color_viridis_c( option="inferno", name="Average annual\npairwise province Pearson's\ncorrelation coefficients", na.value = "black") +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.background = element_rect(fill="gray50"),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=.8) 

#and the annual distribution

dist.dat <- ddply(pearsons.df, .(year), summarise, median_corr = quantile(corr, na.rm=T)["50%"],  min_corr = quantile(corr, na.rm=T)["25%"], max_corr = quantile(corr, na.rm=T)["75%"])
head(dist.dat)
dist.dat$year <- as.numeric(as.character(dist.dat$year))

#and plot
FigS14a <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=year, ymin=min_corr, ymax=max_corr), alpha=.3) + ylab("annual distribution\npairwise province Pearson's\ncorrelation coefficients") +
  geom_line(aes(x=year, y=median_corr), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS14 <- cowplot::plot_grid(FigS14a, FigS14b, ncol=1, nrow=2, rel_heights = c(.25,1))


#now add some GAM plots


ggsave(file = paste0(homewd, "/final-figures/FigS14.png"),
       plot= FigS14,
       units="mm",  
       width=65, 
       height=70, 
       scale=3, 
       dpi=300)

p3 <-  ggplot(pearsons.df) + 
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_point(aes(x=year, y=corr, color=comp_prov)) + facet_wrap(~provname)

p3

# add pop size, mean temp, mean precip, and elevation as covariates in this
# and calculate predictors of the correlation coefficient
centroid.add <- dplyr::select(centroid.prov, provname, mean_elevation_m)
pearsons.df <- merge(pearsons.df, centroid.add, by="provname")
tsir.dat$provname[tsir.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

tsir.add <- ddply(tsir.dat, .(provname, year), summarise, pop=mean(pop))
#tsir.add <- ddply(tsir.dat, .(provname), summarise, pop=mean(pop))
pearsons.df <- merge(pearsons.df, tsir.add, by= c("provname", "year"))

#pearsons.df <- merge(pearsons.df, tsir.add, by= c("provname"))

#and temp and precip
mean.temp <- ddply(temp.dat, .(provname, year), summarise, mean_temp_C = mean(temp_C))

#mean.temp <- ddply(temp.dat, .(provname), summarise, mean_temp_C = mean(temp_C))
setdiff(unique(mean.temp$year), unique(pearsons.df$year))
setdiff(unique(pearsons.df$year), unique(mean.temp$year)) #no 2020 in climate data
setdiff(unique(mean.temp$provname), unique(pearsons.df$provname))

pearsons.df <- merge(pearsons.df, mean.temp, by=c("provname", "year"))

#pearsons.df <- merge(pearsons.df, mean.temp, by=c("provname"))

unique(precip.dat$provname)
precip.sum <- ddply(precip.dat, .(provname, year), summarise, sum_precip = sum(precip_mm))

#precip.sum <- ddply(precip.dat, .(provname), summarise, sum_precip = sum(precip_mm))

pearsons.df <- merge(pearsons.df, precip.sum, by=c("provname", "year"))

#pearsons.df <- merge(pearsons.df, precip.sum, by=c("provname"))

head(pearsons.df)
pearsons.df$log10_elevation_m <- log10(pearsons.df$mean_elevation_m)
pearsons.df$log10_temp_C <- log10(pearsons.df$mean_temp_C)
pearsons.df$log10_pop <- log10(pearsons.df$pop)
pearsons.df$log10_precip <- log10(pearsons.df$sum_precip)
pearsons.df$provname <- as.factor(pearsons.df$provname)
pearsons.df$year <- as.factor(pearsons.df$year)

gam1 <- gam(corr ~ dist_m:provname +
                  s(log10_temp_C, bs="tp") +
                  #s(sum_precip, bs="tp") +
                  #s(log10_elevation_m, bs="tp") +
                  s(log10_pop, bs="tp") +
                  s(log10_precip, bs="tp") +
                  #s(pop, bs="tp") +
                  s(year, bs="re") +
                  s(provname, bs="re"), #different y-intercept for each province
            data = pearsons.df)

summary(gam1) 
#strong negative correlation with increasing distance and correlation, sig for most provinces
source(paste0(homewd, "/figure-development/Fig1/mollentze-streicker-2020-functions.R"))

year.df <- get_partial_effects(fit=gam1, var="year")
prov.df <- get_partial_effects(fit=gam1, var="provname")
temp.df <- get_partial_effects_continuous(gamFit=gam1, var="log10_temp_C")
precip.df <- get_partial_effects_continuous(gamFit=gam1, var="log10_precip")
pop.df <- get_partial_effects_continuous(gamFit=gam1, var="log10_pop")

plot.partial(year.df, var="year", response_var = "corr") #all epi years significantly positive
plot.partial(prov.df, var="provname", response_var = "corr") #some provs more or less correlated
plot.partial.cont(df=temp.df, log=T, var="log10_temp_C", response_var = "corr", alt_var = "log10(temp (*C))")
plot.partial.cont(df=precip.df, log=T, var="log10_precip", response_var = "corr", alt_var ="sum precip (mm)")
plot.partial.cont(df=pop.df, log=T, var="log10_pop", response_var = "corr", alt_var = "pop")

#and 

#   
# 
# 
# p2 <- ggplot(pearsons.sum) + theme_bw() + ylab("full time series pearson's correlation between provinces") +
#   xlab("distance from Phnom Penh (km)") + coord_cartesian(ylim=c(0,1))+
#   geom_point(aes(x=dist_from_PP_km, y=full_ts_corr, color=comp_prov), size=3) +
#   scale_color_manual(values=colz) +theme(panel.grid = element_blank(), legend.title = element_blank(),
#                                          legend.position = "bottom") +facet_wrap(~provname)
