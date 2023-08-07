rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

#load the province-level transmission
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
head(beta.df)


#also load and prep the climate data at the province level
katwd = "/Users/carabrook/Developer/cambodia-dengue-province"
temp.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )

names(temp.dat)[names(temp.dat)=="adm1_name"] <- "provname"
names(precip.dat)[names(precip.dat)=="adm1_name"] <- "provname"

temp.dat = subset(temp.dat, provname!="Administrative unit not available")
precip.dat = subset(precip.dat, provname!="Administrative unit not available")
head(temp.dat) # year and biweek
head(precip.dat) # year and biweek

setdiff(unique(temp.dat$provname), unique(beta.df$provname)) #"Mondul Kiri"    "Oddar Meanchey" "Ratanak Kiri"   "Siemreap"       "Tboung Khmum"  
temp.dat$provname[temp.dat$provname=="Siemreap"] <- "Siem Reap"
temp.dat$provname[temp.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"
precip.dat$provname[precip.dat$provname=="Siemreap"] <- "Siem Reap"
precip.dat$provname[precip.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"

#and plot to check
temp.dat$year <-as.factor(temp.dat$year)

ggplot(temp.dat) + geom_line(aes(x=biwk, y=temp_C, color=year)) + facet_wrap(provname~.)

precip.dat$year <-as.factor(precip.dat$year)
ggplot(precip.dat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) + facet_wrap(provname~.)

#need to correct precip.dat for 2020 and 2021 (even though will be dropped in merge)
precip.dat$year <- as.numeric(as.character(precip.dat$year))
precip.dat$precip_mm[precip.dat$year>=2020] <- precip.dat$precip_mm[precip.dat$year>=2020]/1000
precip.dat$precip_m[precip.dat$year>=2020] <- precip.dat$precip_m[precip.dat$year>=2020]/1000

precip.dat$year <-as.factor(precip.dat$year)
ggplot(precip.dat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) + facet_wrap(provname~.)

head(precip.dat)
head(temp.dat)
# remove year 2020 since we don't use it anyway
precip.dat$year <- as.numeric(as.character(precip.dat$year))
temp.dat$year <- as.numeric(as.character(temp.dat$year))
precip.dat = subset(precip.dat, year>2001 & year <2020)
temp.dat = subset(temp.dat, year>2001 & year <2020)

# and test whether precip and temp are increasing through time
precip.dat$provname <- as.factor(precip.dat$provname)
temp.dat$provname <- as.factor(temp.dat$provname)

#remove Tboung Khmum becuase it was created too recently
precip.dat = subset(precip.dat, provname!="Tboung Khmum")
temp.dat = subset(temp.dat, provname!="Tboung Khmum")

# Here, we look at precipitation through time, allowing for 
# a different intercept for each province
gam1 <- gam(precip_mm~ year + 
                       s(biwk, bs="cc", k=7) + 
                       s(provname, bs="re"), data=precip.dat)
summary(gam1) #precip not changing over the longer term


gam2 <- gam(precip_mm~ 
              s(year, by=provname, k=3)+ 
              s(biwk, bs="cc", k=7) + 
              s(provname, bs="re"), data=precip.dat)
summary(gam2) #precip not changing over the longer term for any province

AIC(gam1,gam2) #slightly better fit when allowing for different slopes by province

precip.predict <- cbind.data.frame(provname = rep(unique(precip.dat$provname), each=18), year = rep(unique(precip.dat$year), 24))
precip.predict$biwk = 10
precip.predict$precip_mm <- predict.gam(gam2, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit
precip.predict$precip_mm_lci <- predict.gam(gam2, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit - 1.96*predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se
precip.predict$precip_mm_uci <- predict.gam(gam2, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit + 1.96*predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se


# Here, we look at temperature through time, allowing for
# a different intercept for each province
gam3 <- gam(temp_C ~ year + 
              s(biwk, bs="cc", k=7) + 
              s(provname, bs="re"),
              data=temp.dat)

summary(gam3) 
# temp is increasing slightly, with different intercepts for province

# Here, a version of the GAM, allowing for both different annual slopes and intercepts by province
gam4 <- gam(temp_C ~ s(year, k=3, by=provname) + 
                     s(biwk, bs="cc", k=7) +
                     s(provname, bs="re"),
                     data=temp.dat)

summary(gam4)
# temp is still increasing, with some variation among provinces - but generally all positive

AIC(gam3, gam4) # version with distinct slope and intercept is the best

#extract partial effects of province on temperature
source(paste0(homewd, "/figure-development/Fig1-S1/mollentze-streicker-2020-functions.R"))
prov.df <- get_partial_effects(fit=gam4, var="provname")$effects

#and plot it 

prov.df$color <- "NotSig"
prov.df$color[prov.df$IsSignificant=="Yes" & prov.df$y>0] <- "Pos"
prov.df$color[prov.df$IsSignificant=="Yes" & prov.df$y<0] <- "Neg"

colz = c('Pos' = "red", 'Neg'="blue", "NotSig"="gray65")
prov.df$provname <- as.character(prov.df$provname)
prov.df$provname[prov.df$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
prov.df$provname <- as.factor(prov.df$provname)

p1 <- ggplot(data=prov.df) +
      geom_point(aes(x=provname, y=y, color=color), show.legend = F, size=3) +
      geom_linerange(aes(x=provname, ymin=ylower, ymax=yupper, color=color), show.legend = F, size=.8) +
      scale_color_manual(values=colz) + theme_bw() +
      theme(panel.grid = element_blank(), axis.title.x = element_blank(),
            axis.text.x = element_text(angle=90, size=13),
            plot.margin = unit(c(.2,.1,.4,.1), "cm"),
            axis.text.y = element_text(size=13), axis.title.y = element_text(size=16)) +
      ylab("partial effect of province on temperature") + geom_hline(aes(yintercept=0))





# here are predictions for best fit model
temp.predict <- cbind.data.frame(provname = rep(unique(temp.dat$provname), each=18), year = rep(unique(temp.dat$year), 24))
temp.predict$biwk = 10
temp.predict$temp_C <- predict.gam(gam4, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit
temp.predict$temp_C_lci <- predict.gam(gam4, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit - 1.96*predict.gam(gam1, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se
temp.predict$temp_C_uci <- predict.gam(gam4, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit + 1.96*predict.gam(gam1, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se


# and merge with case data
temp.dat <- dplyr::select(temp.dat, -(X), -(temp_K))
names(temp.dat)[names(temp.dat)=="biwk"] <- "biweek"

precip.dat <- dplyr::select(precip.dat, -(X))
names(precip.dat)[names(precip.dat)=="biwk"] <- "biweek"

head(beta.df)
unique(beta.df$year)

climate.merge <- merge(beta.df, temp.dat, by = c("provname", "year", "biweek"), all.x = T)
head(climate.merge)
climate.merge <- merge(climate.merge, precip.dat, by = c("provname", "year", "biweek"), all.x = T)
head(climate.merge)

climate.merge$year <- as.factor(climate.merge$year)

#change provname for plotting
climate.merge$provname[climate.merge$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

pS2 <- ggplot(climate.merge) + geom_line (aes(x=biweek, y=precip_mm, color=year)) + facet_wrap(provname~.) + 
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        ylab("biweekly sum precipitation (mm)") + xlab("biweek of year")

ggsave(file = paste0(homewd,"/final-figures/FigS2.png"),
       plot = pS2,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


pS3 <- ggplot(climate.merge) + geom_line(aes(x=biweek, y=temp_C, color=year)) + facet_wrap(provname~.) + 
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
        ylab(bquote('biweekly mean temperature ('^0~'C)')) + xlab("biweek of year")

ggsave(file = paste0(homewd,"/final-figures/FigS3.png"),
       plot = pS3,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


# and plot the time series foward
head(climate.merge)
head(precip.predict)
precip.predict$time <- precip.predict$year +.5
temp.predict$time <- temp.predict$year +.5


precip.predict$provname <- as.character(precip.predict$provname)
temp.predict$provname <- as.character(temp.predict$provname)
precip.predict$provname[precip.predict$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
temp.predict$provname[temp.predict$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

precip.predict$provname <- as.factor(precip.predict$provname)
temp.predict$provname <- as.factor(temp.predict$provname)


pS4 <- ggplot(subset(climate.merge, time<2020)) +
       geom_line(aes(x=time, y=precip_mm, color=provname), show.legend = F) + facet_wrap(provname~.) + 
       theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), axis.title.x = element_blank())+
       geom_ribbon(data=precip.predict, aes(time, ymin=precip_mm_lci, ymax=precip_mm_uci),alpha=.4) +
       geom_line(data=precip.predict, aes(time, precip_mm), size=.8) + ylab("biweekly sum precipitation (mm)")
  
ggsave(file = paste0(homewd,"/final-figures/FigS4.png"),
       plot = pS4,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


pS5 <- ggplot(subset(climate.merge, time<2020)) +
    geom_line(aes(x=time, y=temp_C, color=provname), show.legend = F) + facet_wrap(provname~.) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), axis.title.x = element_blank())+
    geom_ribbon(data=temp.predict, aes(time, ymin=temp_C_lci, ymax=temp_C_uci),alpha=.4) +
    geom_line(data=temp.predict, aes(time, temp_C), size=.8) + ylab(bquote('biweekly mean temperature ('^0~'C)'))

ggsave(file = paste0(homewd,"/final-figures/FigS5.png"),
       plot = pS5,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)

#climate predictors by province
p2 <- ggplot(temp.predict) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(),
                     axis.title.y = element_text(size=16), axis.text = element_text(size=13),
                     plot.margin = unit(c(.2,.3,.1,.3), "cm"), legend.position = "bottom")+
  geom_ribbon(data=temp.predict, aes(time, ymin=temp_C_lci, ymax=temp_C_uci, fill=provname),alpha=.4) +
  geom_line(data=temp.predict, aes(time, temp_C, color=provname), size=.8) + 
  ylab(bquote('annual mean temperature prediction ('^0~'C)')) #+ 
  #guides(fill=guide_legend(ncol=1), color= guide_legend(ncol=1))


pS6 <- cowplot::plot_grid(p1,p2, ncol=2, nrow=1, rel_widths = c(1,1.2))

ggsave(file = paste0(homewd,"/final-figures/FigS6.png"),
       plot = pS6,
       units="mm",  
       width=110, 
       height=60, 
       scale=3, 
       dpi=300)


#drop the sus reconstruction info
climate.merge <- dplyr::select(climate.merge, -(rsquared), -(sus_reconstruction))

#save this

write.csv(climate.merge, file = paste0(homewd, "/data/climate_beta_prov.csv"), row.names = F)

# and test for lags now - split by province first

climate.split <- dlply(climate.merge, .(provname, epiyr))

# Look for cross-correlations between 
# Look for cross-correlations between 
find.lags <- function(dat){
  
  #first, take out the actual epiyear that was not involved in the caculation
  dat = subset(dat, year!=unique(dat$epiyr))
  
  #precip
  dat.lag <- cbind.data.frame(lag = print(ccf(dat$precip_mm, dat$beta))$lag, acf=print(ccf(dat$precip_mm, dat$beta))$acf)
  dat.lag$variable <- "precip_mm"
  dat.lag$lag[dat.lag$acf==max(dat.lag$acf[dat.lag$lag<0 & dat.lag$lag>-27])] 
  
  
  
  dat2 = cbind.data.frame(lag = print(ccf(dat$temp_C, dat$beta))$lag, acf=print(ccf(dat$temp_C, dat$beta))$acf)
  dat2$variable <- "temp_C"
  dat2$lag[dat2$acf==max(dat2$acf[dat2$lag<0 & dat2$lag>-27])]
  dat.lag <- rbind(dat.lag, dat2)
  
  max.lag <- dlply(dat.lag, .(variable))
  get.lag <- function(df){
    lag = df$lag[df$acf==max(df$acf[df$lag<0 &df$lag>-27])]
    if(length(lag)>1){
      lag = max(lag[lag<0])
    }
    df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
    return(df.out)
  }
  max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))
  
  max.lag$epiyr <- unique(dat$epiyr)
  max.lag$provname <- unique(dat$provname)
  
  return(max.lag)
  
  
}


lag.climate.out <- lapply(climate.split, find.lags)
lag.climate.df <- data.table::rbindlist(lag.climate.out)

ddply(lag.climate.df, .(variable), summarise, mean_lag=mean(lag), median_lag = median(lag)) 
#mean precip = -5, temp = -6
#median precip = -2, temp = -7

#and save lag data
write.csv(lag.climate.df, file=paste0(homewd,"/data/lags_climate_prov.csv"), row.names = F)


#make lagged data based on optimal lags

#make shifted dataset and save
head(climate.merge)

climate.shift <- climate.merge[8:nrow(climate.merge),] #starts in year 8 after 7 timestep shift
climate.shift$temp_C_lag <- climate.merge$temp_C[1:(length(climate.merge$temp_C)-7)]
head(climate.shift)

#and the shifted precip.
climate.shift$precip_mm_lag <- climate.merge$precip_mm[6:(length(climate.merge$precip_mm)-2)]

#write data and do lagged regression in another script
write.csv(climate.shift, file = paste0(homewd, "/data/lagged-prov-clim-beta.csv"), row.names = F)

