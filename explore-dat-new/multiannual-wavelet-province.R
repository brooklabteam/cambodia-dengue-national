rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
  
head(dat) 
min(dat$year)

#plot time series of each type by province by year

unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date)
dat$epimonth <- cut.Date(dat$date, breaks="months", start.on.monday = T)
dat$epimonth <- as.Date(as.character(dat$epimonth))

head(dat)

class(dat$epimonth)

dat$case=1

#total cases by week by year by province
dat.prov <- ddply(dat, .(epimonth, provname), summarise, cases = sum(case))
head(dat.prov)

#and for the whole dataset

dat.nat  <- ddply(dat, .(epimonth), summarise, cases = sum(case))
head(dat.nat)
dat.nat$epimonth_num <- 1:length(unique(dat.nat$epimonth))

p4 <- ggplot(data=dat.nat) + geom_line(aes(x=epimonth, y=cases)) 
   
p4


#now, see if you can get the wavelet periodicity by province
library(WaveletComp)

#try the national level
dat.nat <- arrange(dat.nat, epimonth)

#now try it for above 1 year!
anal.dat.nat.multi <-  analyze.wavelet(dat.nat,
                                 my.series = 2,
                                 #loess.span = 0,
                                 dt = 1/12,#this allows for annual timestep
                                 dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                 lowerPeriod = 2,#shortest possible period (2)
                                 upperPeriod = 20, #largest possible period (in weeks; here, 10 years)
                                 make.pval = TRUE, n.sim = 1000)

dat.nat.anual <-  analyze.wavelet(dat.nat,
                                       my.series = 2,
                                       #loess.span = 0,
                                       dt = 1/12,#this allows for annual timestep
                                       dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                       lowerPeriod = 1/12,#shortest possible period (2)
                                       upperPeriod = 2, #largest possible period (in weeks; here, 10 years)
                                       make.pval = TRUE, n.sim = 1000)

#multi
wt.image(anal.dat.nat.multi, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         plot.coi = F,
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))


#annual
wt.image(dat.nat.anual, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         plot.coi = F,
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

reconstruct(dat.nat.anual, only.ridge = T, period=1)
wt.sel.phases(dat.nat.anual,    spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                  labels = 2002:2019))


avg.powr <- colSums(anal.dat.nat.multi$Power)/333
plot(avg.powr) #multi-annual

avg.powr.annual <- colSums(dat.nat.anual$Power)/333
plot(avg.powr.annual) #multi-annual
#save... this is the national multi-annual periodicity - both cycles getting longer

maximum.level = 1.001*max(anal.dat.nat.multi$Power.avg, anal.dat.nat.multi$Power.avg)
wt.avg(anal.dat.nat.multi, maximum.level = maximum.level)

avg.dat <- cbind.data.frame(power = anal.dat.nat.multi$Power.avg, log2period=log2(anal.dat.nat.multi$Period), period = anal.dat.nat.multi$Period)
head(avg.dat)
ggplot(avg.dat) + geom_point(aes(x=power, y=log2period))
ggplot(avg.dat) + geom_point(aes(y=power, x=period)) + theme_bw()


avg.dat$period[avg.dat$power== max(avg.dat$power)]#3.11
avg.dat$period[avg.dat$period>5 & avg.dat$power== max(avg.dat$power[avg.dat$period>5])]#5.81589


reconstruct(anal.dat.nat.multi, only.ridge = T)
reconstruct(anal.dat.nat.multi, sel.period = 3,  spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                                       labels = 2002:2019))
reconstruct(anal.dat.nat.multi, sel.period = 6,  spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                                       labels = 2002:2019))

#was 2019 the compilation of two different offset multiannual cycles????


#now split by province, deconvolute and plot them all

prov.split <- dlply(dat.prov, .(provname))

#and apply
wavelet.app <- function(dat, my.series, n.sim, dt, lower.period, upper.period){
  anal.dat <-  analyze.wavelet(dat,
                              my.series = my.series,
                              #loess.span = 0,
                              dt = dt,#this allows for annual
                              dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                              lowerPeriod = lower.period,#shortest possible period (in weeks; here, one year)
                              upperPeriod = upper.period, #largest possible period (in weeks; here, 10 years)
                              make.pval = TRUE, n.sim = n.sim)
  return (anal.dat)
}

prov.out <- lapply(prov.split, wavelet.app, dt = 1/12, my.series=3, lower.period=2, upper.period=30, n.sim=100)

#and save each and populate plots

for (i in 1:length(prov.split)){
  png(paste0(homewd, "/test-figs/wave/multi_prov_", names(prov.split[i]),".png"))
  wt.image(prov.out[[i]], color.key = "quantile", n.levels = 250,
                             legend.params = list(lab = "wavelet power levels", mar = 4.7),
                             label.time.axis = TRUE,
                             plot.coi = F,
                             main= names(prov.split[i]),
           spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)",
           spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                 labels = 2002:2019))
  
  
  dev.off()

  #save then import them back in as ggplot images
}

#and load them and arrange
library(png)

p.list <- list()
for (i in 1:25){
  tmp <- readPNG(paste0(homewd, "/test-figs/wave/multi_prov_", names(prov.split[i]),".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}

pall <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)



ggsave(file = paste0(homewd,"/test-figs/MultiAllProvinceWavelet.pdf"),
       plot=pall,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)

#arrange in reverse order of mean age of infection?
#lowest to highest


load(paste0(homewd, "/data/meanagebyprov.Rdata"))
mean.age.df <- arrange(mean.age.df, avg_age)
p.list <- list()
for (i in 1:25){
  #tmp <- readPNG(paste0(homewd, "/test-figs/wave/prov_", names(prov.split[i]),".png"))
  tmp <- readPNG(paste0(homewd, "/test-figs/wave/multi_prov_", mean.age.df$ADM1_EN[i],".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}


pall2 <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)


ggsave(file = paste0(homewd,"/test-figs/AllProvinceWaveletByMeanAgeMulti.pdf"),
       plot=pall2,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


#now look for mult-annual oscillations in temperature and precip
#load in temperature data from Katie
katwd = "/Users/carabrook/Developer/cambodia-dengue-province"
temp.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )


#get mean temp across study period for precip and temp
temp.mean <- ddply(temp.dat, .(adm1_name), summarise, mean_temp_C = mean(temp_C), median_temp_C = median(temp_C), quarter_quant_temp_C=quantile(temp_C)[2], threequarter_quant_temp=quantile(temp_C)[4])
#and precip
precip.mean <- ddply(precip.dat, .(adm1_name), summarise, mean_precip_mm = mean(precip_mm), median_precip_mm = median(precip_mm), quarter_quant_precip_mm=quantile(precip_mm)[2], threequarter_quant_precip_mm=quantile(precip_mm)[4])


#and merge
climate.mean <-  temp.mean
climate.mean <- merge(climate.mean, precip.mean, by="adm1_name")
climate.mean <- climate.mean[2:nrow(climate.mean),]

#save 
#write.csv(climate.mean, file = paste0(homewd, "/data/climate_mean_dat.csv"), row.names = F)
head(temp.dat)
temp.dat <- dplyr::select(temp.dat, -(X))

#make a time series with date
unique(temp.dat$biwk)
months <- rep(seq(1,12,1), each=2)
biwks <- seq(1,26,1)
days <- seq(1,365,14)[-1]
dates <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), 1)
dates <- paste(sapply(strsplit(as.character(dates[days]), "-"), "[",2), sapply(strsplit(as.character(dates[days]), "-"), "[",3), sep="-")
years <- unique(temp.dat$year)


#link dates to the data
date.bind <- cbind.data.frame(year = rep(years, each=26), dates = rep(dates, time=length(years)), biwk= rep(biwks, time=length(years)) )
date.bind$dates <- paste0(date.bind$year, "-", date.bind$dates)
head(date.bind)

temp.dat <- merge(temp.dat, date.bind, by=c("year", "biwk"))
head(temp.dat)

temp.dat$dates <- as.Date(temp.dat$dates)

ggplot(data=temp.dat) + geom_line(aes(x=dates, y=temp_C,color=adm1_name)) +
      facet_wrap(~adm1_name, ncol=5)


#get national avg. - is this useful?
temp.dat$month <- month(temp.dat$dates)
temp.dat$month_date <- paste0(temp.dat$year, "-", temp.dat$month, "-01")
temp.dat$month_date <- as.Date(as.character(temp.dat$month_date))

temp.dat <- ddply(temp.dat, .(adm1_name, month_date, year, month), summarise, temp_K = mean(temp_K), temp_C=mean(temp_C) )
head(temp.dat)
nat.temp <- ddply(temp.dat, .(month_date, year, month), summarise, temp_K = mean(temp_K), temp_C=mean(temp_C) )
ggplot(data=nat.temp) + geom_line(aes(x=month_date, y=temp_C)) 

head(nat.temp)

nat.temp <- arrange(nat.temp, month_date)

anal.temp <-  analyze.wavelet(nat.temp,
                                my.series = 5,
                                #loess.span = 0,
                                dt = 1/12,#these are annual cycles
                                dj = 1/100,
                                lowerPeriod = 2,
                                upperPeriod = 20,
                                make.pval = TRUE, n.sim = 100)

wt.image(anal.temp, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         plot.coi = F,
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))



maximum.level = 1.001*max(anal.temp$Power.avg, anal.temp$Power.avg)
wt.avg(anal.temp, maximum.level = maximum.level)
reconstruct(anal.temp,  spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))


avg.dat <- cbind.data.frame(power = anal.temp$Power.avg, log2period=log2(anal.temp$Period), period = anal.temp$Period)
head(avg.dat)
ggplot(avg.dat) + geom_point(aes(x=power, y=log2period))
ggplot(avg.dat) + geom_point(aes(y=power, x=period)) + theme_bw()


avg.dat$period[avg.dat$power== max(avg.dat$power)]#2.94


        
#and also precip

precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )

head(precip.dat)
precip.dat <- dplyr::select(precip.dat, -(X))

#make a time series with date
unique(precip.dat$biwk)
months <- rep(seq(1,12,1), each=2)
biwks <- seq(1,26,1)
days <- seq(1,365,14)[-1]
dates <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), 1)
dates <- paste(sapply(strsplit(as.character(dates[days]), "-"), "[",2), sapply(strsplit(as.character(dates[days]), "-"), "[",3), sep="-")
years <- unique(temp.dat$year)


#link dates to the data
date.bind <- cbind.data.frame(year = rep(years, each=26), dates = rep(dates, time=length(years)), biwk= rep(biwks, time=length(years)) )
date.bind$dates <- paste0(date.bind$year, "-", date.bind$dates)
head(date.bind)

precip.dat <- merge(precip.dat, date.bind, by=c("year", "biwk"))
head(precip.dat)

precip.dat$month <- month(precip.dat$dates)
precip.dat$month_date <- as.Date(as.character(paste0(precip.dat$year, "-", precip.dat$month, "-01")))


precip.dat <- ddply( precip.dat, .(adm1_name, month_date, year, month), summarise, precip_mm = max(precip_mm), precip_m = max(precip_m) )

precip.nat <- ddply( precip.dat, .(month_date, year, month), summarise, precip_mm = mean(precip_mm), precip_m = mean(precip_m) )

ggplot(data=precip.dat) + geom_line(aes(x=month_date, y=precip_mm)) 

head(precip.nat)

precip.dat <- arrange(precip.dat, month_date)

anal.precip <-  analyze.wavelet(precip.nat,
                              my.series = 4,
                              #loess.span = 0,
                              dt = 1/12,#these are annual cycles
                              dj = 1/100,
                              lowerPeriod = 2,
                              upperPeriod = 20,
                              make.pval = TRUE, n.sim = 100)

wt.image(anal.precip, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         plot.coi = F,
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

reconstruct(anal.precip,  spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                labels = 2002:2019))

maximum.level = 1.001*max(anal.precip$Power.avg, anal.precip$Power.avg)
wt.avg(anal.precip, maximum.level = maximum.level)
reconstruct(anal.precip,  spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))


avg.dat <- cbind.data.frame(power = anal.precip$Power.avg, log2period=log2(anal.precip$Period), period = anal.precip$Period)
head(avg.dat)
ggplot(avg.dat) + geom_point(aes(x=power, y=log2period))
ggplot(avg.dat) + geom_point(aes(y=power, x=period)) + theme_bw()


avg.dat$period[avg.dat$power== max(avg.dat$power)]#4.72





#now get province-level for temp and precip.
temp.dat = subset(temp.dat, adm1_name!="Administrative unit not available")
prov.temp <- dlply(temp.dat, .(adm1_name))

prov.out <- lapply(prov.temp, wavelet.app, dt = 1/12, my.series=6, lower.period=2, upper.period=30, n.sim=100)


for (i in 1:length(prov.temp)){
  png(paste0(homewd, "/test-figs/tempwave/multi_prov_", names(prov.out[i]),".png"))
  wt.image(prov.out[[i]], color.key = "quantile", n.levels = 250,
           legend.params = list(lab = "wavelet power levels", mar = 4.7),
           label.time.axis = TRUE,
           plot.coi = F,
           main= names(prov.split[i]),
           spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)",
           spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                 labels = 2002:2019))
  
  
  dev.off()
  
  #save then import them back in as ggplot images
}

#and load them and arrange
library(png)

p.list <- list()
for (i in 1:25){
  tmp <- readPNG(paste0(homewd, "/test-figs/tempwave/multi_prov_", names(prov.out[i]),".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}

pall <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)



ggsave(file = paste0(homewd,"/test-figs/MultiAllProvinceTemp.pdf"),
       plot=pall,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


#and precip
precip.dat = subset(precip.dat, adm1_name!="adm1_name")
prov.precip <- dlply(precip.dat, .(adm1_name))

precip.out <- lapply(prov.precip, wavelet.app, dt = 1/12, my.series=5, lower.period=2, upper.period=30, n.sim=100)


for (i in 1:length(prov.precip)){
  png(paste0(homewd, "/test-figs/precipwave/multi_prov_", names(prov.precip[i]),".png"))
  wt.image(precip.out[[i]], color.key = "quantile", n.levels = 250,
           legend.params = list(lab = "wavelet power levels", mar = 4.7),
           label.time.axis = TRUE,
           plot.coi = F,
           main= names(prov.split[i]),
           spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)",
           spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                 labels = 2002:2019))
  
  
  dev.off()
  
  #save then import them back in as ggplot images
}

#and load them and arrange
library(png)

p.list <- list()
for (i in 1:25){
  tmp <- readPNG(paste0(homewd, "/test-figs/precipwave/multi_prov_", names(prov.out[i]),".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}

pall <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)



ggsave(file = paste0(homewd,"/test-figs/MultiAllProvincePrecip.pdf"),
       plot=pall,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)














