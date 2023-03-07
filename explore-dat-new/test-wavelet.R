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
dat$epiwk <- cut.Date(dat$date, breaks="weeks", start.on.monday = T)
dat$epiwk <- as.Date(as.character(dat$epiwk))
dat$epiwk[dat$epiwk< "2002-01-01"] <- "2002-01-01"
head(dat)

class(dat$epiwk)

#now attach biwks (2 biweeks per week)
wks = 1:52
biwks <- rep(1:26, each=2)
wk.biwk <- cbind.data.frame(week_report=wks, biwk=biwks)

sort(unique(dat$week_report))
unique(dat$provname)

dat <- merge(dat, wk.biwk, by="week_report")
head(dat) #now plot as the lowest epiweek date per biweek
dat$case=1

#total cases by week by year by province
dat.prov <- ddply(dat, .(epiwk, provname), summarise, cases = sum(case))
head(dat.prov)

#and for the whole dataset

dat.nat  <- ddply(dat, .(epiwk), summarise, cases = sum(case))
head(dat.nat)
dat.nat$epiwk_num <- 1:length(unique(dat.nat$epiwk))

p4 <- ggplot(data=dat.nat) + geom_line(aes(x=epiwk, y=cases)) 
   
p4


#now, see if you can get the wavelet periodicity by province
library(WaveletComp)

#try the national level
dat.nat <- arrange(dat.nat, epiwk)
anal.dat.nat <-  analyze.wavelet(dat.nat,
                                my.series = 2,
                                #loess.span = 0,
                                dt = 1/52,#this allows for annual
                                dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                lowerPeriod = 1,#shortest possible period (in weeks; here, one year)
                                upperPeriod = 30, #largest possible period (in weeks; here, 10 years)
                                make.pval = TRUE, n.sim = 100)



wt.image(anal.dat.nat, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         spec.period.axis = list(at = c(1:10,12,14,16), labels = c(1:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1, (18*52), 52), #converting weeks to years
                               labels = 2002:2019))


#extract the dates and corresponding period where the power spectrum is maximized

hist(anal.dat.nat$Power) #let's say above 80
which(anal.dat.nat$Power>80)
anal.dat.nat$Power[442883]

power.dat <- melt(anal.dat.nat$Power)
power.p.dat <- melt(anal.dat.nat$Power.pval)
head(power.dat)
head(power.p.dat)

epitime <- cbind.data.frame(Var2=dat.nat$epiwk_num, epiwk= dat.nat$epiwk)

power.dat <- merge(power.dat, epitime, by="Var2")
power.dat <- dplyr::select(power.dat, -(Var2))
head(power.dat)
names(power.dat)[names(power.dat)=="value"] <- "power"
#var 1 is period
#var 2 is time
ggplot(power.dat) + geom_tile(aes(x=epiwk, y=Var1, fill=power))


period.dat <- cbind.data.frame(period = unique(anal.dat.nat$Period))
period.dat$Var1 = 1:nrow(period.dat)
head(period.dat)


power.dat <- merge(power.dat, period.dat, all.x = T, by="Var1")
power.dat <- dplyr::select(power.dat, -(Var1))

power.dat$pval <- power.p.dat$value
head(power.dat)

power.dat$period_trunc <- trunc(power.dat$period)
subset(power.dat, pval<0.01 & period_trunc>1)


# ggplot(power.dat) + geom_tile(aes(x=epiwk, y=period, fill=log100(power))) 
# 
# subset(power.dat, power>80)
# test = subset(power.dat, power>60 & power<80)
# 
# power.dat =arrange(power.dat, desc(power))
# 
# #remove the ones
# power.sub = subset(power.dat, period_trunc>1)
# 
# head(power.sub)

#annual cycles with some synchrony at 3 year cycles too that emerges later in the time series
#there could be a demographic model added in here to show if this is the direction in which this goes


#now split by province, deconvolute and plot them all

prov.split <- dlply(dat.prov, .(provname))

#and apply
wavelet.app <- function(dat, my.series, n.sim){
  anal.dat <-  analyze.wavelet(dat,
                              my.series = my.series,
                              #loess.span = 0,
                              dt = 1/52,#this allows for annual
                              dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                              lowerPeriod = 1,#shortest possible period (in weeks; here, one year)
                              upperPeriod = 30, #largest possible period (in weeks; here, 10 years)
                              make.pval = TRUE, n.sim = n.sim)
  return (anal.dat)
}

prov.out <- lapply(prov.split, wavelet.app, my.series=3, n.sim=100)

#and save each and populate plots



for (i in 1:length(prov.split)){
  png(paste0(homewd, "/test-figs/wave/prov_", names(prov.split[i]),".png"))
  wt.image(prov.out[[i]], color.key = "quantile", n.levels = 250,
                             legend.params = list(lab = "wavelet power levels", mar = 4.7),
                             label.time.axis = TRUE,
                             plot.coi = F,
                             spec.period.axis = list(at = c(1:10,12,14,16), labels = c(1:10,12,14,16)),
                             periodlab= "period (in years)",
                             main= names(prov.split[i]),
                             spec.time.axis = list(at = seq(1, (18*52), 52), #converting weeks to years
                                                   labels = 2002:2019))
  
  dev.off()

  #save then import them back in as ggplot images
}

#and load them and arrange
library(png)

p.list <- list()
for (i in 1:25){
  tmp <- readPNG(paste0(homewd, "/test-figs/wave/prov_", names(prov.split[i]),".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}

pall <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)



ggsave(file = paste0(homewd,"/test-figs/AllProvinceWavelet.pdf"),
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
  tmp <- readPNG(paste0(homewd, "/test-figs/wave/prov_", mean.age.df$ADM1_EN[i],".png"))
  p.list[[i]] <- rasterGrob(tmp, interpolate=TRUE)
}


pall2 <- cowplot::plot_grid(plotlist = p.list, ncol=5, nrow=5)


ggsave(file = paste0(homewd,"/test-figs/AllProvinceWaveletByMeanAge.pdf"),
       plot=pall2,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


#and order by birth rate

#and plot 

#now try looking at periodicity by province and latitude???





#pull one province and try
dat.SiemReap = subset(dat.prov, provname=="Siem Reap")

head(dat.SiemReap)
dat.SiemReap <- arrange(dat.SiemReap, epiwk)
anal.dat.SR <-  analyze.wavelet(dat.SiemReap,
                              my.series = 4,
                              #loess.span = 0,
                              dt = 1,#these are biwks
                              dj = 1/250,
                              lowerPeriod = 1,
                              upperPeriod = 494,
                              make.pval = TRUE, n.sim = 10)

wt.image(anal.dat.SR, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),
         label.time.axis = TRUE,
         spec.time.axis = list(at = seq(1, (17*26+26), 26),
                               labels = 2001:2018))
#annual cycles with some synchrony at 3 year cycles too

dat.Kandal = subset(dat.prov, provname=="Kandal")
anal.dat.Kandal <-  analyze.wavelet(dat.Kandal,
                                my.series = 3,
                                #loess.span = 0,
                                dt = 1/52,#these are annual cycles
                                dj = 1/20,
                                make.pval = TRUE, n.sim = 10)

wt.image(anal.dat.Kandal, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7)) #same

#
dat.PP = subset(dat.prov, provname=="Phnom Penh")
anal.dat.PP <-  analyze.wavelet(dat.PP,
                                    my.series = 3,
                                    #loess.span = 0,
                                    dt = 1/52,#these are annual cycles
                                    dj = 1/20,
                                    make.pval = TRUE, n.sim = 10)

wt.image(anal.dat.PP, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7)) #same


#compare PP with SR
SR.merge <- dplyr::select(dat.SiemReap, epiwk, cases)
dat.PPSR <- merge(dat.PP, SR.merge, by="epiwk")
head(dat.PPSR)
coherent.dat <- analyze.coherency(dat.PPSR, my.pair = c(3,4),
                                  #loess.span = 0,
                                  dt = 1/52,#these are annual cycles
                                  dj = 1/20,
                                  make.pval = TRUE, n.sim = 10)
wt.image(coherent.dat, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7)) #same


# then, later, do the same for climate and compare temp cycling with cases

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
write.csv(climate.mean, file = paste0(homewd, "/data/climate_mean_dat.csv"), row.names = F)
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

temp.dat <- arrange(temp.dat, adm1_name, dates)
test.dat = subset(temp.dat, adm1_name=="Siemreap")
anal.temp.SR <-  analyze.wavelet(test.dat,
                                my.series = 4,
                                #loess.span = 0,
                                dt = 1,#these are annual cycles
                                dj = 1/250,
                                lowerPeriod = 1,
                                upperPeriod = 494,
                                make.pval = TRUE, n.sim = 10)

wt.image(anal.temp.SR, color.key = "quantile", n.levels = 250, periodlab = "period (biweeks)",
         legend.params = list(lab = "wavelet power levels"),
         label.time.axis = TRUE) #same

wt.image(anal.temp.SR, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"),
         label.time.axis = TRUE,
         spec.time.axis = list(at = seq(1, (17*26+26), 26),
                               labels = 2001:2018))
         

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

precip.dat$dates <- as.Date(precip.dat$dates)

p1 <- ggplot(data=precip.dat) + 
  geom_line(aes(x=dates, y=precip_mm,color=adm1_name)) +
  facet_wrap(~adm1_name, ncol=5) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


#and zero in
ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(precip.dat, adm1_name=="Koh Kong" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)





ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(precip.dat, adm1_name=="Koh Kong" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)




ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(precip.dat, adm1_name=="Koh Kong" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=precip_mm)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)








#and zero in
ggplot(data=subset(temp.dat, adm1_name=="Phnom Penh" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(temp.dat, adm1_name=="Koh Kong" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(temp.dat, adm1_name=="Mondul Kiri" & dates > "2018-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)





ggplot(data=subset(temp.dat, adm1_name=="Phnom Penh" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(temp.dat, adm1_name=="Koh Kong" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(temp.dat, adm1_name=="Mondul Kiri" & dates < "2014-01-01" & dates>"2011-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)




ggplot(data=subset(temp.dat, adm1_name=="Phnom Penh" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


ggplot(data=subset(temp.dat, adm1_name=="Koh Kong" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

ggplot(data=subset(temp.dat, adm1_name=="Mondul Kiri" & dates < "2009-01-01" & dates>"2006-01-01")) + geom_line(aes(x=dates, y=temp_C)) +
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)
















precip.dat <- arrange(precip.dat, adm1_name, dates)
test.precip = subset(precip.dat, adm1_name=="Siemreap")
anal.precip.SR <-  analyze.wavelet(test.precip,
                                 my.series = 4,
                                 #loess.span = 0,
                                 dt = 1,#these are annual cycles
                                 dj = 1/250,
                                 lowerPeriod = 1,
                                 upperPeriod = 494,
                                 make.pval = TRUE, n.sim = 10)

wt.image(anal.precip.SR, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"),
         label.time.axis = TRUE,
         spec.time.axis = list(at = seq(1, (17*26+26), 26),
                               labels = 2001:2018))


head(test.precip)

