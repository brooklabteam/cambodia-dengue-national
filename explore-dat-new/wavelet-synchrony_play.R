rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

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

#synchrony of national dengue ts and climate data
#then try synchrony in the cases by province


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


#now look at synchrony in national cases and prcip/temp
head(nat.temp)
head(dat.nat)

names(dat.nat)[1] <- "month_date"
temp.case <- merge(nat.temp, dat.nat, by="month_date")
head(temp.case)

ggplot(temp.case) + geom_line(aes(x=month_date, y=cases)) +
  geom_line(aes(x=month_date, y=temp_C*500), color="red")

my.wc <- analyze.coherency(temp.case, my.pair = c("temp_C", "cases"),
                           loess.span = 0,
                           dt = 1/12, dj = 1/100,
                           lowerPeriod = 2, #shortes possible period in years
                           make.pval = TRUE, n.sim = 100)

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)", plot.ridge = T,
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#synchrony intensifies from 2008 onwards (again that 3-4 year period!)
#arrows to the right signify that these are in phase. 
#arrows up signify that x is leading
#arrows down signify that y is leading...
#temp leads cases until 219 when they largely converge
my.wc$Angle #give lead (in years of x over y per location)

wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)")

wc.sel.phases(my.wc,    spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))

#and try precip
head(temp.case)
head(precip.nat)
precip.nat <- dplyr::select(precip.nat, -(year), -(month))
temp.case <- merge(temp.case, precip.nat, by="month_date")

my.precipcase <- analyze.coherency(temp.case, my.pair = c("precip_mm","cases"),
                           loess.span = 0,
                           dt = 1/12, dj = 1/100,
                           lowerPeriod = 2, #shortes possible period in years
                           make.pval = TRUE, n.sim = 100)

wc.avg(my.precipcase , siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #not significant

wc.image(my.precipcase, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#arrows to the left suggest that these are out of phase... x follows y
wc.sel.phases(my.precipcase, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))




#3-4 years just later on
precip.temp <- analyze.coherency(temp.case, my.pair = c("temp_C","precip_mm"),
                                   loess.span = 0,
                                   dt = 1/12, dj = 1/100,
                                   lowerPeriod = 2, #shortes possible period in years
                                   make.pval = TRUE, n.sim = 100)


wc.image(precip.temp, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))
#only coherency in precip and temp later in the time series
#out of phase

wc.avg(precip.temp, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #not significant

wc.sel.phases(precip.temp, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                   labels = 2002:2019))


#what about coherency between provinces???
head(dat.nat)
head(dat.prov)
dat.pp = subset(dat.prov, provname=="Phnom Penh")
dat.pp <- dplyr::select(dat.pp, epimonth, cases)
names(dat.pp)[names(dat.pp)=="cases"] <- "cases_PP"
dat.prov.merge <- merge(dat.prov, dat.pp, by="epimonth")
dat.prov.merge <- arrange(dat.prov.merge, provname, cases)


Kandal.match <- analyze.coherency(subset(dat.prov.merge, provname=="Kandal"), my.pair = c("cases_PP", "cases"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 lowerPeriod = 2, #shortes possible period in years
                                 make.pval = TRUE, n.sim = 100)

wc.image(Kandal.match, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))
#only coherency later in the time series: x leading y (PP leading kandal)


Takeo.match <- analyze.coherency(subset(dat.prov.merge, provname=="Takeo"), my.pair = c("cases_PP", "cases"),
                                  loess.span = 0,
                                  dt = 1/12, dj = 1/100,
                                  lowerPeriod = 2, #shortes possible period in years
                                  make.pval = TRUE, n.sim = 100)

wc.image(Takeo.match, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))



KP.match <- analyze.coherency(subset(dat.prov.merge, provname=="Kampong Speu"), my.pair = c("cases_PP", "cases"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 lowerPeriod = 2, #shortes possible period in years
                                 make.pval = TRUE, n.sim = 100)

wc.image(KP.match, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))


wc.sel.phases(KP.match, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                 labels = 2002:2019))


#what about temp in other regions
dat.pp <- subset(dat.prov, provname=="Phnom Penh")
temp.pp <- subset(temp.dat, adm1_name=="Phnom Penh")
precip.pp <- subset(precip.dat, adm1_name=="Phnom Penh")
head(precip.pp)
head(temp.pp)
temp.pp <- dplyr::select(temp.pp, month_date, temp_C)
precip.pp <- dplyr::select(precip.pp, month_date, precip_mm)
names(temp.pp)[1] <- "epimonth"
names(precip.pp)[1] <- "epimonth"

pp.temp <- merge(dat.pp, temp.pp, by="epimonth")
pp.temp <- merge(pp.temp, precip.pp, by="epimonth")

tmp.pp <- analyze.coherency(pp.temp, my.pair = c("temp_C", "cases"),
                              loess.span = 0,
                              dt = 1/12, dj = 1/100,
                              lowerPeriod = 1/12, #shortes possible period in years
                              make.pval = TRUE, n.sim = 100)


wc.image(tmp.pp, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))


wc.sel.phases(tmp.pp, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))


precip.cor.pp <- analyze.coherency(pp.temp, my.pair = c("precip_mm", "cases"),
                            loess.span = 0,
                            dt = 1/12, dj = 1/100,
                            lowerPeriod = 1/12, #shortes possible period in years
                            make.pval = TRUE, n.sim = 100)


wc.image(precip.cor.pp, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))


wc.sel.phases(precip.cor.pp, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                            labels = 2002:2019))




#what about temp in other regions
dat.tk <- subset(dat.prov, provname=="Takeo")
temp.tk <- subset(temp.dat, adm1_name=="Takeo")
precip.tk <- subset(precip.dat, adm1_name=="Takeo")
head(temp.tk)
temp.tk <- dplyr::select(temp.tk, month_date, temp_C)
precip.tk <- dplyr::select(precip.tk, month_date, precip_mm)
names(temp.tk)[1] <- "epimonth"
names(precip.tk)[1] <- "epimonth"

tk.temp <- merge(dat.tk, temp.tk, by="epimonth")
tk.temp <- merge(tk.temp, precip.tk, by="epimonth")

tmp.tk <- analyze.coherency(tk.temp, my.pair = c("temp_C", "cases"),
                            loess.span = 0,
                            dt = 1/12, dj = 1/100,
                            lowerPeriod = 1/12, #shortes possible period in years
                            make.pval = TRUE, n.sim = 100)



wc.image(tmp.tk, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#here, temp lines up only early in time series
wc.sel.phases(tmp.tk, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                            labels = 2002:2019))


precip.cor.tk <- analyze.coherency(tk.temp, my.pair = c("precip_mm", "cases"),
                            loess.span = 0,
                            dt = 1/12, dj = 1/100,
                            lowerPeriod = 1/12, #shortes possible period in years
                            make.pval = TRUE, n.sim = 100)

wc.image(precip.cor.tk, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))





#what about temp in other regions
dat.tk <- subset(dat.prov, provname=="Siem Reap")
temp.tk <- subset(temp.dat, adm1_name=="Siemreap")
precip.tk <- subset(precip.dat, adm1_name=="Siemreap")
head(temp.tk)
temp.tk <- dplyr::select(temp.tk, month_date, temp_C)
precip.tk <- dplyr::select(precip.tk, month_date, precip_mm)
names(temp.tk)[1] <- "epimonth"
names(precip.tk)[1] <- "epimonth"

tk.temp <- merge(dat.tk, temp.tk, by="epimonth")
tk.temp <- merge(tk.temp, precip.tk, by="epimonth")

tmp.tk <- analyze.coherency(tk.temp, my.pair = c("temp_C", "cases"),
                            loess.span = 0,
                            dt = 1/12, dj = 1/100,
                            lowerPeriod = 2, #shortes possible period in years
                            make.pval = TRUE, n.sim = 100)



wc.image(tmp.tk, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#here, temp lines up only early in time series
wc.sel.phases(tmp.tk, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                            labels = 2002:2019))


precip.cor.tk <- analyze.coherency(tk.temp, my.pair = c("precip_mm", "cases"),
                                   loess.span = 0,
                                   dt = 1/12, dj = 1/100,
                                   lowerPeriod = 2, #shortes possible period in years
                                   make.pval = TRUE, n.sim = 100)

wc.image(precip.cor.tk, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#and load oni and mei
oni <- read.csv(file = paste0(homewd, "/data/oni_long.csv"), header=T, stringsAsFactors = F)
head(oni)
oni$month_date <- as.Date(oni$month_date)
ggplot(oni) + geom_point(aes(month_date, oni))

#and mei
mei <- read.csv(file = paste0(homewd, "/data/mei_long.csv"), header=T, stringsAsFactors = F)
head(mei)
mei$month_date <- as.Date(mei$month_date)
ggplot(mei) + geom_point(aes(month_date, mei))
