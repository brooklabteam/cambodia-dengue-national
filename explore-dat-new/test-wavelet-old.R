rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
  
head(dat) 

#plot time series of each type by province by year

unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date)
dat$epiwk <- cut.Date(dat$date, breaks="weeks", start.on.monday = T)
dat$epiwk <- as.Date(as.character(dat$epiwk))
dat$epiwk[dat$epiwk< "2002-01-01"] <- "2002-01-01"
head(dat)

#now attach biwks (2 biweeks per week)
wks = 1:52
biwks <- rep(1:26, each=2)
wk.biwk <- cbind.data.frame(week_report=wks, biwk=biwks)

sort(unique(dat$week_report))
unique(dat$provname)

dat <- merge(dat, wk.biwk, by="week_report")
head(dat) #now plot as the lowest epiweek date per biweek
dat$case=1


#sum by epidemic biweek, type, province
dat.sum <- ddply(dat, .(year, biwk, provname, dianostic), summarise, cases = sum(case), epiwk = min(epiwk))

head(dat.sum)

#plot infection type by province by year
p1 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=dianostic)) +
      facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2018-01-01"), as.Date("2021-01-01")))
p1


p2 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=dianostic)) +
  facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2011-01-01"), as.Date("2014-01-01")))
p2

p3 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=dianostic)) +
  facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2006-01-01"), as.Date("2009-01-01")))
p3

#total cases by year by province
dat.prov <- ddply(dat, .(year, biwk, provname), summarise, cases = sum(case), epiwk=min(epiwk))
head(dat.prov)



#now, see if you can get the wavelet periodicity by province
library(WaveletComp)

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

