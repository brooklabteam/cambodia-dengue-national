rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)


dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
head(dat) 

#plot time series of each type by province by year
unique(dat$diagnostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
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
#head(dat) #now plot as the lowest epiweek date per biweek

dat$age <- ceiling(dat$age)


#head(dat)

dat <- arrange(dat, procode, date)
class(dat$epimonth)

dat$case=1


#make a time series with date for biweek
sort(unique(dat$biwk))
months <- rep(seq(1,12,1), each=2)
biwks <- seq(1,26,1)
days <- seq(1,365,14)[-1]
dates <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), 1)
dates <- paste(sapply(strsplit(as.character(dates[days]), "-"), "[",2), sapply(strsplit(as.character(dates[days]), "-"), "[",3), sep="-")
years <- unique(dat$year)


#link dates to the data
date.bind <- cbind.data.frame(year = rep(years, each=26), dates = rep(dates, time=length(years)), biwk= rep(biwks, time=length(years)) )
date.bind$dates <- paste0(date.bind$year, "-", date.bind$dates)
head(date.bind)

dat <- merge(dat, date.bind, by=c("year", "biwk"))
head(dat)

#and make a national version
dat.nat <- ddply(dat, .(dates), summarise,  cases = sum(case) )
head(dat.nat)
dat.nat$dates <- as.Date(dat.nat$dates)
ggplot(dat.nat) + geom_line(aes(x=dates, y=cases), show.legend = F) 




#looks good! 
#now load climate data to combine


#now look for mult-annual oscillations in temperature and precip
#load in temperature data from Katie
katwd = "/Users/carabrook/Developer/cambodia-dengue-province"
temp.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )


head(temp.dat)
head(precip.dat)

#compile by year and biweek for national data - then merge with case data
temp.nat <- ddply(temp.dat, .(year, biwk), summarise, temp_C = mean(temp_C))
precip.nat <- ddply(precip.dat, .(year, biwk), summarise, precip_mm = sum(precip_mm) )


#precip in 2020 and 2021 needs to be corrected
precip.nat$year <- as.factor(precip.nat$year)
ggplot(precip.nat) + geom_line(aes(x=biwk, y=precip_mm, color=year, group=year))

precip.nat$precip_mm[precip.nat$year==2020 | precip.nat$year==2021] <- precip.nat$precip_mm[precip.nat$year==2020 | precip.nat$year==2021]/1000

ggplot(precip.nat) + geom_line(aes(x=biwk, y=precip_mm, color=year, group=year))

#now merge
clim.nat <- merge(temp.nat, precip.nat, by =c("year", "biwk"))
head(clim.nat)

clim.nat <- arrange(clim.nat, year, biwk)



write.csv(clim.nat, file = paste0(homewd, "/data/clim_biwk_nat.csv"), row.names = F)
#now merge with case data by biweek


#make a time series with date
unique(temp.dat$biwk)
months <- rep(seq(1,12,1), each=2)
biwks <- seq(1,26,1)
days <- seq(1,365,14)[-1]
days[4] <- days[4]-1
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
head(temp.dat)
max(temp.dat$dates) #end of year 2021



#precip
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )

head(precip.dat)
precip.dat <- dplyr::select(precip.dat, -(X))

#make a time series with date
unique(precip.dat$biwk)
months <- rep(seq(1,12,1), each=2)
biwks <- seq(1,26,1)
days <- seq(1,365,14)[-1]
days[4] <- days[4]-1
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
precip.dat$precip_m <- as.numeric(precip.dat$precip_m)
precip.dat$precip_mm <- as.numeric(precip.dat$precip_mm)


ggplot(data=precip.dat) + 
  geom_line(aes(x=dates, y=precip_mm, color=adm1_name), show.legend = F) +
  facet_wrap(~adm1_name, ncol=5) + 
  coord_cartesian(xlim=c(as.Date("2002-01-01"), as.Date("2020-12-31")))



# #2020 and 2021 are messed up!
# 
# 
# #rescale 
# precip.dat$year <- as.factor(precip.dat$year)
# ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_m,color=year)) + facet_wrap(~adm1_name)
#   
# max(precip.dat$dates) #end of year 2021
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates<as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri"& dates<as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates>=as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri"& dates>=as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates<as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_m,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh" & dates>=as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_m,color=year)) + facet_wrap(~adm1_name)
# 
# 
 precip.dat$precip_mm[precip.dat$dates>=as.Date("2020-01-01")] <- precip.dat$precip_mm[precip.dat$dates>=as.Date("2020-01-01")]/1000
 precip.dat$precip_m[precip.dat$dates>=as.Date("2020-01-01")] <- precip.dat$precip_m[precip.dat$dates>=as.Date("2020-01-01")]/1000
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
# ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)


#okay, precip and temp are corrected - 
#now, get national average/sum

temp.nat <- ddply(temp.dat, .(dates, biwk), summarise, temp_K = mean(temp_K), temp_C=mean(temp_C))
ggplot(temp.nat) + geom_line(aes(x=dates, y=temp_C))

precip.nat <- ddply(precip.dat, .(dates, biwk), summarise, precip_mm = sum(precip_mm), precip_m=sum(precip_m))
ggplot(precip.nat) + geom_line(aes(x=dates, y=precip_mm))

#for now, do a mock-correct of the precip data
head(precip.nat)
precip.nat$year <- year(precip.nat$dates)
precip.nat$year <- as.factor(precip.nat$year )
ggplot(precip.nat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) 

ggplot(data=subset(precip.nat, dates<as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_m,color=year)) 
ggplot(data=subset(precip.nat, dates>=as.Date("2020-01-01"))) + geom_line(aes(x=biwk, y=precip_m,color=year))

# precip.nat$precip_mm[precip.nat$dates>=as.Date("2020-01-01")] <- precip.nat$precip_mm[precip.nat$dates>=as.Date("2020-01-01")]/1000
# precip.nat$precip_m[precip.nat$dates>=as.Date("2020-01-01")] <- precip.nat$precip_m[precip.nat$dates>=as.Date("2020-01-01")]/1000

ggplot(precip.nat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) 
ggplot(precip.nat) + geom_line(aes(x=dates, y=precip_mm))

min(precip.nat$dates) #2001
min(temp.nat$dates) #2001


#now merge climate data by biweek and save
head(precip.nat)
head(temp.nat)

precip.nat <- dplyr::select(precip.nat, dates, precip_mm)
temp.nat <- dplyr::select(temp.nat, dates, biwk, temp_C)

clim.nat <- merge(temp.nat, precip.nat, by ="dates")
head(clim.nat)





write.csv(clim.nat, file = paste0(homewd, "/data/clim_biwk.csv"), row.names = F)
#now merge with case data by biweek

#national first
head(dat.nat)
head(precip.nat)
dat.all.nat <- merge(precip.nat, dat.nat, by="dates")
head(dat.all.nat)
head(temp.nat)
temp.merge <- dplyr::select(temp.nat, -(biwk))
dat.all.nat <- merge(dat.all.nat, temp.merge, by="dates")
head(dat.all.nat)


#plot cases and climate

merge.melt <- melt(dat.all.nat,id.vars=c("dates", "biwk", "year"))
case.dat1 = subset(merge.melt,variable=="cases")
case.dat1$variable <- "temp_C"
case.dat2 = subset(merge.melt,variable=="cases")
case.dat2$variable <- "precip_mm"
case.dat <- rbind(case.dat1,case.dat2)

head(case.dat)



# Plot cases and climate variables
case.dat$value_correct <- case.dat$value/100
ggplot(data=subset(merge.melt, variable=="precip_mm")) + 
  geom_line(data=subset(case.dat, variable=="precip_mm"), aes(x=dates, y=value_correct),  size=1, alpha=.2) +
  geom_point(data=subset(case.dat, variable=="precip_mm"), aes(x=dates, y=value_correct), size=3, alpha=.2) +
  geom_point(aes(x=dates, y=value, color=variable),  size=3, show.legend = F) +
  geom_line(aes(x=dates, y=value, color=variable),  size=1, show.legend = F) +
  facet_grid(variable~., scales = "free") + ylim(c(0, NA)) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     axis.title = element_blank(), 
                     strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"),
                     legend.text = element_text(size=12),
                     plot.margin = unit(c(.2,.1,1.3,1.1), "lines"),
                     axis.text = element_text(size=14))

ggplot(data=subset(merge.melt, variable=="temp_C" )) + 
  geom_line(data=subset(case.dat, variable=="temp_C"), aes(x=dates, y=value_correct),  size=1, alpha=.2) +
  geom_point(data=subset(case.dat, variable=="temp_C"), aes(x=dates, y=value_correct), size=3, alpha=.2) +
  geom_point(aes(x=dates, y=value, color=variable),  size=3, show.legend = F) +
  geom_line(aes(x=dates, y=value, color=variable),  size=1, show.legend = F) +
  facet_grid(variable~., scales = "free") + ylim(c(0, NA)) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     axis.title = element_blank(), 
                     strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"),
                     legend.text = element_text(size=12),
                     plot.margin = unit(c(.2,.1,1.3,1.1), "lines"),
                     axis.text = element_text(size=14))




#now calculate lags - here with national data
# Look for cross-correlations
dat.lag <- cbind.data.frame(lag = print(ccf(dat.all.nat$precip_mm, dat.all.nat$cases))$lag, acf=print(ccf(dat.all.nat$precip_mm, dat.all.nat$cases))$acf)
dat.lag$variable <- "precip_mm"
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)] #cases preced ppt by 1 biweek??? Sum ppt precedes cases by 1months
dat2 = cbind.data.frame(lag = print(ccf(dat.all.nat$temp_C, dat.all.nat$cases))$lag, acf=print(ccf(dat.all.nat$temp_C, dat.all.nat$cases))$acf)
dat2$variable <- "temp_C"
dat2$lag[dat2$acf==max(dat2$acf)] #Mean temp precedes cases by 7 months
dat.lag <- rbind(dat.lag, dat2)

# Plot acf
# include the optimal lag on plot
max.lag <- dlply(dat.lag, .(variable))
get.lag <- function(df){
  lag = df$lag[df$acf==max(df$acf)]
  df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
  return(df.out)
}
max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))
max.lag$label = paste0("lag=", max.lag$lag, "epimonth")


ggplot(dat.lag) + geom_label(data=max.lag, aes(x=18,y=.4, label=label), label.size = 0) +
  geom_bar(aes(x=lag, y=acf), stat = "identity") + ylim(c(NA,.45)) +
  geom_hline(aes(yintercept=0.09), color="blue", linetype=2) +
  geom_hline(aes(yintercept=-0.09), color="blue", linetype=2) +
  facet_grid(variable~.) + theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                                              legend.title = element_blank(),
                                              axis.title = element_text(size=16),
                                              strip.background = element_rect(fill="white"),
                                              strip.text = element_text(size=14),
                                              legend.text = element_text(size=12),
                                              plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                                              axis.text = element_text(size=14))

#now make shifted dataset....


# Make a dataframe with lagged climate variables
##Temp lagged 7 biweeks, Precipitation lagged not lagged at all - or do we do 25 biweeks?
merge.shift <- dat.all.nat[8:length(dat.all.nat$dates),] #starts in year 2
head(merge.shift)
#merge.shift$precip_lag <- dat.all.nat$precip_mm[1:(length(dat.all.nat$precip_mm)-11)] #no shift for precip
merge.shift$meantempLag <- dat.all.nat$temp_C[1:length(dat.all.nat$temp_C[1:(length(dat.all.nat$temp_C)-7)])]

head(merge.shift)  # here is your lagged dataset for regression / tsir

#no longer starts at the beginning of the year. But we can recover those values using the 2001 climate data
merge.start = dat.all.nat[1:7,] 
temp.nat.2001 = subset(temp.nat, dates < as.Date("2002-01-01"))
tail(temp.nat.2001)
temp.merge = temp.nat.2001[20:26,]
merge.start$meantempLag <- temp.merge$temp_C



merge.shift <- rbind(merge.start, merge.shift)

head(merge.shift)

#write data
write.csv(merge.shift, file = paste0(homewd, "/data/lagged-nat-clim.csv"), row.names = F)







