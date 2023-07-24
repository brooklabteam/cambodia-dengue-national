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

#and sum cases by biweek
dat.prov <- ddply(dat, .(provname, dates), summarise,  cases = sum(case) )
head(dat.prov)
unique(dat.prov$provname)

dat.prov$dates <- as.Date(dat.prov$dates)

#and plot 
ggplot(dat.prov) + geom_line(aes(x=dates, y=cases, color=provname), show.legend = F) + facet_wrap(~provname)



#looks good! 
#now load climate data to combine


#now look for mult-annual oscillations in temperature and precip
#load in temperature data from Katie
katwd = "/Users/carabrook/Developer/cambodia-dengue-province"
temp.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )


# #get mean temp across study period for precip and temp
# temp.mean <- ddply(temp.dat, .(adm1_name), summarise, mean_temp_C = mean(temp_C), median_temp_C = median(temp_C), quarter_quant_temp_C=quantile(temp_C)[2], threequarter_quant_temp=quantile(temp_C)[4])
# #and precip
# precip.mean <- ddply(precip.dat, .(adm1_name), summarise, mean_precip_mm = mean(precip_mm), median_precip_mm = median(precip_mm), quarter_quant_precip_mm=quantile(precip_mm)[2], threequarter_quant_precip_mm=quantile(precip_mm)[4])
# 
# 
# #and merge
# climate.mean <-  temp.mean
# climate.mean <- merge(climate.mean, precip.mean, by="adm1_name")
# climate.mean <- climate.mean[2:nrow(climate.mean),]
# 
# #save 
# #write.csv(climate.mean, file = paste0(homewd, "/data/climate_mean_dat.csv"), row.names = F)
# head(temp.dat)
# temp.dat <- dplyr::select(temp.dat, -(X))

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

ggplot(data=precip.dat) + geom_line(aes(x=dates, y=precip_mm,color=adm1_name), show.legend = F) +
  facet_wrap(~adm1_name, ncol=5) + coord_cartesian(xlim=c(as.Date("2002-01-01"), as.Date("2020-12-31")))

# #2020 and 2021 are messed up!
# 
# 
# #rescale 
 precip.dat$year <- as.factor(precip.dat$year)
 ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
 ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_m,color=year)) + facet_wrap(~adm1_name)
#   
# max(precip.dat$dates) #end of year 2021
# 
# ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
 ggplot(data=subset(precip.dat, adm1_name=="Mondul Kiri")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
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
precip.dat$precip_m[precip.dat$dates>=as.Date("2020-01-01")] <- precip.dat$precip_m[precip.dat$dates>=as.Date("2020-01-01")]/100
# 

 ggplot(data=subset(precip.dat, adm1_name=="Phnom Penh")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)
 ggplot(data=subset(precip.dat, adm1_name=="Kep")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name) # o data for 2020 or 2021
 ggplot(data=subset(precip.dat, adm1_name=="Preah Sihanouk")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name) # o data for 2020 or 2021
 ggplot(data=subset(precip.dat, adm1_name=="Pailin")) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name) # o data for 2020 or 2021
 ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)

 #for these 4 provinces, take the average precip across the entire dataset by biweek and replace values with this for 2020 and 2021
 #it does not actually matter for our TSIR which stops in 2019 anyway!
avg.phnom.penh <- ddply(subset(precip.dat, adm1_name=="Phnom Penh" & precip.dat$dates<as.Date("2020-01-01")), .(biwk), summarise, precip_m=mean(precip_m), precip_mm=mean(precip_mm))
avg.kep <- ddply(subset(precip.dat, adm1_name=="Kep" & precip.dat$dates<as.Date("2020-01-01")), .(biwk), summarise, precip_m=mean(precip_m), precip_mm=mean(precip_mm))
avg.PS <- ddply(subset(precip.dat, adm1_name=="Preah Sihanouk" & precip.dat$dates<as.Date("2020-01-01")), .(biwk), summarise, precip_m=mean(precip_m), precip_mm=mean(precip_mm))
avg.pailin <- ddply(subset(precip.dat, adm1_name=="Pailin" & precip.dat$dates<as.Date("2020-01-01")), .(biwk), summarise, precip_m=mean(precip_m), precip_mm=mean(precip_mm))

precip.dat <- arrange(precip.dat, adm1_name, year, biwk)

precip.dat$precip_mm[precip.dat$year==2020 & precip.dat$adm1_name=="Phnom Penh"] <- avg.phnom.penh$precip_mm
precip.dat$precip_mm[precip.dat$year==2021 & precip.dat$adm1_name=="Phnom Penh"] <- avg.phnom.penh$precip_mm
precip.dat$precip_m[precip.dat$year==2020 & precip.dat$adm1_name=="Phnom Penh"] <- avg.phnom.penh$precip_m
precip.dat$precip_m[precip.dat$year==2021 & precip.dat$adm1_name=="Phnom Penh"] <- avg.phnom.penh$precip_m


precip.dat$precip_mm[precip.dat$year==2020 & precip.dat$adm1_name=="Kep"] <- avg.kep$precip_mm
precip.dat$precip_mm[precip.dat$year==2021 & precip.dat$adm1_name=="Kep"] <- avg.kep$precip_mm
precip.dat$precip_m[precip.dat$year==2020 & precip.dat$adm1_name=="Kep"] <- avg.kep$precip_m
precip.dat$precip_m[precip.dat$year==2021 & precip.dat$adm1_name=="Kep"] <- avg.kep$precip_m

precip.dat$precip_mm[precip.dat$year==2020 & precip.dat$adm1_name=="Preah Sihanouk"] <- avg.PS$precip_mm
precip.dat$precip_mm[precip.dat$year==2021 & precip.dat$adm1_name=="Preah Sihanouk"] <- avg.PS$precip_mm
precip.dat$precip_m[precip.dat$year==2020 & precip.dat$adm1_name=="Preah Sihanouk"] <- avg.PS$precip_m
precip.dat$precip_m[precip.dat$year==2021 & precip.dat$adm1_name=="Preah Sihanouk"] <- avg.PS$precip_m

precip.dat$precip_mm[precip.dat$year==2020 & precip.dat$adm1_name=="Pailin"] <- avg.pailin$precip_mm
precip.dat$precip_mm[precip.dat$year==2021 & precip.dat$adm1_name=="Pailin"] <- avg.pailin$precip_mm
precip.dat$precip_m[precip.dat$year==2020 & precip.dat$adm1_name=="Pailin"] <- avg.pailin$precip_m
precip.dat$precip_m[precip.dat$year==2021 & precip.dat$adm1_name=="Pailin"] <- avg.pailin$precip_m


# and plot 
ggplot(data=precip.dat) + geom_line(aes(x=biwk, y=precip_mm,color=year)) + facet_wrap(~adm1_name)

#okay, precip is corrected - 
#now, get national average/sum

#now merge with case data by biweek
head(dat.prov)
head(precip.dat)
head(temp.dat)
names(temp.dat)[names(temp.dat)=="adm1_name"] <- "provname"
names(precip.dat)[names(precip.dat)=="adm1_name"] <- "provname"

setdiff(unique(temp.dat$provname), unique(dat.prov$provname))
setdiff(unique(dat.prov$provname), unique(temp.dat$provname))

setdiff(unique(precip.dat$provname), unique(dat.prov$provname))
setdiff(unique(dat.prov$provname), unique(precip.dat$provname))


temp.dat$provname[temp.dat$provname=="Siemreap"] <- "Siem Reap"
precip.dat$provname[precip.dat$provname=="Siemreap"] <- "Siem Reap"

temp.dat$provname[temp.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"
precip.dat$provname[precip.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"

temp.dat = subset(temp.dat, provname!="Administrative unit not available")
precip.dat = subset(precip.dat, provname!="Administrative unit not available")

temp.merge <- dplyr::select(temp.dat, -(year), -(biwk))

all.dat.prov  <- merge(temp.merge, dat.prov, by = c("provname", "dates"))
head(all.dat.prov)

all.dat.prov  <- merge(all.dat.prov, precip.dat, by = c("provname", "dates"))
head(all.dat.prov)


write.csv(all.dat.prov, file = paste0(homewd, "/data/climate_cases_prov.csv"), row.names = F)


# save this climate/cases dataset
# estimate beta from the TSIR 
# find province-level lags
# also calculate national lag
# do panel regression of beta predicted by climate variables and 
# use to predict 2019 data


#now calculate lags - need to split by province first

# Look for cross-correlations
dat.lag <- cbind.data.frame(lag = print(ccf(all.dat.prov$precip_mm, all.dat.prov$cases))$lag, acf=print(ccf(all.dat.prov$precip_mm, all.dat.prov$cases))$acf)
dat.lag$variable <- "precip_mm"
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)] #cases preced ppt by 1 biweek??? Sum ppt precedes cases by 1months
dat2 = cbind.data.frame(lag = print(ccf(all.dat.prov$temp_C, all.dat.prov$cases))$lag, acf=print(ccf(all.dat.prov$temp_C, all.dat.prov$cases))$acf)
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
##Temp lagged 3 months, Precipitation lagged 11 months
merge.shift <- merge.dat[12: length(merge.dat$epimonth),]
merge.shift$precip_lag <- merge.dat$sumppt[1:(length(merge.dat$sumppt)-11)]
merge.shift$meantempLag <- merge.dat$meant[9:(length(merge.dat$meant)-3)]



#then, do same for national
#then, calc lags, 
#then do at national level too...
#then fit tsir

