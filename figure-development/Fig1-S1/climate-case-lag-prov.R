rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)


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

#still some funky stuff in Kep, Pailin, Phnom Penh, Preah Sihanouk

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
ggplot(climate.merge) + geom_line(aes(x=biweek, y=precip_mm, color=year)) + facet_wrap(provname~.)
ggplot(climate.merge) + geom_line(aes(x=biweek, y=temp_C, color=year)) + facet_wrap(provname~.)

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

