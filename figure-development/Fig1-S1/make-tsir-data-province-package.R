rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)
library(tsiR)

#make national population and birth data

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

#feed to tsir data
# time = time series of timesteps by 1 week
# cases = series of weekly cases (same length as time)
# births = annual time series of births
# population = annual time series of population

#load data 
#and get national pop data
pop.dat <- read.csv(file=paste0(homewd, "data/world_bank_cambodia.csv"), header = T, stringsAsFactors = F)
head(pop.dat)

#get population vector
pop.vec <- pop.dat[2,5:ncol(pop.dat)]
names(pop.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
#so one timestep before gives you the population at the beginning of the year
pop.vec <- c(unlist(pop.vec[which(names(pop.vec)=="2001"):which(names(pop.vec)=="2020")]))
#pop.vec[length(pop.vec)] <- pop.vec[length(pop.vec)-1]

#do the same for births - these are births per 1000 people
#get total births
birth.vec <- pop.dat[1,5:ncol(pop.dat)]
names(birth.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="2002"):which(names(birth.vec)=="2020")]
birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year

#now scale up by population size to get total births per year
birth.vec = c(unlist(birth.vec*(pop.vec/1000)))
#time = c(2002:2021)

#and cases

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
head(dat) 

#plot time series of each type by province by year
unique(dat$diagnostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$epiwk <- cut.Date(dat$date, breaks="weeks", start.on.monday = T)
dat$epiwk <- as.Date(as.character(dat$epiwk))
dat$epiwk[dat$epiwk< "2002-01-01"] <- "2002-01-01"
head(dat)

#now sum by week and split by province 
dat.sum <- ddply(dat,.(provname,epiwk), summarise, cases=sum(case))
head(dat.sum)

#and split by the province level
dat.split <- dlply(dat.sum, .(provname))

# now write a function that multiplies pop and births by province proportion,
# creates a time column, and generates tsir data by province

load(paste0(homewd, "/data/cambodia_province_proportions.Rdata"))
head(prop.prov)

build.prov.tsir <- function(prov.case.dat, nat.births, nat.pop, prop.prov){
  #first, scale the pop and birth vectors by prov
  prov.proportion = prop.prov$pop_prop[prop.prov==unique(prov.case.dat$provname)]
  prop.pop <- nat.pop*prov.proportion
  prop.births <- nat.births*prov.proportion
  
  #add time column
  prov.case.dat$time <- year(prov.case.dat$epiwk) + yday(prov.case.dat$epiwk)/365
  
  #and use tsir
  
  out.prov <- tsiRdata(time = prov.case.dat$time, cases = prov.case.dat$cases, births = prop.births, pop = prop.births)
  out.prov$provname <- unique(prov.case.dat$provname)
  
  return(out.prov)
  
  
}


# and build tSIR dataset by province
tsir.prov <- lapply(dat.split, build.prov.tsir, nat.births=birth.vec, nat.pop=pop.vec, prop.prov=prop.prov)


#and save - no climate, but this is fine
tsir.prov <- data.table::rbindlist(tsir.prov)

head(tsir.prov)


write.csv(tsir.prov, file = paste0(homewd, "/data/tsir_dat_province.csv"), row.names = F)



