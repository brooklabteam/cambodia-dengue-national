rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(geosphere)

#calculate pearson's cross correlation coefficient 
#between time series of two provinces. First, do it annually with the annual data
#then do it in overlapping 5-year windows for the multi-annual data

#first, load the time series data
homewd= "/Users/carabrook/Developer/cambodia-dengue-national"

#load tsir data
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
tail(tsir.dat) #goes from beginning of 2002 to end of 2020

#link to centorid data so we can also calculate geographic
#distance between provinces


#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)
centroid.merge <- dplyr::select(centroid.prov, -(mean_elevation_m))

setdiff(unique(tsir.dat$provname), unique(centroid.merge$provname))
setdiff(unique(centroid.merge$provname), unique(tsir.dat$provname))
tsir.dat$provname[tsir.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

#merge
tsir.dat <- merge(tsir.dat, centroid.merge, by = "provname")
tsir.dat$cases[tsir.dat$time<2016 & tsir.dat$provname=="Tboung Khmum"] <- NA
tsir.dat <- tsir.dat[complete.cases(tsir.dat),]

cross.corr <- function(df){
  out.df <- cbind.data.frame(year=unique(df$year), corr=cor(df$cases_this_prov, df$cases_other_prov))
  return(out.df)
  
}
assess.corr.annual <- function(df2, df1){
  #set aside prior data
  df1.hold = df1
  df2.hold = df2
  
  #first, join the two series for the full length that they match up
  df1 <- dplyr::select(df1, time, year, cases)
  df2 <- dplyr::select(df2, time, year, cases)
  names(df1) <- c("time", "year", "cases_this_prov")
  names(df2) <- c("time", "year", "cases_other_prov")
  df.join <- merge(df1,df2, by = c("time", "year"))
  
  #now get the full time series correlation
  full.corr <- cor(df.join$cases_this_prov, df.join$cases_other_prov)
  
  #then, split both by year and look within
  df.year.split <- dlply(df.join, .(year))
  
  df.year.out <- lapply(df.year.split, cross.corr)
  df.year.out <- data.table::rbindlist(df.year.out)
  #and assess cross corr within each year
  
  # now, add in the identifying details
  df.year.out$provname = unique(df1.hold$provname)
  df.year.out$comp_prov <- unique(df2.hold$provname)
  df.year.out$full_ts_corr <- full.corr
  #and add in the phylogenetic distance between these two
  df.year.out$dist_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(unique(df2.hold$longitude), unique(df2.hold$latitude)), fun = distHaversine)
  df.year.out$dist_km <- df.year.out$dist_m/1000
  df.year.out$dist_from_PP_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(104.8397, 11.58035), fun = distHaversine)
  df.year.out$dist_from_PP_km <- df.year.out$dist_from_PP_m/1000
  
  return(df.year.out)
}
assess.cross.corr <- function(provname1, df.all){
  df.other <- subset(df.all, provname!=provname1)
  
  #now split other by province 
  df.other <- arrange(df.other, provname, time)
  df.split <- dlply(df.other, .(provname))
  
  df.now = subset(df.all, provname==provname1)
  df.now <- arrange(df.now, provname, time)
  
  #and apply other function across this list of provinces to compare
  out.provs.list <- lapply(df.split, assess.corr.annual, df1 = df.now)
  out.provs.df <- data.table::rbindlist(out.provs.list)
  
  #ggplot(out.provs.df) + geom_line(aes(x=year, y = corr, color = comp_prov)) + geom_line(aes(x=year, y=full_ts_corr, color=comp_prov),size=1) + facet_wrap(~comp_prov)
   return(out.provs.df) 
}

#now, split the provinces and run for all
provname.list <- as.list(unique(tsir.dat$provname))

out.pearsons <- lapply(provname.list, assess.cross.corr, df.all=tsir.dat)
pearsons.df <- data.table::rbindlist(out.pearsons)
head(pearsons.df)

#save these data
write.csv(pearsons.df, file = paste0(homewd, "/data/pearsons_correlations_provinces.csv"),row.names = F)
