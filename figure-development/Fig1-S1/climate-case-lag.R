rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

#load the province-level transmission
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
head(beta.df)


#also load the climate data at the province level
clim.dat <- read.csv(file = paste0(homewd, "/data/clim_biwk.csv"), header = T, stringsAsFactors = F)
head(clim.dat)

#cut to length of the time series

#ggplot(clim.dat) + 


#join together, and calculate lags