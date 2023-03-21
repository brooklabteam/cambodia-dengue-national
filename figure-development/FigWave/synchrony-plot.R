rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(WaveletComp)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)


dat <- read.csv(file = paste0(homewd, "/data/all_case_climate.csv"), header=T, stringsAsFactors = F)

head(dat) 
dat$month_date <- as.Date(dat$month_date)


#first, for each time series per province, run a function to collect:
#(a) the reconstructed period for each timestep, both annual
#(b) and multi-annual
#(c) the average wavelet power per timestep, for annual
#(d) and multi-annual
#(e) the proportion of other provinces with which it shares a statistically significant coherency
#(f) the average wavelet coherency with ONI
#(g) the average wavelet coherency with temperature for that province
#(h) the average wavelet cohernecy with precipitation for that province