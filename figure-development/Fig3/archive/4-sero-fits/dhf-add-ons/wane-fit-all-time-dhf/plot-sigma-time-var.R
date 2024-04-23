

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)


# at the province level, compare each run with and without age data
# and with and without waning immunity

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/4-sero-fits/dhf-add-ons/wane-fit-all-time-dhf/"))
load("sigma-fit-all-time-var.Rdata")


p1 <- ggplot(sigma.fit) + geom_line(aes(x=year, y=sigma)) + 
                          #geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma), alpha=.3) + 
                          geom_point(aes(x=year, y=sigma)) #+ scale_y_log10()
p1

p2 <- ggplot(sigma.fit) + geom_line(aes(x=year, y=dur_immunity)) + 
#  geom_ribbon(aes(x=year, ymin=lci_dur_imm, ymax=uci_dur_imm), alpha=.3) + 
  geom_point(aes(x=year, y=dur_immunity)) + scale_y_log10()
p2




p3 <- ggplot(subset(sigma.fit, year>2018)) + geom_line(aes(x=year, y=sigma)) + 
#  geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma), alpha=.3) + 
  geom_point(aes(x=year, y=sigma)) 
p3

p4 <- ggplot(subset(sigma.fit, year>2018))+ geom_line(aes(x=year, y=dur_immunity)) + 
#  geom_ribbon(aes(x=year, ymin=lci_dur_imm, ymax=uci_dur_imm), alpha=.3) + 
  geom_point(aes(x=year, y=dur_immunity)) + scale_y_log10()
p4
