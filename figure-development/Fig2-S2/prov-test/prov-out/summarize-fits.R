rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig2-S2/prov-test/prov-out/"))

#load and combine the fits and plot them together (except not the ones with Inf)
load("fit-prov-Banteay-Meanchey.Rdata")
fit.dat <- out

#load("fit-prov-Battambang.Rdata")
#fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Speu.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampot.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kep.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kratie.Rdata")
fit.dat <- rbind(fit.dat, out)


#load("fit-prov-Kandal.Rdata")
#fit.dat <- rbind(fit.dat, out)

load("fit-prov-Koh-Kong.Rdata")
fit.dat <- rbind(fit.dat, out)


#load("fit-prov-Otdar-Meanchey.Rdata")
#fit.dat <- rbind(fit.dat, out)

load("fit-prov-Phnom-Penh.Rdata")
fit.dat <- rbind(fit.dat, out)

#load("fit-prov-Preah-Sihanouk.Rdata")
#fit.dat <- rbind(fit.dat, out)

load("fit-prov-Prey-Veng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Siem-Reap.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Stung-Treng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Takeo.Rdata")
fit.dat <- rbind(fit.dat, out)


head(fit.dat)
fit.dat <- subset(fit.dat, lambda<10)
 
p1 <- ggplot (fit.dat) + geom_point(aes(x=year, y=lambda, color=provname)) +
      geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + 
      theme(panel.grid = element_blank()) + coord_cartesian(ylim = c(0,1)) +
      geom_vline(xintercept = 2007, linetype=2) +
      geom_vline(xintercept = 2012, linetype=2) +
      geom_vline(xintercept = 2019, linetype=2)

p2 <- ggplot(subset(fit.dat, year>=1995)) + geom_point(aes(x=year, y=lambda, color=provname)) +
  geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + 
  theme(panel.grid = element_blank()) + 
  geom_vline(xintercept = 2007, linetype=2) +
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2)


head(fit.dat)

#now, take the mean of these and submit it to feed the others
mean.fit <- ddply(fit.dat, .(year), summarise, mean_lambda = mean(lambda), median_lambda = median(lambda))
mean.fit

p3 <- ggplot(mean.fit) + geom_point(aes(x=year, y=mean_lambda)) +
  geom_point(aes(x=year, y=median_lambda), color="red") +
  geom_line(aes(x=year, y=mean_lambda)) +
  geom_line(aes(x=year, y=median_lambda), color="red") +
  theme_bw() + 
  theme(panel.grid = element_blank())  +
  geom_vline(xintercept = 2007, linetype=2) +
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2)

