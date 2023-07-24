rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig2-S2/refit-new/refit-out/"))

#load and combine the fits and plot them together (except not the ones with Inf)

load("fit-prov-Banteay-Meanchey.Rdata")
fit.dat <- out

load("fit-prov-Battambang.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Cham.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Kampong-Chhnang.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Speu.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Thom.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampot.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kandal.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kep.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Koh-Kong.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Kratie.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Mondul-Kiri.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Otdar-Meanchey.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Pailin.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Phnom-Penh.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Preah-Sihanouk.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Preah-Vihear.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Prey-Veng.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Pursat.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Ratanak-Kiri.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Siem-Reap.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Stung-Treng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Svay-Rieng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Takeo.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Tboung-Khmum.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-national.Rdata")
fit.dat <- rbind(fit.dat, out)


head(fit.dat)
#fit.dat <- subset(fit.dat, lambda<10)

p1 <- ggplot(fit.dat) + geom_point(aes(x=year, y=lambda, color=provname)) +
      geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + #facet_wrap(~provname, scales = "free_y", ncol=2)+
      geom_line(data=subset(fit.dat, provname=="national"), aes(x=year, y=lambda), linewidth=1, color="black") +
      theme(panel.grid = element_blank()) + #coord_cartesian(ylim = c(0,1)) +
      geom_vline(xintercept = 2007, linetype=2) +
      geom_vline(xintercept = 2012, linetype=2) +
      geom_vline(xintercept = 2019, linetype=2)

p2 <-  ggplot(fit.dat) + geom_point(aes(x=year, y=lambda, color=provname)) +
  geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + facet_wrap(~provname, ncol=4)+
  theme(panel.grid = element_blank()) + #coord_cartesian(ylim = c(0,1)) +
  geom_vline(xintercept = 2007, linetype=2) +
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2)


head(fit.dat)


min.all <- ddply(fit.dat, .(provname), summarise, min_year = min(year))
sort(unique(min.all$min_year)) #1966 forTboung Khmum and 1981 after that

p3 <-  ggplot(fit.dat) + geom_point(aes(x=year, y=lambda, color=provname)) +
  geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + facet_wrap(~provname, ncol=4)+
  theme(panel.grid = element_blank()) + coord_cartesian(xlim = c(1981,2020)) +
  geom_vline(xintercept = 2007, linetype=2) +
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2)

# 
# p4 <-  ggplot(fit.dat) + geom_point(aes(x=year, y=lambda, color=provname)) +
#   geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + facet_wrap(~provname, ncol=4)+
#   theme(panel.grid = element_blank()) + coord_cartesian(xlim = c(2002,2020), ylim=c(0,0.5)) +
#   geom_vline(xintercept = 2007, linetype=2) +
#   geom_vline(xintercept = 2012, linetype=2) +
#   geom_vline(xintercept = 2019, linetype=2)

unique(fit.dat$provname[fit.dat$convergence==1]) #still lots that did not converge!

head(fit.dat)
write.csv(fit.dat, file = "fitted_province_par.csv", row.names = F)
#was the refit any different than the original fit???

rm(list=ls())

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig2-S2/refit-new/refit-out/"))


fit.dat.new <- read.csv(file = "fitted_province_par.csv", header = T, stringsAsFactors = F)
head(fit.dat.new)
unique(fit.dat.new$provname[fit.dat.new$convergence==0])

#now load the old...
fit.dat.old <- read.csv(file = paste0(homewd, "/figure-development/Fig2-S2/fit-new/fit-out/fitted_province_par_original.csv"), header = T, stringsAsFactors = F)
head(fit.dat.old)
unique(fit.dat.old$provname[fit.dat.old$convergence==0])

fit.dat.old <-dplyr::select(fit.dat.old, year, lambda, provname)
names(fit.dat.old)[names(fit.dat.old)=="lambda"]<- "lambda_old"

fit.dat <- merge(fit.dat.new, fit.dat.old, by=c("year", "provname"))
head(fit.dat)
fit.dat$lambda <- round(fit.dat$lambda, 3)
fit.dat$lambda_old <- round(fit.dat$lambda_old, 3)
fit.dat$lambda_diff = fit.dat$lambda_old-fit.dat$lambda
unique(fit.dat$lambda_diff) #0 = all the same


