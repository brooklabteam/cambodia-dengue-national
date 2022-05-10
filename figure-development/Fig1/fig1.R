# This file is for figure 1. 

######################################################################################
#### Fig A | TSIR with increased S,time vs case ######################################
######################################################################################


rm(list = ls())


library(ggplot2)
library(reshape2)
library(pracma)
library(astsa)
library(sp)
library(GADMTools)
library(plotrix)
library(tsiR)
library(lubridate)
library(plm)
library(interp)
library(foreign)
library(gplots)
library(lfe)
library(zoo)
library(e1071)
library(xts)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(graphics)
library(stargazer)
library(grid)
library(gridExtra)
library(png)
library(kernlab)
require("knitr")

# homewd = "/home/mae/0BrookLab/TSIR/0main_r_tsir/6_fig4/data"
 homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
#homewd <- "/home/rstudio"
# setwd(homewd)



# set up size and color
xy_size<-13
my_red="tomato"
my_blue<-"#2980B9"
my_orange<-"#D35400"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
color_2007<-"purple"
color_2012<-"turquoise3"
color_2019<-"tomato"


#epidemic years 2007, 2012, 2019

my_blue<-"#2980B9"
my_orange<-"#D35400"
my_red<-"tomato"
my_green<-"seagreen"
my_yellow<-"#F1C40F"





sim.com.epi <- function(time.start, dat, epiyr, family){
  #
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    ## The type of regression used in susceptible reconstruction. 
    # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
    
    ## The family in the GLM regression. 
    # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
    
  } 
  #saveRDS(fittedpars, paste0(homewd,"/sim.com.epi_estpars_",epiyr,".RDS") )
  #saveRDS(simfitted, paste0(homewd,"/sim.com.epi_",epiyr,".RDS") )
  
  
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  fracincvec<-seq(1,3,0.01)
  #fracincdist<-matrix(0,1,1)
  ratioImax<-matrix(0,1,length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  
  for(i in 1:length(fracincvec)){
    
    Sbegin=simfitted$sbar+simfitted$Z
    Sepi<-Sbegin[length(Sbegin)]
    SepiInc<-fracincvec[i]*Sbegin[length(Sbegin)]
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time<epiyr)
    
    predict_ts <- predicttsir(times=dat.fit$time,
                              births = dat.fit$births,
                              beta = fittedpars$beta,
                              alpha = fittedpars$alpha,
                              S0 = Sbegin[1],
                              I0 = dat.fit$cases[1],
                              nsim=100,
                              stochastic = T)
    
    
    
    #and the epi year-no increase
    dat.pred= subset(dat, time >= epiyr & time<=(epiyr+1))
    predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                     births = dat.pred$births,
                                     beta = fittedpars$beta,
                                     alpha = fittedpars$alpha,
                                     S0 = Sepi,
                                     I0 = Ifinal,
                                     nsim=100,
                                     stochastic = T)
    
    #saveRDS(predict_epi_noInc, paste0(homewd,"/sim.com.epi_predict_epi_noInc_",epiyr,".RDS") )
    
    #and including an increase
    predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = SepiInc,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
    
    
    
    ##saveRDS(predict_epi_Inc, paste0(homewd,"/sim.com.epi_predict_epi_Inc_",epiyr,".RDS") )
    
    IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
    ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
    
    Iall <- rbind(IPredEpi, ISimTS)
    Iall$variable <- as.character(Iall$variable)
    Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
    Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
    
    all.dat <- subset(dat, time<=(epiyr+1))
    all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
    ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
      geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
      geom_line(aes(x=time, y=value, color=variable)) 
    
    ggplot(data=Iall) + 
      geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
      geom_line(aes(x=time, y=value, color=variable)) 
    
    
    
    ratioImax[i]<-max(Iall$value[Iall$variable=="prediction TSIR-increased S"])/max(all.dat$cases)
  }
  spot<-which(abs(ratioImax-1)==min(abs(ratioImax-1)))
  fracincout <- fracincvec[spot]
  ratioImaxout<-ratioImax[spot]
  S_noinc<-Sepi
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr, epiSnoic= S_noinc, frac_incS = fracincout)
  
  return(dat.out)
}

sim.com.noS <- function(time.start, dat, epiyr, family){
  
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
  }
  
  
  
  
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time<=(epiyr+1))
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_noInc$I$time, mean_noInc=predict_epi_noInc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  #Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable)) 
  
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases/mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall)=="value"] <- "cases"
  Iall$reported_cases <- Iall$cases/mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data=all.out) + 
    geom_line(aes(x=time, y=cases, color=variable)) 
  
  return(all.out)
  
  
  
}

#here, run fitted version with and without increased S
sim.with.increaseS <- function(fracincS,time.start, dat, epiyr, family){
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if (family=="gaussian"){
    # 
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    # 
    # #now simulate tsir
    
    # 
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois', 
                              nsim=100)
  }
  #saveRDS(fittedpars, paste0(homewd,"/sim.with.increaseS_estpars_",epiyr,".RDS") )
  #saveRDS(simfitted, paste0(homewd,"/sim.with.increaseS_simfitted_",epiyr,".RDS") )
  #
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  #and sim with and without increased S
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  
  Sepi<-Sbegin[length(Sbegin)]
  SepiInc<-fracincS*Sepi
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time< (epiyr+1))
  
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  #saveRDS(predict_epi_noInc, paste0(homewd,"/sim.with.increaseS_predict_epi_noInc_",epiyr,".RDS") )
  
  #and including an increase
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #saveRDS(predict_epi_Inc, paste0(homewd,"/sim.with.increaseS_predict_epi_Inc_",epiyr,".RDS") )
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable)) 
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases/mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall)=="value"] <- "cases"
  Iall$reported_cases <- Iall$cases/mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data=all.out) + 
    geom_line(aes(x=time, y=cases, color=variable)) 
  
  return(all.out)
  
}
#here, plot the output
plot.comp <- function(dat, filename){
  dat$variable = factor(dat$variable, levels=c("data",
                                               "TSIR fit",
                                               "prediction TSIR",
                                               "prediction TSIR-increased S"))
  colz = c('data' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR-increased S' =my_red)
  typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1)
  p1 <- ggplot(data=dat) + 
    geom_line(aes(x=time, y=reported_cases, color=variable, linetype=variable), size=1) +
    scale_color_manual(values=colz) +
    scale_linetype_manual(values=typez) +theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=18),
          legend.text = element_text(size=10),
          legend.position = c(.25,.82),
          axis.text = element_text(size=14), legend.title = element_blank()) +
    ylab("reported cases")
  
  print(p1)
  
  
}




get.data.tsir <- function(dat, births, N, natl, prop.prov) {
  #make a vector of weeks, years, months
  week <- week(dat$date)
  week[week == 53] <- 52 #now have a vector of 52 weeks in a year
  year <- year(dat$date) #here we do the same for years
  month <- month(dat$date) #and here the same for months
  month.index <- rep(NA, 52) #here make a vector of NA 52 times
  for (j in 1:52) {
    month.index[j] <-
      mean(month[week == j])
  }#here fill in that vector to assign the 52 weeks in the year to a given month (1-12)
  month.index[1] <- 1
  month.index[52] <- 12
  
  
  u.yr <-
    unique(year)
  u.yr <-
    u.yr[order(u.yr)]  #here make a vector of the unique years in increasing order over which the time series spans
  yr.index <-
    rep(u.yr, each = 52) #here make the same vector that has a year assigned to each of those 52 entries. you now have 2 vectors, for the weeks and years the dataset spans
  #make vector of number of positive cases by biweek
  #across the entire  time series for dengue (bwn)
  
  wks <- 1:52
  biwks <- rep(1:26, each = 2)
  
  dat$biwk <- dat$week
  for (i in 1:52) {
    dat$biwk[dat$biwk == wks[i]] <- biwks[i]
  }
  
  #summarize case data by biweek
  date.sum <-
    ddply(dat, .(year, biwk), summarise, cases = length(year))
  #head(date.sum)
  index.dat <-
    cbind.data.frame(biwk = rep(1:26, length(u.yr)),
                     year = rep(u.yr, each = 26))
  index.dat <-
    merge(
      index.dat,
      date.sum,
      by = c("year", "biwk"),
      all.x = T,
      sort = F
    )
  index.dat$cases[is.na(index.dat$cases)] <- 0
  
  #vector of cases by biweek across the timeseries
  bwn <- index.dat$cases
  yrs <-
    rep(u.yr, each = 26) #make a vector that is the same length and assigns a year to each of these 2 week periods
  weeks <- rep(1:26, length(u.yr)) #do the same for weeks
  months <-
    round(rep(0.5 * (month.index[seq(1, 52, by = 2)] + month.index[seq(1, 52, by =
                                                                         2) + 1]), length(u.yr))) #do the same for months
  
  
  ## get at births in biweeks annually across the time series
  ##input vector taken from World Bank databank which reports births per 1000 ppl
  ## this is at the national level, so convert based on total population
  #now spread across each biweek of each year
  
  #now scale births (here as births per 1000 ppl across the time series) to population,
  #where pop is either a vector the same length as births or a single integer
  
  if (natl == TRUE) {
    births <- births * (N / 1000)
  } else{
    #first, scale N by the province proportion
    N = N * (prop.prov$pop_prop[prop.prov$provname == unique(dat$provname)])
    births <- births * (N / 1000)
  }
  
  
  
  births <- rep(births, each = 26) / 26
  #now have a timeseries of births (per 1000 ppl), and you need population to scale the births
  
  #and make a vector of population to match the length of the timeseries
  pop = rep(N, each = 26)
  
  
  ## plot it out for dengue- this is raw reported, not estimated number of cases
  #plot(rep(u.yr,each=26)+rep(1:26,length(u.yr))/26,bwn, type="l", xlab="", ylab="Reported dengue incidence", xlim=c(2002,2021))
  #abline(v=2000:2020, col="grey")
  
  # list all the components of your TSIR
  #list(bwn=bwn,births=births, yrs=yrs,weeks=weeks, months=months)
  #then, prep data for return
  if (natl == TRUE) {
    out.dat <-
      cbind.data.frame(
        yr = yrs,
        month = months,
        week = weeks,
        cases = bwn,
        births = births,
        population = pop
      )
    out.dat$province <- "national"
  } else{
    out.dat <-
      cbind.data.frame(
        yr = yrs,
        month = months,
        week = weeks,
        cases = bwn,
        births = births,
        population = N,
        province = unique(dat$provname)
      )
    
    
  }
  ## return things needed for TSIR
  return(out.dat)
  
  
}





# tsir_data.csv is generated by make-tsir-data.R
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)




# provname<-unique(tsir.dat.prepare$province)
# print(provname)
# tsir.dat <- tsir.dat.prepare %>% select(time,cases,births,pop,births_per_1000)

sim.2007 <- sim.com.epi(
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)




#now return and simulate with results
comp.2007 <- sim.with.increaseS(
  fracincS = sim.2007$frac_incS,
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)
# #and plot
# plot.comp(dat=comp.2007, filename = paste0(homewd,
#                                            "/Fig4-2007-fit-gaussian.png"))


#now for 2012
sim.2012 <- sim.com.epi(
  dat = tsir.dat,
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)


#now return and simulate with results
comp.2012 <- sim.with.increaseS(
  fracincS = sim.2012$frac_incS,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2008]) & time < 2013),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)
#and plot
# plot.comp(dat=comp.2012, filename = paste0(homewd,
#                                            "/Fig4-2012-fit-gaussian.png"))


#and the results for 2019
sim.2019 <- sim.com.epi(
  dat = tsir.dat,
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2013]),
  family = "gaussian",
  epiyr = 2019
)


#now return and simulate with results
comp.2019 <- sim.with.increaseS(
  fracincS = sim.2019$frac_incS,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2013]) & time < 2020),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2013]),
  family = "gaussian",
  epiyr = 2019
)
# #and plot
# plot.comp(dat=comp.2019, filename = paste0(homewd,
#                                            "/Fig4-2019-fit-gaussian.png"))


#now plot the whole thing
all.comp <- rbind(comp.2007, comp.2012, comp.2019)

#and plot all
dat=all.comp
dat$variable = factor(dat$variable, levels=c("data",
                                             "TSIR fit",
                                             "prediction TSIR",
                                             "prediction TSIR-increased S"))
colz = c('data' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR-increased S' =my_red)
typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1)


pl_a <- ggplot(data=dat, aes(x=time, y=reported_cases, color=variable, linetype=variable)) + 
  geom_line(data=dplyr::filter(dat, variable=="data"), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time<2007), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2008 & time<2012), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2019 & time<2020), size=1) +
  scale_color_manual(values=colz) + # coord_cartesian(ylim = c(0,1500))+
  scale_linetype_manual(values=typez) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=10),
        legend.position = c(.1,.82),
        axis.text = element_text(size=14), legend.title = element_blank()) +
  ylab(paste0("reported cases"))


plot(pl_a)








######################################################################################
#### Fig B | log ######################################
######################################################################################



sim.fit.return.beta <- function(time.start, dat, epiyr, family){
  
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                                                 IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                                                 regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                                                         IP = 2,
                                                         parms=fittedpars,
                                                         #epidemics='break', threshold=3,
                                                         #method='pois',
                                                         nsim=100)
                           
    
  }
  
   
   
   df.out = fittedpars$contact
   df.out$epiyr <- epiyr
   
  
  
  return(df.out)
}

#load data into tsir form
# tsir_data.csv is generated by make-tsir-data.R
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

#first try to predict 2007-epidemic
#here we fit to 2002:2007

## @Yimei, I know you made some changes to the fitting script for TSIR
## Below assumes regtype = "lm" and family = "gaussian".
## Please make sure it matches

sim.2007 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)


#now for 2012
sim.2012 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                        family="gaussian",
                        epiyr = 2012)


#and the results for 2019
sim.2019 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                        family="gaussian",
                        epiyr = 2019)

#and combine
fit.all <- rbind(sim.2007, sim.2012, sim.2019)
fit.all$epiyr[fit.all$epiyr==2007] <- "2002-2006"
fit.all$epiyr[fit.all$epiyr==2012] <- "2008-2011"
fit.all$epiyr[fit.all$epiyr==2019] <- "2013-2018"
fit.all$epiyr <- as.factor(fit.all$epiyr)


p1 <- ggplot(data=fit.all) +#facet_grid(epiyr~., scales = "free_y") +
      geom_line(aes(x=time, y=beta, color=epiyr)) +
      geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)




#or as lines and points
pl_b <- ggplot(data=fit.all) + #facet_grid(epiyr~., scales = "free_y") +
  geom_point(aes(x=time, y=beta, color=epiyr)) +
  geom_linerange(aes(x=time, ymin=betalow, ymax=betahigh, color=epiyr), alpha=.3)+
  geom_line(aes(x=time, y=beta, color=epiyr)) +
  #geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)+
  theme_bw() + scale_y_log10()+
  scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                     name="fitting period") +
  scale_fill_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                       name="fitting period") + 
  ylab(bquote(log[10](beta)~',transmission rate')) +
  #xlab("week of year")+
  scale_x_continuous(breaks=c(1*2,3*2,5*2,7*2,9*2,11*2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  theme(panel.grid = element_blank(), 
        legend.position = c(.128,.885),
        legend.text = element_text(size=14),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.text = element_text(size=14), legend.title = element_blank(),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))









######################################################################################
#### Fig C | time & sus & color birth rate ######################################
######################################################################################



# Sdatcombined is generated by "tsir-increase-return-S-sdatcombined.R"
Sdatcombined <- readRDS(paste0(homewd, "/data/Sdatcombined.rds"))
colz = c("other years"="black", "epidemic years"="red")
Sdatcombined



Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
    
}
# Sdatcombined




################################################################
# Fig 3C, time vs sus, birth change to birth/1000 ppl
################################################################


# this is from the original birth data with added  births_per_1000 column for each year. (The birth rate is actually declining while the population size is increasing.)
find_births<-read.csv(paste0(homewd,"/data/tsir_data_birth_updated.csv"))


# choose a number every 26 rows
births_update<- cbind.data.frame(seq(2002,2020,1), sapply(split(find_births$births_per_1000, rep(1:(nrow(find_births)/26), each=26)), mean))
colnames(births_update)<-c("time","births_per1000")
births_update


Sdatcombined_update<-merge(x = Sdatcombined, y =births_update, by = "time", all.x = TRUE)

Sdatcombined = subset(Sdatcombined_update, select = -c(births) )
names(Sdatcombined)[names(Sdatcombined) == 'births_per1000'] <- 'births per 1000 ppl'
Sdatcombined



colz = c("other years"="black", "epidemic years"="red")

Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
    
}

# remove TSIR and only show TSIR-increased S

Sdatcombined<-Sdatcombined[!(Sdatcombined$sim %in% "TSIR"),]



Sdatcombined$variable<-sub("-", "- \n ",Sdatcombined$variable)  


a<-c("black","red","red","black","red")

my_blue<-"#2980B9"
my_green<-"green3"
colz = c("prediction TSIR"=my_green, "prediction TSIR- \n increased S"="tomato","TSIR fit"=my_blue)
options(scipen=10000) # we don't want to plot with 1e+00


pl_c<-ggplot(data=Sdatcombined) + 
  geom_point(aes(x=year, y=sus_mean, 
                 color=variable, shape=variable), size=5, stroke = 1) + labs(shape=NULL, colour=NULL) +# only remove two legend
  scale_color_manual(values=colz) +

  ggnewscale::new_scale_color() +
  geom_point(aes(x=year, y=sus_mean, 
                 color=`births per 1000 ppl`, shape=variable), size=3) + 
  
  theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white")) +
  labs(x ="Year", y = "Susceptible Count")+
  theme(
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14))+
  scale_y_log10() + scale_color_viridis_c()+ geom_vline(xintercept = c(2007,2012,2019) , linetype=2, 
                                                                                        color = "red", size=0.5)+
  scale_x_continuous(limits=c(2002, 2020), breaks=c(2002,2007,2012,2015,2019))+
  theme(axis.text.x = element_text( hjust = 0.5, colour = a),strip.text = element_text(colour = 'black',size = 14))  +
  
  theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),legend.title=element_text(size=14),
        legend.position = c(.847,.723),legend.box.background = element_rect(colour = "black"),legend.background = element_blank(),
        axis.text = element_text(size=14))+# rremove("xy.title")+# theme(legend.title = element_blank()) + 
  theme(legend.spacing.y = unit(0, "cm"))+ # space between legend rows
  theme(legend.title.align = 1.08)




######################################################################################
#### Fig D | reported case & predicted case TSIR & increasedS ########################
######################################################################################




#load data into tsir form

#homewd<-"/home/rstudio"
# tsir_data.csv is generated by make-tsir-data.R
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)




# this part is for increased S

sim.com.epi <- function(time.start, dat, epiyr, family) {
  if (family == "poisson") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'poisson',
        link = 'log'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        method = 'pois',
        nsim = 100
      )
    
  } else if (family == "gaussian") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'gaussian'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        #method='pois',
        nsim = 100
      )
    
    
  }
  
  
  
  
  
  
  Its = simfitted$res$mean * fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
  Ifinal <- Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the
  #closest match of model to observed cases
  
  fracincvec <- seq(1, 3, 0.01)
  #fracincdist<-matrix(0,1,1)
  ratioImax <- matrix(0, 1, length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  
  for (i in 1:length(fracincvec)) {
    Sbegin = simfitted$sbar + simfitted$Z
    Sepi <- Sbegin[length(Sbegin)]
    SepiInc <- fracincvec[i] * Sbegin[length(Sbegin)]
    
    
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time < epiyr)
    
    predict_ts <- predicttsir(
      times = dat.fit$time,
      births = dat.fit$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sbegin[1],
      I0 = dat.fit$cases[1],
      nsim = 100,
      stochastic = T
    )
    
    
    
    #and the epi year-no increase
    dat.pred = subset(dat, time >= epiyr & time <= (epiyr + 1))
    predict_epi_noInc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sepi,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    #and including an increase
    predict_epi_Inc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = SepiInc,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    # print("predict_epi_Inc:")
    # print(predict_epi_Inc)
    
    
    IPredEpi = melt(
      cbind.data.frame(
        time = predict_epi_Inc$I$time,
        mean_noInc = predict_epi_noInc$I$mean,
        mean_Inc = predict_epi_Inc$I$mean
      ),
      id.vars = "time"
    )
    ISimTS = cbind.data.frame(
      time = predict_ts$I$time,
      variable = rep("TSIR fit", length(predict_ts$I$time)),
      value = predict_ts$I$mean
    )
    
    # print("IPredEpi:")
    # print(IPredEpi)
    # print("ISimTS:")
    # print(ISimTS)
    
    Iall <- rbind(IPredEpi, ISimTS)
    Iall$variable <- as.character(Iall$variable)
    Iall$variable[Iall$variable == "mean_Inc"] <-
      "prediction TSIR-increased S"
    Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
    
    all.dat <- subset(dat, time <= (epiyr + 1))
    all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
    ggplot(data = subset(Iall, variable != "prediction TSIR-increased S")) +
      geom_line(data = all.dat,
                aes(x = time, y = cases),
                linetype = 2) +
      geom_line(aes(x = time, y = value, color = variable))
    
    ggplot(data = Iall) +
      geom_line(data = all.dat,
                aes(x = time, y = cases),
                linetype = 2) +
      geom_line(aes(x = time, y = value, color = variable))
    
    
    
    ratioImax[i] <-
      max(Iall$value[Iall$variable == "prediction TSIR-increased S"]) / max(all.dat$cases)
  }
  
  
  spot <- which(abs(ratioImax - 1) == min(abs(ratioImax - 1)))
  fracincout <- fracincvec[spot]
  ratioImaxout <- ratioImax[spot]
  S_noinc <- Sepi
  
  dat.out <-
    cbind.data.frame(epidemic_year = epiyr,
                     epiSnoic = S_noinc,
                     frac_incS = fracincout)
  
  
  return(dat.out)
}


# no increase part & plot

sim.com.noS <- function(time.start, dat, epiyr, family) {
  if (family == "poisson") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'poisson',
        link = 'log'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        method = 'pois',
        nsim = 100
      )
    
  } else if (family == "gaussian") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'gaussian'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        #method='pois',
        nsim = 100
      )
    
    
  }
  
  
  
  
  
  
  Its = simfitted$res$mean * fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
  Ifinal <- Its[length(Its)]
  
  
  
  Sbegin = simfitted$sbar + simfitted$Z
  Sepi <- Sbegin[length(Sbegin)]
  
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time < epiyr)
  
  predict_ts <- predicttsir(
    times = dat.fit$time,
    births = dat.fit$births,
    beta = fittedpars$beta,
    alpha = fittedpars$alpha,
    S0 = Sbegin[1],
    I0 = dat.fit$cases[1],
    nsim = 100,
    stochastic = T
  )
  
  
  
  #and the epi year-no increase
  dat.pred = subset(dat, time >= epiyr & time <= (epiyr + 1))
  predict_epi_noInc <- predicttsir(
    times = dat.pred$time,
    births = dat.pred$births,
    beta = fittedpars$beta,
    alpha = fittedpars$alpha,
    S0 = Sepi,
    I0 = Ifinal,
    nsim = 100,
    stochastic = T
  )
  
  
  
  IPredEpi = melt(
    cbind.data.frame(time = predict_epi_noInc$I$time, mean_noInc = predict_epi_noInc$I$mean),
    id.vars = "time"
  )
  ISimTS = cbind.data.frame(
    time = predict_ts$I$time,
    variable = rep("TSIR fit", length(predict_ts$I$time)),
    value = predict_ts$I$mean
  )
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  #Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time <= (epiyr + 1))
  all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) +
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable))
  
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases / mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall) == "value"] <- "cases"
  Iall$reported_cases <- Iall$cases / mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data = all.out) +
    geom_line(aes(x = time, y = cases, color = variable))
  
  return(all.out)
  
  
  
}
#here, run fitted version with and without increased S
sim.with.increaseS <-
  function(fracincS, time.start, dat, epiyr, family) {
    if (family == "poisson") {
      fittedpars <-
        estpars(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          alpha = NULL,
          sbar = NULL,
          xreg = "cumcases",
          regtype = 'lm',
          family = 'poisson',
          link = 'log'
        )
      
      simfitted <-
        simulatetsir(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          parms = fittedpars,
          #epidemics='break', threshold=3,
          method = 'pois',
          nsim = 100
        )
      
    } else if (family == "gaussian") {
      #
      fittedpars <-
        estpars(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          alpha = NULL,
          sbar = NULL,
          xreg = "cumcases",
          regtype = 'lm',
          family = 'gaussian'
        )
      #
      # #now simulate tsir
      
      #
      simfitted <-
        simulatetsir(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          parms = fittedpars,
          #epidemics='break', threshold=3,
          #method='pois',
          nsim = 100
        )
    }
    #
    Its = simfitted$res$mean * fittedpars$rho #account for underreporting
    dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
    Ifinal <- Its[length(Its)]
    
    #and sim with and without increased S
    
    
    #now try to predict the epidemic both with and without the increased S
    #compute the most likely model fit with the proper increase which gives the
    #closest match of model to observed cases
    
    
    Sbegin = simfitted$sbar + simfitted$Z
    Sepi <- Sbegin[length(Sbegin)]
    SepiInc <- fracincS * Sepi
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time < epiyr)
    
    predict_ts <- predicttsir(
      times = dat.fit$time,
      births = dat.fit$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sbegin[1],
      I0 = dat.fit$cases[1],
      nsim = 100,
      stochastic = T
    )
    
    
    
    #and the epi year-no increase
    dat.pred = subset(dat, time >= epiyr & time < (epiyr + 1))
    predict_epi_noInc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sepi,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    #and including an increase
    predict_epi_Inc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = SepiInc,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    
    
    
    IPredEpi = melt(
      cbind.data.frame(
        time = predict_epi_Inc$I$time,
        mean_noInc = predict_epi_noInc$I$mean,
        mean_Inc = predict_epi_Inc$I$mean
      ),
      id.vars = "time"
    )
    ISimTS = cbind.data.frame(
      time = predict_ts$I$time,
      variable = rep("TSIR fit", length(predict_ts$I$time)),
      value = predict_ts$I$mean
    )
    
    Iall <- rbind(IPredEpi, ISimTS)
    Iall$variable <- as.character(Iall$variable)
    Iall$variable[Iall$variable == "mean_Inc"] <-
      "prediction TSIR-increased S"
    Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
    
    all.dat <- subset(dat, time <= (epiyr + 1))
    all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
    # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) +
    #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    #   geom_line(aes(x=time, y=value, color=variable))
    
    all.dat <- dplyr::select(all.dat, time, cases)
    all.dat$reported_cases <- all.dat$cases / mean(fittedpars$rho)
    all.dat$variable <- "data"
    
    names(Iall)[names(Iall) == "value"] <- "cases"
    Iall$reported_cases <- Iall$cases / mean(fittedpars$rho)
    Iall <- dplyr::select(Iall, names(all.dat))
    
    all.out <- rbind(Iall, all.dat)
    
    
    p2 <- ggplot(data = all.out) +
      geom_line(aes(x = time, y = cases, color = variable))
    
    return(all.out)
    
  }
#here, plot the output
plot.comp <- function(dat, filename) {
  dat$variable = factor(
    dat$variable,
    levels = c(
      "data",
      "TSIR fit",
      "prediction TSIR",
      "prediction TSIR-increased S"
    )
  )
  colz = c(
    'data' = "black",
    'TSIR fit' = "blue",
    'prediction TSIR' = "red",
    'prediction TSIR-increased S' = "seagreen"
  )
  typez = c(
    'data' = 2,
    'TSIR fit' = 1,
    'prediction TSIR' = 1,
    'prediction TSIR-increased S' = 1
  )
  p1 <- ggplot(data = dat) +
    geom_line(aes(
      x = time,
      y = reported_cases,
      color = variable,
      linetype = variable
    ),
    size = 1) +
    scale_color_manual(values = colz) +
    scale_linetype_manual(values = typez) + theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size = 10),
      legend.position = c(.25, .82),
      axis.text = element_text(size = 14),
      legend.title = element_blank()
    ) +
    ylab("reported cases")
  
  print(p1)
  
  ggsave(
    file = filename,
    plot = p1,
    units = "mm",
    width = 50,
    height = 40,
    scale = 3,
    dpi = 300
  )
  
  
}




# gaussian

#######now try as gaussian

#first try to predict 2007-epidemic
sim.2007 <- sim.com.epi(dat = tsir.dat,
                    time.start =  min(tsir.dat$time),
                    family="gaussian",
                    epiyr = 2007)


#now return and simulate with results
comp.2007 <- sim.with.increaseS(fracincS = sim.2007$frac_incS,
                            dat=tsir.dat,
                            time.start =  min(tsir.dat$time),
                            family="gaussian",
                            epiyr = 2007)

#now for 2012
sim.2012 <- sim.com.epi(dat = tsir.dat,
                    time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                    family="gaussian",
                    epiyr = 2012)


#now return and simulate with results
comp.2012 <- sim.with.increaseS(fracincS = sim.2012$frac_incS,
                            dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2008]) & time<2013),
                            time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                            family="gaussian",
                            epiyr = 2012)

#and the results for 2019
sim.2019 <- sim.com.epi(dat = tsir.dat,
                    time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                    family="gaussian",
                    epiyr = 2019)


#now return and simulate with results
comp.2019 <- sim.with.increaseS(fracincS = sim.2019$frac_incS,
                            dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2013]) &time<2020),
                            time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                            family="gaussian",
                            epiyr = 2019)




# prediction TSIR is green
# data is black
# prediction TSIR-increased S is red
# TSIR fit is blue

dat<-rbind(comp.2007,comp.2012,comp.2019)
dat_true<-dat[dat$variable=="data",]

dat_pred_update<-merge(x = dat, y =dat_true, by = "time", all.x = TRUE)

dat_df<-cbind.data.frame(dat_pred_update$time, dat_pred_update$reported_cases.x, dat_pred_update$reported_cases.y, dat_pred_update$variable.x,dat_pred_update$variable.y)

colnames(dat_df)<-c("time","reported_cases","true_cases","pred_type","data")




# make every year different shape


dat_pred_3<-dat_df
dat_pred_3$year<-floor(dat_pred_3$time)

non_pei_years<- as.data.frame(dplyr::filter(dat_pred_3,  (year>=2002 & year <2007) | (year>=2008 & year <2012) | (year>=2013 & year <2019) | (year>2019) ) )


dat_pred_3[dat_pred_3$time %in% non_pei_years$time,]$year<-"other_year"




dat <- dat_pred_3

colz = c(
  "prediction TSIR" = "seagreen",
  "prediction TSIR-increased S" = "tomato",
  "TSIR fit" = my_blue,
  "data" = "black"
)


circle_size <- 3
pl_d <-
  ggplot(data = dat, aes(x = true_cases, y = reported_cases)) + geom_point(
    aes(
      x = true_cases,
      y = reported_cases,
      color = pred_type,
      shape = year
    ),
    size = circle_size,
    stroke = 2
  ) +
  labs(shape = NULL, colour = NULL) +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # aspect.ratio = 1,
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14)
  ) +theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white"))+ 
# +geom_point(aes(x=true_cases, y=reported_cases,
#                fill=year), size=2, stroke = 0.05)
  ylim(0, 9000) + xlim(0, 9000) + scale_color_manual(values = colz) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(.218, .765),
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    axis.text = element_text(size = 14)
  ) + labs(x = "Reported Cases", y = "Predicted Cases")+theme(legend.spacing.y = unit(-0.1, "cm"))# + # space between legend rows
  # theme(legend.title.align = 1.08)














######################################################################################
######################################################################################
######################################################################################
############################# subplots together ######################################
######################################################################################
######################################################################################



pl<-plot_grid(pl_a, NULL, pl_b, pl_c, NULL, pl_d, ncol = 3, nrow = 2,align = "vh",labels = c("a","", "b","c","","d"), label_size =xy_size+20, hjust=-0.5, vjust = c(1.3,1.3,1.4,1.4),rel_widths = c(1, 0.08,1))+ theme(plot.margin = unit(c(1,1,2,1), "cm")) 




ggsave("fig1.jpg",
  plot = pl,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 17,
  height = 12,
  bg='white',
  units = c("in"))









