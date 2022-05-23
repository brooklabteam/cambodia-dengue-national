


rm(list=ls())
## tsiR examples

## load the package and dependencies
## kernlab is for the gaussian process 
## the rest are for plotting 
require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
require(plyr)

#epidemic years 2007, 2012, 2019


#homewd = "/home/rstudio/"
homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd))




#load the functions
#here, fit the model in an "endemic" period to predict an epidemic year
#fits with an increased S
sim.com.epi <- function(time.start, dat, epiyr, family){
  
  
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
    
    #and including an increase
    predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = SepiInc,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
    
    
    
    
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
sim.return.S.incS <- function(fracincS,time.start, dat, epiyr, family){
  
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
  
  
  #now, return the data frame of susceptibles vs births
  
  dat.S <- cbind.data.frame(time = dat.fit$time, births=dat.fit$births, sus_mean = predict_ts$S$mean, sus_low= predict_ts$S$low, sus_high=predict_ts$S$high)
  dat.S$variable <- "TSIR fit"
  #head(dat.S)
  
  
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
  
  #and including an increase
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #now combine the S data
  dat.S.noinc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_noInc$S$mean, sus_low= predict_epi_noInc$S$low, sus_high=predict_epi_noInc$S$high)
  dat.S.noinc$variable <- "prediction TSIR"
  
  dat.S.inc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_Inc$S$mean, sus_low= predict_epi_Inc$S$low, sus_high=predict_epi_Inc$S$high)
  dat.S.inc$variable <- "prediction TSIR-increased S"
  
  dat.S.all <- rbind(dat.S, dat.S.noinc, dat.S.inc)
  dat.S.all$year <- trunc(dat.S.all$time)
  dat.S.yr <- dlply(dat.S.all, .(year))
  
  take.min <- function(df){
    df1 <- subset(df, time == min(df$time))
    return(df1)
  }
  
  #then, take only the even year data
  dat.S.yr <- data.table::rbindlist(lapply(dat.S.yr, take.min))
  
  #and return this
  return(dat.S.yr)
}
sim.return.S.noinc <- function(fracincS,time.start, dat, epiyr, family){
  
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
  #SepiInc<-fracincS*Sepi
  
  
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
  
  
  #now, return the data frame of susceptibles vs births
  
  dat.S <- cbind.data.frame(time = dat.fit$time, births=dat.fit$births, sus_mean = predict_ts$S$mean, sus_low= predict_ts$S$low, sus_high=predict_ts$S$high)
  dat.S$variable <- "TSIR fit"
  #head(dat.S)
  
  
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
  # 
  # #and including an increase
  # predict_epi_Inc <- predicttsir(times=dat.pred$time,
  #                                births = dat.pred$births,
  #                                beta = fittedpars$beta,
  #                                alpha = fittedpars$alpha,
  #                                S0 = SepiInc,
  #                                I0 = Ifinal,
  #                                nsim=100,
  #                                stochastic = T)
  # 
  #now combine the S data
  dat.S.noinc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_noInc$S$mean, sus_low= predict_epi_noInc$S$low, sus_high=predict_epi_noInc$S$high)
  dat.S.noinc$variable <- "prediction TSIR"
  # 
  # dat.S.inc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_Inc$S$mean, sus_low= predict_epi_Inc$S$low, sus_high=predict_epi_Inc$S$high)
  # dat.S.inc$variable <- "prediction TSIR-increased S"
  
  dat.S.all <- rbind(dat.S, dat.S.noinc)#, dat.S.inc)
  dat.S.all$year <- trunc(dat.S.all$time)
  dat.S.yr <- dlply(dat.S.all, .(year))
  
  take.min <- function(df){
    df1 <- subset(df, time == min(df$time))
    return(df1)
  }
  
  #then, take only the even year data
  dat.S.yr <- data.table::rbindlist(lapply(dat.S.yr, take.min))
  
  #and return this
  return(dat.S.yr)
}


#here, plot the output
plot.comp <- function(dat, filename){
  dat$variable = factor(dat$variable, levels=c("data",
                                               "TSIR fit",
                                               "prediction TSIR",
                                               "prediction TSIR-increased S"))
  colz = c('data' = "black", 'TSIR fit' = "blue", 'prediction TSIR' = "red", 'prediction TSIR-increased S' = "seagreen")
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
  
  ggsave(file =filename,
         plot=p1,
         units="mm",  
         width=50, 
         height=40, 
         scale=3, 
         dpi=300)
  
  
}




#load data into tsir form
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

#first try to predict 2007-epidemic
sim.2007 <- sim.com.epi(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)


#now return and recover Sdat
Sdat2007 <- sim.return.S.incS(fracincS = sim.2007$frac_incS,
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
Sdat2012 <- sim.return.S.incS(fracincS = sim.2012$frac_incS,
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
Sdat2019 <- sim.return.S.incS(fracincS = sim.2019$frac_incS,
                              dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2013]) &time<2020),
                              time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                              family="gaussian",
                              epiyr = 2019)

#and combine
SdatIncSim <- rbind(Sdat2007, Sdat2012, Sdat2019)

SdatIncSim$sim <- "TSIR-increased-S"


#and without the breaks
###and the whole time series
#and the results for 2019

sim.all <- sim.com.epi(dat = tsir.dat,
                       time.start =  min(tsir.dat$time[tsir.dat$time>=2002]),
                       family="gaussian",
                       epiyr = 2019)

Sdat2019fullsim <- sim.return.S.noinc(fracincS = sim.all$frac_incS,
                                      dat=subset(tsir.dat, time>=2002),
                                      time.start =  min(tsir.dat$time[tsir.dat$time>=2002]),
                                      family="gaussian",
                                      epiyr = 2019)

Sdat2019fullsim$sim <- "TSIR"

#and combine and plot

Sdatcombined <- rbind(SdatIncSim, Sdat2019fullsim)
Sdatcombined$epiyr <- 0
Sdatcombined$epiyr[Sdatcombined$year==2007|Sdatcombined$year==2012|Sdatcombined$year==2019] <- 1
Sdatcombined$epiyr <- factor(Sdatcombined$epiyr)




#and write data file
#write.csv(Sdatcombined, file = paste0(homewd, "/SdatCara.csv"), row.names = F)
saveRDS(Sdatcombined,file = paste0(homewd, "/data/Sdatcombined.rds"))



colz = c('0'="black", '1'="red")




# 
# #and plot these possible figures
# p1 <- ggplot(data=Sdatcombined) + 
#   geom_point(aes(x=births, y=sus_mean, 
#                  color=epiyr, shape=variable), size=5) +
#   scale_color_manual(values=colz) +
#   ggnewscale::new_scale_color() + theme_bw() +
#   theme(panel.grid = element_blank()) +
#   geom_point(aes(x=births, y=sus_mean, 
#                  color=year, shape=variable), size=3) +
#   facet_grid(~sim) + scale_y_log10() + scale_color_viridis_c()
# 
# 
# 
# ggsave(file =paste0(homewd,"/figures/final-figures/SusFig-draft-births.png"),
#        plot=p1,
#        units=c("mm"),  
#        width=80, 
#        height=40, 
#        scale=3, 
#        dpi=300)
# 
# 
# 
# 
# p2 <- ggplot(data=Sdatcombined) + 
#   geom_point(aes(x=year, y=sus_mean, 
#                  color=epiyr, shape=variable), size=5) +
#   scale_color_manual(values=colz) +
#   ggnewscale::new_scale_color() + theme_bw() +
#   theme(panel.grid = element_blank()) +
#   geom_point(aes(x=year, y=sus_mean, 
#                  color=year, shape=variable), size=3) +
#   facet_grid(~sim) + scale_y_log10() + scale_color_viridis_c()
# 
# 
# ggsave(file =paste0(homewd,"/figures/final-figures/SusFig-draft-time.png"),
#        plot=p2,
#        units=c("mm"),  
#        width=80, 
#        height=40, 
#        scale=3, 
#        dpi=300)


