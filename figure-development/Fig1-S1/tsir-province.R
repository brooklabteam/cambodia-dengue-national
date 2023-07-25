rm(list=ls())

require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
require(plyr)

#epidemic years 2007, 2012, 2019


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd))

# generate tsir estimates of the biweekly transmission rate for each 
# province for the three inter-epidemic years
# then, calculate lags for this transmission rate with the climate 
# variables (temp and precip) by province

# then, shift the climate variables by that lag (or the national lag)
# and fit a panel regression model to predict log-beta from climate
# then project beta using climate and try to recover outcomes for the 
# epidemic year. calculate the additional susceptibles needed to recover
# the epidemic year prediction both with and without the role of climate

# first, just fit the tsir data by province


#load the functions
#here, fit the model in an "endemic" period to predict an epidemic year
#fits TSIR
fit.TSIR <- function(dat, epiyr, family){
  #
  print(unique(dat$provname))
  time.start =  min(dat$time)
  
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, 
                          alpha=NULL, 
                          sbar=NULL, 
                          xreg = "cumcases",
                          regtype='lm',
                          family='poisson',
                          link='log')
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
                          IP=2, alpha=NULL,
                          sbar=NULL, xreg = "cumcases",
                          regtype='gaussian',
                          family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
    
  } 
  
  return(simfitted)
}
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
  #ratioImax<-matrix(0,1,length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  sm.sq <- rep(NA, length(fracincvec))
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
    
    #now look at data and model predictions for the epidemic year only
    IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
    names(IPredEpi) <- c("time", "model_predicted_absolute_cases")
    IPredEpi$model_predicted_reported_cases <- IPredEpi$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
    
    #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
    #head(Iall)
    #head(all.dat)
    
    #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
    #and merge on time
    all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
    
    
    #and get sum of sq differences
    all.dat.merge$sq_diff <- (all.dat.merge$cases - all.dat.merge$model_predicted_reported_cases)^2
    
    sm.sq[i] = sum(all.dat.merge$sq_diff)
    
    # ##saveRDS(predict_epi_Inc, paste0(homewd,"/sim.com.epi_predict_epi_Inc_",epiyr,".RDS") )
    # 
    # IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
    # ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
    # 
    # Iall <- rbind(IPredEpi, ISimTS)
    # Iall$variable <- as.character(Iall$variable)
    # Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
    # Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
    # 
    # all.dat <- subset(dat, time<=(epiyr+1))
    # all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
    # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
    #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    #   geom_line(aes(x=time, y=value, color=variable)) 
    # 
    # ggplot(data=Iall) + 
    #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    #   geom_line(aes(x=time, y=value, color=variable)) 
    # 
    # 
    # #capture the increased S that best recovers the peak
    # ratioImax[i]<-max(Iall$value[Iall$variable=="prediction TSIR-increased S"])/max(all.dat$cases)
  }
  #spot<-which(abs(ratioImax-1)==min(abs(ratioImax-1)))
  spot<-which(sm.sq==min(sm.sq))
  fracincout <- fracincvec[spot]
  #ratioImaxout<-ratioImax[spot]
  S_noinc<-Sepi
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr, epiSnoic= S_noinc, frac_incS = fracincout)
  
  return(dat.out)
}
#fits with an increased Beta
sim.com.beta <- function(time.start, dat, epiyr, family){
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
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  fracincvec<-seq(1,5,0.01)
  #fracincdist<-matrix(0,1,1)
  #ratioImax<-matrix(0,1,length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  sm.sq <- rep(NA, length(fracincvec))
  for(i in 1:length(fracincvec)){
    
    Sbegin=simfitted$sbar+simfitted$Z
    Sepi<-Sbegin[length(Sbegin)]
    #SepiInc<-fracincvec[i]*Sbegin[length(Sbegin)]
    
    
    #and look for the beta to multiply
    incBeta <- fracincvec[i]*fittedpars$beta
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
    
    
    
    #and including an increase on the beta
    predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = incBeta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
    
    
    #now look at data and model predictions for the epidemic year only
    IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
    names(IPredEpi) <- c("time", "model_predicted_absolute_cases")
    IPredEpi$model_predicted_reported_cases <- IPredEpi$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]

    #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
    #head(Iall)
    #head(all.dat)
    
    #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
    #and merge on time
    all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
    
    
    #and get sum of sq differences
    all.dat.merge$sq_diff <- (all.dat.merge$cases - all.dat.merge$model_predicted_reported_cases)^2
    
    sm.sq[i] = sum(all.dat.merge$sq_diff)
    
  }
  spot<-which(sm.sq==min(sm.sq))
  fracincout <- fracincvec[spot]
  #ratioImaxout<-ratioImax[spot]
  S_noinc<-Sepi
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr, epiSnoic= S_noinc, frac_incBeta = fracincout)
  
  return(dat.out)
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
sim.with.increaseBeta <- function(fracincBeta,time.start, dat, epiyr, family){
  
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
  

  
  Sbegin=simfitted$sbar+simfitted$Z
  
  Sepi<-Sbegin[length(Sbegin)]
  #SepiInc<-fracincS*Sepi
  BetaInc <- fittedpars$beta*fracincBeta
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
  
  #and including an increase in beta
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = BetaInc,
                                 alpha = fittedpars$alpha,
                                 S0 = Sepi,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #saveRDS(predict_epi_Inc, paste0(homewd,"/sim.with.increaseS_predict_epi_Inc_",epiyr,".RDS") )
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased Beta"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
   # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
   #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
   #   geom_line(aes(x=time, y=value, color=variable)) 
   # 
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
  #betabegin=simfitted$beta
  #betaepi<-betabegin[length(betabegin)]
  #betaepiInc<-fracincS*betaepi
  
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


estpars <- function (data, xreg = "cumcases", IP = 2, seasonality = "standard", 
             regtype = "gaussian", sigmamax = 3, family = "gaussian", 
             link = "identity", userYhat = numeric(), alpha = NULL, sbar = NULL, 
             printon = F) 
  {
    datacheck <- c("time", "cases", "pop", "births")
    if (sum(datacheck %in% names(data)) < length(datacheck)) {
      stop("data frame must contain \"time\", \"cases\", \"pop\", and \"births\" columns")
    }
    na.casescheck <- sum(is.na(data$cases))
    if (na.casescheck > 0) {
      stop("there cannot be any NAs in the cases vector -- please correct")
    }
    na.birthscheck <- sum(is.na(data$births))
    if (na.casescheck > 0) {
      stop("there cannot be any NAs in the births vector -- please correct")
    }
    xregcheck <- c("cumcases", "cumbirths")
    if (xreg %in% xregcheck == F) {
      stop("xreg must be either \"cumcases\" or \"cumbirths\"")
    }
    regtypecheck <- c("gaussian", "lm", "spline", "lowess", "loess", 
                      "user")
    if (regtype %in% regtypecheck == F) {
      stop("regtype must be one of 'gaussian','lm','spline','lowess','loess','user'")
    }
    if (length(sbar) == 1) {
      if (sbar > 1 || sbar < 0) {
        stop("sbar must be a percentage of the population, i.e. between zero and one.")
      }
    }
    linkcheck <- c("log", "identity")
    if (link %in% linkcheck == F) {
      stop("link must be either 'log' or 'identity'")
    }
    seasonalitycheck <- c("standard", "schoolterm", "none")
    if (seasonality %in% seasonalitycheck == F) {
      stop("seasonality must be either 'standard' or 'schoolterm' or 'none'")
    }
    input.alpha <- alpha
    input.sbar <- sbar
    cumbirths <- cumsum(data$births)
    cumcases <- cumsum(data$cases)
    if (xreg == "cumcases") {
      X <- cumcases
      Y <- cumbirths
    }
    if (xreg == "cumbirths") {
      X <- cumbirths
      Y <- cumcases
    }
    x <- seq(X[1], X[length(X)], length = length(X))
    y <- approxfun(X, Y)(x)
    y[1] <- y[2] - (y[3] - y[2])
    if (regtype == "lm") {
      Yhat <- predict(lm(Y ~ X))
    }
    if (regtype == "lowess") {
      Yhat <- lowess(X, Y, f = 2/3, iter = 1)$y
    }
    if (regtype == "loess") {
      Yhat <- predict(loess(y ~ x, se = T, family = "gaussian", 
                            degree = 1, model = T), X)
    }
    if (regtype == "spline") {
      Yhat <- predict(smooth.spline(x, y, df = 2.5), X)$y
    }
    if (regtype == "gaussian") {
      sigvec <- seq(sigmamax, 0, -0.1)
      for (it in 1:length(sigvec)) {
        if (printon == T) {
          print(sprintf("gaussian regression attempt number %d", 
                        it))
        }
        Yhat <- predict(gausspr(x, y, variance.model = T, 
                                fit = T, tol = 1e-07, var = 0.01, kernel = "rbfdot", 
                                kpar = list(sigma = sigvec[it])), X)
        if (sigvec[it] <= min(sigvec)) {
          print("gaussian regressian failed -- switching to loess regression")
          Yhat <- predict(loess(y ~ x, se = T, family = "gaussian", 
                                degree = 1, model = T), X)
        }
        if (xreg == "cumcases") {
          Z <- residual.cases(Yhat, Y)
          rho <- derivative(X, Yhat)
          if (length(which(rho <= 1)) == 0) {
            (break)()
          }
        }
        if (xreg == "cumbirths") {
          rho <- derivative(X, Yhat)
          Z <- residual.births(rho, Yhat, Y)
          if (length(which(rho >= 1)) == 0 && length(which(rho < 
                                                           0)) == 0) {
            (break)()
          }
        }
      }
    }
    if (regtype == "user") {
      Yhat <- userYhat
      if (length(Yhat) == 0) {
        stop("Yhat returns numeric(0) -- make sure to input a userYhat under regtype=user")
      }
    }
    rho <- derivative(X, Yhat)
    if (xreg == "cumcases") {
      Z <- residual.cases(Yhat, Y)
    }
    if (xreg == "cumbirths") {
      Z <- residual.births(rho, Yhat, Y)
    }
    if (xreg == "cumcases") {
      adj.rho <- rho
    }
    if (xreg == "cumbirths") {
      adj.rho <- 1/rho
    }
    if (regtype == "lm") {
      adj.rho <- signif(adj.rho, 3)
    }
    if (length(which(adj.rho < 1)) > 1) {
      warning("Reporting exceeds 100% -- use different regression")
    }
    Iadjusted <- data$cases * adj.rho
    datacopy <- data
    if (seasonality == "standard") {
      period <- rep(1:(52/IP), round(nrow(data) + 1))[1:(nrow(data) - 
                                                           1)]
      if (IP == 1) {
        period <- rep(1:(52/2), each = 2, round(nrow(data) + 
                                                  1))[1:(nrow(data) - 1)]
      }
    }
    if (seasonality == "schoolterm") {
      term <- rep(1, 26)
      term[c(1, 8, 15, 16, 17, 18, 19, 23, 26)] <- 2
      iterm <- round(approx(term, n = 52/IP)$y)
      period <- rep(iterm, round(nrow(data) + 1))[1:(nrow(data) - 
                                                       1)]
    }
    if (seasonality == "none") {
      period <- rep(1, nrow(data) - 1)
      period[nrow(data) - 1] <- 2
    }
    Inew <- tail(Iadjusted, -1) + 1
    lIminus <- log(head(Iadjusted, -1) + 1)
    Zminus <- head(Z, -1)
    pop <- data$pop
    minSmean <- max(0.01 * pop, -(min(Z) - 1))
    Smean <- seq(minSmean,  min(pop), length = 500)
    loglik <- rep(NA, length(Smean))
    if (link == "identity") {
      Inew <- log(Inew)
    }
    if (family %in% c("poisson", "quasipoisson")) {
      Inew <- round(Inew)
    }
    if (length(input.alpha) == 0 && length(input.sbar) == 0) {
      for (i in 1:length(Smean)) {
        lSminus <- log(Smean[i] + Zminus)
        glmfit <- glm(Inew ~ -1 + as.factor(period) + (lIminus) + 
                        offset(lSminus), family = eval(parse(text = family))(link = link))
        loglik[i] <- glmfit$deviance
      }
      sbar <- Smean[which.min(loglik)]
      lSminus <- log(sbar + Zminus)
      glmfit <- glm(Inew ~ -1 + as.factor(period) + (lIminus) + 
                      offset(lSminus), family = eval(parse(text = family))(link = link))
      beta <- exp(head(coef(glmfit), -1))
      alpha <- tail(coef(glmfit), 1)
    }
    if (length(input.alpha) == 1 && length(input.sbar) == 0) {
      for (i in 1:length(Smean)) {
        lSminus <- log(Smean[i] + Zminus)
        glmfit <- glm(Inew ~ -1 + as.factor(period) + offset(alpha * 
                                                               lIminus) + offset(lSminus), family = eval(parse(text = family))(link = link))
        loglik[i] <- glmfit$deviance
      }
      sbar <- Smean[which.min(loglik)]
      lSminus <- log(sbar + Zminus)
      glmfit <- glm(Inew ~ -1 + as.factor(period) + offset(alpha * 
                                                             lIminus) + offset(lSminus), family = eval(parse(text = family))(link = link))
      beta <- exp(coef(glmfit))
    }
    if (length(input.alpha) == 0 && length(input.sbar) == 1) {
      sbar <- sbar * mean(pop)
      lSminus <- log(sbar + Zminus)
      glmfit <- glm(Inew ~ -1 + as.factor(period) + (lIminus) + 
                      offset(lSminus), family = eval(parse(text = family))(link = link))
      beta <- exp(head(coef(glmfit), -1))
      alpha <- tail(coef(glmfit), 1)
    }
    if (length(input.alpha) == 1 && length(input.sbar) == 1) {
      sbar <- sbar * mean(pop)
      lSminus <- log(sbar + Zminus)
      glmfit <- glm(Inew ~ -1 + as.factor(period) + offset(alpha * 
                                                             lIminus) + offset(lSminus), family = eval(parse(text = family))(link = link))
      beta <- exp(coef(glmfit))
    }
    if (seasonality == "none") {
      beta[2] <- beta[1]
      beta <- mean(beta)
      period <- rep(1, nrow(data) - 1)
    }
    confinterval <- suppressMessages(confint(glmfit))
    continterval <- confinterval[1:length(unique(period)), ]
    betalow <- exp(confinterval[, 1])
    betahigh <- exp(confinterval[, 2])
    glmAIC <- AIC(glmfit)
    contact <- as.data.frame(cbind(time = seq(1, length(beta[period]), 
                                              1), betalow[period], beta[period], betahigh[period]), 
                             row.names = F)
    names(contact) <- c("time", "betalow", "beta", "betahigh")
    contact <- head(contact, 52/IP)
    return(list(X = X, Y = Y, Yhat = Yhat, Smean = Smean, contact = contact, 
                period = period, IP = IP, beta = beta, rho = adj.rho, 
                Z = Z, pop = pop, time = data$time, AIC = glmAIC, sbar = sbar, 
                alpha = alpha, loglik = loglik))
  }





###################################################################################
###################################################################################
###################################################################################
###################################################################################


#load province data into tsir form
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

length(unique(tsir.dat$provname)) #25 unique provinces



# # #remove 8 provinces with few cases and sporadic dynamics
# tsir.dat.split <- subset(tsir.dat, provname!="Kep")
# tsir.dat.split <- subset(tsir.dat.split, provname!="Banteay Meanchey")
# tsir.dat.split <- subset(tsir.dat.split, provname!= "Mondul Kiri")
# tsir.dat.split <- subset(tsir.dat.split, provname!="Pailin")
# tsir.dat.split <- subset(tsir.dat.split, provname!="Koh Kong")
# tsir.dat.split <- subset(tsir.dat.split, provname!="Ratanak Kiri")
# tsir.dat.split <- subset(tsir.dat.split, provname!="Stung Treng")
 tsir.dat <- subset(tsir.dat, provname!="Tboung Khmum")
# 
#split by province and pre-epidemic period
tsir.split.2007 <- dlply(tsir.dat,.(provname))

tsir.dat.2012 <- subset(tsir.dat, time > 2007.9999)
tsir.split.2012 <- dlply(tsir.dat.2012,.(provname))

tsir.dat.2019 <- subset(tsir.dat, time > 2012.9999)
tsir.split.2019 <- dlply(tsir.dat.2019,.(provname))



plot.test.tsir <- function(df, epiyr, sbar1){ 
  
  #get the suffix to all the filenames
  
  suffix = unique(df$provname)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  print(suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  
  
  
  p1 <- plotdata(dplyr::select(dat, -(provname), -(biweek), -(year)))
  p1 <- p1 + theme(plot.margin = unit(c(1,.2,.2,.2), "cm"))
  # 
  # ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-data-plots/data_", epiyr,"_", suffix, ".png"),
  #        units="mm",  
  #        width=100, 
  #        height=110, 
  #        scale=3, 
  #        dpi=300)
  
  if(is.null(sbar1)){
    
    
    
    fittedpars1 <- estpars(data=dat,
                           IP=2, 
                           alpha=NULL, 
                           sbar=NULL, 
                           xreg = "cumcases",
                           regtype='lm',
                           family='poisson',
                           link='log')
    
    
    p2 <- plotsbar(fittedpars1)
    
    # ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/lm_", epiyr,"_", suffix, ".png"),
    #        units="mm",  
    #        width=100, 
    #        height=110, 
    #        scale=3, 
    #        dpi=300)
    # 
    
    # if ((fittedpars1$sbar + fittedpars1$Z[1])>0){
    #   
    #   
    simfitted1 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars1,
                               #epidemics='break', threshold=3,
                               method='pois', nsim=100)
    #   
    # }
    
    p3 <- plotrho(fittedpars1)
    
    fittedpars2 <- estpars(data=dat,
                           IP=2, 
                           alpha=NULL, 
                           sbar=NULL, 
                           xreg = "cumcases",
                           regtype='gaussian',
                           family='poisson',
                           link='log')
    
    
    p4 <- plotsbar(fittedpars2)
    
    simfitted2 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars2,
                               #epidemics='break', threshold=3,
                               method='pois', nsim=100)
    
    p5 <- plotrho(fittedpars2)
    #   
    
  }else{
    
    
    fittedpars1 <- estpars(data=dat,
                           IP=2, 
                           alpha=NULL, 
                           sbar=sbar1, 
                           xreg = "cumcases",
                           regtype='lm',
                           family='poisson',
                           link='log')
    
    
    p2 <- plotsbar(fittedpars1)
    
    # ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/lm_", epiyr,"_", suffix, ".png"),
    #        units="mm",  
    #        width=100, 
    #        height=110, 
    #        scale=3, 
    #        dpi=300)
    # 
    
    # if ((fittedpars1$sbar + fittedpars1$Z[1])>0){
    #   
    #   
    # simfitted1 <- simulatetsir(data=dat,
    #                            IP = 2,
    #                            parms=fittedpars1,
    #                            #epidemics='break', threshold=3,
    #                            method='pois', nsim=100)
    # #   
    # }
    
    p3 <- plotrho(fittedpars1)
    
    fittedpars2 <- estpars(data=dat,
                           IP=2, 
                           alpha=NULL, 
                           sbar=sbar1, 
                           xreg = "cumcases",
                           regtype='gaussian',
                           family='poisson',
                           link='log')
    
    
    p4 <- plotsbar(fittedpars2)
    
    # simfitted2 <- simulatetsir(data=dat,
    #                            IP = 2,
    #                            parms=fittedpars2,
    #                            #epidemics='break', threshold=3,
    #                            method='pois', nsim=100)
    # 
    p5 <- plotrho(fittedpars2)
  }
  
  
  pB <- cowplot::plot_grid(p2, p3, ncol = 1, nrow=2) + theme(plot.margin = unit(c(1,.2,.2,.2), "cm"))
  pC <- cowplot::plot_grid(p4, p5, ncol = 1, nrow=2) + theme(plot.margin = unit(c(1,.2,.2,.2), "cm"))
  
  out.plot <- cowplot::plot_grid(p1, pB, pC, ncol = 3, nrow=1, 
                                 rel_widths = c(1,1.1,1.1), 
                                 labels = c("A. Data", "B. Linear regression fit", "C. Gaussian regression fit"), 
                                 label_size = 18, label_x = c(0,-.2,-.2)) + theme(plot.background = element_rect(fill="white"))
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/comp-plots/", suffix, "_", epiyr, ".png"),
         units="mm",  
         width=110, 
         height=60, 
         scale=3, 
         dpi=300)
  
  # if ((fittedpars2$sbar + fittedpars2$Z[1])>0){
  #   
  
  # }
  # 
  # ## The type of regression used in susceptible reconstruction. 
  # # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
  # 
  # ## The family in the GLM regression. 
  # # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
  # if ( ((fittedpars1$sbar + fittedpars1$Z[1])>0) & ((fittedpars2$sbar + fittedpars2$Z[1])>0)) {
  #   
  out.df <- cbind.data.frame(provname = c(rep(unique(dat$provname),2)),
                             regtype=c("lm", "gaussian"), 
                             AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(simfitted1$rsquared, simfitted2$rsquared), 
                             Prop_Susceptible = c((signif(fittedpars1$sbar/mean(dat$pop) * 100, 2)), (signif(fittedpars2$sbar/mean(dat$pop) * 100, 2))))
  
  # } else if ( ((fittedpars1$sbar + fittedpars1$Z[1])<0) & ((fittedpars2$sbar + fittedpars2$Z[1])>0)){
  #   
  #   out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(NA, simfitted2$rsquared))
  #   
  # }else if ( ((fittedpars1$sbar + fittedpars1$Z[1])>0) & ((fittedpars2$sbar + fittedpars2$Z[1])<0)){
  #   
  #   out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(simfitted1$rsquared, NA))
  # }
  # out.df$provname <- suffix
  # 
  return(out.df)
  
}



#Run this to test all the time series
fit.2007.plot <- lapply(tsir.split.2007, plot.test.tsir, epiyr = 2007, sbar=NULL)
fit.2007.plot <- data.table::rbindlist(fit.2007.plot)
fit.2007.plot$epiyr = 2007
# choose.reg <- function(df){
#   df1 <- df[df$rsquared ==min(df$rsquared),]
#   return(df1)
# }
# fit.2007.plot <- dlply(fit.2007.plot, .(provname))
# fit.2007.slim <-  lapply(fit.2007.plot, choose.reg)
# fit.2007.slim <- data.table::rbindlist(fit.2007.slim ) #10 lm vs. 7 gaussian

fit.2012.plot <- lapply(tsir.split.2012, plot.test.tsir, epiyr = 2012, sbar=NULL)
fit.2012.plot <- data.table::rbindlist(fit.2012.plot)
fit.2012.plot$epiyr = 2012
# fit.2012.plot <- dlply(fit.2012.plot, .(provname))
# fit.2012.slim <-  lapply(fit.2012.plot, choose.reg)
# fit.2012.slim <- data.table::rbindlist(fit.2012.slim ) #6 lm vs. 11 gaussian

fit.2019.plot <- lapply(tsir.split.2019, plot.test.tsir, epiyr = 2019, sbar=NULL)
fit.2019.plot <- data.table::rbindlist(fit.2019.plot)
fit.2019.plot$epiyr = 2019
# fit.2019.plot <- dlply(fit.2019.plot, .(provname))
# fit.2019.slim <-  lapply(fit.2019.plot, choose.reg)
# fit.2019.slim <- data.table::rbindlist(fit.2019.slim ) #11 lm vs. 6 gaussian

fit.comp <- rbind(fit.2007.plot, fit.2012.plot, fit.2019.plot)

#combine all the fits
# 
# #and write over the few that need help
# plot.test.tsir(tsir.split.2012$Kandal, epiyr = 2012, sbar1 = median(c(.57, .62)))
# plot.test.tsir(df = tsir.split.2019$`Phnom Penh`, epiyr = 2019, sbar1 = median(c(.43, .46)))
# plot.test.tsir(df = tsir.split.2019$`Kampong Thom`, epiyr = 2019, sbar1 = median(c(.51, .59)))
# plot.test.tsir(df = tsir.split.2019$`Kampong Cham`, epiyr = 2019, sbar1 = median(c(.44, .84)))
# plot.test.tsir(df = tsir.split.2019$`Siem Reap`, epiyr = 2019, sbar1 = .36)
# plot.test.tsir(df = tsir.split$`Siem Reap`, epiyr = 2007, sbar1 = .36)
# plot.test.tsir(df = tsir.split$`Preah Vihear`, epiyr = 2007, sbar1 =.45)


runFulltSIR <- function(dat, fit.comp.df, epiyr1, sbar1){
  #first, subselect the comp data to the right province and year
  prov.comp= subset(fit.comp.df, provname==unique(dat$provname) & epiyr==epiyr1)
  
  prov.choose = prov.comp[prov.comp$rsquared==max(prov.comp$rsquared),]
  if(nrow(prov.choose)>1){
    prov.choose = subset(prov.choose, regtype=="lm") #choose the simpler model
   }
  
  if(prov.choose$regtype=="lm"){
    
   out.df <-  fit_tsir(df=dat, epiyr = epiyr1, sbar=sbar1)
   out.df$sus_reconstruction = "lm"
   
  }else if (prov.choose$regtype=="gaussian"){
    
   out.df <-  fit_tsir_gaussian(df=dat, epiyr = epiyr1, sbar=sbar1)  
   out.df$sus_reconstruction = "gaussian"
  }
  
  
  
  return(out.df)
}

plotres <- function (dat, max.plot = 10) {
  if (is.null(dat$SIRS) == TRUE) {
    dat$SIRS <- FALSE
  }
  if (dat$SIRS == TRUE) {
    ll.melt <- dat$ll.melt
    p1 <- ggplot(ll.melt, aes_string(x = "X1", y = "X2", 
                                     z = "loglik")) + geom_tile(aes_string(fill = "loglik")) + 
      scale_fill_gradient(low = "white", high = "red") + 
      theme_bw() + geom_contour(col = "black") + geom_point(aes(x = dat$k, 
                                                                y = dat$m), col = "black") + xlab("strength of immunity") + 
      ylab("duration of immunity") + ggtitle(sprintf("duration = %g, strength = %g", 
                                                     round(dat$m), signif(dat$k, 2)))
    update.ll <- ll.melt[, which(!apply(ll.melt == 0, 2, 
                                        all))]
    if (nrow(update.ll) * ncol(update.ll) != nrow(ll.melt) * 
        ncol(ll.melt)) {
      lldf <- NULL
      lldf$loglik <- update.ll$loglik
      lldf$Var <- update.ll[names(update.ll)[which(names(update.ll) != 
                                                     "loglik")]]
      lldf <- as.data.frame(lldf)
      names(lldf) <- c("loglik", "Var")
      p1 <- ggplot(lldf, aes_string("Var", "loglik")) + 
        geom_line(size = 2) + theme_bw() + ggtitle(sprintf("duration = %g, strength = %g", 
                                                           round(dat$m), signif(dat$k, 2)))
    }
    Sdf <- NULL
    Sdf$time <- head(dat$res$time, -1)
    Sdf$S <- dat$S
    Sdf <- as.data.frame(Sdf)
    p2 <- ggplot(Sdf, aes_string("time", "S")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(rho) == .(signif(mean(1/dat$rho), 
                                                       2))))
    p3 <- ggplot(data = dat$contact, aes_string("time", "beta")) + 
      geom_line(size = 2) + geom_ribbon(aes_(ymin = ~betalow, 
                                             ymax = ~betahigh), alpha = 0.5, col = "dodgerblue", 
                                        fill = "dodgerblue") + ylim(c(min(dat$contact$betalow), 
                                                                      max(dat$contact$betahigh))) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                                                                                                  .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == .(signif(dat$alpha, 
                                                                                                                                                                                         3)))) + ylab(bquote(beta))
    if (sum(sum(is.na(dat$contact))) > 0) {
      p3 <- ggplot(betadf, aes_string("time", "beta")) + 
        geom_line(size = 2) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                            .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                            .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
    }
    p4 <- logcorr(dat) + geom_abline(slope = 1, colour = "dodgerblue")
    drops <- c("mean", "sd", "error", "cases", "time")
    sim.only <- dat$res[, !(names(dat$res) %in% drops)]
    n <- ncol(sim.only)
    error <- qt(0.975, df = n - 1) * dat$res$sd/sqrt(n)
    dat$res$error <- error
    eb <- aes(ymax = mean + error, ymin = mean - error)
    p6 <- ggplot(data = dat$res, aes_string("time")) + theme(legend.position = "none") + 
      geom_line(aes_string(y = "cases"), colour = "dodgerblue", 
                size = 1) + xlab("year") + ylab("cases") + geom_line(aes_string(y = "mean"), 
                                                                     colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                                                                    alpha = 0.3) + theme_bw()
    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases
    p7 <- ggplot(data = inversecases, aes_string("time")) + 
      theme(legend.position = "none") + geom_line(aes_string(y = "cases"), 
                                                  colour = "dodgerblue", size = 1) + xlab("time") + 
      ylab("cases") + geom_line(aes_string(y = "mean"), 
                                colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                               alpha = 0.3) + theme_bw()
    if (dat$nsim > max.plot) {
      sampledat <- sample(sim.only, max.plot)
      sampledat$time <- dat$res$time
    }
    else {
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }
    meltdf <- melt(sampledat, id = "time")
    p8 <- ggplot(meltdf, aes_string(x = "time", y = "value", 
                                    fill = "variable")) + geom_line(alpha = 0.6, colour = "orangered4") + 
      xlab("time") + ylab("cases") + geom_line(data = dat$res, 
                                               aes_string(x = "time", y = "cases", fill = NA), colour = "dodgerblue", 
                                               size = 1) + theme_bw()
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 2)))
    print(p1, vp = vplayout(1, 1))
    print(p2, vp = vplayout(1, 2))
    print(p3, vp = vplayout(2, 1))
    print(p4, vp = vplayout(2, 2))
    print(p8, vp = vplayout(3, 1:2))
    print(p7, vp = vplayout(4, 1:2))
  }
  else {
    regdf <- NULL
    regdf$X <- dat$X
    regdf$Y <- dat$Y
    regdf$Yhat <- dat$Yhat
    regdf$time <- dat$res$time
    regdf <- as.data.frame(regdf)
    meltdf <- melt(regdf, id.vars = "time")
    p1 <- ggplot(meltdf, aes_string("time", "value", col = "variable")) + 
      geom_line(size = 2) + theme_bw()
    rhodf <- NULL
    rhodf$time <- dat$res$time
    rhodf$rho <- 1/dat$rho
    rhodf <- as.data.frame(rhodf)
    p2 <- ggplot(rhodf, aes_string("time", "rho")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(rho) == .(signif(mean(1/dat$rho), 
                                                       2)))) + ylab(bquote(1/rho))
    resdf <- NULL
    resdf$time <- dat$res$time
    resdf$Z <- dat$Z
    resdf$S <- dat$Z + dat$sbar
    resdf <- as.data.frame(resdf)
    meltres <- melt(resdf, id.vars = "time")
    p3 <- ggplot(meltres, aes_string("time", "value", col = "variable")) + 
      geom_line(size = 2) + theme_bw()
    loglikdf <- NULL
    loglikdf$sbar <- dat$Smean
    loglikdf$loglik <- dat$loglik
    loglikdf <- as.data.frame(loglikdf)
    p9 <- ggplot(loglikdf, aes_string("sbar", "loglik")) + 
      geom_line() + geom_point() + theme_bw() + geom_vline(xintercept = dat$sbar, 
                                                           linetype = "longdash") + ggtitle(bquote(bar(S) == 
                                                                                                     .(signif(dat$sbar, 2)) ~ "," ~ .(signif(dat$sbar/mean(dat$pop) * 
                                                                                                                                               100, 2)) ~ "%")) + xlab(bquote(bar(S)))
    betadf <- NULL
    betadf <- dat$contact
    betadf <- as.data.frame(betadf)
    p4 <- ggplot(betadf, aes_string("time", "beta")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(beta) == .(signif(mean(dat$beta), 
                                                        2)) ~ "," ~ alpha == .(signif(dat$alpha, 3)))) + 
      ylab(bquote(beta))
    if ("contact" %in% names(dat)) {
      p4 <- ggplot(data = dat$contact, aes_string("time", 
                                                  "beta")) + geom_line(size = 2) + geom_ribbon(aes_(ymin = ~betalow, 
                                                                                                    ymax = ~betahigh), alpha = 0.5, col = "dodgerblue", 
                                                                                               fill = "dodgerblue") + ylim(c(min(dat$contact$betalow), 
                                                                                                                             max(dat$contact$betahigh))) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                                                                                                                                                         .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                                                                                                                                                         .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
      if (sum(sum(is.na(dat$contact))) > 0) {
        p4 <- ggplot(betadf, aes_string("time", "beta")) + 
          geom_line(size = 2) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                              .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                              .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
      }
    }
    p4 <- p4 + xlab(sprintf("time mod %g", nrow(dat$contact)))
    if (dat$inits.fit == TRUE) {
      inits.grid <- dat$inits.grid
      p5 <- ggplot(inits.grid, aes_string(x = "S0", y = "I0", 
                                          z = "log10LS")) + geom_tile(aes_string(fill = "log10LS")) + 
        scale_fill_gradient(low = "white", high = "black") + 
        theme_bw() + geom_contour(col = "black") + geom_point(aes(x = dat$inits[1]/mean(dat$pop), 
                                                                  y = dat$inits[2]/mean(dat$pop)), col = "red") + 
        xlab("prop. init. sus.") + ylab("prop. init. inf.")
    }
    drops <- c("mean", "sd", "error", "cases", "time")
    sim.only <- dat$res[, !(names(dat$res) %in% drops)]
    n <- ncol(sim.only)
    error <- qt(0.975, df = n - 1) * dat$res$sd/sqrt(n)
    dat$res$error <- error
    eb <- aes(ymax = mean + error, ymin = mean - error)
    p6 <- ggplot(data = dat$res, aes_string("time")) + theme(legend.position = "none") + 
      geom_line(aes_string(y = "cases"), colour = "dodgerblue", 
                size = 1) + xlab("year") + ylab("cases") + geom_line(aes_string(y = "mean"), 
                                                                     colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                                                                    alpha = 0.3) + theme_bw()
    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases
    p7 <- ggplot(data = inversecases, aes_string("time")) + 
      theme(legend.position = "none") + geom_line(aes_string(y = "cases"), 
                                                  colour = "dodgerblue", size = 1) + xlab("time") + 
      ylab("cases") + geom_line(aes_string(y = "mean"), 
                                colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                               alpha = 0.3) + theme_bw()
    if (dat$nsim > max.plot) {
      sampledat <- sample(sim.only, max.plot)
      sampledat$time <- dat$res$time
    }
    else {
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }
    meltdf <- melt(sampledat, id = "time")
    p8 <- ggplot(meltdf, aes_string(x = "time", y = "value", 
                                    fill = "variable")) + geom_line(alpha = 0.6, colour = "orangered4") + 
      xlab("time") + ylab("cases") + geom_line(data = dat$res, 
                                               aes_string(x = "time", y = "cases", fill = NA), colour = "dodgerblue", 
                                               size = 1) + theme_bw()
    #grid.newpage()
    #pushViewport(viewport(layout = grid.layout(5, 2)))
    if (dat$inits.fit == TRUE) {
      if (all(is.na(dat$loglik)) == T) {
        return(list(p1,p2,p3,p4,p5,p8,p7))
        # print(p1, vp = vplayout(1, 1))
        # print(p2, vp = vplayout(1, 2))
        # print(p3, vp = vplayout(2, 1:2))
        # print(p4, vp = vplayout(3, 1))
        # print(p5, vp = vplayout(3, 2))
        # print(p8, vp = vplayout(4, 1:2))
        # print(p7, vp = vplayout(5, 1:2))
      }
      else {
        return(list(p1,p2,p3,p9,p4,p5,p8,p7))
        # print(p1, vp = vplayout(1, 1))
        # print(p2, vp = vplayout(1, 2))
        # print(p3, vp = vplayout(2, 1))
        # print(p9, vp = vplayout(2, 2))
        # print(p4, vp = vplayout(3, 1))
        # print(p5, vp = vplayout(3, 2))
        # print(p8, vp = vplayout(4, 1:2))
        # print(p7, vp = vplayout(5, 1:2))
      }
    }
    else {
      if (all(is.na(dat$loglik)) == T) {
        
        return(list(p1,p2,p3, p9,p4,p8,p7))
        # print(p1, vp = vplayout(1, 1))
        # print(p2, vp = vplayout(1, 2))
        # print(p3, vp = vplayout(2, 1:2))
        # print(p4, vp = vplayout(3, 1:2))
        # print(p8, vp = vplayout(4, 1:2))
        # print(p7, vp = vplayout(5, 1:2))
      }
      else {
        
        return(list(p1,p2,p3,p9,p4,p8,p7))
        # print(p1, vp = vplayout(1, 1))
        # print(p2, vp = vplayout(1, 2))
        # print(p3, vp = vplayout(2, 1))
        # print(p9, vp = vplayout(2, 2))
        # print(p4, vp = vplayout(3, 1:2))
        # print(p8, vp = vplayout(4, 1:2))
        # print(p7, vp = vplayout(5, 1:2))
      }
    }
  }
}
fit_tsir_gaussian <- function(df, epiyr, sbar){
  
  #first, fit using the gaussian assumption
  
  suffix = unique(df$provname)
  print(suffix)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  
  
  fittedpars <- estpars(data=dat,
                        IP=2, 
                        alpha=NULL, 
                        sbar=sbar, 
                        xreg = "cumcases",
                        regtype='gaussian',
                        family='poisson',
                        link='log')
  
  
  #p1 <- plotsbar(fittedpars2)
  
  simfitted <- simulatetsir(data=dat,
                            IP = 2,
                            parms=fittedpars,
                            #epidemics='break', threshold=3,
                            method='pois', nsim=100)
  
  
  out.plot <- plotres(dat=simfitted)
  
  
  pAB <- cowplot::plot_grid(out.plot[[1]], out.plot[[2]], nrow=1, ncol=2, align = "vh", labels = c("A", "B"))
  pCD <- cowplot::plot_grid(out.plot[[3]], out.plot[[4]],  nrow=1, ncol=2, align = "vh", labels = c("C", "D"))
  pEFG <- cowplot::plot_grid(out.plot[[5]], out.plot[[6]], out.plot[[7]],  nrow=3, ncol=1, align = "vh", labels = c("E", "F", "G"))
  
  
  
  pS2top <- cowplot::plot_grid(pAB, pCD, nrow = 2, ncol=1,align = "vh")
  FigOut <- cowplot::plot_grid(pS2top, pEFG, nrow=2, ncol=1, align = "vh", rel_heights = c(1,1.3))
  
  
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/TSIR-prov-fits/", suffix, "_", epiyr, ".png"),
         plot = FigOut,
         units="mm",  
         width=50, 
         height=60, 
         scale=3, 
         dpi=300)
  
  
  #and return beta
  beta.df <- cbind.data.frame(biwk=1:26, beta=simfitted$beta, betalow = simfitted$contact$betalow, betahigh = simfitted$contact$betahigh)
  beta.df$provname <- unique(dat$provname)
  rownames(beta.df) <- c()
  beta.df$epiyr <- epiyr
  beta.df$rsquared <- simfitted$rsquared
  
  return(beta.df)
  
}
fit_tsir <- function(df, epiyr, sbar){
  
  #first, fit using the gaussian assumption
  
  suffix = unique(df$provname)
  print(suffix)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  
  
  fittedpars <- estpars(data=dat,
                        IP=2, 
                        alpha=NULL, 
                        sbar=sbar, 
                        xreg = "cumcases",
                        regtype='lm',
                        family='poisson',
                        link='log')
  
  
  #p1 <- plotsbar(fittedpars)
  
  simfitted <- simulatetsir(data=dat,
                            IP = 2,
                            parms=fittedpars,
                            #epidemics='break', threshold=3,
                            method='pois', nsim=100)
  
  
  out.plot <- plotres(dat=simfitted)
  
  
  pAB <- cowplot::plot_grid(out.plot[[1]], out.plot[[2]], nrow=1, ncol=2, align = "vh", labels = c("A", "B"))
  pCD <- cowplot::plot_grid(out.plot[[3]], out.plot[[4]],  nrow=1, ncol=2, align = "vh", labels = c("C", "D"))
  pEFG <- cowplot::plot_grid(out.plot[[5]], out.plot[[6]], out.plot[[7]],  nrow=3, ncol=1, align = "vh", labels = c("E", "F", "G"))
  
  
  
  pS2top <- cowplot::plot_grid(pAB, pCD, nrow = 2, ncol=1,align = "vh")
  FigOut <- cowplot::plot_grid(pS2top, pEFG, nrow=2, ncol=1, align = "vh", rel_heights = c(1,1.3))
  
  
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/TSIR-prov-fits/", suffix, "_", epiyr, ".png"),
         plot = FigOut,
         units="mm",  
         width=50, 
         height=60, 
         scale=3, 
         dpi=300)
  
  
  #and return beta
  beta.df <- cbind.data.frame(biwk=1:26, beta=simfitted$beta, betalow = simfitted$contact$betalow, betahigh = simfitted$contact$betahigh)
  beta.df$provname <- unique(dat$provname)
  rownames(beta.df) <- c()
  beta.df$epiyr <- epiyr
  beta.df$rsquared <- simfitted$rsquared
  
  return(beta.df)
  
}


beta.df.2007 <- lapply(tsir.split.2007, runFulltSIR, epiyr1 = 2007, sbar1=NULL, fit.comp.df=fit.comp)
beta.df.2007 <- data.table::rbindlist(beta.df.2007)
head(beta.df.2007)
beta.df.2012 <- lapply(tsir.split.2012, runFulltSIR, epiyr1 = 2012, sbar1=NULL, fit.comp.df=fit.comp)
beta.df.2012 <- data.table::rbindlist(beta.df.2012)
head(beta.df.2012)
beta.df.2019 <- lapply(tsir.split.2019, runFulltSIR, epiyr1 = 2019, sbar1=NULL, fit.comp.df=fit.comp)
beta.df.2019 <- data.table::rbindlist(beta.df.2019)
head(beta.df.2019)

beta.all <- rbind(beta.df.2007, beta.df.2012, beta.df.2019)
head(beta.all) #phnom penh 2019 amd preah sihanouk 2019 are bad (below .2)

# pick those time series that are the most reliable (e.g. those with rsquared >.2) to do the climate lags
# everything else gets rejected
low.rsq <- subset(beta.all,rsquared<.4)
low.rsq <- ddply(low.rsq , .(provname, epiyr), summarise, rsquared = unique(rsquared))
low.rsq <- subset(beta.all,rsquared<.2)
low.rsq <- ddply(low.rsq , .(provname, epiyr), summarise, rsquared = unique(rsquared))

#maybe not worth trusting Mondul Kiri or Ratanak Kiri



#beta.all <- subset(beta.all,rsquared>=.2)

# Now, supplementary plot of Beta fitted by district
beta.all$epiyr <- as.factor(beta.all$epiyr)

pSupp <- ggplot(data=subset(beta.all, provname!="Mondul Kiri" & provname!="Ratanak Kiri")) + theme_bw() +
         geom_ribbon(aes(x=biwk, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.3) + 
         geom_line(aes(x=biwk, y= beta, color=epiyr), size=1) + scale_fill_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
         scale_color_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
         facet_wrap(provname~., scales = "free_y")+ylab(bquote(beta~', biweekly transmission')) +
         xlab("biweek of year") +
         theme(panel.grid = element_blank(), legend.position = c(.94,.06),
               strip.background = element_rect(fill="white"),
               axis.title = element_text(size=18), axis.text = element_text(size=13))


ggsave(file = paste0(homewd, "/final-figures/SuppFigBetaProv.png"),
       plot = pSupp,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


#write.csv(beta.all, file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), row.names = F)


# Now attach as an output on the case data

beta.merge = subset(beta.all, provname!="Mondul Kiri" & provname!="Ratanak Kiri")
beta.merge <- dplyr::select(beta.merge, provname, biwk, beta, betahigh, betalow, epiyr)
head(beta.merge)
unique(beta.merge$provname)
names(beta.merge)[names(beta.merge)=="biwk"] <- "biweek"

head(tsir.dat)
tsir.merge <- dplyr::select(tsir.dat, time, cases, year, biweek, provname)
tsir.merge$epiyr = 2007
tsir.merge$epiyr[tsir.merge$time>=2008] <- 2012
tsir.merge$epiyr[tsir.merge$time>=2013] <- 2019
tsir.merge = subset(tsir.merge, time<2020)
head(tsir.merge)
tail(tsir.merge)
unique(tsir.merge$provname)

setdiff(unique(tsir.merge$provname), unique(beta.merge$provname)) #"Mondul Kiri"  "Ratanak Kiri"
# and merge 

out.merge <- merge(tsir.merge, beta.merge, by = c("provname", "epiyr", "biweek"))
unique(out.merge$provname)
out.merge <- arrange(out.merge, provname, time)
head(out.merge) # pSupp above could be made from this dataset

# Now, write to data folder

write.csv(out.merge, file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), row.names = F)

# And in another attached this beta to the climate data (another script)

# This should be a supplementary plot those


#now do climate regression with beta

#and finally go to predictions with susceptibles

#and remake figure 1 - just a few examples???

