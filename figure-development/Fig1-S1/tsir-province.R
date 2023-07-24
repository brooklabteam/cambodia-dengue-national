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



# #remove 8 provinces with few cases and sporadic dynamics
tsir.dat.split <- subset(tsir.dat, provname!="Kep")
tsir.dat.split <- subset(tsir.dat.split, provname!="Banteay Meanchey")
tsir.dat.split <- subset(tsir.dat.split, provname!= "Mondul Kiri")
tsir.dat.split <- subset(tsir.dat.split, provname!="Pailin")
tsir.dat.split <- subset(tsir.dat.split, provname!="Koh Kong")
tsir.dat.split <- subset(tsir.dat.split, provname!="Ratanak Kiri")
tsir.dat.split <- subset(tsir.dat.split, provname!="Stung Treng")
tsir.dat.split <- subset(tsir.dat.split, provname!="Tboung Khmum")

#split by province and pre-epidemic period
tsir.split.2007 <- dlply(tsir.dat.split,.(provname))

tsir.dat.2012 <- subset(tsir.dat.split, time > 2007.9999)
tsir.split.2012 <- dlply(tsir.dat.2012,.(provname))

tsir.dat.2019 <- subset(tsir.dat.split, time > 2012.9999)
tsir.split.2019 <- dlply(tsir.dat.2019,.(provname))


plot.test.tsir <- function(df, epiyr, sbar1){ 
  
  #get the suffix to all the filenames
  
  suffix = unique(df$provname)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  
  
  
  p1 <- plotdata(dplyr::select(dat, -(provname)))
  
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
 
  
  pB <- cowplot::plot_grid(p2, p3, ncol = 1, nrow=2)
  pC <- cowplot::plot_grid(p4, p5, ncol = 1, nrow=2)
  
  out.plot <- cowplot::plot_grid(p1, pB, pC, ncol = 3, nrow=1, rel_widths = c(1,1.1,1.1), labels = c("A. Data", "B. lm", "C. gaussian"), label_size = 22)
  
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
    # out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(simfitted1$rsquared, simfitted2$rsquared), Prop_Susceptible = c((signif(fittedpars1$sbar/mean(dat$pop) * 100, 2)), (signif(fittedpars2$sbar/mean(dat$pop) * 100, 2))))
     
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
  # return(out.df)
  
}


#Run this to test all the time series
fit.2007.plot <- lapply(tsir.split.2007, plot.test.tsir, epiyr = 2007, sbar=NULL)
fit.2012.plot <- lapply(tsir.split.2012, plot.test.tsir, epiyr = 2012, sbar=NULL)
fit.2019.plot <- lapply(tsir.split.2019, plot.test.tsir, epiyr = 2019, sbar=NULL)


#and write over the few that need help
plot.test.tsir(tsir.split.2012$Kandal, epiyr = 2012, sbar1 = median(c(.57, .62)))
plot.test.tsir(df = tsir.split.2019$`Phnom Penh`, epiyr = 2019, sbar1 = median(c(.43, .46)))
plot.test.tsir(df = tsir.split.2019$`Kampong Thom`, epiyr = 2019, sbar1 = median(c(.51, .59)))
plot.test.tsir(df = tsir.split.2019$`Kampong Cham`, epiyr = 2019, sbar1 = median(c(.44, .84)))
plot.test.tsir(df = tsir.split.2019$`Siem Reap`, epiyr = 2019, sbar1 = .36)
plot.test.tsir(df = tsir.split$`Siem Reap`, epiyr = 2007, sbar1 = .36)
plot.test.tsir(df = tsir.split$`Preah Vihear`, epiyr = 2007, sbar1 =.45)


# Now rerun gaussian and fit and save data
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

beta.df.2007 <- lapply(tsir.split.2007, fit_tsir, epiyr = 2007, sbar=NULL)
beta.df.2007 <- data.table::rbindlist(beta.df.2007)
head(beta.df.2007)
beta.df.2012 <- lapply(tsir.split.2012, fit_tsir, epiyr = 2012, sbar=NULL)
beta.df.2012 <- data.table::rbindlist(beta.df.2012)
head(beta.df.2012)
beta.df.2019 <- lapply(tsir.split.2019, fit_tsir, epiyr = 2019, sbar=NULL)
beta.df.2019 <- data.table::rbindlist(beta.df.2019)
head(beta.df.2019)

beta.all <- rbind(beta.df.2007, beta.df.2012, beta.df.2019)
head(beta.all)


#now, re-do a few for which you decided to fix the mean susceptible pop
#check each first before you overwrite it
kandall.2012 <- fit_tsir(df = tsir.split.2012$Kandal, epiyr = 2012, sbar = median(c(.57, .62)))
phnom.penh.2019 <- fit_tsir(df = tsir.split.2019$`Phnom Penh`, epiyr = 2019, sbar = .7)
kampong.thom.2019 <- fit_tsir(df = tsir.split.2019$`Kampong Thom`, epiyr = 2019, sbar = .9)
kampong.cham.2019 <- fit_tsir(df = tsir.split.2019$`Kampong Cham`, epiyr = 2019, sbar = .95)
siem.reap.2019 <- fit_tsir(df = tsir.split.2019$`Siem Reap`, epiyr = 2019, sbar = .9)
siem.reap.2007 <- fit_tsir(df = tsir.split.2007$`Siem Reap`, epiyr = 2007, sbar = .7)
preah.vihear.2007 <- fit_tsir(df = tsir.split.2007$`Preah Vihear`, epiyr = 2007, sbar =.45)

#now write over this in the overall beta fits
beta.all[beta.all$provname=="Kandal" & beta.all$epiyr==2012,] <- kandall.2012 
beta.all[beta.all$provname=="Phnom Penh" & beta.all$epiyr==2019,] <- phnom.penh.2019
beta.all[beta.all$provname=="Kampong Thom" & beta.all$epiyr==2019,] <- kampong.thom.2019
beta.all[beta.all$provname=="Kampong Cham" & beta.all$epiyr==2019,] <-kampong.cham.2019
beta.all[beta.all$provname=="Siem Reap" & beta.all$epiyr==2019,] <- siem.reap.2019
beta.all[beta.all$provname=="Siem Reap" & beta.all$epiyr==2007,] <- siem.reap.2007
beta.all[beta.all$provname=="Preah Vihear" & beta.all$epiyr==2007,] <- preah.vihear.2007


# pick those time series that are the most reliable (e.g. those with rsquared >.2)
# everything else gets rejected
low.rsq <- subset(beta.all,rsquared<.2)
low.rsq <- ddply(low.rsq , .(provname, epiyr), summarise, rsquared = unique(rsquared))


beta.all <- subset(beta.all,rsquared>=.2)

# Now, supplementary plot of Beta fitted by district
beta.all$epiyr <- as.factor(beta.all$epiyr)

pSupp <- ggplot(data=beta.all)+ theme_bw() +
         geom_ribbon(aes(x=biwk, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.3) + 
         geom_line(aes(x=biwk, y= beta, color=epiyr), size=1) + scale_fill_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
         scale_color_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
         facet_wrap(~provname, scales = "free_y")+ylab(bquote(beta~', biweekly transmission')) +
         xlab("biweek of year") +
         theme(panel.grid = element_blank(), legend.position = c(.9,.1),
               strip.background = element_rect(fill="white"),
               axis.title = element_text(size=18), axis.text = element_text(size=13))


ggsave(file = paste0(homewd, "/final-figures/SuppFigBeta.png"),
       plot = pSupp,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


# Now, attached this beta to the climate data (another script)
# This should be a supplementary plot those


#now do climate regression with beta

#and finally go to predictions with susceptibles

#and remake figure 1 - just a few examples???

#and bind all
beta.all$epiyr <- as.factor(beta.all$epiyr)
ggplot(data=beta.all) + geom_line(aes(x=biwk, y= beta, color=provname), size=1) + facet_grid(~epiyr)
ggplot(data=subset(beta.all, rsquared>.2)) + geom_line(aes(x=biwk, y= beta, color=provname), size=1) + facet_grid(~epiyr)
ggplot(data=subset(beta.all, rsquared>.2))+
  geom_line(aes(x=biwk, y= beta, color=epiyr), size=1) + facet_wrap(~provname, scales = "free_y")


tmp <- subset(beta.all, rsquared<0.5)
head(tmp)
tmp <- ddply(tmp,.(provname, epiyr), summarise, rsquared = unique(rsquared))
# Then, plot all beta together and with climate parameters
# Then calc all lags at province and national level
# Then fit panel regression and repull climate data




#fit gaussian
fit.2007.gauss <- lapply(tsir.split, fit.gauss.tsir, epiyr = 2007)
fit.2007.gauss <- data.table::rbindlist(fit.2007.gauss)



plot.test.tsir(tsir.dat, epiyr = 2007, sbar1 =.5)



#and apply fits on the 2007 data
fit.2007 <- lapply(tsir.split, fit.test.tsir, epiyr = 2007)
fit.2007.df <- data.table::rbindlist(fit.2007)

fit.2007.gauss <- lapply(tsir.split, fit.gauss.tsir, epiyr = 2007)




fit.test.tsir <- function(df, epiyr){
  
  #get the suffix to all the filenames
  
  suffix = unique(df$provname)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  plotdata(dplyr::select(dat, -(provname)))
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-data-plots/data_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  fittedpars1 <- estpars(data=dat,
                         IP=2, 
                         alpha=NULL, 
                         sbar=NULL, 
                         xreg = "cumcases",
                         regtype='lm',
                         family='poisson',
                         link='log')
  
  
  
  plotsbar(fittedpars1)
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/lm_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  
  if ((fittedpars1$sbar + fittedpars1$Z[1])>0){
    
    
    simfitted1 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars1,
                               epidemics='break', threshold=3,
                               method='pois', nsim=100)
    
  }
  
  fittedpars2 <- estpars(data=dat,
                         IP=2, 
                         alpha=NULL, 
                         sbar=NULL, 
                         xreg = "cumcases",
                         regtype='gaussian',
                         family='poisson',
                         link='log')
  
  
  plotsbar(fittedpars2)
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/gaussian_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  if ((fittedpars2$sbar + fittedpars2$Z[1])>0){
    
    simfitted2 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars2,
                               #epidemics='break', threshold=3,
                               method='pois', nsim=100)
    
  }
  
  ## The type of regression used in susceptible reconstruction. 
  # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
  
  ## The family in the GLM regression. 
  # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
  if ( ((fittedpars1$sbar + fittedpars1$Z[1])>0) & ((fittedpars2$sbar + fittedpars2$Z[1])>0)) {
    
    out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(simfitted1$rsquared, simfitted2$rsquared))
  } else if ( ((fittedpars1$sbar + fittedpars1$Z[1])<0) & ((fittedpars2$sbar + fittedpars2$Z[1])>0)){
    
    out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(NA, simfitted2$rsquared))
    
  }else if ( ((fittedpars1$sbar + fittedpars1$Z[1])>0) & ((fittedpars2$sbar + fittedpars2$Z[1])<0)){
    
    out.df <- cbind.data.frame(regtype=c("lm", "gaussian"), AIC = c(fittedpars1$AIC, fittedpars2$AIC), rsquared=c(simfitted1$rsquared, NA))
  }
  out.df$provname <- suffix
  
  return(out.df)
  
}
fit.gauss.tsir <- function(df, epiyr){
  
  #get the suffix to all the filenames
  
  suffix = unique(df$provname)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  plotdata(dplyr::select(dat, -(provname)))
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-data-plots/data_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  
  fittedpars2 <- estpars(data=dat,
                         IP=2, 
                         alpha=NULL, 
                         sbar=NULL, 
                         xreg = "cumcases",
                         regtype='gaussian',
                         family='poisson',
                         link='log')
  
  
  plotsbar(fittedpars2)
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/gaussian_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  if ((fittedpars2$sbar + fittedpars2$Z[1])>0){
    
    simfitted2 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars2,
                               #epidemics='break', threshold=3,
                               method='pois', nsim=100)
    
  }
  
  
  
  out.df <- cbind.data.frame(regtype=c("gaussian"), AIC = c( fittedpars2$AIC), rsquared=c(simfitted2$rsquared))
  
  
  
  out.df$provname <- suffix
  
  out.df$Sbar = fittedpars2$sbar
  out.df$SbarProp <- out.df$Sbar/mean(df$pop)
  
  return(out.df)
  
}
fit.loess.tsir <- function(df, epiyr){
  
  #get the suffix to all the filenames
  
  suffix = unique(df$provname)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  dat = subset(df, time >= time.start & time<epiyr)
  # plotdata(dplyr::select(dat, -(provname)))
  # 
  # ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-data-plots/data_", epiyr,"_", suffix, ".png"),
  #        units="mm",  
  #        width=100, 
  #        height=110, 
  #        scale=3, 
  #        dpi=300)
  # 
  
  fittedpars2 <- estpars(data=dat,
                         IP=2, 
                         alpha=NULL, 
                         sbar=NULL, 
                         xreg = "cumcases",
                         regtype='loess',
                         family='poisson',
                         link='log')
  
  
  plotsbar(fittedpars2)
  
  ggsave(file = paste0(homewd, "/figure-development/Fig1-S1/prov-Sbar-plots/loess_", epiyr,"_", suffix, ".png"),
         units="mm",  
         width=100, 
         height=110, 
         scale=3, 
         dpi=300)
  
  if ((fittedpars2$sbar + fittedpars2$Z[1])>0){
    
    simfitted2 <- simulatetsir(data=dat,
                               IP = 2,
                               parms=fittedpars2,
                               #epidemics='break', threshold=3,
                               method='pois', nsim=100)
    
  }
  
  
  
  out.df <- cbind.data.frame(regtype=c("loess"), AIC = c( fittedpars2$AIC), rsquared=c(simfitted2$rsquared))
  
  
  
  out.df$provname <- suffix
  
  out.df$Sbar = fittedpars2$sbar
  out.df$SbarProp <- out.df$Sbar/mean(df$pop)
  
  return(out.df)
  
}

fit.2012.gauss <- lapply(tsir.split.2012, fit.gauss.tsir, epiyr = 2012)
fit.2012.gauss <- data.table::rbindlist(fit.2012.gauss)
fit.2019.gauss <- lapply(tsir.split.2019, fit.gauss.tsir, epiyr = 2019)
fit.2019.gauss <- data.table::rbindlist(fit.2019.gauss)



fit.2007.loess <- lapply(tsir.split, fit.loess.tsir, epiyr = 2007)
fit.2007.loess <- data.table::rbindlist(fit.2007.loess)
















#split by province
tsir.split <- dlply(tsir.dat,.(provname))



fit.2012.gauss <- data.table::rbindlist(fit.2007.gauss)



  
#and apply fits on the 2007 data
fit.2007 <- lapply(tsir.split, fit.TSIR, family="poisson", epiyr = 2007)

dat = tsir.split$`Phnom Penh`

fit.all <- lapply( tsir.split, fit.TSIR, family="gaussian", epiyr = 2019)

#and apply fits on the 2007 data
fit.2007 <- lapply(tsir.split, fit.TSIR, family="gaussian", epiyr = 2007)



simulatetsir <- function (data, nsim = 100, IP = 2, parms, method = "deterministic", 
                                    epidemics = "cont", pred = "forward", threshold = 1, inits.fit = FALSE, 
                                    add.noise.sd = 0, mul.noise.sd = 0) 
{
  nzeros <- length(which(data$cases == 0))
  ltot <- length(data$cases)
  if (nzeros > 0.3 * ltot && epidemics == "cont") {
    print(sprintf("time series is %.0f%% zeros, consider using break method", 
                  100 * nzeros/ltot))
  }
  Smean <- parms$Smean
  beta <- parms$beta
  adj.rho <- parms$rho
  Z <- parms$Z
  sbar <- parms$sbar
  alpha <- parms$alpha
  X <- parms$X
  Y <- parms$Y
  Yhat <- parms$Yhat
  contact <- parms$contact
  alphalow <- parms$alphalow
  alphahigh <- parms$alphahigh
  loglik <- parms$loglik
  pop <- data$pop
  datacopy <- data
  period <- rep(1:(52/IP), round(nrow(data) + 1))[1:(nrow(data) - 
                                                       1)]
  if (IP == 1) {
    period <- rep(1:(52/2), each = 2, round(nrow(data) + 
                                              1))[1:(nrow(data) - 1)]
  }
  S <- rep(0, length(data$cases))
  I <- rep(0, length(data$cases))
  nsample <- 30
  inits.grid <- expand.grid(S0 = seq(0.01 * mean(pop), 0.1 * 
                                       mean(pop), length = nsample), I0 = seq(0.01 * 0.001 * 
                                                                                mean(pop), 1 * 0.001 * mean(pop), length = nsample))
  if (inits.fit == TRUE) {
    inits.res <- rep(NA, nsample * nsample)
    for (it in 1:nrow(inits.grid)) {
      S0 <- inits.grid[it, 1]
      I0 <- inits.grid[it, 2]
      S[1] <- S0
      I[1] <- I0
      for (t in 2:(nrow(data))) {
        lambda <- min(S[t - 1], unname(beta[period[t - 
                                                     1]] * S[t - 1] * (I[t - 1])^alpha))
        if (is.nan(lambda) == T) {
          lambda <- 0
        }
        I[t] <- lambda
        if (epidemics == "cont") {
          I[t] <- I[t]
        }
        if (epidemics == "break") {
          t0s <- epitimes(data, threshold)$start
          if (t %in% t0s) {
            I[t] <- adj.rho[t] * data$cases[t]
          }
        }
        S[t] <- max(S[t - 1] + data$births[t - 1] - I[t], 
                    0)
      }
      inits.res[it] <- sum((I - data$cases * adj.rho)^2)
    }
    inits <- inits.grid[which.min(inits.res), ]
    inits.grid$S0 <- inits.grid$S0/mean(pop)
    inits.grid$I0 <- inits.grid$I0/mean(pop)
    inits.grid$log10LS <- log10(inits.res)
    S_start <- inits[[1]]
    I_start <- inits[[2]]
  }else {
    S_start <- sbar + Z[1]
    I_start <- adj.rho[1] * datacopy$cases[1]
  }
  IC <- c(S_start, I_start)
  print(c(alpha = unname(signif(alpha, 2)), `mean beta` = unname(signif(mean(beta), 
                                                                        3)), `mean rho` = unname(signif(mean(1/adj.rho), 3)), 
          `mean sus` = unname(signif(sbar, 3)), `prop. init. sus.` = unname(signif(S_start/mean(pop), 
                                                                                   3)), `prop. init. inf.` = unname(signif(I_start/mean(pop), 
                                                                                                                           3))))
  nsim <- nsim
  res <- matrix(0, length(data$cases), nsim)
  Sres <- matrix(0, length(data$cases), nsim)
  for (ct in 1:nsim) {
    S <- rep(0, length(data$cases))
    I <- rep(0, length(data$cases))
    S[1] <- S_start
    I[1] <- I_start
    for (t in 2:(nrow(data))) {
      if (pred == "step-ahead") {
        lambda <- min(S[t - 1], unname(beta[period[t - 
                                                     1]] * S[t - 1] * (adj.rho[t - 1] * data$cases[t - 
                                                                                                     1])^alpha))
      }
      if (pred == "forward") {
        I <- I
        lambda <- min(S[t - 1], unname(beta[period[t - 
                                                     1]] * S[t - 1] * (I[t - 1])^alpha))
      }
      if (is.nan(lambda) == T) {
        lambda <- 0
      }
      if (method == "deterministic") {
        I[t] <- lambda * rnorm(n = 1, mean = 1, sd = mul.noise.sd)
        if (I[t] < 0 && lambda >= 0) {
          warning("infected overflow  -- reduce multiplicative noise sd")
        }
      }
      if (method == "negbin") {
        I[t] <- rnbinom(n = 1, mu = lambda, size = I[t - 
                                                       1] + 1e-10)
      }
      if (method == "pois") {
        I[t] <- rpois(n = 1, lambda = lambda)
      }
      if (epidemics == "cont") {
        I[t] <- I[t]
      }
      if (epidemics == "break") {
        t0s <- epitimes(data, threshold)$start
        if (t %in% t0s) {
          I[t] <- adj.rho[t] * data$cases[t]
        }
      }
      S[t] <- max(S[t - 1] + data$births[t - 1] - I[t] + 
                    rnorm(n = 1, mean = 0, sd = add.noise.sd), 0)
      if (S[t] < 0 && (S[t - 1] + data$births[t - 1] - 
                       I[t]) > 0) {
        warning("susceptible overflow  -- reduce additive noise sd")
      }
    }
    res[, ct] <- I/adj.rho
    Sres[, ct] <- S
  }
  res[is.nan(res)] <- 0
  res[res < 1] <- 0
  res <- as.data.frame(res)
  Sres <- as.data.frame(Sres)
  Sres$mean <- rowMeans(Sres, na.rm = T)
  Sres$sd <- apply(Sres, 1, function(row) sd(row[-1], na.rm = T))
  Sres$time <- data$time
  res$mean <- rowMeans(res, na.rm = T)
  res$sd <- apply(res, 1, function(row) sd(row[-1], na.rm = T))
  res$time <- data$time
  res$cases <- data$cases
  obs <- res$cases
  pred <- res$mean
  fit <- lm(pred ~ obs)
  rsquared <- signif(summary(fit)$adj.r.squared, 2)
  return(list(X = X, Y = Y, Yhat = Yhat, pop = pop, Smean = Smean, 
              IP = IP, beta = beta, rho = adj.rho, Z = Z, sbar = sbar, 
              alpha = alpha, pop = pop, alphalow = alphalow, alphahigh = alphahigh, 
              res = res, simS = Sres, loglik = loglik, nsim = nsim, 
              contact = contact, rsquared = rsquared, inits.fit = inits.fit, 
              inits.grid = inits.grid, inits = IC))
}








#and apply fits on the 2007 data
fit.2007 <- lapply( tsir.split, fit.TSIR, family="gaussian", epiyr = 2007)

extract.rsq <- function(df){
  rsquared = df$rsquared
  return(rsquared)
}

fitR2.df <- cbind.data.frame(names = names(fit.2007), year=2007)
fitR2.df$rsquared <-  c(unlist(lapply(fit.2007, extract.rsq)))
fitR2.df$flag = 0
fitR2.df$flag[fitR2.df$rsquared>.5] <- 1


tsir.dat.2012 <- subset(tsir.dat, year>2007)

#split by province
tsir.split.2012 <- dlply(tsir.dat.2012,.(provname))


#and apply fits on the 2012 data
fit.2012 <- lapply( tsir.split.2012, fit.TSIR, family="gaussian", epiyr = 2012)


fitR2.df2 <- cbind.data.frame(names = names(fit.2012), year=2012)
fitR2.df2$rsquared <-  c(unlist(lapply(fit.2012, extract.rsq)))
fitR2.df2$flag = 0
fitR2.df2$flag[fitR2.df2$rsquared>.5] <- 1


tsir.dat.2019 <- subset(tsir.dat, year>2012)

#split by province
tsir.split.2019 <- dlply(tsir.dat.2019,.(provname))

fit.2019 <- lapply( tsir.split.2019, fit.TSIR, family="gaussian", epiyr = 2019)


fitR2.df3 <- cbind.data.frame(names = names(fit.2019), year=2019)
fitR2.df3$rsquared <-  c(unlist(lapply(fit.2019, extract.rsq)))
fitR2.df3$flag = 0
fitR2.df3$flag[fitR2.df3$rsquared>.5] <- 1


#extract and plot beta
extract.beta <- function(df, name.df){
  beta.df <- cbind.data.frame(biwk = 1:26, beta=df$beta)
  beta.df$provname <- name.df
  beta.df$rsquared <- df$rsquared
  return(beta.df)
}
beta.df <- data.table::rbindlist(mapply(extract.beta, fit.2019, as.list(names(fit.2019)), SIMPLIFY = F))

head(beta.df)

ggplot(beta.df) + geom_point(aes(x=biwk, y=beta, color=provname)) +
    geom_line(aes(x=biwk, y=beta, color=provname)) + 
    facet_wrap(~provname, scales = "free_y")

ggplot(data=subset(beta.df, rsquared>.3)) + geom_point(aes(x=biwk, y=beta, color=provname)) +
  geom_line(aes(x=biwk, y=beta, color=provname)) + 
  facet_wrap(~provname, scales = "free_y")


#fit to pre epidemic period with increased S
sim.2007 <- sim.com.epi(dat = tsir.df,#at,
                        time.start =  min(tsir.df$time),
                        family="gaussian",
                        epiyr = 2007)

#and fit with increased beta
sim.beta.2007 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)

#now return and recover Sdat for increased S
Sdat2007 <- sim.return.S.incS(fracincS = sim.2007$frac_incS,
                              dat=tsir.dat,
                              time.start =  min(tsir.dat$time),
                              family="gaussian",
                              epiyr = 2007)


#and sim with increased S
comp.2007 <- sim.with.increaseS(
  fracincS = sim.2007$frac_incS,
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)


#and sim with increased beta
comp.beta.2007 <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2007$frac_incBeta,
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)


#now for 2012 - fit with increased S
sim.2012 <- sim.com.epi(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                        family="gaussian",
                        epiyr = 2012)

#and with increased beta
sim.beta.2012 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                        family="gaussian",
                        epiyr = 2012)


#now return S from increased S
Sdat2012 <- sim.return.S.incS(fracincS = sim.2012$frac_incS,
                              dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2008]) & time<2013),
                              time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                              family="gaussian",
                              epiyr = 2012)


#now sim with increased S
comp.2012 <- sim.with.increaseS(
  fracincS = sim.2012$frac_incS,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2008]) & time < 2013),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)

#and sim with increased beta
comp.beta.2012  <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2012$frac_incBeta,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2008]) & time < 2013),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)


#and fit 2019 with increased S
sim.2019 <- sim.com.epi(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                        family="gaussian",
                        epiyr = 2019)

#and with increased beta
sim.beta.2019 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                        family="gaussian",
                        epiyr = 2019)


#now return S from increased S
Sdat2019 <- sim.return.S.incS(fracincS = sim.2019$frac_incS,
                              dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2013]) &time<2020),
                              time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                              family="gaussian",
                              epiyr = 2019)


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


#and sim with increased beta
comp.beta.2019  <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2019$frac_incBeta,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2013]) & time < 2020),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2013]),
  family = "gaussian",
  epiyr = 2019
)





#and combine for the susceptible plot
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

head(Sdatcombined)
unique(Sdatcombined$sim)

#and write data file
#write.csv(Sdatcombined, file = paste0(homewd, "/SdatCara.csv"), row.names = F)
saveRDS(Sdatcombined,file = paste0(homewd, "/data/Sdatcombined.rds"))



##################################################################
####################################################################
#####################################################################
######################################################################

#and save the sims of increased S and not for the joint dataset
all.comp <- rbind(comp.2007, comp.2012, comp.2019)

#and combine beta sims and add to dataset
comp.beta <- rbind(comp.beta.2007, comp.beta.2012, comp.beta.2019)
comp.beta = subset(comp.beta, variable == "prediction TSIR-increased Beta")

#and combine
all.comp <- rbind(all.comp, comp.beta)

#and save all comp as data
write.csv(all.comp, file = paste0(homewd, "/data/TSIR_fitted_timeseries_estimates.csv"), row.names = F)

my_blue<-"#2980B9"
my_orange<-"#D35400"
my_red<-"tomato"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
my_purple <- "purple"



#now plot the whole thing
# all.IncBeta <- rbind(comp.2007, comp.2012, comp.2019)
# all.IncBeta <- subset(all.IncBeta, variable == "prediction TSIR-increased Beta")
# all.IncBeta$epiyr <- 1
# 
# #and join
# Sdatcombined <- rbind(Sdatcombined, all.IncBeta)
#and plot all
dat=all.comp
dat$variable = factor(dat$variable, levels=c("data",
                                             "TSIR fit",
                                             "prediction TSIR",
                                             "prediction TSIR-increased S",
                                             "prediction TSIR-increased Beta"))
colz = c('data' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR-increased S' =my_red, 'prediction TSIR-increased Beta' = my_purple)
typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1, 'prediction TSIR-increased Beta' = 1)


pl_a <- ggplot(data=dat, aes(x=time, y=reported_cases, color=variable, linetype=variable)) + 
  geom_line(data=dplyr::filter(dat, variable=="data"), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time<2007), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2008 & time<2012), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2019 & time<2020), size=1) +
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















