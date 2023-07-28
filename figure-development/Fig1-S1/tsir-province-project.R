rm(list=ls())

require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
require(plyr)



# as one script, split by province, pull in and reconstruct susceptibles by best method,
# fit TSIR, and run future year with climate projected beta that you have already saved elsewhere


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd))

#load the functions
###################################################################################

#this function finds the optimal fraction of increased susceptibles both
#with and without the climate-driven beta
find.frac.incS <- function(simfitted1, fittedpars1, dat1, epiyr1, sbar1, clim.beta){
  
  
  
  Its=simfitted1$res$mean*fittedpars1$rho #account for underreporting
  Ifinal<-Its[length(Its)]
  
  dat.sim <- cbind.data.frame(time=simfitted1$res$time, I=Its)
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  #make a vector of the increases to tested
  
  fracincvec <- seq(1,3,0.01)
  fracincvec.list <- as.list(fracincvec)
  
  #and apply this function over that vector to get deviations from the actual data
  out.inc.S <- lapply(fracincvec.list, run.increased.S, data.series=dat1, epiyr1=epiyr1,  simfitted1=simfitted1, fittedpars1=fittedpars1, beta.clim=clim.beta$beta)
  
  #now, fit the increase that has the lowest sum of squares for both the standard beta and the climate-driven beta
  sm.sq.noclim <- c(unlist(sapply(out.inc.S, '[[', 1)))
  sm.sq.clim <- c(unlist(sapply(out.inc.S, '[[', 2)))
  
  spot.noclim<-which(sm.sq.noclim==min(sm.sq.noclim))
  spot.clim<-which(sm.sq.clim==min(sm.sq.clim))
  fracincout.noclim <- fracincvec[spot.noclim]
  fracincout.clim <- fracincvec[spot.clim]
  
  Sbegin=simfitted1$sbar+simfitted1$Z
  S_noinc<-Sbegin[length(Sbegin)]
  
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr1, epiSnoic= S_noinc, frac_incS_noclim = fracincout.noclim, frac_incS_clim=fracincout.clim)
  
  return(dat.out)
  
  
  
  
}
run.increased.S <- function(fracincS, data.series, epiyr1, simfitted1, fittedpars1, beta.clim){
  
  Its=simfitted1$res$mean*fittedpars1$rho #account for underreporting
  Ifinal<-Its[length(Its)]
  
  time.start =  min(data.series$time)
  
  #the fitting data
  dat.fit = subset(data.series, time >= time.start & time<epiyr1)
  #the prediction data
  dat.pred= subset(data.series, time >= epiyr1 & time<=(epiyr1+1))
  
  
  
  Sbegin=simfitted1$sbar+simfitted1$Z
  Sepi<-Sbegin[length(Sbegin)]
  SepiInc<-fracincS*Sbegin[length(Sbegin)]
  
  #this is for no increased S
  
  
  # predict_ts <- predicttsir(times=dat.fit$time,
  #                           births = dat.fit$births,
  #                           beta = fittedpars1$beta,
  #                           alpha = fittedpars1$alpha,
  #                           S0 = Sbegin[1],
  #                           I0 = dat.fit$cases[1],
  #                           nsim=100,
  #                           stochastic = T)
  # 
  
  
  # and the epi prediction year-no increased S 
  # and using the same beta as the rest of the year
  
  
  
  # predict_epi_noInc <- predicttsir(times=dat.pred$time,
  #                                  births = dat.pred$births,
  #                                  beta = fittedpars1$beta,
  #                                  alpha = fittedpars1$alpha,
  #                                  S0 = Sepi,
  #                                  I0 = Ifinal,
  #                                  nsim=100,
  #                                  stochastic = T)
  # 
  
  # and epi-year prediction including increased S
  # and using the same beta as the rest of the year
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars1$beta,
                                 alpha = fittedpars1$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  
  # # and epi-year prediction with climate-driven beta but no increased S
  # predict_epi_noInc_clim <- predicttsir(times=dat.pred$time,
  #                                       births = dat.pred$births,
  #                                       beta = beta.clim,
  #                                       alpha = fittedpars1$alpha,
  #                                       S0 = Sepi,
  #                                       I0 = Ifinal,
  #                                       nsim=100,
  #                                       stochastic = T)
  # 
  
  
  #and the epiyear prediction with climate-driven beta AND increased S
  predict_epi_Inc_clim <- predicttsir(times=dat.pred$time,
                                      births = dat.pred$births,
                                      beta = beta.clim,
                                      alpha = fittedpars1$alpha,
                                      S0 = SepiInc,
                                      I0 = Ifinal,
                                      nsim=100,
                                      stochastic = T)
  
  
  
  
  #now look at data and model predictions for the epidemic year only across all the possibilities
  
  #here just TSIR
  
  #  IPredEpi1 = cbind.data.frame(time=predict_epi_noInc$I$time, mean_Inc=predict_epi_noInc$I$mean)
  #  names(IPredEpi1) <- c("time", "model_predicted_absolute_cases")
  #  IPredEpi1$model_predicted_reported_cases <- IPredEpi1$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  #  IPredEpi1$sim_type <- "TSIR fit"
  # # 
  # # here with increased S
  IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
  names(IPredEpi) <- c("time", "no_clim_model_predicted_absolute_cases")
  IPredEpi$no_clim_model_predicted_reported_cases <- IPredEpi$no_clim_model_predicted_absolute_cases/fittedpars1$rho[length(fittedpars1$rho)]
  
  
  #here with increased S and climate-driven beta
  IPredEpi2 = cbind.data.frame(time=predict_epi_Inc_clim$I$time, mean_Inc=predict_epi_Inc_clim$I$mean)
  names(IPredEpi2) <- c("time", "clim_model_predicted_absolute_cases")
  IPredEpi2$clim_model_predicted_reported_cases <- IPredEpi2$clim_model_predicted_absolute_cases/fittedpars1$rho[length(fittedpars1$rho)]
  
  
  
  
  #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
  #head(Iall)
  #head(all.dat)
  
  #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
  #and merge on time
  all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
  all.dat.merge <- merge(all.dat.merge, IPredEpi2, by="time", all.x = T, sort=F)
  
  # ggplot(data=all.dat.merge) + 
  #   geom_line(aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=model_predicted_reported_cases), color="tomato") 
  # # 
  # 
  
  #and get sum of sq differences for increased S
  all.dat.merge$sq_diff_noclim <- (all.dat.merge$cases - all.dat.merge$no_clim_model_predicted_reported_cases)^2
  all.dat.merge$sq_diff_clim <- (all.dat.merge$cases - all.dat.merge$clim_model_predicted_reported_cases)^2
  
  #and the sum of squared differences for increased S from the climate driven beta
  
  sm.sq.noclim = sum(all.dat.merge$sq_diff_noclim)
  sm.sq.clim = sum(all.dat.merge$sq_diff_clim)
  
  return(list(sm.sq.noclim, sm.sq.clim))
  
  # # ##saveRDS(predict_epi_Inc, paste0(homewd,"/sim.com.epi_predict_epi_Inc_",epiyr,".RDS") )
  # # 
  #  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
  #  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  # # 
  #  Iall <- rbind(IPredEpi, ISimTS)
  #  Iall$variable <- as.character(Iall$variable)
  #  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  #  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  # # 
  #  all.dat <- subset(dat, time<=(epiyr+1))
  #  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
  #  ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
  #    geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #    geom_line(aes(x=time, y=value, color=variable)) 
  #  
  #  ggplot(data=Iall) + 
  #    geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #    geom_line(aes(x=time, y=value, color=variable)) 
  # # 
  # 
  # #capture the increased S that best recovers the peak
  # ratioImax[i]<-max(Iall$value[Iall$variable=="prediction TSIR-increased S"])/max(all.dat$cases)
}
wrap.pipeline.TSIR <- function(df, epiyr1, sbar1, sus.dat, epi.beta.df){
  
  #first, identify the right susceptible reconstruction
  prov.choose = subset(sus.dat, provname == unique(df$provname) & epiyr == epiyr1)
  
  #and pick the right epiyear beta
  beta.epi = subset(epi.beta.df, provname == unique(df$provname) & epiyr == epiyr1)
  beta.epi <- arrange(beta.epi, time)
  
  suffix = unique(df$provname)
  #print(suffix)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  
  dat.fit = subset(df, time >= time.start & time<epiyr1)
  
  
  if(prov.choose$sus_reconstruction=="lm"){
    fittedpars <- estpars(data=dat.fit,
                          IP=2, 
                          alpha=NULL, 
                          sbar=sbar1,
                          xreg = "cumcases",
                          regtype='lm',
                          family='poisson',
                          link='log')
    
    
    simfitted <- simulatetsir(data=dat.fit,
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(prov.choose$sus_reconstruction=="gaussian"){
    
    fittedpars <- estpars(data=dat.fit,
                          IP=2, 
                          alpha=NULL, 
                          sbar=sbar1, 
                          xreg = "cumcases",
                          regtype='gaussian',
                          family='poisson',
                          link='log')
    
    
    #p1 <- plotsbar(fittedpars2)
    
    simfitted <- simulatetsir(data=dat.fit,
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  } 
  
  
  #now, take above and find the optimal increased S both with and without climate-driven beta
  incS.df <- find.frac.incS(simfitted1 = simfitted, fittedpars1 = fittedpars, dat1 = df, epiyr1 = epiyr1, sbar1 = sbar1, clim.beta =  epi.beta.df)
  
  
 
  # Run the fitted data series
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  Ifinal<-Its[length(Its)]
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  SepiIncNoClim<-Sbegin[length(Sbegin)]*incS.df$frac_incS_noclim
  SepiIncClim<-Sbegin[length(Sbegin)]*incS.df$frac_incS_clim
  
  
  
   predict_ts <- predicttsir(times=dat.fit$time,
                             births = dat.fit$births,
                             beta = fittedpars$beta,
                             alpha = fittedpars$alpha,
                             S0 = Sbegin[1],
                             I0 = dat.fit$cases[1],
                             nsim=100,
                             stochastic = T)
   
  
   #now run  all 4 time series for the epidemic year
   
   #prediction data
   dat.pred= subset(df, time >= epiyr1 & time<=(epiyr1+1))
   
   #1. project TSIR fit (no increased S + normal beta)
   
   predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                    births = dat.pred$births,
                                    beta = fittedpars$beta,
                                    alpha = fittedpars$alpha,
                                    S0 = Sepi,
                                    I0 = Ifinal,
                                    nsim=100,
                                    stochastic = T)
  
   
  #2. project epi year with increased S and normal beta
   
  # and epi-year prediction including increased S
  # and using the same beta as the rest of the year
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiIncNoClim,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  
  
  #3. project epi year with no increased S but climate-driven beta
  
  
  
   predict_epi_noInc_clim <- predicttsir(times=dat.pred$time,
                                         births = dat.pred$births,
                                         beta = epi.beta.df$beta,
                                         alpha = fittedpars$alpha,
                                         S0 = Sepi,
                                         I0 = Ifinal,
                                         nsim=100,
                                         stochastic = T)
  
  
  #4. project epi year with increased S and climate-driven beta
   
  predict_epi_Inc_clim <- predicttsir(times=dat.pred$time,
                                      births = dat.pred$births,
                                      beta = epi.beta.df$beta,
                                      alpha = fittedpars$alpha,
                                      S0 = SepiIncClim,
                                      I0 = Ifinal,
                                      nsim=100,
                                      stochastic = T)
  
  
  
  
  #now look at data and model predictions for the epidemic year only across all the possibilities
  
  
  #rho will depend on gaussian vs. linear
  
  #1. here just TSIR
  
  IPredEpi1 = cbind.data.frame(time=predict_epi_noInc$I$time, mean_Inc=predict_epi_noInc$I$mean)
  names(IPredEpi1) <- c("time", "model_predicted_absolute_cases")
  IPredEpi1$model_predicted_reported_cases <- IPredEpi1$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  IPredEpi1$sim_type <- "TSIR-prediction"
   
  
  #2.  here with increased S and normal beta
  IPredEpi2 = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
  names(IPredEpi2) <- c("time", "model_predicted_absolute_cases")
  IPredEpi2$model_predicted_reported_cases <- IPredEpi2$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  IPredEpi2$sim_type <- "increased-S-standard-beta"
  
  
  #3. here with no increase in S but a climate-driven beta
  IPredEpi3 = cbind.data.frame(time=predict_epi_noInc_clim$I$time, mean_Inc=predict_epi_noInc_clim$I$mean)
  names(IPredEpi3) <- c("time", "model_predicted_absolute_cases")
  IPredEpi3$model_predicted_reported_cases <- IPredEpi3$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  IPredEpi3$sim_type <- "no-increase-S-climate-beta"
  
  
  #4. here with both increased S and climate driven beta
  
  IPredEpi4 = cbind.data.frame(time=predict_epi_Inc_clim$I$time, mean_Inc=predict_epi_Inc_clim$I$mean)
  names(IPredEpi4) <- c("time", "model_predicted_absolute_cases")
  IPredEpi4$model_predicted_reported_cases <- IPredEpi4$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  IPredEpi4$sim_type <- "increased-S-climate-beta"
  
  
  # fitted TSIR
  IFit = cbind.data.frame(time=predict_ts$I$time, mean_Inc=predict_ts$I$mean)
  names(IFit) <- c("time", "model_predicted_absolute_cases")
  IFit$model_predicted_reported_cases <- IFit$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  IFit$sim_type <- "TSIR-fit"
  
  #bind all the predictions
  
  IPred <- rbind(IPredEpi1, IPredEpi2, IPredEpi3, IPredEpi4, IFit)
  
  
  

  #ggplot(data = IPred) + geom_line(aes(x=time, y= model_predicted_reported_cases, color=sim_type))
  
  
  #and merge with data on time
  
  df.merge = subset(df, time<=(epiyr1+1))
  all.dat.merge <- merge(df.merge, IPred, by="time", all.x = T, sort=F)
  
   # ggplot(data=all.dat.merge) + 
   #   geom_line(aes(x=time, y=cases),linetype=2) +
   #   geom_line(aes(x=time, y=model_predicted_reported_cases, color=sim_type)) 
   # 
   
   #and return both sets of data
   
   return(list(incS.df, all.dat.merge))
}
estpars <- function (data, xreg = "cumcases", IP = 2, seasonality = "standard", 
             regtype = "gaussian", sigmamax = 3, family = "gaussian", 
             link = "identity", userYhat = numeric(), alpha = NULL, sbar = NULL, 
             printon = F) {
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
  } #edited from TSIR package

###################################################################################

#load province data in tsir form - with associated beta info and climate (including climate projections)
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
 
length(unique(tsir.dat$provname)) #25 unique provinces
 
 #remove province with uneven time series
tsir.dat <- subset(tsir.dat, provname!="Tboung Khmum" & provname!="Mondul Kiri" & provname!= "Ratanak Kiri")

length(unique(tsir.dat$provname)) #22 used for climate projections

#load beta info with susceptible reconstruction
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
head(beta.df)
sus.merge <- ddply(beta.df, .(provname, epiyr), summarise, rsquared=unique(rsquared), sus_reconstruction=unique(sus_reconstruction))
head(sus.merge)

#and load the epiyear beta predictions using climate
clim.dat <- read.csv(file= paste0(homewd, "/data/tsir_dat_beta_climate_province.csv"), header = T, stringsAsFactors = F)
clim.dat <- read.csv(file= paste0(homewd, "/data/tsir_dat_beta_climate_province_linear.csv"), header = T, stringsAsFactors = F)
clim.dat <- read.csv(file= paste0(homewd, "/data/tsir_dat_beta_climate_province_2007.csv"), header = T, stringsAsFactors = F)
clim.dat <- read.csv(file= paste0(homewd, "/data/tsir_dat_beta_climate_province_linear_2007.csv"), header = T, stringsAsFactors = F)
#slim to epi year and beta

clim.dat <- arrange(clim.dat, provname, time)
test = subset(clim.dat, provname=="Battambang" & year == 2006)

ggplot(test) + geom_line(aes(x=time, y=beta))


clim.epi.dat = subset(clim.dat, year == 2007 | year == 2012 | year == 2019)

test.epi = subset(clim.dat, provname=="Battambang" & year == 2007)

ggplot(test.epi) + geom_line(aes(x=biweek, y=beta), color="red") + geom_line(data=test, aes(x=biweek, y=beta))

# 
#split by province and pre-epidemic period
tsir.split.2007 <- dlply(tsir.dat,.(provname))

tsir.dat.2012 <- subset(tsir.dat, time > 2007.9999)
tsir.split.2012 <- dlply(tsir.dat.2012,.(provname))

tsir.dat.2019 <- subset(tsir.dat, time > 2012.9999)
tsir.split.2019 <- dlply(tsir.dat.2019,.(provname))


#and try it on one subset
out.Battambang <- wrap.pipeline.TSIR(df = tsir.split.2007[[2]], epiyr1 = 2007, sbar1 = NULL, sus.dat = sus.merge, epi.beta.df = clim.epi.dat)

out.Battambang[[1]]
 ggplot(data=out.Battambang[[2]]) + 
   geom_line(aes(x=time, y=cases),linetype=2) +
   geom_line(aes(x=time, y=model_predicted_reported_cases, color=sim_type)) 


out.Battambang.2012 <- wrap.pipeline.TSIR(df = tsir.split.2012[[2]], epiyr1 = 2012, sbar1 = NULL, sus.dat = sus.merge, epi.beta.df = clim.epi.dat)
 
out.Battambang.2012[[1]]
 ggplot(data=out.Battambang.2012[[2]]) + 
   geom_line(aes(x=time, y=cases),linetype=2) +
   geom_line(aes(x=time, y=model_predicted_reported_cases, color=sim_type)) 
 
 

#Run this to test all the time series
fit.2007.plot <- lapply(tsir.split.2007, plot.test.tsir, epiyr = 2007, sbar=NULL)
fit.2007.plot <- data.table::rbindlist(fit.2007.plot)
fit.2007.plot$epiyr = 2007

fit.2012.plot <- lapply(tsir.split.2012, plot.test.tsir, epiyr = 2012, sbar=NULL)
fit.2012.plot <- data.table::rbindlist(fit.2012.plot)
fit.2012.plot$epiyr = 2012

fit.2019.plot <- lapply(tsir.split.2019, plot.test.tsir, epiyr = 2019, sbar=NULL)
fit.2019.plot <- data.table::rbindlist(fit.2019.plot)
fit.2019.plot$epiyr = 2019

fit.comp <- rbind(fit.2007.plot, fit.2012.plot, fit.2019.plot)

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

#remove Mondul Kiri and Ratanak Kiri

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

# And in another script attach this beta to the climate data,
# do a climate regression with beta, and use this regression to 
# predict beta for the epidemic years
