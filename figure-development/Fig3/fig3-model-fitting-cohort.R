library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(sjPlot)
library(MuMIn)
library(relaimpo)
library(mgcv)
library(stringr)
library(tidyr)


rm(list=ls())

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

#load the master list
dat <- read.csv(file = paste0(homewd, "/data/IDseq_PAGODAS_ALL_metadata_through_2020_CLEAN.csv"), header=T, stringsAsFactors = F)
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
head(dat)
max(dat$date) #2020-12-31
min(dat$date) # 2018-07-16



#select the positives
dat.pos = subset(dat, dengue_result==1) #106
#there are 4 'positives' with no corresponding sequence: "109-0051" "109-0211" "109-0342" "100-0424"
max(dat.pos$date) #"2020-09-23"
min(dat.pos$date) # "2019-03-15"


#there is just one missing igg data - 
subset(dat.pos, is.na(igg_res)) #100-0052, 100-0632, 109-0211

#for this analysis, we only care about age, infection status, and igg status
dat.pos$age <- as.numeric(dat.pos$age)
dat.pos$age <- round(dat.pos$age, 0)
dat.pos$year <- year(dat.pos$date)
dat.pos$igg_res[dat.pos$igg_res=="pos" & !is.na(dat.pos$igg_res)] <- 1
dat.pos$igg_res[dat.pos$igg_res=="neg"  & !is.na(dat.pos$igg_res)] <- 0


#just 2019
dat.2019 = subset(dat.pos, year==2019) #80
dat.2020 = subset(dat.pos, year==2020) #26
dat.2019 = subset(dat.2019, !is.na(igg_res)) #80
dat.2020 = subset(dat.2020, !is.na(igg_res)) #25. removing one denv 2 idseq: 109-0211


dat.2019 = subset(dat.2019, !is.na(age)) #80.
dat.2020 = subset(dat.2020, !is.na(age)) #25.

#plot age-sero of naive vs. multi
dat.2019$age <- ceiling(dat.2019$age)


sum.yr.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
sum.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out[is.na(df.out)]<- 0
  df.out <- rbind(c(0,0), df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.out$Nage)
  df.out$cum_prop_cases <- cumsum(df.out$Nage)/sum(df.out$Nage)
  df.out$year = unique(df$year)
  return(df.out)
  
}
model.age.incidence.series <- function(par.dat, age_vect){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  #N_sero = rep(N.sero.guess, length=lts)
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  #pexposed= rep(NA, (max(age_vect)))
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  #pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  #pnaive= rep(NA, (max(age_vect)))
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #pnaiveone= rep(NA, (max(age_vect)))
  pnaiveone= rep(NA, (length(age_vect)-1))
  pnaiveone[1] <- 0
  pnaiveone = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  #pprim=as.list(rep(NA, (length(age_vect)*lts)))
  #pprim= rep(NA, (max(age_vect)))
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  #pmulti= rep(NA, (max(age_vect)))
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker[1] <- 1980
  year_tracker = rep(list(year_tracker),lts)
  
  year.start= min(par.dat$year)
  
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ #for-loops over all possible ages in our data range across all the years
      
      
      
      #remake the age vector for each year
      #print(a)
      age = age_vect_tmp[a]
      age_trunc = floor(age)
      age_current = age-age_trunc
      age_current[age_current==0] <- 1 #if you are at the end of the year, you spent the whole year here
      #age_class = ceiling(age)
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero#[2:nrow(par.dat)]
      #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
      lambda =  year.par$lambda#[2:nrow(par.dat)] #one per year
      
      #make a vector of durations for the lambdas across the time series
      if (age>1){
        dur = rep(1, length(lambda))  
      }else if(age<=1){
        dur = rep(0, length(lambda))  
      }
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>2){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      dur[length(dur)] <- age_current
      
      
      # #age should never exceed i under our new regime
      # dur1 <- age - i
      # if(dur1 <=0){
      #   dur[1] <- 0
      # }else{
      #   dur[1] <- dur[1] + dur1
      # }
      #dur1[dur1<0] <- 0
      
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      
      #if you are < 1 year in age, we still want it to be possible for you
      #to experience a multitypic infection (we know that cases are reported <1 year in
      #our natl data, and we make the assumption that all reported data signifies secondary
      #infections)
      
      #therefore, if you are <1 year, we have to account for this by breaking down the
      #foi across each quarter
      
      #if (age<1){
      #this just means, we split the current year into time before this timestep and now
      lambda <- c(lambda, lambda[i])
      N_sero <- c(N_sero, N_sero[i])
      dur <- c(dur, dur[i])
      
      #and we account for the current timestep accordingly
      dur[length(dur)-1] <-   dur[length(dur)]-.25
      dur[length(dur)] <- .25 #always just a quarter year in here
      
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      #}else if (age>=1){
      #still get the bonus of more time in this timstep
      # dur[i] <-  dur[i]
      #}
      
      #
      
      #now, we integrate the hazard of becoming infected over all the years in each
      
      #and, finally, because this is dengue, we also need to track the 
      #rate of getting exposed in this timestep, assuming you were naive at 
      #all timesteps before so we separate out the current and past timesteps 
      #in our integration
      
      
      
      #now we distinguish between exposures in the past and the present
      
      # #if this is the first step in the time series, 
      # we can't tell the difference, so we set everything to inte_pre 
      # #assuming that everything reported is a multitypic infection
      
      inte_pre = sum(dur[1:(length(dur)-1)]*lambda[1:(length(lambda)-1)]*N_sero[1:(length(N_sero)-1)])
      inte_now = sum(dur[length(dur)]*lambda[length(dur)]*N_sero[length(dur)])  
      
      #now sum them
      inte_all <-  inte_pre + inte_now
      
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_all))
      
      #this is the probability of being naive to all serotypes one year prior
      #(its the inverse of being exposed at one timestep prior)
      pnaiveone[[i]][[a]] <- 1-(1-exp(-inte_pre))
      
      
      # this is the probability of being naive to all serotypes up to one year prior AND 
      # getting infected with any serotype in this year
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting exposed this year AND the probability of 
      # NOT GETTING any serotypes at all up to now
      pprim[[i]][[a]] =  pnaiveone[[i]][[a]]*(1-exp(-inte_now))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive = c(unlist(pnaive)),
                             all_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  # #head(p.out)
  # p.add <- cbind.data.frame(year = seq((year.start), (year.start+lts-1), 1),
  #                           age=rep(0, length=lts), 
  #                           exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
  #                           multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  # p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases)) +facet_wrap(~year)
  p.out <- p.out[complete.cases(p.out),]
  #and split by age category
  p.out$age_trunc = ceiling(p.out$age)
  
  #p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  p.add <- cbind.data.frame(year = unique(p.sum$year), age = rep(0, length(unique(p.sum$year))), exposed=0, naive=1, all_prim=0, multi=0, sum_exp_naive=2, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  p.sum <- rbind(p.add, p.sum)
  p.sum <- arrange(p.sum, year, age)
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.last.yr <- function(par, par.dat, age_vect, dat){
  
  par.dat$lambda[length(par.dat$lambda)] <- exp(par)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  
  
  #now, select only the years for fitting for which there are data
  out.mod = subset(out.mod, year ==unique(dat$year))
  
  
  #plot(out.mod$cum_prop_cases, type="b")
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
fit.one.year.BFGS <- function(dat, par.dat, lambda.guess, N.sero.fix, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  dist.back <- max(dat$age[dat$year==min(dat$year)])#16
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  # ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #       geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # 
  
  
  
  
  #and fit it cumulatively
  
  #now test the next year with all 4 serotype assumptions
  log.lambda.guess <- log(lambda.guess)
  
  out.NS <- optim(par = log.lambda.guess, 
                  fn=log.lik.last.yr, 
                  method = "BFGS",
                  par.dat=par.dat, 
                  age_vect=age_vect, 
                  dat=df.out)
  
  
  par.dat$lambda_KP <- NA
  par.dat$lambda_KP[length(par.dat$lambda_KP)] <- exp(out.NS$par)
  par.dat$llik_KP <- NA
  par.dat$llik_KP <- out.NS$value
  par.dat$convergence_KP <- NA
  par.dat$convergence_KP <- out.NS$convergence
  
  
  
  #and return
  
  return(par.dat)
  
}

sum.dat.igg <- function(dat, age_vect){
df.sum <- ddply(dat, .(age, igg_res), summarise, Nage = length(age))
  df.out = cbind.data.frame(age=1:max(age_vect), igg_res = rep(c(0,1), each=max(age_vect)))
  df.out <- merge(df.out, df.sum, by=c("age", "igg_res"), all.x = T, sort = F)
  df.out[is.na(df.out)]<- 0
  df.out <- rbind(c(0,0,0), df.out)
  df.out <- rbind(c(0,1,0), df.out)
  df.out$igg_res <- as.factor(df.out$igg_res)
  df.out <- arrange(df.out, igg_res, age)
  #ggplot(dat = df.out) + geom_point(aes(x=age, y=Nage, color=igg_res))
  df.out$n_tot <- nrow(dat) #sum(df.out$Nage)
  df.out$proportion_by_age <- df.out$Nage/df.out$n_tot
  #ggplot(dat = df.out) + geom_point(aes(x=age, y=proportion_by_age*100, color=igg_res), size=3) + ylab("percent all cases by age")
  
  #then get the cumulative observed cases
  df.out.0 = df.out[df.out$igg_res==0,]
  df.out.1 = df.out[df.out$igg_res==1,]
  df.out.1$cum_cases <- cumsum(df.out.1$Nage)
  df.out.1$n <- sum(df.out.1$Nage)
  df.out.1$cum_prop_cases = df.out.1$cum_cases/df.out.1$n
  #plot(df.out.1$cum_prop_cases, type="b")
  
  #and add in the proportion of primary
  df.out.1$n_prim_cases <- df.out.0$Nage
  df.out.1$ntot <- df.out.0$n_tot
  
  df.out.1$prop_primary_by_age <- df.out.0$proportion_by_age
  
  df.out <- dplyr::select(df.out.1, age, n, cum_cases, cum_prop_cases, n_prim_cases, ntot, prop_primary_by_age)
  return(df.out)
}
sum.yr.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.out$Nage)
  df.out$cum_prop_cases <- cumsum(df.out$Nage)/sum(df.out$Nage)
  
  return(df.out)
  
}
fit.one.yr.only.BFGS <- function(dat, par.dat, lambda.guess, N.sero.fix, age_vect, fit.CI){
  
  
  #first, prep the data
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  df.out = sum.dat.igg(dat = dat, age_vect = age_vect_year)
  df.out$year <- max(dat$year)
  
  #head(df.out)
  #ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #      geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  #we are only guessing the one lambda!
  
  #now test the next year with all 4 serotype assumptions
  
  
  out.NS <- optim(par = log(lambda.guess), 
                  fn=log.lik.fit.two.ts, 
                  method = "BFGS",
                  par.dat=par.dat, 
                  age_vect=age_vect, 
                  dat=df.out, hessian = T)
  
  par.dat$lambda_KP <- NA
  par.dat$lambda_KP[length(par.dat$lambda_KP)] <- exp(out.NS$par)
  par.dat$llik_KP <- NA
  par.dat$llik_KP <- out.NS$value
  par.dat$convergence_KP <- NA
  par.dat$convergence_KP <- out.NS$convergence
  
  #and profile likelihood for CI
  
  if(fit.CI==TRUE){
    lambda.test = log(seq(0.0001, 1, length = 100))  
    llik.list <- list()
      for (i in 1:length(lambda.test)){
        
        llik.list[[i]] <- log.lik.fit.one.yr(par=lambda.test[i], par.dat=par.dat, age_vect=age_vect, dat=df.out)
      }
     
    dat.CI <- cbind.data.frame(lambda = exp(lambda.test), neg_llik=c(unlist(llik.list)))
    
    
    dat.CI <- subset(dat.CI, neg_llik<Inf)
    tmp=smooth.spline(dat.CI$lambda, dat.CI$neg_llik) 
    new=seq(0.0001, 1, by=.001) 
    interp = predict(tmp, new)$y
    #and then use the fact that the profile likelihood is χ2-distributed to erect 95% confidence intervals:
    mle1=new[which.min(interp)]
    tmp3=(predict(tmp, new)$y-min(predict(tmp, new)$y))-qchisq(0.95,1)
    CI <- range(new[tmp3<0])
    dat.CI$lambda[dat.CI$neg_llik==min(dat.CI$neg_llik)]
     
    par.dat$lambda_KP_lci <- NA
    par.dat$lambda_KP_uci <- NA
    par.dat$lambda_KP_lci[length(par.dat$lambda_KP)] <-CI[1]
    par.dat$lambda_KP_uci[length(par.dat$lambda_KP)] <-CI[2]
  }
  
  ret.par <- par.dat[complete.cases(par.dat),]
  
  #and return
  
  return(ret.par)
  
}
fit.two.ts.one.yr.BFGS <- function(dat, par.dat, lambda.guess, N.sero.fix, age_vect, fit.CI){
  
  
  #first, prep the data
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  df.out = sum.dat.igg(dat = dat, age_vect = age_vect_year)
  df.out$year <- max(dat$year)
  
  # #head(df.out)
  #  ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #        geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # # 
  # #we are only guessing the one lambda!
  
  #now test the next year with all 4 serotype assumptions
  par.dat.hold = par.dat
  
  out.NS <- optim(par = log(lambda.guess), 
                  fn=log.lik.fit.two.ts, 
                  method = "BFGS",
                  par.dat=par.dat, 
                  age_vect=age_vect, 
                  dat=df.out, hessian = T)
  
  par.dat$lambda_KP <- NA
  par.dat$lambda_KP[length(par.dat$lambda_KP)] <- exp(out.NS$par)
  par.dat$llik_KP <- NA
  par.dat$llik_KP <- out.NS$value
  par.dat$convergence_KP <- NA
  par.dat$convergence_KP <- out.NS$convergence
  
  #and profile likelihood for CI
  
  if(fit.CI==TRUE){
    lambda.test = log(seq(0.0001, 1, length = 100))  
    llik.list <- list()
    for (i in 1:length(lambda.test)){
      
      llik.list[[i]] <- log.lik.fit.two.ts(par=lambda.test[i], par.dat=par.dat, age_vect=age_vect, dat=df.out)
    }
    
    dat.CI <- cbind.data.frame(lambda = exp(lambda.test), neg_llik=c(unlist(llik.list)))
    
    
    dat.CI <- subset(dat.CI, neg_llik<Inf)
    tmp=smooth.spline(dat.CI$lambda, dat.CI$neg_llik) 
    new=seq(0.0001, 1, by=.001) 
    interp = predict(tmp, new)$y
    #and then use the fact that the profile likelihood is χ2-distributed to erect 95% confidence intervals:
    mle1=new[which.min(interp)]
    tmp3=(predict(tmp, new)$y-min(predict(tmp, new)$y))-qchisq(0.95,1)
    CI <- range(new[tmp3<0])
    dat.CI$lambda[dat.CI$neg_llik==min(dat.CI$neg_llik)]
    
    par.dat$lambda_KP_lci <- NA
    par.dat$lambda_KP_uci <- NA
    par.dat$lambda_KP_lci[length(par.dat$lambda_KP)] <-CI[1]
    par.dat$lambda_KP_uci[length(par.dat$lambda_KP)] <-CI[2]
  }
  
  par.dat$lambda[length(par.dat$lambda)] <- par.dat.hold$lambda[length(par.dat.hold$lambda)] 
  par.dat$lambda_min[length(par.dat$lambda_min)] <- par.dat.hold$lambda_min[length(par.dat.hold$lambda_min)] 
  par.dat$lambda_max[length(par.dat$lambda_max)] <- par.dat.hold$lambda_max[length(par.dat.hold$lambda_max)] 
  
  ret.par <- par.dat[complete.cases(par.dat),]
  
  #and return
  
  return(ret.par)
  
}
model.age.incidence.series <- function(par.dat, age_vect){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  #N_sero = rep(N.sero.guess, length=lts)
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  #pexposed= rep(NA, (max(age_vect)))
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  #pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  #pnaive= rep(NA, (max(age_vect)))
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #pnaiveone= rep(NA, (max(age_vect)))
  pnaiveone= rep(NA, (length(age_vect)-1))
  pnaiveone[1] <- 0
  pnaiveone = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  #pprim=as.list(rep(NA, (length(age_vect)*lts)))
  #pprim= rep(NA, (max(age_vect)))
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  #pmulti= rep(NA, (max(age_vect)))
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker[1] <- 1980
  year_tracker = rep(list(year_tracker),lts)
  
  year.start= min(par.dat$year)
  
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ #for-loops over all possible ages in our data range across all the years
      
      
      
      #remake the age vector for each year
      #print(a)
      age = age_vect_tmp[a]
      age_trunc = floor(age)
      age_current = age-age_trunc
      age_current[age_current==0] <- 1 #if you are at the end of the year, you spent the whole year here
      #age_class = ceiling(age)
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero#[2:nrow(par.dat)]
      #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
      lambda =  year.par$lambda#[2:nrow(par.dat)] #one per year
      
      #make a vector of durations for the lambdas across the time series
      if (age>1){
        dur = rep(1, length(lambda))  
      }else if(age<=1){
        dur = rep(0, length(lambda))  
      }
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>2){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      dur[length(dur)] <- age_current
      
      
      # #age should never exceed i under our new regime
      # dur1 <- age - i
      # if(dur1 <=0){
      #   dur[1] <- 0
      # }else{
      #   dur[1] <- dur[1] + dur1
      # }
      #dur1[dur1<0] <- 0
      
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      
      #if you are < 1 year in age, we still want it to be possible for you
      #to experience a multitypic infection (we know that cases are reported <1 year in
      #our natl data, and we make the assumption that all reported data signifies secondary
      #infections)
      
      #therefore, if you are <1 year, we have to account for this by breaking down the
      #foi across each quarter
      
      #if (age<1){
      #this just means, we split the current year into time before this timestep and now
      lambda <- c(lambda, lambda[i])
      N_sero <- c(N_sero, N_sero[i])
      dur <- c(dur, dur[i])
      
      #and we account for the current timestep accordingly
      dur[length(dur)-1] <-   dur[length(dur)]-.25
      dur[length(dur)] <- .25 #always just a quarter year in here
      
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      #}else if (age>=1){
      #still get the bonus of more time in this timstep
      # dur[i] <-  dur[i]
      #}
      
      #
      
      #now, we integrate the hazard of becoming infected over all the years in each
      
      #and, finally, because this is dengue, we also need to track the 
      #rate of getting exposed in this timestep, assuming you were naive at 
      #all timesteps before so we separate out the current and past timesteps 
      #in our integration
      
      
      
      #now we distinguish between exposures in the past and the present
      
      # #if this is the first step in the time series, 
      # we can't tell the difference, so we set everything to inte_pre 
      # #assuming that everything reported is a multitypic infection
      
      inte_pre = sum(dur[1:(length(dur)-1)]*lambda[1:(length(lambda)-1)]*N_sero[1:(length(N_sero)-1)])
      inte_now = sum(dur[length(dur)]*lambda[length(dur)]*N_sero[length(dur)])  
      
      #now sum them
      inte_all <-  inte_pre + inte_now
      
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_all))
      
      #this is the probability of being naive to all serotypes one year prior
      #(its the inverse of being exposed at one timestep prior)
      pnaiveone[[i]][[a]] <- 1-(1-exp(-inte_pre))
      
      
      # this is the probability of being naive to all serotypes up to one year prior AND 
      # getting infected with any serotype in this year
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting exposed this year AND the probability of 
      # NOT GETTING any serotypes at all up to now
      pprim[[i]][[a]] =  pnaiveone[[i]][[a]]*(1-exp(-inte_now))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive = c(unlist(pnaive)),
                             all_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  # #head(p.out)
  # p.add <- cbind.data.frame(year = seq((year.start), (year.start+lts-1), 1),
  #                           age=rep(0, length=lts), 
  #                           exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
  #                           multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  # p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases)) +facet_wrap(~year)
  p.out <- p.out[complete.cases(p.out),]
  #and split by age category
  p.out$age_trunc = ceiling(p.out$age)
  
  #p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  p.add <- cbind.data.frame(year = unique(p.sum$year), age = rep(0, length(unique(p.sum$year))), exposed=0, naive=1, all_prim=0, multi=0, sum_exp_naive=2, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  p.sum <- rbind(p.add, p.sum)
  p.sum <- arrange(p.sum, year, age)
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.fit.two.ts <- function(par, par.dat, age_vect, dat){
  
  
  par.dat$lambda[length(par.dat$lambda)] <- exp(par)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= all_prim)) + facet_wrap(~year)
  
  #now, select only the years for fitting for which there are data
  out.mod = subset(out.mod, year ==min(dat$year))
  
  
  #plot(out.mod$cum_prop_cases, type="b")
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
    ll=ll+dbinom(dat$n_prim_cases[i],dat$ntot[i],p=out.mod$all_prim[i],log=T) #and prim
  }
  
  return(-ll)
}
log.lik.fit.one.yr <- function(par, par.dat, age_vect, dat){
  
  
  par.dat$lambda[length(par.dat$lambda)] <- exp(par)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  #now, select only the years for fitting for which there are data
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  #plot(out.mod$cum_prop_cases, type="b")
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
plot.model.series.data <- function(par.dat, age_vect, dat){
  
  
  #first, prep the data
  dat = subset(dat, year<=max(par.dat$year))
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
    geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, age_vect = age_vect)
  head(out.mod)
  
  out.mod <- subset(out.mod, year >= min(df.out$year))
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  
  df.all<- rbind(df.mod, df.dat)
  
  
  p1 <- ggplot(data=df.all) + geom_point(aes(x=age, y=cum_prop_cases, color=type))  + 
    geom_line(aes(x=age, y=cum_prop_cases, color=type)) +facet_wrap(~year)
  print(p1)
  
  
  # 
  # p1 <- ggplot(data=subset(df.all, year==year.in)) + geom_point(aes(x=age, y=cum_prop_cases, color=type))  + 
  #   geom_line(aes(x=age, y=cum_prop_cases, color=type)) +facet_wrap(~year)
  # print(p1)
  
} 


#load the fitted data from the national level and estimate just the 2019 results
par.dat <- read.csv(file = paste0(homewd, "/data/foi-fit-national.csv"), stringsAsFactors = F, header = T)
head(par.dat)


fit.KP.2019 <- fit.two.ts.one.yr.BFGS(dat=dat.2019,
                                    par.dat = subset(par.dat, year<=2019),
                                    lambda.guess=c(0.01),
                                    N.sero.fix=2,
                                    age_vect=seq(0,22, by=1/4),
                                    fit.CI = T)
#now write data
write.csv(fit.KP.2019, file = paste0(homewd, "/data/foi-kampong-speu-2019.csv"), row.names = F)

