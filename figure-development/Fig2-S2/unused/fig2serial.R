

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#compare N serotype hypothesis across all years to cumulative case data
#functions
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


#load the age structured national data
dat <- read.csv(file = paste0(homewd, "/data/DENV-Nat-Aged.csv") , header = T, stringsAsFactors = F)
# dat
dat %>% filter(epiwk=='2001-12-31')
dat$epiwk <- gsub('2001-12-31','2002-01-01',dat$epiwk)
dat %>% filter(epiwk=='2001-12-31')
dat <- arrange(dat, date, age)

unique(dat$age) #round to years
#dat$age <- round(dat$age, 0)
dat$age <- ceiling(dat$age)

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
  
  year.start= min(par.dat$year)
  
  
  
  for (i in 1:lts){
    for(a in 2:(length(age_vect))){ #for-loops over all possible ages in our data range across all the years
      

      age = age_vect[a]
      age_trunc = floor(age)
      age_current = age-age_trunc
      age_current[age_current==0] <- 1 #if you are at the end of the year, you spent the whole year here
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero#[2:nrow(par.dat)]
      #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
      lambda =  year.par$lambda#[2:nrow(par.dat)] #one per year
      
      #print(i)
      #print(age)
      
      
      #then, ask, was the kid alive before the time series started?
      #so this is the duration of time prior to the first instance of data
      #we assume that 2002
      dur1 <- age - i
      dur1 <- ceiling(dur1)
      
      #if negative, this kid was born after the time series started and spent no time in 
      #the pre-data classes
      if(dur1<=0){
      #if the kid spent no time in prior age classes
      #then make your duration vector reflect the time series only
        dur = rep(1, length(lambda))
        dur[length(dur)] <- age_current #you may not have completed the whole year here yet
      
      }else if(dur1>0){
        dur = rep(1, ((length(lambda)) +dur1))
        dur[length(dur)] <- age_current
      }
      #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      
      #now, make sure there is a lambda for each of the 
      
      #if you are < 1 year in age, we still want it to be possible for you
      #to experience a multitypic infection (we know that cases are reported <1 year in
      #our natl data, and we make the assumption that all reported data signifies secondary
      #infections)
      
      #therefore, if you are <1 year, we have to account for this by breaking down the
      #foi across each quarter
      
      if (age<1){
        #this just means, we split the current year into time before this timestep and now
        lambda <- c(lambda, lambda[i])
        N_sero <- c(N_sero, N_sero[i])
        dur <- c(dur, dur[i])
        
        #and we account for the current timestep accordingly
        dur[length(dur)] <- .25 #always just a quarter year in here
        dur[length(dur)-1] <-   age_current-.25
      }else if (age>=1){
        #still get the bonus of more time in this timstep
        dur[i] <- age_current + dur[i]
      }
      
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
      
    }
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
                           age=rep(age_vect, lts), 
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
  
  #and split by age category
  p.out$age_trunc = ceiling(p.out$age)
  
  #p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
   geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==2002)) + geom_point(aes(x=age, y=multi))  + 
   #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year
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
plot.model.series.data.oneyr <- function(par.dat, age_vect, dat, year.in){
  
  
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
  
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  
  df.all<- rbind(df.mod, df.dat)
  
  
  # p1 <- ggplot(data=df.all) + geom_point(aes(x=age, y=cum_prop_cases, color=type))  + 
  #   geom_line(aes(x=age, y=cum_prop_cases, color=type)) +facet_wrap(~year)
  # print(p1)
  
  
  
  p1 <- ggplot(data=subset(df.all, year==year.in)) + geom_point(aes(x=age, y=cum_prop_cases, color=type))  + 
    geom_line(aes(x=age, y=cum_prop_cases, color=type)) +facet_wrap(~year)
  print(p1)
  
} 

#this fits one year of model output to one year of data:
log.lik.prop.yearly <- function(par, par.dat, age_vect, dat, year.now){
  
  par.dat$lambda[par.dat$year==year.now] <- exp(par)
  
  par.dat = subset(par.dat, year <=year.now)
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  #with(out.mod, plot(age, cum_prop_cases, type="b"))
  
  #and subset to just this one year for comparison
  
  dat.fit <- subset(dat, year==year.now)
  #with(dat.fit, plot(age, cum_prop_cases, type="b"))
  
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat.fit$age))
  
  out.mod <- subset(out.mod, year==year.now)
  
  #plot(out.mod$cum_prop_cases, type="b")
  # # # # 
  #how likely are the data, given the model as truth?
  
  ll=0
  #dat.fit$cum_fake <- round(dat.fit$n*out.mod$cum_prop_cases,0)
  for (i in 1:length(dat.fit$age)){
    ll=ll+dbinom(dat.fit$cum_cases[i],size=dat.fit$n[i],prob=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
yearly.RSS <- function(par, par.dat, age_vect, dat, year.now){
  
  par.dat$lambda[par.dat$year==year.now] <- exp(par)
  
  par.dat = subset(par.dat, year <=year.now)
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  #with(out.mod, plot(age, cum_prop_cases, type="b"))
  
  #and subset to just this one year for comparison
  
  dat.fit <- subset(dat, year==year.now)
  #with(dat.fit, plot(age, cum_prop_cases, type="b"))
  
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat.fit$age))
  
  out.mod <- subset(out.mod, year==year.now)
  
  #and get the RSS
  dat.fit$diff_sq <-  (dat.fit$cum_prop_cases - out.mod$cum_prop_cases)^2
  
  RSS <- sum(dat.fit$diff_sq)
  
  return(RSS)
}
#this fits across all the data all together - but does so very poorly
log.lik.prop.cumulative <- function(par, par.dat, age_vect, dat, year.now){
  
  par.dat$lambda[par.dat$year==year.now] <- exp(par)
  
  par.dat = subset(par.dat, year <=year.now)
  out.mod <- model.age.incidence.series(par.dat = par.dat, 
                                        age_vect=age_vect)  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  dat.fit <- subset(dat, year<=year.now)
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat.fit$age))
  #plot(out.mod$cum_prop_cases, type="b")
  out.mod <- arrange(out.mod, year, age)
  dat.fit <- arrange(dat.fit, year, age)
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat.fit$age)){
    ll=ll+pbinom(dat.fit$cum_cases[i],dat.fit$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}

#this one fits each year in sequence and optimizes the current year foi 
#to the current year data only
fit.all.yrs.seq.yr <- function(dat, lambda.guess, N.sero.vect, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  dist.back <- max(dat$age[dat$year==min(dat$year)])#22
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
   #ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
    #       geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  
  
  #make your guess parameters
  par.dat <- cbind.data.frame(year=min(dat$year):max(dat$year),
                              lambda = rep(lambda.guess, length(unique(dat$year))),
                              N_sero = N.sero.vect)

  #then also make a list to hold the comparisons
  fit.par <- list(rep(NA, length(unique(dat$year))))
  
  for (i in 1:nrow(par.dat)){
  
  #now test the next year with all 4 serotype assumptions
  log.lambda.guess <- log(par.dat$lambda[i])
  par.dat$N_sero[i] <- 1
  out.NS1 <- optim(par = log.lambda.guess, 
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  log.lambda.guess <- log(par.dat$lambda[i]/2)
  par.dat$N_sero[i] <- 2
  out.NS2 <- optim(par = log.lambda.guess, #guess the same as last time
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  
  log.lambda.guess <- log(par.dat$lambda[i]/3)
  par.dat$N_sero[i] <- 3
  out.NS3 <- optim(par = log.lambda.guess, #guess the same as last time
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  log.lambda.guess <- log(par.dat$lambda[i]/4)
  par.dat$N_sero[i] <- 4
  out.NS4 <- optim(par = log.lambda.guess, #guess the same as last time
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])

  
  #and combine them and compare
  fit.dat <- cbind.data.frame(year = rep(par.dat$year[i], 4),
                              lambda = exp(c(out.NS1$par,out.NS2$par, out.NS3$par, out.NS4$par)),
                              neg_loglik = c(out.NS1$value,out.NS2$value, out.NS3$value, out.NS4$value),
                              convergence = c(out.NS1$convergence,out.NS2$convergence, out.NS3$convergence, out.NS4$convergence),
                              N_sero = c(1,2,3,4))
                              
  fit.dat = subset(fit.dat, convergence ==0)
  fit.dat$AIC <- 2*fit.dat$neg_loglik+2*fit.dat$N_sero
  fit.dat$delta_AIC <- fit.dat$AIC-min(fit.dat$AIC)
  
  #now select the best fit and feed into next timestep
  best.fit = subset(fit.dat, AIC==min(fit.dat$AIC))
  
  #and store it in the par.dat
  par.dat$lambda[i] <- best.fit$lambda
  par.dat$N_sero[i] <- best.fit$N_sero
  
  fit.par[[i]] <- fit.dat
  
  
  }
  
  #and return
  
  return(list(par.dat, fit.par))

}
fit.all.yrs.seq.yr.extend <- function(dat, lambda.guess, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  dist.back <- max(dat$age[dat$year==min(dat$year)])#22
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  #ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #       geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  
  
  #make your guess parameters
  par.dat <- cbind.data.frame(year=min(dat$year):max(dat$year),
                              lambda = rep(lambda.guess, length(unique(dat$year))),
                              N_sero = N.sero.vect)
  
  #then also make a list to hold the comparisons
  fit.par <- list(rep(NA, length(unique(dat$year))))
  
  for (i in 1:nrow(par.dat)){
    
    #now test the next year with all 4 serotype assumptions
    log.lambda.guess <- log(par.dat$lambda[i])
    par.dat$N_sero[i] <- 1
    out.NS1 <- optim(par = log.lambda.guess, 
                     fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/2)
    par.dat$N_sero[i] <- 2
    out.NS2 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    log.lambda.guess <- log(par.dat$lambda[i]/3)
    par.dat$N_sero[i] <- 3
    out.NS3 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/4)
    par.dat$N_sero[i] <- 4
    out.NS4 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    #and combine them and compare
    fit.dat <- cbind.data.frame(year = rep(par.dat$year[i], 4),
                                lambda = exp(c(out.NS1$par,out.NS2$par, out.NS3$par, out.NS4$par)),
                                neg_loglik = c(out.NS1$value,out.NS2$value, out.NS3$value, out.NS4$value),
                                convergence = c(out.NS1$convergence,out.NS2$convergence, out.NS3$convergence, out.NS4$convergence),
                                N_sero = c(1,2,3,4))
    
    fit.dat = subset(fit.dat, convergence ==0)
    fit.dat$AIC <- 2*fit.dat$neg_loglik+2*fit.dat$N_sero
    fit.dat$delta_AIC <- fit.dat$AIC-min(fit.dat$AIC)
    
    #now select the best fit and feed into next timestep
    best.fit = subset(fit.dat, AIC==min(fit.dat$AIC))
    
    #and store it in the par.dat
    par.dat$lambda[i] <- best.fit$lambda
    par.dat$N_sero[i] <- best.fit$N_sero
    
    fit.par[[i]] <- fit.dat
    
    
  }
  
  #and return
  
  return(list(par.dat, fit.par))
  
}
fit.all.yrs.seq.fixsero <- function(dat, lambda.guess, N.sero.vect, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  #ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
   #      geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  #make your guess parameters
  par.dat <- cbind.data.frame(year=min(dat$year):max(dat$year),
                              lambda = rep(lambda.guess, length(unique(dat$year))),
                              N_sero = N.sero.vect,
                              llik = rep(NA, length(unique(dat$year))))
  
  #then also make a list to hold the comparisons
  lik.vect <- list(rep(NA, length(unique(dat$year))))
  
  for (i in 1:nrow(par.dat)){
    
    #now test the next year with all 4 serotype assumptions
    log.lambda.guess <- log(par.dat$lambda[i])
    
    out.NS <- optim(par = log.lambda.guess, 
                     fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    #and store it in the par.dat
    par.dat$lambda[i] <- exp(out.NS$par)
    
    
    par.dat$llik[i]  <- out.NS$value
    
    
  }
  
  #and compile
  
  
  return(par.dat)
  
}
fit.all.yrs.seq.yr.rss <- function(dat, lambda.guess, N.sero.vect, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  #ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #       geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  #make your guess parameters
  par.dat <- cbind.data.frame(year=min(dat$year):max(dat$year),
                              lambda = rep(lambda.guess, length(unique(dat$year))),
                              N_sero = N.sero.vect)
  
  #then also make a list to hold the comparisons
  fit.par <- list(rep(NA, length(unique(dat$year))))
  
  for (i in 1:nrow(par.dat)){
    
    #now test the next year with all 4 serotype assumptions
    log.lambda.guess <- log(par.dat$lambda[i])
    par.dat$N_sero[i] <- 1
    out.NS1 <- optim(par = log.lambda.guess, 
                     fn=yearly.RSS, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/2)
    par.dat$N_sero[i] <- 2
    out.NS2 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=yearly.RSS, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    log.lambda.guess <- log(par.dat$lambda[i]/3)
    par.dat$N_sero[i] <- 3
    out.NS3 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=yearly.RSS, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/4)
    par.dat$N_sero[i] <- 4
    out.NS4 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=yearly.RSS, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    #and combine them and compare
    fit.dat <- cbind.data.frame(year = rep(par.dat$year[i], 4),
                                lambda = exp(c(out.NS1$par,out.NS2$par, out.NS3$par, out.NS4$par)),
                                RSS = c(out.NS1$value,out.NS2$value, out.NS3$value, out.NS4$value),
                                convergence = c(out.NS1$convergence,out.NS2$convergence, out.NS3$convergence, out.NS4$convergence),
                                N_sero = c(1,2,3,4))
    
    fit.dat = subset(fit.dat, convergence ==0)
    #fit.dat$AIC <- 2*fit.dat$neg_loglik+2*fit.dat$N_sero
    #fit.dat$delta_AIC <- fit.dat$AIC-min(fit.dat$AIC)
    
    #now select the best fit and feed into next timestep
    #best.fit = subset(fit.dat, AIC==min(fit.dat$AIC))
    best.fit = subset(fit.dat, RSS==min(fit.dat$RSS))
    
    #and store it in the par.dat
    par.dat$lambda[i] <- best.fit$lambda
    par.dat$N_sero[i] <- best.fit$N_sero
    
    fit.par[[i]] <- fit.dat
    
    
  }
  
  #and return
  
  return(list(par.dat, fit.par))
  
}
fit.all.yrs.seq.cumulative <- function(dat, lambda.guess, N.sero.vect, age_vect){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
    geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  #make your guess parameters
  par.dat <- cbind.data.frame(year=min(dat$year):max(dat$year),
                              lambda = rep(lambda.guess, length(unique(dat$year))),
                              N_sero = N.sero.vect)
  
  #then also make a list to hold the comparions
  fit.par <- list(rep(NA, length(unique(dat$year))))
  
  for (i in 1:nrow(par.dat)){
    
    #now test the next year with all 4 serotype assumptions
    log.lambda.guess <- log(par.dat$lambda[i])
    par.dat$N_sero[i] <- 1
    out.NS1 <- optim(par = log.lambda.guess, 
                     fn=log.lik.prop.cumulative, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/2)
    par.dat$N_sero[i] <- 2
    out.NS2 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.cumulative, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    log.lambda.guess <- log(par.dat$lambda[i]/3)
    par.dat$N_sero[i] <- 3
    out.NS3 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.cumulative, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    log.lambda.guess <- log(par.dat$lambda[i]/4)
    par.dat$N_sero[i] <- 4
    out.NS4 <- optim(par = log.lambda.guess, #guess the same as last time
                     fn=log.lik.prop.cumulative, #not cumulative... we for-loop through and test year by year
                     method = "Brent",
                     lower = -100, upper=20,
                     par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
    
    
    #and combine them and compare
    fit.dat <- cbind.data.frame(year = rep(par.dat$year[i], 4),
                                lambda = exp(c(out.NS1$par,out.NS2$par, out.NS3$par, out.NS4$par)),
                                neg_loglik = c(out.NS1$value,out.NS2$value, out.NS3$value, out.NS4$value),
                                convergence = c(out.NS1$convergence,out.NS2$convergence, out.NS3$convergence, out.NS4$convergence),
                                N_sero = c(1,2,3,4))
    
    fit.dat = subset(fit.dat, convergence ==0)
    fit.dat$AIC <- 2*fit.dat$neg_loglik+2*fit.dat$N_sero
    fit.dat$delta_AIC <- fit.dat$AIC-min(fit.dat$AIC)
    
    #now select the best fit and feed into next timestep
    best.fit = subset(fit.dat, AIC==min(fit.dat$AIC))
    
    #and store it in the par.dat
    par.dat$lambda[i] <- best.fit$lambda
    par.dat$N_sero[i] <- best.fit$N_sero
    
    fit.par[[i]] <- fit.dat
    
    
  }
  
  #and return
  
  return(list(par.dat, fit.par))
  
}

par.dat.seq.yr <- fit.all.yrs.seq.yr(dat=dat,
                      lambda.guess = 0.09, 
                      #N.sero.vect=rep(2, length(unique(dat$year))),
                      age_vect = seq(0,22, by=1/4))
#and save values 
#save(par.dat.seq.yr, file =paste0(homewd, "/figure-development/Fig2-S1-S2/tmp-dat/par.dat.seq.yr.Rdata"))


par.dat.seq.fixsero <- fit.all.yrs.seq.fixsero(dat=dat,
                                     lambda.guess = 0.09, 
                                     N.sero.vect=rep(2, length(unique(dat$year))),
                                     age_vect = seq(0,35, by=1/4))
with(par.dat.seq.fixsero, plot(year, lambda, type="b"))

c(rep(1, 5), 2, rep(1,4), 2,rep(1,6), 2, 1)
par.dat.epi.invasion <- fit.all.yrs.seq.fixsero(dat=dat,
                                               lambda.guess = 0.09, 
                                               N.sero.vect=c(rep(1, 5), 2, rep(1,4), 2,rep(1,6), 2, 1),
                                               age_vect = seq(0,35, by=1/4))

with(par.dat.epi.invasion, plot(year, lambda, type="b"))

par.dat.seq.cumulative <- fit.all.yrs.seq.cumulative(dat=dat,
                                     lambda.guess = 0.09, 
                                     N.sero.vect=rep(2, length(unique(dat$year))),
                                     age_vect = seq(0,35, by=1/4))
#save(par.dat.seq.cumulative , file =paste0(homewd, "/figure-development/Fig2-S1-S2/tmp-dat/par.dat.seq.cumulative.Rdata"))



#and plot both with the data
plot.model.series.data(par.dat = par.dat.seq.yr[[1]], 
                       age_vect = seq(0,35, by=1/4), dat=dat)

plot.model.series.data.oneyr(par.dat = par.dat.seq.yr[[1]], 
                       age_vect = seq(0,35, by=1/4), dat=dat, year.in = 2016)

plot.model.series.data(par.dat = par.dat.seq.fixsero, 
                       age_vect = seq(0,35, by=1/4), dat=dat)

plot.model.series.data(par.dat = par.dat.epi.invasion, 
                       age_vect = seq(0,35, by=1/4), dat=dat)

plot.model.series.data(par.dat = par.dat.seq.cumulative[[1]],
                       age_vect = seq(0,35, by=1/4), dat=dat)

plot.model.series.data(par.dat = par.dat.seq.fixsero[[1]],
                       age_vect = seq(0,35, by=1/4), dat=dat)




tmp.dat <- par.dat.seq.yr[[1]]
tmp.dat$lambda <- tmp.dat$lambda[1]

plot.model.series.data(par.dat = tmp.dat,
                       age_vect = seq(0,35, by=1/4), dat=dat)



par.dat.RSS <- fit.all.yrs.seq.yr.rss(dat=dat,
                                     lambda.guess = 0.09, 
                                     N.sero.vect=rep(2, length(unique(dat$year))),
                                     age_vect = seq(0,35, by=1/4))


plot.model.series.data(par.dat = par.dat.RSS[[1]],
                       age_vect = seq(0,35, by=1/4), dat=dat)



