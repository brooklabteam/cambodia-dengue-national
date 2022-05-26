

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
year.split <- dlply(dat, .(year))

#get age frequency for all years
dat.sum <- lapply(year.split, sum.yr, age_vect = 1:max(dat$age))
dat.sum <- data.table::rbindlist(dat.sum)
dat.sum$year <- rep(unique(dat$year), each=(max(dat$age)+1))
dat.sum$year <- as.factor(dat.sum$year)

dat.sum$prop_prev <- dat.sum$Nage/dat.sum$n

#now pull out the first year
dat.start = subset(dat.sum, year==min(as.numeric(as.character(dat.sum$year))))

dat.sum = subset(dat.sum, year!=min(as.numeric(as.character(dat.sum$year))))


#get the best fit foi and N-sero for year one
model.age.incidence <- function(log.lambda, N_sero, age_vect){
  
  lambda= exp(as.numeric(log.lambda))
  pexposed=as.list(rep(NA,length(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,length(age_vect))) #naive
  pprim=as.list(rep(NA,length(age_vect))) #first infection with a single serotype
  pmulti=as.list(rep(NA,length(age_vect))) #secondary infection
  #prob_naive_one=as.list(rep(NA,max(age_vect))) #prob of being naive to all strains at one timestep prior to the current
  
  
  for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range
    
    age = age_vect[a]
    #index = which(age_vect==age)
    #age_trunc = trunc(age)
    
    #i gives you years since the beginning of the time series.
    #so this is the duration of time under the first foi (pre-2003) in years
    #dur1=(age_trunc)
    dur12=age #you assume foi before this year was the same as now
    #all years in between year 1 of the time series and current time had a duration of 1 year
    
    #then the duration of time in the current foi, also in years
    #dur2 = age - age_trunc
    
    #and the instantaneous hazard of infection, since cases are reported in weeks
    dur3 = 1/52 #here also in years
    
    
    #here is the rate of exposure before the beginning of our time series
    #for some, this will be 0 if they were born during the time series
    inte_exposed_all_pre = dur12*lambda*N_sero
    
    
    
    #rate of getting exposed to any given serotype in the current timestep
    #basically instantaneous but the case reports are in weeks
    inte_now_immediate = lambda*N_sero*dur3
    
    ########## now to sum ########## 
    
    #and this is the rate of being exposed to all serotypes at this point in time
    inte_exposed_all <- inte_exposed_all_pre +  inte_now_immediate
    
    
    #this is the probability of being exposed to any serotype at the current time
    pexposed[[a]] <- (1-exp(-inte_exposed_all))
    
    #this is the probability of being naive to all serotypes at this point in time
    pnaive[[a]] <- 1-(1-exp(-inte_exposed_all))
    
    
    # this is the probability of being naive to all serotypes up to this timepoint AND 
    # getting infected with a single serotype (any serotype) now in the immediate timestep
    # (e.g. this is the probability of a primary infection)
    # this will be the combined probability of getting serotype X NOW AND the probability of 
    # NOT GETTING any serotypes at all up to now
    
    pprim[[a]] =  pnaive[[a]]*(1-exp(-inte_now_immediate))
    
    
    #if not primarily infected or naive, this should be a multitypic infection
    #but remember it could be a primary infection with any of the four serotypes - that
    # is already captured above by multiplying the rate of infection by the serotype number
    pmulti[[a]] = 1 - pnaive[[a]] - (pprim[[a]])
    
    
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(age=age_vect, 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_prim=c(unlist(pprim)*N_sero), 
                           multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  p.add <- cbind.data.frame(age=0, exposed=0, naive=1, all_prim=0, 
                            multi=0, sum_exp_naive=1, sum_naive_prim_multi=1)
  p.out <- rbind(p.add, p.out)
  
  
  
  ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + geom_line(aes(x=age, y=multi))
  # #need to cap when the proportion reaches 1
  #  if(max(p.out$multitypic_inf>=1)){
  #    first.one <- min(which(p.out$multitypic_inf>=1))
  #    p.out = p.out[1:first.one,]
  #  }
  # 
  #plot(ecdf(p.out$multitypic_inf))
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  
  p.out$cum_prop_cases <- p.out$multi
  
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=all_primary_inf))  + geom_line(aes(x=age, y=all_primary_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multitypic_inf))  + geom_line(aes(x=age, y=multitypic_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases))
  
  return(p.out) #returns prevalence by age for the length of all possible ages in our dataset
}
model.age.one <- function(log.lambda, age_vect){
  
  lambda= exp(as.numeric(log.lambda))
  pexposed=as.list(rep(NA,length(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,length(age_vect))) #naive
  #pprim=as.list(rep(NA,max(age_vect))) #first infection with a single serotype
  #pmulti=as.list(rep(NA,max(age_vect))) #secondary infection
  #prob_naive_one=as.list(rep(NA,max(age_vect))) #prob of being naive to all strains at one timestep prior to the current
  for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range
    
    age = age_vect[a]
    #index = which(age_vect==age)
    #age_trunc = trunc(age)
    
    #i gives you years since the beginning of the time series.
    #so this is the duration of time under the first foi (pre-2003) in years
    #dur1=(age_trunc)
    dur12=age #you assume foi before this year was the same as now
    #all years in between year 1 of the time series and current time had a duration of 1 year
    
    #then the duration of time in the current foi, also in years
    #dur2 = age - age_trunc
    
    #and the instantaneous hazard of infection, since cases are reported in weeks
    dur3 = 1/52 #here also in years
    
    
    #here is the rate of exposure before the beginning of our time series
    #for some, this will be 0 if they were born during the time series
    inte_exposed_all_pre = dur12*lambda
    
    inte_now_immediate = lambda*dur3
    
    ########## now to sum ########## 
    
    #and this is the rate of being exposed to all serotypes at this point in time
    inte_exposed_all <- inte_exposed_all_pre +  inte_now_immediate
    
    
    
    pexposed[[a]] <- (1-exp(-inte_exposed_all))
    pnaive[[a]] <- 1-(1-exp(-inte_exposed_all))
    
    #if not primarily infected or naive, you should be a multitypic infection
    #we could fit to that generally but we have serotype-specific info that can 
    #fill in instead
    
    #for it to be a multitypic infection with D1, you need to multiply
    #prob exposure to D2 up to a-1 x  prob naive to D1 up to a-1 x prob to D1 now
    # pmulti[[a]] = 1 - pnaive[[a]] - (N_sero*pprim[[a]])
    
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(age=age_vect, 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)))
  #all_prim=c(unlist(pprim)*N_sero), 
  #multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  #p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  p.add <- cbind.data.frame(age=0, exposed=0, naive=1, 
                            sum_exp_naive=1)
  p.out <- rbind(p.add, p.out)
  
  
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + geom_line(aes(x=age, y=multi))
  # #need to cap when the proportion reaches 1
  #  if(max(p.out$multitypic_inf>=1)){
  #    first.one <- min(which(p.out$multitypic_inf>=1))
  #    p.out = p.out[1:first.one,]
  #  }
  # 
  #plot(ecdf(p.out$multitypic_inf))
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  
  p.out$cum_prop_cases <- p.out$exposed
  
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=all_primary_inf))  + geom_line(aes(x=age, y=all_primary_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multitypic_inf))  + geom_line(aes(x=age, y=multitypic_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases))
  
  return(p.out) #returns prevalence by age for the length of all possible ages in our dataset
}
log.lik.prop <- function(par, age_vect, dat, N_sero){
  
  #first run the model with the specified par
  if(N_sero==1){
    
    out.mod <- model.age.one(log.lambda = par, age_vect = age_vect)  
    
  }else{
    out.mod <- model.age.incidence(log.lambda = par, age_vect = age_vect,  N_sero=N_sero)  
  }
  
  
  
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat$age))
  #plot(out.mod$cum_prop_cases, type="b")
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
wrap.all.llik.multi <- function(lambda.guess, age_vect, dat, N_sero){
  
  #first, prep the data
  df.out <- sum.yr(dat, age_vect = age_vect)
  #head(df.out)
  
  #  
  out.fit <- optim(log(lambda.guess), 
                   fn=log.lik.prop, method = "Brent",
                   lower = -20, upper=20,
                   age_vect=age_vect,
                   dat=df.out,
                   N_sero=N_sero, hessian = T)
  
  
  
  fisher_info <- solve(out.fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  upper <- exp(out.fit$par+1.96*prop_sigma)
  lower <- exp(out.fit$par-1.96*prop_sigma)
  
  
  out.fit <- cbind.data.frame(lambda_fit = exp(out.fit$par), 
                              lambda_lower=lower,
                              lambda_upper=upper,
                              neg_llik=out.fit$value,
                              convergence=out.fit$convergence)#, 
  
  
  return(out.fit)
}
wrap.all.llik <- function(lambda.guess, age_vect, dat, N_sero){
  
  #first, prep the data
  df.out <- sum.yr(dat, age_vect = age_vect)
  
  
  #  
  out.fit <- optim(log(lambda.guess), fn=log.lik.prop, method = "Brent",
                   lower = -20, upper=20,
                   age_vect=age_vect, 
                   year=unique(dat$year),
                   dat=df.out,
                   N_sero=N_sero, hessian = T)
  
  
  
  fisher_info <- solve(out.fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  upper <- exp(out.fit$par+1.96*prop_sigma)
  lower <- exp(out.fit$par-1.96*prop_sigma)
  
  
  out.fit <- cbind.data.frame(year=unique(dat$year),
                              lambda_fit = exp(out.fit$par), 
                              lambda_lower=lower,
                              lambda_upper=upper,
                              sum_sq=out.fit$value,
                              convergence=out.fit$convergence)#, 
  
  
  return(out.fit)
}
run.plot.single <- function(mod.dat, dat, age_vect, N_sero){
  
  #and run over many years
  out.run <- model.age.incidence(log.lambda = log(mod.dat$lambda_fit), 
                                 age_vect=age_vect, N_sero=N_sero)
  
  out.run.lci <- model.age.incidence(log.lambda = log(mod.dat$lambda_lower), 
                                     age_vect=age_vect, N_sero=N_sero)
  
  out.run.uci <- model.age.incidence(log.lambda = log(mod.dat$lambda_upper), 
                                     age_vect=age_vect, N_sero=N_sero)
  
  #out.run <- data.table::rbindlist(out.run)
  
  dat.mod <- dplyr::select(out.run,  age,  cum_prop_cases)
  dat.mod$lci = out.run.lci$cum_prop_cases
  dat.mod$uci = out.run.uci$cum_prop_cases
  dat.mod$type<- "model"
  names(dat.mod) <- c("age", "prop", "lci", "uci", "type")
  
  
  
  #now join with data and plot
  
  dat.out <- sum.yr(df=dat,age_vect = age_vect)
  #dat.out$year <- unique(dat$year)
  dat.dat <- dplyr::select(dat.out,age, cum_prop_cases)
  dat.dat$lci <- dat.dat$uci <- NA
  dat.dat$type<- "data"
  names(dat.dat) <-  c( "age", "prop", "lci", "uci", "type")
  dat.all <- rbind(dat.dat, dat.mod)
  
  p1 <- ggplot(data=dat.all) + 
    geom_point(aes(x=age, y=prop, color=type)) +
    geom_line(aes(x=age, y=prop, color=type, linetype=type)) +
    geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type), alpha=.3)+
    xlim(0,max(age_vect)) + ylab("cumulative proportion of cases") +theme_bw()+
    theme(panel.grid = element_blank())
  
  print(p1)
  
  return(dat.all)
  
  
}
comp.all.NS.multi <- function(lambda.guess, age_vect, dat){
  
  out.NS1 <- wrap.all.llik.multi(lambda.guess = lambda.guess,
                                 age_vect = age_vect,
                                 dat=dat,
                                 N_sero=1)
  
  out.NS2 <- wrap.all.llik.multi(lambda.guess = lambda.guess,
                                 age_vect = age_vect,
                                 dat=dat,
                                 N_sero=2)
  
  out.NS3 <- wrap.all.llik.multi(lambda.guess = lambda.guess,
                                 age_vect = age_vect,
                                 dat=dat,
                                 N_sero=3)
  
  out.NS4 <- wrap.all.llik.multi(lambda.guess = lambda.guess,
                                 age_vect = age_vect,
                                 dat=dat,
                                 N_sero=4)
  
  
  out.NS1$N_sero=1
  out.NS2$N_sero=2
  out.NS3$N_sero=3
  out.NS4$N_sero=4
  
  out.all <- rbind(out.NS1, out.NS2, out.NS3, out.NS4)
  out.all$year <- unique(dat$year)
  out.all$AIC <- 2*(out.all$neg_llik) + 2*out.all$N_sero
  return(out.all)
  
  
}



#try for one year
out.2002 <- comp.all.NS.multi(lambda.guess = c(.01),
                              age_vect = 1:35,
                               dat=year.split[[1]])
out.2002[out.2002$AIC==min(out.2002$AIC),]
out.2002$deltaAIC = out.2002$AIC-min(out.2002$AIC)
#2 serotypes is best with lambda = 0.08904973  if you include age within the year
#2 serotypes is best with lambda = 0.09263662  if you take ceiling age

#plot model with data to be sure
run.plot.single(mod.dat=out.2002[out.2002$AIC==min(out.2002$AIC),], 
                dat=year.split[[1]], 
                age_vect = 1:35, 
                N_sero=out.2002$N_sero[out.2002$AIC==min(out.2002$AIC)])


model.age.incidence.cumulative <- function(log.lambda.start, log.lambda.guess, N.sero.start, N.sero.vector, age_vect, lts, year.start){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = exp(as.numeric(log.lambda.start))
  #N_sero = rep(N.sero.guess, length=lts)
  N_sero = N.sero.vector
  #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
  lambda = exp(log.lambda.guess)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  pexposed= rep(NA, (length(age_vect)))
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  pnaive= rep(NA, (length(age_vect)))
  pnaive = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  pprim=as.list(rep(NA, (length(age_vect)*lts)))
  pprim= rep(NA, (length(age_vect)))
  pprim = rep(list(pexposed),lts)
  
  pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  pmulti= rep(NA, (length(age_vect)))
  pmulti = rep(list(pexposed),lts)
  
  
  
  for (i in 1:lts){
  for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range across all the years
    
    
    #print(i)
    #print(a)
    
    age = age_vect[a]
    index = which(age_vect==age)
    age_trunc = trunc(age)
    
    #i gives you years since the beginning of the time series.
    #so this is the duration of time under the first foi (pre-2003) in years
    dur1=(age_trunc-i)
    dur1[dur1<0]<- 0
    
    
    #all years in between year 1 of the time series and current time had a duration of 1 year
    
    #then the duration of time in the current foi, also in years
    dur2 = age - age_trunc
    
    #and the instanteous hazard of infection, since cases are reported in weeks
    dur3 = 1/52 #here also in years
    
    #here is the rate of exposure before the beginning of our time series
    #for some, this will be 0 if they were born during the time series
    inte_exposed_all_pre_lts = dur1*lambda.start*N.sero.start 
    
    #here is the rate of exposure to all serotypes up to the current period
    
    lambda.trunc <- lambda[1:(i-1)]
    Nsero.trunc <- N_sero[1:(i-1)]
    
    
    #here is the rate of exposure to any serotype within the time series but before the current year
    trunc.foi =  lambda.trunc*Nsero.trunc
    inte_exposed_lts = sum(trunc.foi)
    
    #and the rate of getting exposed to any of all serotypes in the current year
    inte_exposed_current = dur2*N_sero[i]*lambda[i]
    
    
    #rate of getting exposed to any given serotype in the current timestep
    #basically instantaneous but the case reports are in weeks
    inte_now_immediate = lambda[i]*N_sero[i]*dur3
    
    ########## now to sum ########## 
    
    #and this is the rate of being exposed to all serotypes at this point in time
    inte_exposed_all <- inte_exposed_all_pre_lts + inte_exposed_lts + inte_exposed_current + inte_now_immediate
  
    #this is the probability of being exposed to any serotype at the current time
    pexposed[[i]][[a]] <- (1-exp(-inte_exposed_all))
    
    #this is the probability of being naive to all serotypes at this point in time
    pnaive[[i]][[a]] <- 1-(1-exp(-inte_exposed_all))
    
    
    # this is the probability of being naive to all serotypes up to this timepoint AND 
    # getting infected with a single serotype (any serotype) now in the immediate timestep
    # (e.g. this is the probability of a primary infection)
    # this will be the combined probability of getting serotype X NOW AND the probability of 
    # NOT GETTING any serotypes at all up to now
    
    pprim[[i]][[a]] =  pnaive[[i]][[a]]*(1-exp(-inte_now_immediate))
    
    
    #if not primarily infected or naive, this should be a multitypic infection
    #but remember it could be a primary infection with any of the four serotypes - that
    # is already captured above by multiplying the rate of infection by the serotype number
    pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
    
  }
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(year = rep(seq((year.start+1), (year.start+lts), 1), each=length(age_vect)),
                           age=rep(age_vect, lts), 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_prim=c(unlist(pprim)), 
                           multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  #head(p.out)
  p.add <- cbind.data.frame(year = seq((year.start+1), (year.start+lts), 1),
                            age=rep(0, length=lts), 
                            exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
                            multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
   # geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  
  #and split by age cat and take the min age
  p.out$age_trunc = floor(p.out$age)
  
  p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  #ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
   #geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #ggplot(data=subset(p.sum, year==2005)) + geom_point(aes(x=age, y=multi))  + 
   # geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year
}
model.age.incidence.series <- function(par.dat, age_vect){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  #N_sero = rep(N.sero.guess, length=lts)
  N_sero = par.dat$N_sero#[2:nrow(par.dat)]
  #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
  lambda = par.dat$lambda#[2:nrow(par.dat)]
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  pexposed= rep(NA, (length(age_vect)))
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  pnaive= rep(NA, (length(age_vect)))
  pnaive = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  pprim=as.list(rep(NA, (length(age_vect)*lts)))
  pprim= rep(NA, (length(age_vect)))
  pprim = rep(list(pexposed),lts)
  
  pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  pmulti= rep(NA, (length(age_vect)))
  pmulti = rep(list(pexposed),lts)
  
  year.start= min(par.dat$year)
  
  for (i in 1:lts){
    for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range across all the years
      
      
      #print(i)
      #print(a)
      
      #age = age_vect[a]
      #index = which(age_vect==age)
      #age_trunc = floor(age)
      
      #you want a vector of the duration of time each age kid from each year spends under each foi
      #so, first, the question is, was the kid alive before the time series started?
      #so this is the duration of time under the first foi (pre-2002) in years
      dur1 <- a-i
      dur1[dur1<0]<- 0
      
      #we don't have data on the foi or N_sero pre-ts, so we assume it is the same as 2002
      inte_pre_ts <- dur1*lambda[1]*N_sero[1]
      
      #then, you assume that each kid spends one year under all subsequent fois up to the current year
      dur2 <- 1
      # n.years <- i-1
      # 
      # if(n.years>0){
      #   dur2=1
      # }else{
      #   dur2=0
      # }
      
      
      #here is the rate of exposure to all serotypes up to the current period
      lambda.trunc <- lambda[1:(i)]
      Nsero.trunc <- N_sero[1:(i)]
      
      #here is the rate of exposure to any serotype within the time series but before the current year
      trunc.foi =  lambda.trunc*Nsero.trunc*dur2
      inte_exposed_lts = sum(trunc.foi)
      
      # #and the rate of getting exposed to any of all serotypes in the current year
      # dur3 = age-age_trunc
      # inte_exposed_current = dur3*N_sero[i]*lambda[i]
      # 
      
      #and, finally, the instantaneous rate of getting exposed to any given serotype in the current timestep
      #basically instantaneous but the case reports are in weeks
      inte_now_immediate = lambda[i]*N_sero[i]*(1/52)
      
      ########## now to sum ########## 
      
      #and this is the rate of being exposed to all serotypes at this point in time
      inte_exposed_all <- inte_pre_ts + inte_exposed_lts + inte_now_immediate
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_exposed_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_exposed_all))
      
      
      # this is the probability of being naive to all serotypes up to this timepoint AND 
      # getting infected with a single serotype (any serotype) now in the immediate timestep
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting serotype X NOW AND the probability of 
      # NOT GETTING any serotypes at all up to now
      
      pprim[[i]][[a]] =  pnaive[[i]][[a]]*(1-exp(-inte_now_immediate))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
      
    }
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
                           age=rep(age_vect, lts), 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_prim=c(unlist(pprim)), 
                           multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  #head(p.out)
  p.add <- cbind.data.frame(year = seq((year.start), (year.start+lts-1), 1),
                            age=rep(0, length=lts), 
                            exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
                            multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
   geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  
  #and split by age cat and take the min age
  p.out$age_trunc = floor(p.out$age)
  
  p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
  geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  print(p1)
  
  #ggplot(data=subset(p.sum, year==2002)) + geom_point(aes(x=age, y=multi))  + 
   #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year
}

model.age.incidence.series(par.dat = par.dat, age_vect = 1:35)

par.dat$lambda[par.dat$year>2007] <- 0.001



model.age.incidence.yearly <- function(par.dat, age_vect, year.now){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  par.dat = subset(par.dat, year <= year.now)
  lambda.start = par.dat$lambda[1]
  #N_sero = rep(N.sero.guess, length=lts)
  N_sero = par.dat$N_sero#[2:nrow(par.dat)]
  #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
  lambda = par.dat$lambda#[2:nrow(par.dat)]
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  pexposed= rep(NA, (length(age_vect)))
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  pnaive= rep(NA, (length(age_vect)))
  pnaive = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  pprim=as.list(rep(NA, (length(age_vect)*lts)))
  pprim= rep(NA, (length(age_vect)))
  pprim = rep(list(pexposed),lts)
  
  pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  pmulti= rep(NA, (length(age_vect)))
  pmulti = rep(list(pexposed),lts)
  
  year.start= min(par.dat$year)
  
  for (i in 1:lts){
    for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range across all the years
      
      
      #print(i)
      #print(a)
      
      #age = age_vect[a]
      #index = which(age_vect==age)
      #age_trunc = floor(age)
      
      #you want a vector of the duration of time each age kid from each year spends under each foi
      #so, first, the question is, was the kid alive before the time series started?
      #so this is the duration of time under the first foi (pre-2002) in years
      dur1 <- a-i
      dur1[dur1<0]<- 0
      
      #we don't have data on the foi or N_sero pre-ts, so we assume it is the same as 2002
      inte_pre_ts <- dur1*lambda[1]*N_sero[1]
      
      #then, you assume that each kid spends one year under all subsequent fois up to the current year
      #IF they were alive during that period
      
      
      
     
      #so, here choose which years this kid was alive
      lambda.trunc <- lambda[(length(lambda)-a+1): length(lambda)]
      Nsero.trunc <- N_sero[(length(N_sero)-a+1): length(N_sero)]
      
      #and was for only 1 year each
      dur2 <- 1
      
      #here is the rate of exposure to any serotype within the time series but before the current year
      trunc.foi =  lambda.trunc*Nsero.trunc*dur2
      inte_exposed_lts = sum(trunc.foi)
      
      # #and the rate of getting exposed to any of all serotypes in the current year
      # dur3 = age-age_trunc
      # inte_exposed_current = dur3*N_sero[i]*lambda[i]
      # 
      
      #and, finally, the instantaneous rate of getting exposed to any given serotype in the current timestep
      #basically instantaneous but the case reports are in weeks
      inte_now_immediate = lambda[i]*N_sero[i]*(1/52)
      
      ########## now to sum ########## 
      
      #and this is the rate of being exposed to all serotypes at this point in time
      inte_exposed_all <- inte_pre_ts + inte_exposed_lts + inte_now_immediate
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_exposed_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_exposed_all))
      
      
      # this is the probability of being naive to all serotypes up to this timepoint AND 
      # getting infected with a single serotype (any serotype) now in the immediate timestep
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting serotype X NOW AND the probability of 
      # NOT GETTING any serotypes at all up to now
      
      pprim[[i]][[a]] =  pnaive[[i]][[a]]*(1-exp(-inte_now_immediate))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
      
    }
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
                           age=rep(age_vect, lts), 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_prim=c(unlist(pprim)), 
                           multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  #head(p.out)
  p.add <- cbind.data.frame(year = seq((year.start), (year.start+lts-1), 1),
                            age=rep(0, length=lts), 
                            exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
                            multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  #   geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  # 
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  
  #and split by age cat and take the min age
  p.out$age_trunc = floor(p.out$age)
  
  p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  # ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
  #   geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  # 
  
  # ggplot(data=subset(p.sum, year==2002)) + geom_point(aes(x=age, y=multi))  + 
  #   geom_line(aes(x=age, y=multi)) 
  
  
  #and if you want to report only results for this year, call it here
  p.return = subset(p.sum, year==year.now)
  
  return(p.return) #returns prevalence by age for each year
}
model.age.incidence.yrly <- function(log.lambda.start, log.lambda.vector, timestep, N.sero.start, N.sero.vector, age_vect){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = exp(as.numeric(log.lambda.start))
  #N_sero = rep(N.sero.guess, length=lts)
  N_sero = N.sero.vector[1:timestep]
  
  #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
  lambda = exp(log.lambda.vector[1:timestep])
  
  #proportion exposed in any way - needs to long enough for each age within each year
  pexposed= as.list(rep(NA, length(age_vect))) 
  #pexposed= rep(NA, length(age_vect))
  #pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= as.list(rep(NA, length(age_vect))) 
  #pnaive= rep(NA, length(age_vect))
  #pnaive = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  pprim=as.list(rep(NA, length(age_vect)))
  #pprim= rep(NA, length(age_vect))
  #pprim = rep(list(pexposed),lts)
  
  pmulti=as.list(rep(NA, length(age_vect))) #secondary infection
  #pmulti= rep(NA, length(age_vect))
  #pmulti = rep(list(pexposed),lts)
  current.year = year.start+timestep
  
    for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range across all the years
      
      
      #print(i)
      #print(a)
      
      age = age_vect[a]
      index = which(age_vect==age)
      age_trunc = trunc(age)
      
      #i gives you years since the beginning of the time series.
      #so this is the duration of time under the first foi (pre-2003) in years
      timesub = timestep-1
      dur1=(age_trunc-timesub)
      dur1[dur1<0]<- 0
      
      
      #all years in between year 1 of the time series and current time had a duration of 1 year
      
      #then the duration of time in the current foi, also in years
      dur2 = age - age_trunc
      
      #and the instantanteous hazard of infection, since cases are reported in weeks
      dur3 = 1/52 #here also in years
      
      #here is the rate of exposure before the beginning of our time series
      #for some, this will be 0 if they were born during the time series
      inte_exposed_all_pre_lts = dur1*lambda.start*N.sero.start 
      
      #here is the rate of exposure to all serotypes up to the current period
      
      lambda.trunc <- lambda[1:(timestep-1)]
      Nsero.trunc <- N_sero[1:(timestep-1)]
      
      
      #here is the rate of exposure to any serotype within the time series but before the current year
      trunc.foi =  lambda.trunc*Nsero.trunc
      inte_exposed_lts = sum(trunc.foi)
      
      #and the rate of getting exposed to any of all serotypes in the current year
      inte_exposed_current = dur2*N_sero[timestep]*lambda[timestep]
      
      
      #rate of getting exposed to any given serotype in the current timestep
      #basically instantaneous but the case reports are in weeks
      inte_now_immediate = lambda[timestep]*N_sero[timestep]*dur3
      
      ########## now to sum ########## 
      
      #and this is the rate of being exposed to all serotypes at this point in time
      inte_exposed_all <- inte_exposed_all_pre_lts + inte_exposed_lts + inte_exposed_current + inte_now_immediate
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[a]] <- (1-exp(-inte_exposed_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[a]] <- 1-(1-exp(-inte_exposed_all))
      
      
      # this is the probability of being naive to all serotypes up to this timepoint AND 
      # getting infected with a single serotype (any serotype) now in the immediate timestep
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting serotype X NOW AND the probability of 
      # NOT GETTING any serotypes at all up to now
      
      pprim[[a]] =  pnaive[[a]]*(1-exp(-inte_now_immediate))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      pmulti[[a]] = 1 - pnaive[[a]] - (pprim[[a]])
      
    }

  
  #and get the estimates of each
  p.out = cbind.data.frame(age=rep(age_vect, lts), 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_prim=c(unlist(pprim)), 
                           multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  #head(p.out)
  p.add <- cbind.data.frame(age=rep(0, length=lts), 
                            exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
                            multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out,  age)
  
  ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
   geom_line(aes(x=age, y=multi)) 
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
   # geom_line(aes(x=age, y=exposed)) 
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  
  #and split by age cat and take the min age
  p.out$age_trunc = floor(p.out$age)
  
  p.split <- dlply(p.out, .( age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), all_prim=mean(all_prim), multi=mean(multi), sum_exp_naive= mean(sum_exp_naive), sum_naive_prim_multi = mean(sum_naive_prim_multi),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
  geom_line(aes(x=age, y=multi))
  #ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
   # geom_line(aes(x=age, y=exposed))
  
  #ggplot(data=subset(p.sum, year==2005)) + geom_point(aes(x=age, y=multi))  + 
  # geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year
}
model.age.incidence.one <- function(log.lambda.start, log.lambda.guess, N.sero.start,  age_vect, lts, year.start){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = exp(as.numeric(log.lambda.start))
  #N_sero = rep(N.sero.guess, length=lts)
  #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
  lambda = exp(log.lambda.guess)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  pexposed= rep(NA, (length(age_vect)))
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  pnaive= rep(NA, (length(age_vect)))
  pnaive = rep(list(pexposed),lts) 
  
  
  #first infection with a single serotype
  pprim=as.list(rep(NA, (length(age_vect)*lts)))
  pprim= rep(NA, (length(age_vect)))
  pprim = rep(list(pexposed),lts)
  
  pmulti=as.list(rep(NA, (length(age_vect)*lts))) #secondary infection
  pmulti= rep(NA, (length(age_vect)))
  pmulti = rep(list(pexposed),lts)
  
  
  
  for (i in 1:lts){
    for(a in 1:length(age_vect)){ #for-loops over all possible ages in our data range across all the years
      
      
      #print(i)
      #print(a)
      
      age = age_vect[a]
      index = which(age_vect==age)
      age_trunc = trunc(age)
      
      #i gives you years since the beginning of the time series.
      #so this is the duration of time under the first foi (pre-2003) in years
      dur1=(age_trunc-i)
      dur1[dur1<0]<- 0
      
      
      #all years in between year 1 of the time series and current time had a duration of 1 year
      
      #then the duration of time in the current foi, also in years
      dur2 = age - age_trunc
      
      #and the instanteous hazard of infection, since cases are reported in weeks
      dur3 = 1/52 #here also in years
      
      #here is the rate of exposure before the beginning of our time series
      #for some, this will be 0 if they were born during the time series
      inte_exposed_all_pre_lts = dur1*lambda.start*N.sero.start 
      
      #here is the rate of exposure to all serotypes up to the current period
      
      lambda.trunc <- lambda[1:(i-1)]
      #Nsero.trunc <- N_sero[1:(i-1)]
      
      
      #here is the rate of exposure to any serotype within the time series but before the current year
      #trunc.foi =  lambda.trunc*Nsero.trunc
      inte_exposed_lts = sum(lambda.trunc)
      
      #and the rate of getting exposed to any of all serotypes in the current year
      inte_exposed_current = dur2*lambda[i]
      
      
      #rate of getting exposed to any given serotype in the current timestep
      #basically instantaneous but the case reports are in weeks
      inte_now_immediate = lambda[i]*dur3
      
      ########## now to sum ########## 
      
      #and this is the rate of being exposed to all serotypes at this point in time
      inte_exposed_all <- inte_exposed_all_pre_lts + inte_exposed_lts + inte_exposed_current + inte_now_immediate
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_exposed_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_exposed_all))
      
      
      # this is the probability of being naive to all serotypes up to this timepoint AND 
      # getting infected with a single serotype (any serotype) now in the immediate timestep
      # (e.g. this is the probability of a primary infection)
      # this will be the combined probability of getting serotype X NOW AND the probability of 
      # NOT GETTING any serotypes at all up to now
      
      #pprim[[i]][[a]] =  pnaive[[i]][[a]]*(1-exp(-inte_now_immediate))
      
      
      #if not primarily infected or naive, this should be a multitypic infection
      #but remember it could be a primary infection with any of the four serotypes - that
      # is already captured above by multiplying the rate of infection by the serotype number
      #pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (pprim[[i]][[a]])
      
    }
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(year = rep(seq((year.start+1), (year.start+lts), 1), each=length(age_vect)),
                           age=rep(age_vect, lts), 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)))
                           #all_prim=c(unlist(pprim)), 
                           #multi = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  #p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  #head(p.out)
  p.add <- cbind.data.frame(year = seq((year.start+1), (year.start+lts), 1),
                            age=rep(0, length=lts), 
                            exposed=rep(0, length=lts), naive=rep(1, length=lts), 
                            sum_exp_naive=rep(1, length=lts))
  p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  # geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$exposed
  
  #and split by age cat and take the min age
  p.out$age_trunc = floor(p.out$age)
  
  p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), sum_exp_naive= mean(sum_exp_naive),cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  #ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #ggplot(data=subset(p.sum, year==2005)) + geom_point(aes(x=age, y=multi))  + 
  # geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year
}
wrap.all.llik.cumulative <- function(log.lambda.start, log.lambda.guess, N.sero.start, N.sero.vector, age_vect, lts, year.start, dat){
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  #head(df.out)
  # ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #         geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # 
  #  
  out.fit <- optim(par = log.lambda.guess, 
                   fn=log.lik.prop.cumulative, 
                   method = "Nelder-Mead",
                   log.lambda.start=log.lambda.start,
                   age_vect=age_vect,
                   N.sero.start=N.sero.start,
                   N.sero.vector = N.sero.vector,
                   dat=df.out,
                   lts=lts,
                   year.start = year.start,
                   hessian = T)
  
  
  
  fisher_info <- solve(out.fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  upper <- exp(out.fit$par+1.96*prop_sigma)
  lower <- exp(out.fit$par-1.96*prop_sigma)
  
  
  out.fit <- cbind.data.frame(lambda_fit = exp(out.fit$par), 
                              lambda_lower=lower,
                              lambda_upper=upper,
                              neg_llik=out.fit$value,
                              convergence=out.fit$convergence)#, 
  
  
  return(out.fit)
}
log.lik.prop.cumulative <- function(par, log.lambda.start, log.lambda.guess, N.sero.start, N.sero.vector, age_vect, lts, year.start, dat){
  
  #first run the model with the specified par
  if(unique(N.sero.vector)==1){
    
    out.mod <- model.age.incidence.one(log.lambda.guess=par, 
                                       log.lambda.start=log.lambda.start, 
                                       age_vect=age_vect, 
                                       N.sero.start=N.sero.start, 
                                       lts=lts,
                                       year.start = year.start)  
    
  }else{
    out.mod <- model.age.incidence.cumulative(log.lambda.guess=par, 
                                              log.lambda.start=log.lambda.start, 
                                              age_vect=age_vect, 
                                              N.sero.start=N.sero.start, 
                                              N.sero.vector = N.sero.vector, 
                                              lts=lts,
                                              year.start = year.start)  
  }
  
  
  
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat$age))
  #head(out.mod)
  #ggplot(data = out.mod) + geom_point(aes(x=age, y=cum_prop_cases)) +
   #                       geom_line(aes(x=age, y=cum_prop_cases)) +
    #                       facet_wrap(~year)
  
  #subset dat to just comparison
  dat.fit <- subset(dat, year >= min(out.mod$year) & year <= max(out.mod$year))
  
  #ggplot(data = dat.fit) + geom_point(aes(x=age, y=cum_prop_cases)) +
   # geom_line(aes(x=age, y=cum_prop_cases)) +
    #facet_wrap(~year)
  
  #make sure they line up
  dat.fit <- arrange(dat.fit, year, age) 
  out.mod <- arrange(out.mod, year, age)
  
  #warning if they are not the same
  if(nrow(dat.fit)!=nrow(out.mod)){
    warning("nrow unequal for data and model")
  }
  
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat.fit$cum_cases[i],dat.fit$n[i],p=out.mod$cum_prop_cases[i],log=T) 
    print(ll)
  }
  
  return(-ll)
}
comp.all.yrs.cumulative <- function(log.lambda.start, log.lambda.guess, N.sero.start,  age_vect, lts, year.start, dat){
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  
  out.NS1 <- wrap.all.llik.cumulative(log.lambda.start = log.lambda.start,
                                      log.lambda.guess = log.lambda.guess, #log(exp(log.lambda.guess)*N.sero.start), # because 
                                      N.sero.start =  N.sero.start,
                                      age_vect = age_vect,
                                      lts = lts,
                                      year.start = year.start,
                                      dat=dat,
                                      N.sero.vector=rep(1, length=lts))
  
  out.NS2 <- wrap.all.llik.cumulative(log.lambda.start = log.lambda.start,
                                      log.lambda.guess = log.lambda.guess,
                                      N.sero.start =  N.sero.start,
                                      age_vect = age_vect,
                                      lts = lts,
                                      year.start = year.start,
                                      dat=dat,
                                      N.sero.vector=rep(2, length=lts))
  
  out.NS3 <- wrap.all.llik.cumulative(log.lambda.start = log.lambda.start,
                                      log.lambda.guess = log.lambda.guess,
                                      N.sero.start =  N.sero.start,
                                      age_vect = age_vect,
                                      lts = lts,
                                      year.start = year.start,
                                      dat=dat,
                                      N.sero.vector=rep(3, length=lts))
  
  out.NS4 <- wrap.all.llik.cumulative(log.lambda.start = log.lambda.start,
                                      log.lambda.guess = log.lambda.guess,
                                      N.sero.start =  N.sero.start,
                                      age_vect = age_vect,
                                      lts = lts,
                                      year.start = year.start,
                                      dat=dat,
                                      N.sero.vector=rep(4, length=lts))
  
  
  out.NS1$N_sero=1
  out.NS2$N_sero=2
  out.NS3$N_sero=3
  out.NS4$N_sero=4
  
  out.all <- rbind(out.NS1, out.NS2, out.NS3, out.NS4)
  out.all$year <- min.year
  
  
  
  #now take this and feed it into the other years
  #so this will be the foi of all prior years
  
  
  
  
  
  return(out.all)
  
  
}
log.lik.prop.single <- function(par, log.lambda.start, log.lambda.vector, N.sero.start, N.sero.vector, timestep,  age_vect, year.start, dat){
  
  log.lambda.vector[timestep] <- par
  
  out.mod <- model.age.incidence.yrly(log.lambda.start=log.lambda.start, 
                                      log.lambda.vector = log.lambda.vector,
                                      age_vect=age_vect, 
                                      N.sero.start=N.sero.start, 
                                      N.sero.vector = N.sero.vector, 
                                      timestep = timestep)  
  
  #with(out.mod, plot(age, cum_prop_cases, type="b"))
  current.year = year.start+timestep
  dat.fit <- subset(dat, year==current.year)
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat.fit$age))
  #plot(out.mod$cum_prop_cases, type="b")
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat.fit$age)){
    ll=ll+dbinom(dat.fit$cum_cases[i],dat.fit$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
log.lik.prop.yearly <- function(par, par.dat, age_vect, dat, year.now){
  
  par.dat$lambda[par.dat$year==year.now] <- exp(par)
  
  out.mod <- model.age.incidence.yearly(par.dat = par.dat, 
                                        age_vect=age_vect, 
                                        year.now = year.now)  
  
  #with(out.mod, plot(age, cum_prop_cases, type="b"))
  
  dat.fit <- subset(dat, year==year.now)
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat.fit$age))
  #plot(out.mod$cum_prop_cases, type="b")
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat.fit$age)){
    ll=ll+dbinom(dat.fit$cum_cases[i],dat.fit$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
}
lambda.guess = 0.089
N.sero.vect = rep(2, 19)
age_vect = 1:35
fit.all.yrs.sequential <- function(dat, lambda.guess, N.sero.vect, age_vect){
  
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
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  log.lambda.guess <- log(par.dat$lambda[i]*2)
  par.dat$N_sero[i] <- 2
  out.NS2 <- optim(par = log.lambda.guess, #guess the same as last time
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  
  log.lambda.guess <- log(par.dat$lambda[i]*3)
  par.dat$N_sero[i] <- 3
  out.NS3 <- optim(par = log.lambda.guess, #guess the same as last time
                   fn=log.lik.prop.yearly, #not cumulative... we for-loop through and test year by year
                   method = "Brent",
                   lower = -100, upper=20,
                   par.dat=par.dat, age_vect=age_vect, dat=df.out, year.now=par.dat$year[i])
  
  log.lambda.guess <- log(par.dat$lambda[i]*4)
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
  
  #and combine the data and return
  
  out.dat <- cbind.data.frame(year = year.start:(year.start+lts),
                              N_sero = c(N.sero.start, N.sero.list),
                              lambda = c(exp(log.lambda.start), lambda.list),
                              AIC = c(AIC.start, AIC.list))
  
  p1 <- ggplot(data=out.dat) + geom_point(aes(x=year, y=lambda, fill=log(delta_AIC))) +
        scale_fill_viridis_c(direction=-1)
  return(out.dat)

}


#and run the model with a time series
run.mod <- function(par.dat = out.dat){
  
}





out.fit <- wrap.all.llik.cumulative(log.lambda.start = log(0.08904973),
                         log.lambda.guess = rep(log(0.08904973*2), length =  length(unique(dat.sum$year))), #log(exp(log.lambda.guess)*N.sero.start), # because 
                         N.sero.start =  2,
                         age_vect = seq(1,35, length = 52*34),
                         lts = length(unique(dat.sum$year)),
                         year.start = 2002,
                         dat=dat,
                         N.sero.vector= rep(1, length =  length(unique(dat.sum$year))))



comp.all.NS.2002 <- comp.all.yrs.cumulative(log.lambda.start = log(0.7461632),
                                            log.lambda.guess = rep(log(0.05149362), length =  length(unique(dat.sum$year))),
                                            N.sero.start = 2,
                                            #N.sero.vector = rep(2, length =  length(unique(dat.sum$year))),
                                            age_vect = seq(1,35, length = 52*34),
                                            lts = length(unique(dat.sum$year)),
                                            year.start = 2002,
                                            dat = dat)
                                 
 
seq.start <- fit.all.yrs.sequential(log.lambda.start = log(0.08904973),
                                    N.sero.start = 2,
                                    age_vect = seq(1,35, length = 52*34),
                                    lts = length(unique(dat.sum$year)),
                                    year.start = 2002,
                                    dat = dat)
                             
#and plot


#and another year
out.2003.likprop <- comp.all.NS.multi(lambda.guess = c(.2),
                                      age_vect = seq(1,35, length = 52*34),
                                      dat=year.split[[1]])

#now calculate for all years
out.optim <- lapply(year.split, comp.all.NS.multi, lambda.guess = .2,
                    age_vect = 1:max(dat$age))

out.df <- data.table::rbindlist(out.optim)
out.df$AIC <- 2*(out.df$neg_llik) + 2*out.df$N_sero

out.df <- ddply(out.df, .(year), summarise, lambda_fit=unique(lambda_fit), lambda_lower=unique(lambda_lower), lambda_upper=unique(lambda_upper), neg_llik=unique(neg_llik), convergence=unique(convergence), N_sero=unique(N_sero), AIC=unique(AIC), delta_AIC = min(AIC))
out.df$delta_AIC=out.df$AIC-out.df$delta_AIC

out.df$delta_AIC <- round(out.df$delta_AIC, 2)
out.df$N_sero <- as.factor(out.df$N_sero)

pC <- ggplot(data=subset(out.df, delta_AIC==0)) + geom_point(aes(x=year, y=lambda_fit, color=N_sero), size=3) +
  geom_line(aes(x=year, y=lambda_fit)) +
  geom_ribbon(aes(x=year, ymin=lambda_lower, ymax=lambda_upper), alpha=.3)+
  ylab("annual force of infection") +theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.9), "cm"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))
#print(pC)     
out.df$log_delta_AIC <- log(out.df$delta_AIC)
out.df$log_delta_AIC[out.df$log_delta_AIC==-Inf] <- 0
#or plot all three and color by AIC
pC <- ggplot(data=out.df) + 
  geom_line(aes(x=year, y=lambda_fit, color=N_sero), size=1) +
  geom_linerange(aes(x=year, ymin=lambda_lower, ymax=lambda_upper, color=N_sero), size=1)+
  geom_point(aes(x=year, y=lambda_fit, fill=log_delta_AIC, color=N_sero), size=3, shape=21, stroke=1) +
  scale_color_manual(values=scales::hue_pal()(4), name="number of\nserotypes") +
  scale_fill_viridis_c(direction = -1, name="log\n(delta AIC)")+
  ylab("annual force of infection") +theme_bw()+
  coord_cartesian(ylim=c(0,.25))+
  theme(panel.grid = element_blank(),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.background = element_rect(color="black"),
        legend.position = c(.5,.89),
        plot.margin = unit(c(.1,.1,.1,.9), "cm"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))+
  guides(color=guide_legend(nrow=2))

pC




#compile
top <- cowplot::plot_grid(pA,pB, nrow = 1, rel_widths = c(1,1), labels = c("A", "B"), label_size = 22)
bottom <- cowplot::plot_grid(pC,pD, nrow = 1, rel_widths = c(1,1), labels = c("C", "D"), label_size = 22)

pFig2 <- cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1,1.1))
pFig2



# ggsave(file = paste0(homewd, "/Fig1-uppercase.png"),
#        plot=pFig1,
#        units=c("mm"),  
#        width=100, 
#        height=65, 
#        scale=3, 
#        dpi=300)


# now we want to make black background line on the subplot b and d 3 epi years

#compile
top <- cowplot::plot_grid(pA,pB_blk_lines, nrow = 1, rel_widths = c(1,1), labels = c("a", "b"), label_size = 22)
bottom <- cowplot::plot_grid(pC,pD_blk_lines, nrow = 1, rel_widths = c(1,1), labels = c("c", "d"), label_size = 22)

pFig2_blk_lines <- cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1,1.1))# +cowplot::draw_text('epi years of\n\t2007/2012/2019', x = 0.945, y = 0.1, size = 10, hjust = 0.5, vjust = 0.5)




ggsave(file = paste0(homewd, "/final-figures/fig1.png"),
       plot=pFig2_blk_lines,
       units=c("mm"),  
       width=100, 
       height=65, 
       scale=3, 
       dpi=300)




#and for supplement:
#now run and plot with data
run.plot <- function(mod.dat, dat, age_vect, filename){
  
  #pick best fit by year
  mod.dat <- subset(mod.dat, delta_AIC==0)
  
  #and run over many years
  out.run <- mapply(FUN=model.age.incidence, 
                    log.lambda = as.list(log(mod.dat$lambda_fit)), 
                    #year = as.list(mod.dat$year),
                    N_sero = as.list(as.numeric(as.character(mod.dat$N_sero))),
                    MoreArgs = list(age_vect=age_vect), SIMPLIFY = F)
  
  #and run over many years
  out.run.lci <- mapply(FUN=model.age.incidence, 
                        log.lambda = as.list(log(mod.dat$lambda_lower)), 
                        #year = as.list(mod.dat$year),
                        N_sero = as.list(as.numeric(as.character(mod.dat$N_sero))),
                        MoreArgs = list(age_vect=age_vect), SIMPLIFY = F)
  
  
  out.run.uci <- mapply(FUN=model.age.incidence, 
                        log.lambda = as.list(log(mod.dat$lambda_upper)), 
                        #year = as.list(mod.dat$year),
                        N_sero = as.list(as.numeric(as.character(mod.dat$N_sero))),
                        MoreArgs = list(age_vect=age_vect), SIMPLIFY = F)
  
  
  
  out.run <- data.table::rbindlist(out.run)
  out.run.lci <- data.table::rbindlist(out.run.lci)
  out.run.uci <- data.table::rbindlist(out.run.uci)
  out.run$year <- rep(unique(year(dat$date)), each=(max(age_vect)+1))
  
  dat.mod <- dplyr::select(out.run, age, cum_prop_cases, year)
  
  dat.mod$lci = out.run.lci$cum_prop_cases
  dat.mod$uci = out.run.uci$cum_prop_cases
  dat.mod$type<- "model"
  names(dat.mod) <- c("age", "prop", "year","lci", "uci",  "type")
  
  #now join with data and plot
  year.split <- dlply(dat, .(year))
  dat.out <- lapply(X=year.split, FUN=sum.yr, age_vect = age_vect)
  dat.out <- data.table::rbindlist(dat.out)
  dat.out$year <- rep(unique(dat$year), each=max(age_vect+1))
  
  dat.dat <- dplyr::select(dat.out, age, cum_prop_cases, year)
  dat.dat$lci <- dat.dat$uci <- NA
  dat.dat$type<- "data"
  names(dat.dat) <-  c("age", "prop", "year","lci", "uci",  "type")
  dat.all <- rbind(dat.dat, dat.mod)
  
  p1 <- ggplot(data=dat.all) + 
    geom_point(aes(x=age, y=prop, color=type)) +
    geom_line(aes(x=age, y=prop, color=type, linetype=type)) +
    facet_wrap(~year, ncol=4) + xlim(0,30) + ylab("proportion of cases") +
    theme_bw() + theme(panel.grid = element_blank())
  
  
  
  
  ggsave(file =filename,
         plot=p1,
         units="mm",  
         width=70, 
         height=70, 
         scale=3, 
         dpi=300)
  
  return(dat.all)
  
  
}



run.plot(mod.dat=out.df,
         dat=dat, 
         age_vect = 1:37, 
         #N_sero=2,
         filename = paste0(homewd, "/final-figures/figS2.png"))















