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

homewd= "/Users/carabrook/Developer/cambodia-dengue-distribute/"
setwd(paste0(homewd, "/figures/Fig3/"))

#load the master list
dat <- read.csv(file = paste0(homewd,"data/denv_seq_prior_infection.csv"), header=T, stringsAsFactors = F)
head(dat)

dat = subset(dat, year==2019)
dat = subset(dat, !is.na(igg_res)) #removing two DENV-1 from pagodas with no rapidtest

dat$igg_res[dat$igg_res=="pos"] <- 1
dat$igg_res[dat$igg_res=="neg"] <- 0

dat$age <- as.numeric(dat$age)

dat = subset(dat, !is.na(age)) #85 cases in 2019

#plot age-sero of naive vs. multi
dat$age <- ceiling(dat$age)

sum.dat <- function(dat, age_vect){
df.sum <- ddply(dat, .(age, igg_res), summarise, Nage = length(age))
  df.out = cbind.data.frame(age=1:max(age_vect), igg_res = rep(c(0,1), each=max(age_vect)))
  df.out <- merge(df.out, df.sum, by=c("age", "igg_res"), all.x = T, sort = F)
  df.out[is.na(df.out)]<- 0
  df.out <- rbind(c(0,0,0), df.out)
  df.out <- rbind(c(0,1,0), df.out)
  df.out$igg_res <- as.factor(df.out$igg_res)
  df.out <- arrange(df.out, igg_res, age)
  #ggplot(dat = df.out) + geom_point(aes(x=age, y=Nage, color=igg_res))
  df.out$n_tot <- sum(df.out$Nage)
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
model.age.incidence <- function(log.lambda, N_sero, age_vect){
  
  lambda= exp(as.numeric(log.lambda))
  pexposed=as.list(rep(NA,max(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,max(age_vect))) #naive
  pprim=as.list(rep(NA,max(age_vect))) #first infection with a single serotype
  pmulti=as.list(rep(NA,max(age_vect))) #secondary infection
  prob_naive_one=as.list(rep(NA,max(age_vect))) #prob of being naive to all strains at one timestep prior to the current
  for(a in 1:max(age_vect)){ #for-loops over all possible ages in our data range
    
    
    dur1=a#up to 1 timestep prior--but we assume they are instantaneous
    dur2=1 #this timestep
    
    inte_exposed_all = a*lambda*N_sero
    inte_exposed_prior = dur1*lambda*N_sero #rate of getting exposed to all strains in the timesteps before
    inte_exposed_one = (dur1*lambda) #rate of getting exposed to a single strain in the timesteps before
    inte_now_one = dur2*lambda
    
    
    #prob_exposed_one = (1-exp(-inte_exposed_one))
    
    
    #this is the probability of being naive to all strains at one timestep prior
    prob_naive_one[[a]] <- 1-(1-exp(-(inte_exposed_prior)))
    
    #this is the probability of being naive to all strains at one timestep prior AND 
    #getting infected with a given strain in this timestep
    pprim[[a]] = prob_naive_one[[a]]*(1-exp(-inte_now_one))
    
    
    pexposed[[a]] <- (1-exp(-inte_exposed_all))
    pnaive[[a]] <- 1-(1-exp(-inte_exposed_all))
    
    #if not primarily infected or naive, you should be a multitypic infection
    #we could fit to that generally but we have serotype-specific info that can 
    #fill in instead
    
    #for it to be a multitypic infection with D1, you need to multiply
    #prob exposure to D2 up to a-1 x  prob naive to D1 up to a-1 x prob to D1 now
    pmulti[[a]] = 1 - pnaive[[a]] - (N_sero*pprim[[a]])
    
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
  pexposed=as.list(rep(NA,max(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,max(age_vect))) #naive
  #pprim=as.list(rep(NA,max(age_vect))) #first infection with a single serotype
  #pmulti=as.list(rep(NA,max(age_vect))) #secondary infection
  #prob_naive_one=as.list(rep(NA,max(age_vect))) #prob of being naive to all strains at one timestep prior to the current
  for(a in 1:max(age_vect)){ #for-loops over all possible ages in our data range
    
    
    #dur1=a#up to 1 timestep prior--but we assume they are instantaneous
    #dur2=1 #this timestep
    
    inte_exposed_all = a*lambda
    #inte_exposed_prior = dur1*lambda #rate of getting exposed to all strains in the timesteps before
    #inte_exposed_one = (dur1*lambda) #rate of getting exposed to a single strain in the timesteps before
    #inte_now_one = dur2*lambda
    
    
    #prob_exposed_one = (1-exp(-inte_exposed_one))
    
    
    #this is the probability of being naive to all strains at one timestep prior
    #prob_naive_one[[a]] <- 1-(1-exp(-(inte_exposed_prior)))
    
    #this is the probability of being naive to all strains at one timestep prior AND 
    #getting infected with a given strain in this timestep
    #pprim[[a]] = prob_naive_one[[a]]*(1-exp(-inte_now_one))
    
    
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
log.lik.prop.both <- function(par, age_vect, dat,  N_sero){
  
  #first run the model with the specified par
  if(N_sero==1){
    
    out.mod <- model.age.one(log.lambda = par, age_vect = age_vect)  
    
    ll=0
    for (i in 1:length(dat$age)){
      ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
    }
    
  }else{
    out.mod <- model.age.incidence(log.lambda = par, age_vect = age_vect,  N_sero=N_sero)
  
  
  
  
  #head(out.mod)
  #head(dat)
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat$age))
  #plot(out.mod$cum_prop_cases, type="b")
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) #and prim
    ll=ll+dbinom(dat$n_prim_cases[i],dat$ntot[i],p=out.mod$all_prim[i],log=T) #and prim
  }
  }
  
  return(-ll)
}


wrap.all.llik.both <- function(lambda.guess, age_vect, dat, N_sero){
  
  #first, prep the data
  df.out <- sum.dat(dat, age_vect = age_vect)
  #head(df.out)
  
  #  
  out.fit <- optim(log(lambda.guess), 
                   fn=log.lik.prop.both, method = "Brent",
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
wrap.all.llik.multi <- function(lambda.guess, age_vect, dat, N_sero){
  
  #first, prep the data
  df.out <- sum.dat(dat, age_vect = age_vect)
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
comp.all.NS.both <- function(lambda.guess, age_vect, dat){
  
  out.NS1 <- wrap.all.llik.both(lambda.guess = lambda.guess,
                                      age_vect = age_vect,
                                      dat=dat,
                                      N_sero=1)
  
  out.NS2 <- wrap.all.llik.both(lambda.guess = lambda.guess,
                                age_vect = age_vect,
                                dat=dat,
                                N_sero=2)
  
  out.NS3 <- wrap.all.llik.both(lambda.guess = lambda.guess,
                                age_vect = age_vect,
                                dat=dat,
                                N_sero=3)
  
  out.NS4 <- wrap.all.llik.both(lambda.guess = lambda.guess,
                                age_vect = age_vect,
                                dat=dat,
                                N_sero=4)
  
  
  out.NS1$N_sero=1
  out.NS2$N_sero=2
  out.NS3$N_sero=3
  out.NS4$N_sero=4
  
  out.all <- rbind(out.NS1, out.NS2, out.NS3, out.NS4)
  
  return(out.all)
  
  
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
  
  return(out.all)
  
  
}

out.2019.cohort <- comp.all.NS.both(lambda.guess = c(.01),
                                  age_vect = 1:37,
                                  dat=dat)

out.2019.multi <- comp.all.NS.multi(lambda.guess = c(.01),
                                   age_vect = 1:37,
                                   dat=dat)


#compare
run.one.par <- function(lambda, lambda_lower, lambda_upper, age_vect, N_sero){
  
  if(N_sero==1){
    
    out.run <- model.age.one(log.lambda = log(lambda), 
                                   #year = mod.dat$year,
                                   age_vect=age_vect)
    
    out.run.lci <- model.age.one(log.lambda = log(lambda_lower), 
                                       #year = mod.dat$year,
                                       age_vect=age_vect)
    # 
    out.run.uci <- model.age.one(log.lambda = log(lambda_upper), 
                                       #year = mod.dat$year,
                                       age_vect=age_vect)
    #head(out.run) 
    #out.run <- data.table::rbindlist(out.run)
    
    dat.mod <- dplyr::select(out.run, age,  cum_prop_cases)
    dat.mod$lci_cum_prop = out.run.lci$cum_prop_cases
    #dat.mod$lci_prim = out.run.lci$all_prim
    dat.mod$uci_cum_prop = out.run.uci$cum_prop_cases
    #dat.mod$uci_prim = out.run.uci$all_prim
    dat.mod$type<- "model"
    
    #names(dat.mod) <- c("year", "age", "prop", "lci", "uci", "type")
    names(dat.mod) <- c( "age", "prop",  "prop_lci", "prop_uci",  "type")
    
    #names(dat.mod) <- c( "age", "cum_prop_multi", "prev_prim", "cum_prop_multi_lci","prev_prim_lci", "cum_prop_multi_uci", "prev_prim_uci",  "type")
    
    dat.mod$class <- "multi_cum"
    
    
  }else{
  out.run <- model.age.incidence(log.lambda = log(lambda), 
                                 #year = mod.dat$year,
                                 age_vect=age_vect, N_sero=N_sero)
  
  out.run.lci <- model.age.incidence(log.lambda = log(lambda_lower), 
                                     #year = mod.dat$year,
                                     age_vect=age_vect, N_sero=N_sero)
  # 
  out.run.uci <- model.age.incidence(log.lambda = log(lambda_upper), 
                                     #year = mod.dat$year,
                                     age_vect=age_vect, N_sero=N_sero)
  #head(out.run) 
  #out.run <- data.table::rbindlist(out.run)
  
  dat.mod <- dplyr::select(out.run, age,  cum_prop_cases)
  dat.mod$lci_cum_prop = out.run.lci$cum_prop_cases
  #dat.mod$lci_prim = out.run.lci$all_prim
  dat.mod$uci_cum_prop = out.run.uci$cum_prop_cases
  #dat.mod$uci_prim = out.run.uci$all_prim
  dat.mod$type<- "model"
  
  #names(dat.mod) <- c("year", "age", "prop", "lci", "uci", "type")
  names(dat.mod) <- c( "age", "prop",  "prop_lci", "prop_uci",  "type")
  
  #names(dat.mod) <- c( "age", "cum_prop_multi", "prev_prim", "cum_prop_multi_lci","prev_prim_lci", "cum_prop_multi_uci", "prev_prim_uci",  "type")
  
  
  #out.run$cum_prop_prim <- cumsum(out.run$all_primary_inf)/sum(out.run$all_primary_inf)
  dat.mod2 <- dplyr::select(out.run, age,  all_prim)
  dat.mod2$lci_prim = out.run.lci$all_prim
  dat.mod2$uci_prim = out.run.uci$all_prim
  #dat.mod$lci = out.run.lci$cum_prop_cases
  #dat.mod$uci = out.run.uci$cum_prop_cases
  dat.mod2$type<- "model"
  #names(dat.mod) <- c("year", "age", "prop", "lci", "uci", "type")
  names(dat.mod2) <- c( "age", "prop",  "prop_lci", "prop_uci",  "type")
  dat.mod$class = "multi_cum"
  dat.mod2$class = "prim_prev"
  
  dat.mod <- rbind(dat.mod, dat.mod2)
  }
  
  return(dat.mod)
}
run.plot.all <- function(mod.dat, dat, age_vect){
  
  #and run over many years
  run1 <- run.one.par(lambda = mod.dat$lambda_fit[mod.dat$N_sero==1],
                      lambda_lower = mod.dat$lambda_lower[mod.dat$N_sero==1],
                      lambda_upper = mod.dat$lambda_upper[mod.dat$N_sero==1],
                      age_vect = age_vect,
                      N_sero = 1)
  run1$N_sero = 1
  
  run2 <- run.one.par(lambda = mod.dat$lambda_fit[mod.dat$N_sero==2],
                      lambda_lower = mod.dat$lambda_lower[mod.dat$N_sero==2],
                      lambda_upper = mod.dat$lambda_upper[mod.dat$N_sero==2],
                      age_vect = age_vect,
                      N_sero = 2)
  run2$N_sero = 2
  
  run3 <- run.one.par(lambda = mod.dat$lambda_fit[mod.dat$N_sero==3],
                      lambda_lower = mod.dat$lambda_lower[mod.dat$N_sero==3],
                      lambda_upper = mod.dat$lambda_upper[mod.dat$N_sero==3],
                      age_vect = age_vect,
                      N_sero = 3)
  run3$N_sero = 3
  
  
  run4 <- run.one.par(lambda = mod.dat$lambda_fit[mod.dat$N_sero==4],
                      lambda_lower = mod.dat$lambda_lower[mod.dat$N_sero==4],
                      lambda_upper = mod.dat$lambda_upper[mod.dat$N_sero==4],
                      age_vect = age_vect,
                      N_sero = 4)
  run4$N_sero = 4
  
  dat.out <- sum.dat(dat,age_vect = age_vect)
  #head(dat.out)
  #now join with data and plot
  
  
  dat.dat <- dplyr::select(dat.out, age, cum_prop_cases)
  dat.dat$lci <- dat.dat$uci <- NA
  dat.dat$type<- "data"
  names(dat.dat) <-  c( "age", "prop", "prop_lci", "prop_uci", "type")
  #names(dat.dat) <-  c("age", "prop",  "type")
  dat.dat$class ="multi_cum"
  dat.dat$N_sero = "data"
  
  dat.all <- rbind(dat.dat, run1, run2, run3, run4)
  
  #and also the primary
  
  dat.dat2 <- dplyr::select(dat.out, age, prop_primary_by_age)
  dat.dat2$lci <- dat.dat2$uci <- NA
  dat.dat2$type <- "data"
  names(dat.dat2) <-  c( "age", "prop", "prop_lci", "prop_uci", "type")
  dat.dat2$class ="prim_prev"
  
  dat.dat2$N_sero = "data"
  
  dat.all <- rbind(dat.dat2, dat.all)
  
  dat.all$type[dat.all$N_sero!="data"] <- paste0(dat.all$type[dat.all$N_sero!="data"], ", N-sero=", dat.all$N_sero[dat.all$N_sero!="data"])
  dat.all$label <- dat.all$class
  dat.all$label[dat.all$class=="prim_prev"] <- "b"
  dat.all$label[dat.all$class=="multi_cum"] <- "a"
  
  dat.all$class[dat.all$class=="prim_prev"] <- "proportion\nprimary infections"
  dat.all$class[dat.all$class=="multi_cum"] <- "cumulative proportion\nsecondary infections"
  
  p1 <- ggplot(data=dat.all) + 
    facet_grid(class~.) +
    geom_line(data=subset(dat.all, type!="data"), aes(x=age, y=prop, color=type), position = position_jitter(width=.2), size=1) +
    geom_point(data=subset(dat.all, type=="data"), aes(x=age, y=prop)) +
    geom_line(data=subset(dat.all, type=="data"), aes(x=age, y=prop), size=.7) +
    geom_ribbon(data=subset(dat.all, type!="data"), aes(x=age, ymin=prop_lci, ymax=prop_uci, fill=type), position = position_jitter(width=.2), alpha=.3) +
    geom_label(aes(x=.1,y=.9, label=label), label.size = NA, size=9, fontface="bold")+
    #geom_line(data=subset(dat.all, class=="primary_infection"),aes(x=age, y=prop, color=type)) +
    #geom_line(data=subset(dat.all, class!="primary_infection"),aes(x=age, y=prop, color=type)) +
    #geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type), alpha=.3)+
    xlim(0,max(age_vect)) + ylab("proportion") +theme_bw()+
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill="white"),
          strip.text = element_text(size=16),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.background = element_rect(color="black"),
          legend.position = c(.83,.64),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14)) 
    
  
  print(p1)
  
  return(p1)
  
  
}


ptop <- run.plot.all(mod.dat = out.2019.multi,
                dat=dat,
                age_vect = 1:37)


#functions
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
  return(df.out)
  
}
dat.nat <- read.csv(file = paste0(homewd, "/data/DENV-Nat-Aged.csv") , header = T, stringsAsFactors = F)
head(dat.nat)
dat.nat <- arrange(dat.nat, date, age)
dat.nat$age <- ceiling(dat.nat$age)
dat.nat = subset(dat.nat, year==2019)

#now calculate annual FOI, using the Cummings method
model.age.incidence <- function(log.lambda, N_sero, age_vect, year){
  
  pexposed=as.list(rep(NA,max(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,max(age_vect))) #naive
  pprim=as.list(rep(NA,max(age_vect))) #first infection with a single serotype
  pmulti=as.list(rep(NA,max(age_vect))) #secondary infection
  for(a in 1:max(age_vect)){ #for-loops over all possible ages in our data range
    
    inte_exposed_all = N_sero*a*exp(as.numeric(log.lambda))
    inte_exposed_all_but_one = (N_sero-1)*a*exp(as.numeric(log.lambda)) 
    inte_exposed_one = a*exp(as.numeric(log.lambda))
    
    #prob_exposed_one = (1-exp(-inte_exposed_one))
    
    #this is the probability of escaping infection with 3 strains at age a
    prob_escaped_one = (1-(1-exp(-inte_exposed_all_but_one))) 
    
    #this is the probaility of being infected with strain i at age a
    prob_inf_one = 1-exp(-inte_exposed_one)
    
    
    pexposed[[a]]=1-exp(-inte_exposed_all) #returns probability of having experience infection (any number of infection) at a given age
    pnaive[[a]] = 1-(1-exp(-inte_exposed_all)) #returns the probability of remaining naive for a given age.
    
    #returns the probability of experiencing your first infection with any serotype i at a given age
    #it is the probability of remaining naive up to that point multiplied by the probability of any exposure at that age
    #pprim[[a]] =  pnaive[[a]]*((1-exp(-inte_exposed_one))-1)#probability of being naive to all but one * being infected with 1... that one???
    pprim[[a]] =  prob_escaped_one*prob_inf_one
    pmulti[[a]] = 1-((pnaive[[a]]))-(N_sero*pprim[[a]])
    
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(age=age_vect, 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)),
                           all_primary_inf=c(unlist(pprim)*N_sero), 
                           multitypic_inf = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  p.add <- cbind.data.frame(age=0, exposed=0, naive=1, all_primary_inf=0, 
                            multitypic_inf=0, sum_exp_naive=1, sum_naive_prim_multi=1)
  p.out <- rbind(p.add, p.out)
  p.out$year <- year
  
  
  ggplot(data=p.out) + geom_point(aes(x=age, y=multitypic_inf))  + geom_line(aes(x=age, y=multitypic_inf))
  # #need to cap when the proportion reaches 1
  #  if(max(p.out$multitypic_inf>=1)){
  #    first.one <- min(which(p.out$multitypic_inf>=1))
  #    p.out = p.out[1:first.one,]
  #  }
  # 
  #plot(ecdf(p.out$multitypic_inf))
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  
  p.out$cum_prop_cases <- p.out$multitypic_inf
  
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=all_primary_inf))  + geom_line(aes(x=age, y=all_primary_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multitypic_inf))  + geom_line(aes(x=age, y=multitypic_inf))
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases))
  
  return(p.out) #returns prevalence by age for the length of all possible ages in our dataset
}
model.age.one <- function(log.lambda, age_vect, year){
  
  pexposed=as.list(rep(NA,max(age_vect))) #proportion exposed in any way
  pnaive=as.list(rep(NA,max(age_vect))) #naive
  #pprim=as.list(rep(NA,max(age_vect))) #first infection with a single serotype
  #pmulti=as.list(rep(NA,max(age_vect))) #secondary infection
  for(a in 1:max(age_vect)){ #for-loops over all possible ages in our data range
    
    inte_exposed_all = a*exp(as.numeric(log.lambda))
    #inte_exposed_all_but_one = (N_sero-1)*a*exp(as.numeric(log.lambda)) 
    #inte_exposed_one = a*exp(as.numeric(log.lambda))
    
    #prob_exposed_one = (1-exp(-inte_exposed_one))
    
    #this is the probability of escaping infection with 3 strains at age a
    #prob_escaped_one = (1-(1-exp(-inte_exposed_all_but_one))) 
    
    #this is the probaility of being infected with strain i at age a
    #prob_inf_one = 1-exp(-inte_exposed_one)
    
    
    pexposed[[a]]=1-exp(-inte_exposed_all) #returns probability of having experience infection (any number of infection) at a given age
    pnaive[[a]] = 1-(1-exp(-inte_exposed_all)) #returns the probability of remaining naive for a given age.
    
    #returns the probability of experiencing your first infection with any serotype i at a given age
    #it is the probability of remaining naive up to that point multiplied by the probability of any exposure at that age
    #pprim[[a]] =  pnaive[[a]]*((1-exp(-inte_exposed_one))-1)#probability of being naive to all but one * being infected with 1... that one???
    #pprim[[a]] =  prob_escaped_one*prob_inf_one
    #pmulti[[a]] = 1-((pnaive[[a]]))-(N_sero*pprim[[a]])
    
  }
  
  #and get the estimates of each
  p.out = cbind.data.frame(age=age_vect, 
                           exposed=c(unlist(pexposed)),
                           naive = c(unlist(pnaive)))
                           #all_primary_inf=c(unlist(pprim)*N_sero), 
                           #multitypic_inf = c(unlist(pmulti)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,2:3])
  #p.out$sum_naive_prim_multi <- rowSums(p.out[,3:5])
  
  p.add <- cbind.data.frame(age=0, exposed=0, naive=1,  sum_exp_naive=1)
  p.out <- rbind(p.add, p.out)
  p.out$year <- year
  
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multitypic_inf))  + geom_line(aes(x=age, y=multitypic_inf))
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
log.lik.prop <- function(par, age_vect, dat, year, N_sero){
  
  #first run the model with the specified par
  if(N_sero==1){
    
    out.mod <- model.age.one(log.lambda = par, age_vect = age_vect, year=year)
    
  }else{
    
    out.mod <- model.age.incidence(log.lambda = par, age_vect = age_vect, year=year, N_sero=N_sero)  
  }
  
  
  
  
  #to fit, just cut to the dataset that is equal
  out.mod <- subset(out.mod, age<=max(dat$age))
  
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(dat$age)){
    ll=ll+dbinom(dat$cum_cases[i],dat$n[i],p=out.mod$cum_prop_cases[i],log=T) 
  }
  
  return(-ll)
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
                              neg_llik=out.fit$value,
                              convergence=out.fit$convergence)#, 
  
  
  return(out.fit)
}
run.plot.single <- function(mod.dat, dat, age_vect, N_sero){
  
  #and run over many years
  out.run <- model.age.incidence(log.lambda = log(mod.dat$lambda_fit), 
                                 year = mod.dat$year,
                                 age_vect=age_vect, N_sero=N_sero)
  
  out.run.lci <- model.age.incidence(log.lambda = log(mod.dat$lambda_lower), 
                                     year = mod.dat$year,
                                     age_vect=age_vect, N_sero=N_sero)
  
  out.run.uci <- model.age.incidence(log.lambda = log(mod.dat$lambda_upper), 
                                     year = mod.dat$year,
                                     age_vect=age_vect, N_sero=N_sero)
  
  #out.run <- data.table::rbindlist(out.run)
  
  dat.mod <- dplyr::select(out.run, year, age,  cum_prop_cases)
  dat.mod$lci = out.run.lci$cum_prop_cases
  dat.mod$uci = out.run.uci$cum_prop_cases
  dat.mod$type<- "model"
  names(dat.mod) <- c("year", "age", "prop", "lci", "uci", "type")
  
  
  
  #now join with data and plot
  
  dat.out <- sum.yr(df=dat,age_vect = age_vect)
  dat.out$year <- unique(dat$year)
  dat.dat <- dplyr::select(dat.out, year, age, cum_prop_cases)
  dat.dat$lci <- dat.dat$uci <- NA
  dat.dat$type<- "data"
  names(dat.dat) <-  c("year", "age", "prop", "lci", "uci", "type")
  dat.all <- rbind(dat.dat, dat.mod)
  
  p1 <- ggplot(data=dat.all) + 
    geom_point(aes(x=age, y=prop, color=type)) +
    geom_line(aes(x=age, y=prop, color=type, linetype=type)) +
    geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type), alpha=.3)+
    xlim(0,max(age_vect)) + ylab("proportion") +theme_bw()+
    theme(panel.grid = element_blank())
  
  print(p1)
  
  return(dat.all)
  
  
}

#compare fits by year

comp.all.Nsero <- function(lambda.guess, age_vect, dat){
  
  out.NS1 <- wrap.all.llik(lambda.guess = lambda.guess,
                           age_vect = 1:max(dat$age),
                           dat=dat,
                           N_sero=1)
  
  out.NS2 <- wrap.all.llik(lambda.guess = lambda.guess,
                                        age_vect = 1:max(dat$age),
                                        dat=dat,
                                        N_sero=2)
  
  out.NS3 <- wrap.all.llik(lambda.guess = lambda.guess,
                                        age_vect = 1:max(dat$age),
                                        dat=dat,
                                        N_sero=3)
  
  out.NS4 <- wrap.all.llik(lambda.guess = lambda.guess,
                                        age_vect = 1:max(dat$age),
                                        dat=dat,
                                        N_sero=4)
  
  
  out.NS1$N_sero=1
  out.NS2$N_sero=2
  out.NS3$N_sero=3
  out.NS4$N_sero=4
  
  out.all <- rbind(out.NS1, out.NS2, out.NS3, out.NS4)
  
  out.all <- dplyr::select(out.all, -(year))
  
  return(out.all)
  
}

comp.nat <- comp.all.Nsero(lambda.guess=.02, 
                           age_vect=1:37, 
                           dat=dat.nat)

#and plot lambda together
comp.nat$data_source <- "national"
comp.nat$AIC <- 2*(comp.nat$neg_llik) + 2*comp.nat$N_sero
#comp.nat$stAIC<- comp.nat$AIC/mean(comp.nat$AIC)
comp.nat$delta_AIC <- comp.nat$AIC-comp.nat$AIC[comp.nat$AIC==min(comp.nat$AIC)]
comp.nat$log_delta_AIC <- log(comp.nat$delta_AIC)
comp.nat$log_delta_AIC[comp.nat$log_delta_AIC==-Inf] <- 0
#comp.nat$delta_AIC <- comp.nat$stAIC-comp.nat$stAIC[comp.nat$stAIC==min(comp.nat$stAIC)]

out.2019.cohort$data_source <- "cohort"
out.2019.cohort$AIC <- 2*(out.2019.cohort$neg_llik) + 2*out.2019.cohort$N_sero
#out.2019.cohort$stAIC<- out.2019.cohort$AIC/min(out.2019.cohort$AIC)
out.2019.cohort$delta_AIC <- out.2019.cohort$AIC-out.2019.cohort$AIC[out.2019.cohort$AIC==min(out.2019.cohort$AIC)]
out.2019.cohort$log_delta_AIC <- log(out.2019.cohort$delta_AIC)
out.2019.cohort$log_delta_AIC[out.2019.cohort$log_delta_AIC==-Inf] <- 0
#out.2019.cohort$delta_AIC <- out.2019.cohort$stAIC-out.2019.cohort$stAIC[out.2019.cohort$stAIC==min(out.2019.cohort$stAIC)]
out.2019.multi$data_source <- "cohort"
out.2019.multi$AIC <- 2*(out.2019.multi$neg_llik) + 2*out.2019.multi$N_sero
#out.2019.multi$stAIC<- out.2019.multi$AIC/min(out.2019.multi$AIC)
#out.2019.multi$delta_AIC <- out.2019.multi$stAIC-out.2019.multi$stAIC[out.2019.multi$stAIC==min(out.2019.multi$stAIC)]
out.2019.multi$delta_AIC <- out.2019.multi$AIC-out.2019.multi$AIC[out.2019.multi$AIC==min(out.2019.multi$AIC)]
out.2019.multi$log_delta_AIC <- log(out.2019.multi$delta_AIC)
out.2019.multi$log_delta_AIC[out.2019.multi$log_delta_AIC==-Inf] <- 0
#comp.lambda <- rbind(comp.nat, out.2019.cohort)
comp.lambda <- rbind(comp.nat, out.2019.multi)



comp.lambda$N_sero <- as.factor(comp.lambda$N_sero)
comp.lambda$facet <- "comparing cohort vs.\nnational force of infection"
pC <- ggplot(data=comp.lambda) + 
      facet_grid(facet~.) +
      geom_linerange(aes(x=N_sero, ymin=lambda_lower, ymax=lambda_upper, group=data_source, color=N_sero), size=.7) +
      geom_point(aes(x=N_sero, y=lambda_fit, shape=data_source, fill=log_delta_AIC, color=N_sero), size=3, stroke=1) +
      theme_bw() + 
      scale_color_manual(values=scales::hue_pal()(4), name="number of\nserotypes") +
      scale_fill_viridis_c(direction = -1, name="log\n(delta AIC)")+
      scale_shape_manual(values = c(21,24), name="data source") +
      theme(panel.grid = element_blank(),legend.position = c(.7,.7),
            axis.title = element_text(size = 16),
            legend.box = "horizontal",
            legend.background = element_rect(color="black"),
            strip.background = element_rect(fill="white"),
            strip.text = element_text(size = 16),
            axis.text = element_text(size = 14)) + xlab("number of serotypes") +
      ylab(bquote('force of infection, '~lambda)) +
      geom_label(aes(x=.6,y=.25, label="c"), label.size = NA, size=9, fontface="bold") +
  coord_cartesian(ylim = c(0,.27)) + 
  guides( shape=guide_legend(order=1), color=guide_legend(order=2))



pFig2 <- cowplot::plot_grid(ptop, pC, nrow = 2, ncol = 1, rel_heights = c(2,1.2))

pFig2

ggsave(file = paste0(homewd, "/figures/final-figures/Fig2_hold.png"),
       units=c("mm"),  
       width=55, 
       height=85, 
       scale=3, 
       dpi=300)
