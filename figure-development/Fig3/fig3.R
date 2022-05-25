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

#functions
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
plot.model.prim.multi.data <- function(par.dat, age_vect, dat, panel){
  
  
  
  #first, prep the data
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  df.out = sum.dat.igg(dat = dat, age_vect = age_vect_year)
  df.out$year <- max(dat$year)
  
  # ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #        geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # 
  
  #and run the model
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, age_vect = age_vect)
  head(out.mod)
  
  par.uci <- par.dat
  par.uci$lambda[par.uci$year==2019] <- par.uci$uci[par.uci$year==2019]
  out.mod.uci <- model.age.incidence.series(par.dat = par.uci, age_vect = age_vect)
  
  par.lci <- par.dat
  par.lci$lambda[par.lci$year==2019] <- par.lci$lci[par.lci$year==2019]
  out.mod.lci <- model.age.incidence.series(par.dat = par.lci, age_vect = age_vect)
  
  
  out.mod <- subset(out.mod, year == unique(df.out$year))
  out.mod.lci <- subset(out.mod.lci, year == unique(df.out$year))
  out.mod.uci <- subset(out.mod.uci, year == unique(df.out$year))
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  df.dat$lci <- df.dat$uci <- NA
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  df.mod$lci <- out.mod.lci$cum_prop_cases
  df.mod$uci <- out.mod.uci$cum_prop_cases
  
  df.all<- rbind(df.mod, df.dat)
  df.all$infection <- "multitypic"
  
  
  #and also get primary cases
  df.dat2 <- dplyr::select(df.out, year, age, prop_primary_by_age)
  df.dat2$type="data"
  
  df.mod2 <-   dplyr::select(out.mod, year, age, all_prim)
  
  df.mod2$type="model"
  
  df.mod2$lci <- out.mod.lci$all_prim
  df.mod2$uci <- out.mod.uci$all_prim
  
  df.dat2$lci <- df.dat2$uci <- NA
  
  names(df.mod2) <- names(df.dat2)
  df.all2 <- rbind(df.dat2, df.mod2)
  
  df.all2$infection <- "primary"
  
  names(df.all) <- names(df.all2) <- c("year", "age", "proportion", "type", "lci", "uci", "infection")
  df.all <- rbind(df.all, df.all2)
  
  shapez = c("model" = 24, "data" = 21)
  linez = c("model" = 2, "data" = 1)
  colz = c("model" = "tomato2", "data" = "red")
  
  if(panel=="a"){
    

  pa <- ggplot(data=subset(df.all, infection == "multitypic")) + 
    geom_point(aes(x=age, y=proportion, shape=type, fill=type), show.legend = F)  + 
    geom_line(aes(x=age, y=proportion, linetype=type, color=type), show.legend = F) + facet_grid(~infection) + 
    geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type), alpha=.3, show.legend = F) +
    theme_bw() + scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
    scale_shape_manual(values=shapez) + scale_linetype_manual(values= linez) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
          strip.background = element_rect(fill="white"), strip.text = element_text(size = 16),
          axis.text = element_text(size = 14), legend.title = element_blank(), 
          legend.text = element_text(size=12), plot.tag = element_text(size=22, face="bold"),
          legend.position = c(.85,.85)) +
    labs(tag="a")
  
  return(pa)
  }else if (panel=="b"){
    

  
  pb <- ggplot(data=subset(df.all, infection == "primary")) + 
    geom_point(aes(x=age, y=proportion, shape=type, fill=type))  + 
    geom_line(aes(x=age, y=proportion, linetype=type, color=type)) + facet_grid(~infection) + 
    geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type),  alpha=.3) +
    theme_bw() + scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
    scale_shape_manual(values=shapez) + scale_linetype_manual(values= linez) +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
          strip.background = element_rect(fill="white"), strip.text = element_text(size = 16),
          axis.text = element_text(size = 14), legend.title = element_blank(), 
          legend.text = element_text(size=12), plot.tag = element_text(size=22, face="bold"),
          legend.position = c(.85,.85)) + coord_cartesian(ylim=c(0,1))+
    labs(tag="b")
  
  return(pb)
  
  }
  
} 
plot.model.data <- function(par.dat, age_vect, dat, panel){
  
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  df.out <- sum.yr.yr(dat, age_vect = age_vect_year)
  
  # 
    ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
           geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # # 
  
  #and run the model
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, age_vect = age_vect)
  head(out.mod)
  
  par.uci <- par.dat
  par.uci$lambda[par.uci$year==2019] <- par.uci$uci[par.uci$year==2019]
  out.mod.uci <- model.age.incidence.series(par.dat = par.uci, age_vect = age_vect)
  
  par.lci <- par.dat
  par.lci$lambda[par.lci$year==2019] <- par.lci$lci[par.lci$year==2019]
  out.mod.lci <- model.age.incidence.series(par.dat = par.lci, age_vect = age_vect)
  
  
  out.mod <- subset(out.mod, year == unique(df.out$year))
  out.mod.lci <- subset(out.mod.lci, year == unique(df.out$year))
  out.mod.uci <- subset(out.mod.uci, year == unique(df.out$year))
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  df.dat$lci <- df.dat$uci <- NA
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  df.mod$lci <- out.mod.lci$cum_prop_cases
  df.mod$uci <- out.mod.uci$cum_prop_cases
  
  df.all<- rbind(df.mod, df.dat)
  #df.all$infection <- "multitypic"
  
  
  
  names(df.all) <- c("year", "age", "proportion", "type", "lci", "uci")
  
  
  shapez = c("model" = 24, "data" = 21)
  linez = c("model" = 2, "data" = 1)
  colz = c("model" = "tomato2", "data" = "red")
  
  
    
    pc <- ggplot(data=df.all) + 
      geom_point(aes(x=age, y=proportion, shape=type, fill=type), show.legend = F)  + 
      geom_line(aes(x=age, y=proportion, linetype=type, color=type), show.legend = F) + 
      geom_ribbon(aes(x=age, ymin=lci, ymax=uci, fill=type), alpha=.3, show.legend = F) +
      theme_bw() + scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
      scale_shape_manual(values=shapez) + scale_linetype_manual(values= linez) +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
            strip.background = element_rect(fill="white"), strip.text = element_text(size = 16),
            axis.text = element_text(size = 14), legend.title = element_blank(), 
            legend.text = element_text(size=12), plot.tag = element_text(size=22, face="bold"),
            legend.position = c(.85,.85)) +
      labs(tag="c")
    
    
  
} 


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


#load the parameters - combine 
ks.dat <- read.csv(file = paste0(homewd, "/data/foi-kampong-speu-2019.csv"), header = T, stringsAsFactors = F)
ks.dat <-dplyr::select(ks.dat, lambda_KP, lambda_KP_lci, lambda_KP_uci)
#and add in the national estimates for previous years
nat.dat <- read.csv(file = paste0(homewd, "/data/foi-fit-national.csv"), stringsAsFactors = F, header = T)
par.dat = nat.dat
par.dat$lambda[par.dat$year==2019] <- ks.dat$lambda_KP
par.dat$lci[par.dat$year==2019] <- ks.dat$lambda_KP_lci
par.dat$uci[par.dat$year==2019] <- ks.dat$lambda_KP_uci

#run the model with the data and plot
pA <- plot.model.prim.multi.data(par.dat=  par.dat,
                                 age_vect=seq(0,22, by=1/4),
                                 dat= dat.2019,
                                 panel="a")
pB <- plot.model.prim.multi.data(par.dat=  par.dat,
                                 age_vect=seq(0,22, by=1/4),
                                 dat= dat.2019,
                                 panel="b")


#now look at cases by age for 2019 from the national data for Kampong Speu
ks.nat.foi <- read.csv(file = paste0(homewd, "/data/foi-kampong-speu-2019-national.csv"), header = T, stringsAsFactors = F)
par.dat = nat.dat
par.dat$lambda[par.dat$year==2019] <- ks.nat.foi$lambda_KP
par.dat$lci[par.dat$year==2019] <- ks.nat.foi$lambda_KP_lci
par.dat$uci[par.dat$year==2019] <- ks.nat.foi$lambda_KP_uci

#add in national KS data

#load the age structured national data
dat <- read.csv(file = paste0(homewd, "/data/DENV-KPS-Aged-2019.csv") , header = T, stringsAsFactors = F)
# dat
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$epiwk <- as.Date(dat$epiwk, format = "%m/%d/%y")
dat <- arrange(dat, date, age)

dat$age <- dat$age + dat$month/12
unique(dat$age) #round to years ()
#dat$age <- round(dat$age, 0)
dat$age <- ceiling(dat$age)

pC  <- plot.model.data(par.dat=  par.dat,
                        age_vect=seq(0,22, by=1/4),
                        dat= dat)




#and, finally, compare all the FOI estimates
#and the fitted Kampong Speu National data

ks.nat.foi <- dplyr::select(ks.nat.foi, lambda_KP, lambda_KP_lci, lambda_KP_uci)
nat.dat = subset(nat.dat, year==2019)
nat.dat <- dplyr::select(nat.dat, lambda, lci, uci)
nat.dat$data_source <- "National\nData"
ks.dat$data_source  <- "Kampong Speu\nCohort Data"
ks.nat.foi$data_source  <- "Kampong Speu\nNational Data"
names(ks.nat.foi)  <- names(ks.dat) <- names(nat.dat)

foi.all <- rbind(nat.dat, ks.dat, ks.nat.foi)
foi.all$data_source <- factor(foi.all$data_source, levels = c("National\nData",
                                                              "Kampong Speu\nNational Data",
                                                              "Kampong Speu\nCohort Data"))

#foi.all$facet <- "Kampong Speu vs.\nNational comparison"
pD <- ggplot(data=foi.all) + 
      geom_linerange(aes(x=data_source, ymin=lci, ymax=uci), size=.7) +
      geom_point(aes(x=data_source, y=lambda, shape=data_source), size=3, stroke=1, show.legend = F) +
      theme_bw() + 
      scale_shape_manual(values = c(24,21, 22), name="data source") +
      theme(panel.grid = element_blank(),#legend.position = c(.65,.7),
            axis.title.y = element_text(size = 16),
            #legend.box = "horizontal",
            #legend.background = element_rect(color="black"),
            strip.background = element_rect(fill="white"),
            strip.text = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.title.x = element_blank()) + #xlab("number of serotypes") +
      ylab(bquote('force of infection ('~lambda~'), 2019')) +
      geom_label(aes(x=.6,y=.25, label="d"), label.size = NA, size=9, fontface="bold") +
  coord_cartesian(ylim = c(0,.27)) 
pD


pFig2 <- cowplot::plot_grid(ptop, pC, nrow = 2, ncol = 1, rel_heights = c(2,1.2))

pFig2

ggsave(file = paste0(homewd, "/final-figures/fig3.png"),
       units=c("mm"),  
       width=55, 
       height=85, 
       scale=3, 
       dpi=300)


#and 2020

comp.nat.2020$data_source <- "National"
comp.nat.2020$AIC <- 2*(comp.nat.2020$neg_llik) + 2*comp.nat.2020$N_sero
comp.nat.2020$delta_AIC <- comp.nat.2020$AIC-comp.nat.2020$AIC[comp.nat.2020$AIC==min(comp.nat.2020$AIC)]
comp.nat.2020$log_delta_AIC <- log(comp.nat.2020$delta_AIC)
comp.nat.2020$log_delta_AIC[comp.nat.2020$log_delta_AIC==-Inf] <- 0


out.2020.cohort$data_source <- "Kampong Speu"
out.2020.cohort$AIC <- 2*(out.2020.cohort$neg_llik) + 2*out.2020.cohort$N_sero
out.2020.cohort$delta_AIC <- out.2020.cohort$AIC-out.2020.cohort$AIC[out.2020.cohort$AIC==min(out.2020.cohort$AIC)]
out.2020.cohort$log_delta_AIC <- log(out.2020.cohort$delta_AIC)
out.2020.cohort$log_delta_AIC[out.2020.cohort$log_delta_AIC==-Inf] <- 0
out.2020.multi$data_source <- "Kampong Speu"
out.2020.multi$AIC <- 2*(out.2020.multi$neg_llik) + 2*out.2020.multi$N_sero
out.2020.multi$delta_AIC <- out.2020.multi$AIC-out.2020.multi$AIC[out.2020.multi$AIC==min(out.2020.multi$AIC)]
out.2020.multi$log_delta_AIC <- log(out.2020.multi$delta_AIC)
out.2020.multi$log_delta_AIC[out.2020.multi$log_delta_AIC==-Inf] <- 0
comp.lambda.2020 <- rbind(comp.nat.2020, out.2020.multi)



comp.lambda.2020$N_sero <- as.factor(comp.lambda.2020$N_sero)
comp.lambda.2020$facet <- "Kampong Speu vs.\nNational comparison"
pC2020 <- ggplot(data=comp.lambda.2020) + 
  facet_grid(facet~.) +
  geom_linerange(aes(x=N_sero, ymin=lambda_lower, ymax=lambda_upper, group=data_source, color=N_sero), size=.7) +
  geom_point(aes(x=N_sero, y=lambda_fit, shape=data_source, fill=log_delta_AIC, color=N_sero), size=3, stroke=1) +
  theme_bw() + 
  scale_color_manual(values=scales::hue_pal()(4), name="number of\nserotypes") +
  scale_fill_viridis_c(direction = -1, name="log\n(delta AIC)")+
  scale_shape_manual(values = c(24,21), name="data source") +
  theme(panel.grid = element_blank(),legend.position = c(.65,.7),
        axis.title = element_text(size = 16),
        legend.box = "horizontal",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14)) + xlab("number of serotypes") +
  ylab(bquote('annual force of infection, '~lambda)) +
  geom_label(aes(x=.6,y=.25, label="c"), label.size = NA, size=9, fontface="bold") +
  coord_cartesian(ylim = c(0,.27)) + 
  guides( shape=guide_legend(order=1), color=guide_legend(order=2))




pFig2extra <- cowplot::plot_grid(ptop2020, pC2020, nrow = 2, ncol = 1, rel_heights = c(2,1.2))

pFig2extra # most support for two serotypes in local data... but so few datapoints
