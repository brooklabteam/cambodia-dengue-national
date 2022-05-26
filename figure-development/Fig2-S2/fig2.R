

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(scales)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#compare N serotype hypothesis across all years to cumulative case data
#functions
sum.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out[is.na(df.out)]<- 0
  df.out <- rbind(c(0,0), df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- cumsum(df.out$Nage)/sum(df.out$Nage)
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
dat$age <- ceiling(dat$age)
year.split <- dlply(dat, .(year))

dat.sum <- lapply(year.split, sum.yr, age_vect = 1:max(dat$age))
dat.sum <- data.table::rbindlist(dat.sum)
dat.sum$year <- rep(unique(dat$year), each=(max(dat$age)+1))
dat.sum$year <- as.factor(dat.sum$year)

dat.sum$prop_prev <- dat.sum$Nage/dat.sum$n


#first, plot annual proportion of cases by age by year

# gradiant color generator https://colordesigner.io/gradient-generator 

colz = hue_pal()(length(unique(dat.sum$year)))
names(colz) <- unique(as.character(dat.sum$year))
colz["2002"] <- "#cdcdcb"
colz["2003"] <- "#bebebc"
colz["2004"] <- "#afafad"
colz["2005"] <- "#a0a09e"
colz["2006"] <- "#929290"
colz["2008"] <- "#838482"
colz["2009"] <- "#767674"
colz["2010"] <- "#686866"
colz["2011"] <- "#5b5b59"
colz["2013"] <- "#4e4e4c"
colz["2014"] <- "#41413f"
colz["2015"] <- "#353533"
colz["2016"] <- "#292927"
colz["2017"] <- "#1e1e1c"
colz["2018"] <- "#121210"
colz["2020"] <- "#000000"

colz["2007"] <- "orange"
colz["2012"] <- "tomato2"
colz["2019"] <- "red"



pA <- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) +
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  scale_color_manual(values=colz) +
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.2,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))


#then plot the FOI through time
#load the FOI under assumptions of a 2-seroptype sysetm
#eventually will move to data folder, but leaving here for now until we have confidence intervals to match
dat.foi <- read.csv(file = paste0(homewd, "/data/foi-fit-national.csv"), header = T, stringsAsFactors = F)
dat.foi$lambda_per_1000 <- dat.foi$lambda*1000
dat.foi$lambda_per_1000_lci <- dat.foi$lci*1000
dat.foi$lambda_per_1000_uci <- dat.foi$uci*1000
dat.foi$year_plot <- dat.foi$year
dat.foi$year_plot <- as.factor(dat.foi$year_plot)

# pB2 <- ggplot(data = subset(dat.foi, year>=2002)) +theme_bw() +
#   theme(panel.grid = element_blank(), axis.title.x = element_blank(),
#         axis.title.y = element_text(size=14),
#         plot.margin = unit(c(.1,.1,.1,.1), "lines"),
#         axis.text = element_text(size=12)) + scale_color_manual(values=colz) +
#   geom_line(aes(x=year, y=lambda_per_1000)) + 
#   ylab(bquote(lambda~' (/1000 ppl)')) +
#   geom_point(aes(x=year, y=lambda_per_1000, color=year_plot), size=3, show.legend = F) +
#   geom_ribbon(aes(x=year, ymin=lambda_per_1000_lci, ymax=lambda_per_1000_uci), alpha=.3) +
#   geom_linerange(aes(x=year, ymin=lambda_per_1000_lci, ymax=lambda_per_1000_uci, color=year_plot), size=1, show.legend = F)

pB2 <- ggplot(data = subset(dat.foi, year>=2002)) +theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        plot.margin = unit(c(.1,.1,.1,.1), "lines"),
        axis.text = element_text(size=12)) + scale_color_manual(values=colz) +
  geom_line(aes(x=year, y=lambda)) + 
  ylab(bquote(lambda)) +
  geom_point(aes(x=year, y=lambda, color=year_plot), size=3, show.legend = F) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci), alpha=.3) +
  geom_linerange(aes(x=year, ymin=lci, ymax=uci, color=year_plot), size=1, show.legend = F)




library(ggmap)
colz.long <- c(rep("#cdcdcb", 21), colz)
names(colz.long) <- 1981:2020
# pB1 <- ggplot(data = subset(dat.foi, year>=1980)) +theme_bw() +
#   theme(panel.grid = element_blank(), axis.title.x = element_blank(),
#         axis.title.y = element_text(size=16),
#         axis.text = element_text(size=14)) + scale_color_manual(values=colz.long) +
#   geom_line(aes(x=year, y=lambda_per_1000)) + 
#   geom_vline(aes(xintercept=2001.5), linetype=2, color="red") +
#   ylab(bquote(lambda~',force of infection (/1000 ppl)')) +
#   geom_point(aes(x=year, y=lambda_per_1000, color=year_plot), size=3, show.legend = F) +
#   geom_ribbon(aes(x=year, ymin=lambda_per_1000_lci, ymax=lambda_per_1000_uci), alpha=.3) 
#   #coord_cartesian(ylim=c(0,1000))

pB1 <- ggplot(data = dat.foi) +theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=14)) + scale_color_manual(values=colz.long) +
  geom_line(aes(x=year, y=lambda)) + 
  geom_vline(aes(xintercept=2001.5), linetype=2, color="red") +
  ylab(bquote(lambda~',force of infection (per capita)')) +
  geom_point(aes(x=year, y=lambda, color=year_plot), size=3, show.legend = F) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci), alpha=.3) 
#coord_cartesian(ylim=c(0,1000))


# 
# pB <- pB1 + annotation_custom(ggplotGrob(pB2), xmin = 2002, xmax = 2021, ymin = 200, ymax = 470) + 
#       theme(plot.margin = unit(c(.5,.5,1.8,.5), "lines"),)

pB <- pB1 + annotation_custom(ggplotGrob(pB2), xmin = 2002, xmax = 2021, ymin = .5, ymax = .98) + 
  theme(plot.margin = unit(c(.5,.5,1.8,.5), "lines"),)



#and part C- which is the actual model fit to the data for the three epidemic years
# then plot the FOI with the data for a few choice epidemic years
# the rest can be found in the supplementary material

plot.model.series.data.whole <- function(par.dat, age_vect, dat){
  
  
  #first, prep the data
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  dist.back <- max(dat$age[dat$year==min(dat$year)])#22
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, age_vect = age_vect, year.start = min(par.dat$year))
  head(out.mod)
  
  #and also run it at the uci and lci
  par.dat.lci <- par.dat.uci <- par.dat
  par.dat.lci$lambda <- par.dat.lci$lci
  par.dat.uci$lambda <- par.dat.uci$uci
  
  out.mod.lci <- model.age.incidence.series(par.dat = par.dat.lci, age_vect = age_vect, year.start = min(par.dat$year))
  out.mod.uci <- model.age.incidence.series(par.dat = par.dat.uci, age_vect = age_vect,year.start = min(par.dat$year))
  
  head(out.mod.lci)
  head(out.mod.uci)
  
  out.mod.lci <- subset(out.mod.lci, year >= min(df.out$year))
  out.mod.uci <- subset(out.mod.uci, year >= min(df.out$year))
  out.mod <- subset(out.mod, year >= min(df.out$year))
  
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  
  df.mod$cum_prop_lci <- out.mod.lci$cum_prop_cases
  df.mod$cum_prop_uci <- out.mod.uci$cum_prop_cases
  
  df.dat$cum_prop_lci <- df.dat$cum_prop_uci <- NA
    
  df.all<- rbind(df.mod, df.dat)
  
  
  colz = c("model" = "turquoise3", "data" = "tomato")
  p1 <- ggplot(data=df.all) + geom_point(aes(x=age, y=cum_prop_cases, color=type))  + 
    geom_line(aes(x=age, y=cum_prop_cases, color=type)) + scale_color_manual(values=colz)+ scale_fill_manual(values=colz)+
    geom_ribbon(data=subset(df.all, type=="model"), aes(x=age, ymin=cum_prop_lci, ymax=cum_prop_uci, fill=type), alpha=.3) + facet_wrap(~year) + 
    theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(),
                       legend.position = c(.9,.1), axis.title = element_text(size=14),
                       axis.text = element_text(size=12), legend.text = element_text(size=10)) +
    ylab("cumulative proportion of cases")
  print(p1)
  
  
  
  return(p1)
  
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
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
model.age.incidence.series <- function(par.dat, age_vect, year.start){
  
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
  year_tracker[1] <- year.start
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


#first, the entire time series in the supplement
pS2 <- plot.model.series.data.whole(par.dat=dat.foi,
                             age_vect=seq(0,22, by=1/4), 
                             dat=dat)
ggsave(file = paste0(homewd, "/final-figures/figS2.png"),
       plot=pS2,
       units=c("mm"),  
       width=85, 
       height=65, 
       scale=3, 
       dpi=300)


#then, just our three epidemic years
plot.model.series.data.epi <- function(par.dat, age_vect, dat){
  
  #first, prep the data
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  dist.back <- max(dat$age[dat$year==min(dat$year)])#22
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.yr, age_vect = age_vect_year)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, age_vect = age_vect, year.start = min(par.dat$year))
  head(out.mod)
  
  #and also run it at the uci and lci
  par.dat.lci <- par.dat.uci <- par.dat
  par.dat.lci$lambda <- par.dat.lci$lci
  par.dat.uci$lambda <- par.dat.uci$uci
  
  out.mod.lci <- model.age.incidence.series(par.dat = par.dat.lci, age_vect = age_vect, year.start = min(par.dat$year))
  out.mod.uci <- model.age.incidence.series(par.dat = par.dat.uci, age_vect = age_vect, year.start = min(par.dat$year))
  
  head(out.mod.lci)
  head(out.mod.uci)
  
  out.mod.lci <- subset(out.mod.lci, year >= min(df.out$year))
  out.mod.uci <- subset(out.mod.uci, year >= min(df.out$year))
  out.mod <- subset(out.mod, year >= min(df.out$year))
  
  df.dat <- dplyr::select(df.out, year, age, cum_prop_cases)
  df.dat$type="data"
  
  df.mod <-   dplyr::select(out.mod, year, age, cum_prop_cases)
  
  df.mod$type="model"
  
  df.mod$cum_prop_lci <- out.mod.lci$cum_prop_cases
  df.mod$cum_prop_uci <- out.mod.uci$cum_prop_cases
  
  df.dat$cum_prop_lci <- df.dat$cum_prop_uci <- NA
  
  df.all<- rbind(df.mod, df.dat)
  
  df.epi = subset(df.all, year==2007 | year==2012 | year==2019)
  
  
  
  shapez = c("model" = 24, "data" = 21)
  linez = c("model" = 2, "data" = 1)
  

  
  yearz = c('2007' = "orange", '2012' = "tomato2", '2019'="red")
  df.epi$year <- as.factor(df.epi$year)
  p1 <- ggplot(data=df.epi) + geom_point(aes(x=age, y=cum_prop_cases, shape=type, fill=year, color=year), size=3)  + 
    geom_line(aes(x=age, y=cum_prop_cases, linetype=type, color=year)) + scale_shape_manual(values=shapez)+ 
    scale_linetype_manual(values=linez)+scale_color_manual(values=yearz, guide="none") + scale_fill_manual(values=yearz, guide="none")+
    geom_ribbon(data=subset(df.epi, type=="model"), aes(x=age, ymin=cum_prop_lci, ymax=cum_prop_uci, fill =year), alpha=.3) + facet_grid(year~.) + 
    theme_bw() + theme(panel.grid = element_blank(), legend.title = element_blank(), strip.background = element_rect(fill="white"),
                       legend.position = c(.85,.05), axis.title = element_text(size=16),strip.text = element_text(size=14),
                       axis.text = element_text(size=14), legend.text = element_text(size=12)) +
    ylab("annual cumulative proportion of cases")
  print(p1)
  
  
  
  return(p1)
  
} 

pC <- plot.model.series.data.epi(par.dat=dat.foi,
                                    age_vect=seq(0,22, by=1/4), 
                                    dat=dat)


#compile
left <- cowplot::plot_grid(pA,pB, nrow = 2, ncol=1, rel_heights = c(1,1), labels = c("a", "b"), label_size = 22)


pFig2 <- cowplot::plot_grid(left, pC, ncol = 2, nrow=1, rel_widths = c(1.3,1), labels = c("", "c"), label_size = 22)# +cowplot::draw_text('epi years of\n\t2007/2012/2019', x = 0.945, y = 0.1, size = 10, hjust = 0.5, vjust = 0.5)




ggsave(file = paste0(homewd, "/final-figures/fig2.png"),
       plot=pFig2,
       units=c("mm"),  
       width=100, 
       height=75, 
       scale=3, 
       dpi=300)


