# This file is for figure 3. There are 4 subplots to support the idea of two-strian

######################################################################################
#### Fig A | TSIR with increased S,time vs case ######################################
######################################################################################


rm(list = ls())


library(ggplot2)
library(reshape2)
library(pracma)
library(astsa)
library(sp)
library(GADMTools)
library(plotrix)
library(tsiR)
library(lubridate)
library(plm)
library(interp)
library(foreign)
library(gplots)
library(lfe)
library(zoo)
library(e1071)
library(xts)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(graphics)
library(stargazer)
library(grid)
library(gridExtra)
library(png)
require("knitr")

#set your home directory here
homewd <- "/home/rstudio" #yimei version
homewd <- "/Users/carabrook/Developer/cambodia-dengue-national" #cara version

# setwd(homewd)



# set up size and color
xy_size<-13
my_red="tomato"
my_blue<-"#2980B9"
my_orange<-"#D35400"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
color_2007<-"purple"
color_2012<-"turquoise3"
color_2019<-"tomato"




load(file=paste0(homewd,"/data/sim.all.figa.RData"))

dat=sim.all
dat$variable<-sub("data", "reported case",dat$variable) 
dat$variable = factor(dat$variable, levels=c("reported case",
                                             "TSIR fit",
                                             "prediction TSIR"))
colz = c('reported case' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green)
typez = c('reported case' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1)


all.comp<-readRDS(paste0(homewd,"/data/all.comp.RDS"))

dat=all.comp


dat$variable<-sub("-", "- \n ",dat$variable)  
dat$variable<-sub("data", "reported case",dat$variable) 

dat$variable = factor(dat$variable, levels=c("reported case",
                                             "TSIR fit",
                                             "prediction TSIR",
                                             "prediction TSIR- \n increased S"))
colz = c('reported case' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR- \n increased S' =my_red)
typez = c('reported case' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR- \n increased S' = 1)
# Fig B | TSIR with increased S
pl_a <- ggplot(data=dat, aes(x=time, y=reported_cases, color=variable, linetype=variable)) + 
  geom_line(data=dplyr::filter(dat, variable=="reported case"), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time<2007), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2008 & time<2012), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR- \n increased S" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR- \n increased S" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR- \n increased S" & time>=2019 & time<2020), size=1) +
  scale_color_manual(values=colz) +
  scale_linetype_manual(values=typez) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=12),
        legend.position = c(.15,0.86),
        # legend.key.size = unit(0.7, "cm"),
        axis.text = element_text(size=14), legend.title = element_blank(),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+theme(legend.spacing.y = unit(0, "cm"))+ rremove("x.title")+ labs(y = "Cases")






######################################################################################
#### Fig B | log ######################################
######################################################################################



sim.fit.return.beta <- function(time.start, dat, epiyr, family){
  
  
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
  
  
  
  df.out = fittedpars$contact
  df.out$epiyr <- epiyr
  
  
  
  return(df.out)
}

#load data into tsir form
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

#first try to predict 2007-epidemic
#here we fit to 2002:2007

## @Yimei, I know you made some changes to the fitting script for TSIR
## Below assumes regtype = "lm" and family = "gaussian".
## Please make sure it matches

sim.2007 <- sim.fit.return.beta(dat = tsir.dat,
                                time.start =  min(tsir.dat$time),
                                family="gaussian",
                                epiyr = 2007)


#now for 2012
sim.2012 <- sim.fit.return.beta(dat = tsir.dat,
                                time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                                family="gaussian",
                                epiyr = 2012)


#and the results for 2019
sim.2019 <- sim.fit.return.beta(dat = tsir.dat,
                                time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                                family="gaussian",
                                epiyr = 2019)

#and combine
fit.all <- rbind(sim.2007, sim.2012, sim.2019)
fit.all$epiyr[fit.all$epiyr==2007] <- "2002-2006"
fit.all$epiyr[fit.all$epiyr==2012] <- "2008-2011"
fit.all$epiyr[fit.all$epiyr==2019] <- "2013-2018"
fit.all$epiyr <- as.factor(fit.all$epiyr)


p1 <- ggplot(data=fit.all) +#facet_grid(epiyr~., scales = "free_y") +
  geom_line(aes(x=time, y=beta, color=epiyr)) +
  geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)
p1




#or as lines and points
pl_b <- ggplot(data=fit.all) + #facet_grid(epiyr~., scales = "free_y") +
  geom_point(aes(x=time, y=beta, color=epiyr)) +
  geom_linerange(aes(x=time, ymin=betalow, ymax=betahigh, color=epiyr), alpha=.3)+
  geom_line(aes(x=time, y=beta, color=epiyr)) +
  #geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)+
  theme_bw() + scale_y_log10()+
  scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                     name="fitting period") +
  scale_fill_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                    name="fitting period") + 
  ylab(bquote(log[10](beta)~',transmission rate')) +
  #xlab("week of year")+
  scale_x_continuous(breaks=c(1*2,3*2,5*2,7*2,9*2,11*2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  theme(panel.grid = element_blank(), 
        legend.position = c(.128,.885),
        legend.text = element_text(size=14),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.text = element_text(size=14), legend.title = element_blank(),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

pl_b













######################################################################################
#### Fig C | time & sus & color birth rate ######################################
######################################################################################



Sdatcombined <- readRDS(paste0(homewd, "/data/Sdatcombined.rds"))
colz = c("other years"="black", "epidemic years"="red")
Sdatcombined



Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
  
}
Sdatcombined




################################################################
# Fig 3C, time vs sus, birth change to birth/1000 ppl
################################################################


find_births<-read.csv(paste0(homewd,"/data/tsir_data_birth_updated.csv"))


# choose a number every 26 rows
births_update<- cbind.data.frame(seq(2002,2020,1), sapply(split(find_births$births_per_1000, rep(1:(nrow(find_births)/26), each=26)), mean))
colnames(births_update)<-c("time","births_per1000")

Sdatcombined_update<-merge(x = Sdatcombined, y =births_update, by = "time", all.x = TRUE)

Sdatcombined = subset(Sdatcombined_update, select = -c(births) )
names(Sdatcombined)[names(Sdatcombined) == 'births_per1000'] <- 'births per 1000 ppl'
Sdatcombined


colz = c("other years"="black", "epidemic years"="red")

Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
  
}

# remove TSIR and only show TSIR-increased S

Sdatcombined<-Sdatcombined[!(Sdatcombined$sim %in% "TSIR"),]




Sdatcombined$variable<-sub("-", "- \n ",Sdatcombined$variable)  





a<-c("black","red","red","black","red")

my_blue<-"#2980B9"
my_green<-"green3"
colz = c("prediction TSIR"=my_green, "prediction TSIR- \n increased S"="tomato","TSIR fit"=my_blue)
options(scipen=10000) # we don't want to plot with 1e+00


pl_c<-ggplot(data=Sdatcombined) + 
  geom_point(aes(x=year, y=sus_mean, 
                 color=variable, shape=variable), size=5, stroke = 1) + labs(shape=NULL, colour=NULL) +# only remove two legend
  scale_color_manual(values=colz) +
  
  ggnewscale::new_scale_color() +
  geom_point(aes(x=year, y=sus_mean, 
                 color=`births per 1000 ppl`, shape=variable), size=3) + 
  
  theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white")) +
  labs(x ="Year", y = "Susceptible Count")+
  theme(
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14))+
  scale_y_log10() + scale_color_viridis_c()+ geom_vline(xintercept = c(2007,2012,2019) , linetype=2, 
                                                        color = "red", size=0.5)+
  scale_x_continuous(limits=c(2002, 2020), breaks=c(2002,2007,2012,2015,2019))+
  theme(axis.text.x = element_text( hjust = 0.5, colour = a),strip.text = element_text(colour = 'black',size = 14))  +
  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size=16),
    legend.text = element_text(size=12),legend.title=element_text(size=14),
    legend.position = c(.847,.723),legend.box.background = element_rect(colour = "black"),legend.background = element_blank(),
    axis.text = element_text(size=14))+# rremove("xy.title")+# theme(legend.title = element_blank()) + 
  theme(legend.spacing.y = unit(0, "cm"))+ # space between legend rows
  theme(legend.title.align = 1.08)








######################################################################################
#### Fig D | reported case & predicted case TSIR & increasedS ########################
######################################################################################




#load data into tsir form

#homewd<-"/home/rstudio"
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)


# this part is for increased S

sim.com.epi <- function(time.start, dat, epiyr, family) {
  if (family == "poisson") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'poisson',
        link = 'log'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        method = 'pois',
        nsim = 100
      )
    
  } else if (family == "gaussian") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'gaussian'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        #method='pois',
        nsim = 100
      )
    
    
  }
  
  
  
  
  
  
  Its = simfitted$res$mean * fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
  Ifinal <- Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the
  #closest match of model to observed cases
  
  fracincvec <- seq(1, 3, 0.01)
  #fracincdist<-matrix(0,1,1)
  ratioImax <- matrix(0, 1, length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  
  for (i in 1:length(fracincvec)) {
    Sbegin = simfitted$sbar + simfitted$Z
    Sepi <- Sbegin[length(Sbegin)]
    SepiInc <- fracincvec[i] * Sbegin[length(Sbegin)]
    
    
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time < epiyr)
    
    predict_ts <- predicttsir(
      times = dat.fit$time,
      births = dat.fit$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sbegin[1],
      I0 = dat.fit$cases[1],
      nsim = 100,
      stochastic = T
    )
    
    
    
    #and the epi year-no increase
    dat.pred = subset(dat, time >= epiyr & time <= (epiyr + 1))
    predict_epi_noInc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sepi,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    #and including an increase
    predict_epi_Inc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = SepiInc,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    # print("predict_epi_Inc:")
    # print(predict_epi_Inc)
    
    
    IPredEpi = melt(
      cbind.data.frame(
        time = predict_epi_Inc$I$time,
        mean_noInc = predict_epi_noInc$I$mean,
        mean_Inc = predict_epi_Inc$I$mean
      ),
      id.vars = "time"
    )
    ISimTS = cbind.data.frame(
      time = predict_ts$I$time,
      variable = rep("TSIR fit", length(predict_ts$I$time)),
      value = predict_ts$I$mean
    )
    
    # print("IPredEpi:")
    # print(IPredEpi)
    # print("ISimTS:")
    # print(ISimTS)
    
    Iall <- rbind(IPredEpi, ISimTS)
    Iall$variable <- as.character(Iall$variable)
    Iall$variable[Iall$variable == "mean_Inc"] <-
      "prediction TSIR-increased S"
    Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
    
    all.dat <- subset(dat, time <= (epiyr + 1))
    all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
    ggplot(data = subset(Iall, variable != "prediction TSIR-increased S")) +
      geom_line(data = all.dat,
                aes(x = time, y = cases),
                linetype = 2) +
      geom_line(aes(x = time, y = value, color = variable))
    
    ggplot(data = Iall) +
      geom_line(data = all.dat,
                aes(x = time, y = cases),
                linetype = 2) +
      geom_line(aes(x = time, y = value, color = variable))
    
    
    
    ratioImax[i] <-
      max(Iall$value[Iall$variable == "prediction TSIR-increased S"]) / max(all.dat$cases)
  }
  
  
  spot <- which(abs(ratioImax - 1) == min(abs(ratioImax - 1)))
  fracincout <- fracincvec[spot]
  ratioImaxout <- ratioImax[spot]
  S_noinc <- Sepi
  
  dat.out <-
    cbind.data.frame(epidemic_year = epiyr,
                     epiSnoic = S_noinc,
                     frac_incS = fracincout)
  
  
  return(dat.out)
}


# no increase part & plot

sim.com.noS <- function(time.start, dat, epiyr, family) {
  if (family == "poisson") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'poisson',
        link = 'log'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        method = 'pois',
        nsim = 100
      )
    
  } else if (family == "gaussian") {
    fittedpars <-
      estpars(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        alpha = NULL,
        sbar = NULL,
        xreg = "cumcases",
        regtype = 'lm',
        family = 'gaussian'
      )
    
    simfitted <-
      simulatetsir(
        data = subset(dat, time >= time.start & time < epiyr),
        IP = 2,
        parms = fittedpars,
        #epidemics='break', threshold=3,
        #method='pois',
        nsim = 100
      )
    
    
  }
  
  
  
  
  
  
  Its = simfitted$res$mean * fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
  Ifinal <- Its[length(Its)]
  
  
  
  Sbegin = simfitted$sbar + simfitted$Z
  Sepi <- Sbegin[length(Sbegin)]
  
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time < epiyr)
  
  predict_ts <- predicttsir(
    times = dat.fit$time,
    births = dat.fit$births,
    beta = fittedpars$beta,
    alpha = fittedpars$alpha,
    S0 = Sbegin[1],
    I0 = dat.fit$cases[1],
    nsim = 100,
    stochastic = T
  )
  
  
  
  #and the epi year-no increase
  dat.pred = subset(dat, time >= epiyr & time <= (epiyr + 1))
  predict_epi_noInc <- predicttsir(
    times = dat.pred$time,
    births = dat.pred$births,
    beta = fittedpars$beta,
    alpha = fittedpars$alpha,
    S0 = Sepi,
    I0 = Ifinal,
    nsim = 100,
    stochastic = T
  )
  
  
  
  IPredEpi = melt(
    cbind.data.frame(time = predict_epi_noInc$I$time, mean_noInc = predict_epi_noInc$I$mean),
    id.vars = "time"
  )
  ISimTS = cbind.data.frame(
    time = predict_ts$I$time,
    variable = rep("TSIR fit", length(predict_ts$I$time)),
    value = predict_ts$I$mean
  )
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  #Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time <= (epiyr + 1))
  all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) +
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable))
  
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases / mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall) == "value"] <- "cases"
  Iall$reported_cases <- Iall$cases / mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data = all.out) +
    geom_line(aes(x = time, y = cases, color = variable))
  
  return(all.out)
  
  
  
}
#here, run fitted version with and without increased S
sim.with.increaseS <-
  function(fracincS, time.start, dat, epiyr, family) {
    if (family == "poisson") {
      fittedpars <-
        estpars(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          alpha = NULL,
          sbar = NULL,
          xreg = "cumcases",
          regtype = 'lm',
          family = 'poisson',
          link = 'log'
        )
      
      simfitted <-
        simulatetsir(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          parms = fittedpars,
          #epidemics='break', threshold=3,
          method = 'pois',
          nsim = 100
        )
      
    } else if (family == "gaussian") {
      #
      fittedpars <-
        estpars(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          alpha = NULL,
          sbar = NULL,
          xreg = "cumcases",
          regtype = 'lm',
          family = 'gaussian'
        )
      #
      # #now simulate tsir
      
      #
      simfitted <-
        simulatetsir(
          data = subset(dat, time >= time.start & time < epiyr),
          IP = 2,
          parms = fittedpars,
          #epidemics='break', threshold=3,
          #method='pois',
          nsim = 100
        )
    }
    #
    Its = simfitted$res$mean * fittedpars$rho #account for underreporting
    dat.sim <- cbind.data.frame(time = simfitted$res$time, I = Its)
    Ifinal <- Its[length(Its)]
    
    #and sim with and without increased S
    
    
    #now try to predict the epidemic both with and without the increased S
    #compute the most likely model fit with the proper increase which gives the
    #closest match of model to observed cases
    
    
    Sbegin = simfitted$sbar + simfitted$Z
    Sepi <- Sbegin[length(Sbegin)]
    SepiInc <- fracincS * Sepi
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time < epiyr)
    
    predict_ts <- predicttsir(
      times = dat.fit$time,
      births = dat.fit$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sbegin[1],
      I0 = dat.fit$cases[1],
      nsim = 100,
      stochastic = T
    )
    
    
    
    #and the epi year-no increase
    dat.pred = subset(dat, time >= epiyr & time < (epiyr + 1))
    predict_epi_noInc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = Sepi,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    #and including an increase
    predict_epi_Inc <- predicttsir(
      times = dat.pred$time,
      births = dat.pred$births,
      beta = fittedpars$beta,
      alpha = fittedpars$alpha,
      S0 = SepiInc,
      I0 = Ifinal,
      nsim = 100,
      stochastic = T
    )
    
    
    
    
    IPredEpi = melt(
      cbind.data.frame(
        time = predict_epi_Inc$I$time,
        mean_noInc = predict_epi_noInc$I$mean,
        mean_Inc = predict_epi_Inc$I$mean
      ),
      id.vars = "time"
    )
    ISimTS = cbind.data.frame(
      time = predict_ts$I$time,
      variable = rep("TSIR fit", length(predict_ts$I$time)),
      value = predict_ts$I$mean
    )
    
    Iall <- rbind(IPredEpi, ISimTS)
    Iall$variable <- as.character(Iall$variable)
    Iall$variable[Iall$variable == "mean_Inc"] <-
      "prediction TSIR-increased S"
    Iall$variable[Iall$variable == "mean_noInc"] <- "prediction TSIR"
    
    all.dat <- subset(dat, time <= (epiyr + 1))
    all.dat$cases <- mean(fittedpars$rho) * all.dat$cases
    # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) +
    #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    #   geom_line(aes(x=time, y=value, color=variable))
    
    all.dat <- dplyr::select(all.dat, time, cases)
    all.dat$reported_cases <- all.dat$cases / mean(fittedpars$rho)
    all.dat$variable <- "data"
    
    names(Iall)[names(Iall) == "value"] <- "cases"
    Iall$reported_cases <- Iall$cases / mean(fittedpars$rho)
    Iall <- dplyr::select(Iall, names(all.dat))
    
    all.out <- rbind(Iall, all.dat)
    
    
    p2 <- ggplot(data = all.out) +
      geom_line(aes(x = time, y = cases, color = variable))
    
    return(all.out)
    
  }
#here, plot the output
plot.comp <- function(dat, filename) {
  dat$variable = factor(
    dat$variable,
    levels = c(
      "data",
      "TSIR fit",
      "prediction TSIR",
      "prediction TSIR-increased S"
    )
  )
  colz = c(
    'data' = "black",
    'TSIR fit' = "blue",
    'prediction TSIR' = "red",
    'prediction TSIR-increased S' = "seagreen"
  )
  typez = c(
    'data' = 2,
    'TSIR fit' = 1,
    'prediction TSIR' = 1,
    'prediction TSIR-increased S' = 1
  )
  p1 <- ggplot(data = dat) +
    geom_line(aes(
      x = time,
      y = reported_cases,
      color = variable,
      linetype = variable
    ),
    size = 1) +
    scale_color_manual(values = colz) +
    scale_linetype_manual(values = typez) + theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size = 10),
      legend.position = c(.25, .82),
      axis.text = element_text(size = 14),
      legend.title = element_blank()
    ) +
    ylab("reported cases")
  
  print(p1)
  
  ggsave(
    file = filename,
    plot = p1,
    units = "mm",
    width = 50,
    height = 40,
    scale = 3,
    dpi = 300
  )
  
  
}




# gaussian

#######now try as gaussian

#first try to predict 2007-epidemic
sim.2007 <- sim.com.epi(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)


#now return and simulate with results
comp.2007 <- sim.with.increaseS(fracincS = sim.2007$frac_incS,
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
comp.2012 <- sim.with.increaseS(fracincS = sim.2012$frac_incS,
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
comp.2019 <- sim.with.increaseS(fracincS = sim.2019$frac_incS,
                                dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2013]) &time<2020),
                                time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                                family="gaussian",
                                epiyr = 2019)




# prediction TSIR is green
# data is black
# prediction TSIR-increased S is red
# TSIR fit is blue

dat<-rbind(comp.2007,comp.2012,comp.2019)
dat_true<-dat[dat$variable=="data",]

dat_pred_update<-merge(x = dat, y =dat_true, by = "time", all.x = TRUE)

dat_df<-cbind.data.frame(dat_pred_update$time, dat_pred_update$reported_cases.x, dat_pred_update$reported_cases.y, dat_pred_update$variable.x,dat_pred_update$variable.y)

colnames(dat_df)<-c("time","reported_cases","true_cases","pred_type","data")


dat_pred_3<-dat_df
dat_pred_3$year<-floor(dat_pred_3$time)

non_pei_years<- as.data.frame(dplyr::filter(dat_pred_3,  (year>=2002 & year <2007) | (year>=2008 & year <2012) | (year>=2013 & year <2019) | (year>2019) ) )


dat_pred_3[dat_pred_3$time %in% non_pei_years$time,]$year<-"other_year"




dat <- dat_pred_3

colz = c(
  "prediction TSIR" = "seagreen",
  "prediction TSIR-increased S" = "tomato",
  "TSIR fit" = my_blue,
  "data" = "black"
)


circle_size <- 3
pl_d <-
  ggplot(data = dat, aes(x = true_cases, y = reported_cases)) + geom_point(
    aes(
      x = true_cases,
      y = reported_cases,
      color = pred_type,
      shape = year
    ),
    size = circle_size,
    stroke = 2
  ) +
  labs(shape = NULL, colour = NULL) +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # aspect.ratio = 1,
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14)
  ) +theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white"))+ 
  # +geom_point(aes(x=true_cases, y=reported_cases,
  #                fill=year), size=2, stroke = 0.05)
  ylim(0, 9000) + xlim(0, 9000) + scale_color_manual(values = colz) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(.218, .765),
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    axis.text = element_text(size = 14)
  ) + labs(x = "Reported Cases", y = "Predicted Cases")+theme(legend.spacing.y = unit(-0.1, "cm"))# + # space between legend rows
# theme(legend.title.align = 1.08)



######################################################################################
######################################################################################
######################################################################################
############################# subplots together ######################################
######################################################################################
######################################################################################



pl<-plot_grid(pl_a, NULL, pl_b, pl_c, NULL, pl_d, ncol = 3, nrow = 2,align = "vh",labels = c("a","", "b","c","","d"), label_size =xy_size+20, hjust=-0.5, vjust = c(1.3,1.3,1.4,1.4),rel_widths = c(1, 0.08,1))+ theme(plot.margin = unit(c(1,1,2,1), "cm")) 



ggsave(paste0(homewd, "/final-figures/fig1.jpg"),
       plot = pl,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 17,
       height = 12,
       bg='white',
       units = c("in"))









