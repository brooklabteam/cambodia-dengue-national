rm(list = ls())


library(ggplot2)
library(bobfunctions2)
library(plyr)
library(dplyr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/"))

#load the output from the previous trials
load(paste0(homewd,"/figure-development/Fig5/sim-final/cam-sim-final-10-10.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-final/cam-sim-final-10-10-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-final/cam-sim-final-10-10-uci.Rdata"))


#head(cam.sim)



#helper functions
mean.age <- function(df){
  df$mult <- df$age*df$count
  mean.age <- sum(df$mult)/sum(df$count)
  df2 <- cbind.data.frame(hyp=unique(df$hyp), year=unique(df$year), mean_age=mean.age)
  
  return(df2)
}
replicate.data <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    new.dat$hyp <- unique(df$hyp)
    new.dat$round <- unique(df$round)
    return(new.dat)
  }
  
  
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  df.sum$hyp <- unique(df$hyp)
  return(df.sum)
}

#first, make the first colum
#column 1- distribution of serotypes
column.1 <- function(dat, year.start,perc.obs){
  #and get total by time
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  
  denv.case <- rbind(denv.case, denv.case.tert, denv.case.tert.4)
  
  denv.case$serotype <- NA
  denv.case$serotype[denv.case$class=="I12" | denv.case$class=="I32"| denv.case$class=="I132" | denv.case$class=="I312"] <- "I2"
  denv.case$serotype[denv.case$class=="I13" | denv.case$class=="I23" | denv.case$class=="I123" | denv.case$class=="I213"] <- "I3"
  denv.case$serotype[denv.case$class=="I31" | denv.case$class=="I21"| denv.case$class=="I231" | denv.case$class=="I321"] <- "I1"
  
  denv.case$keep <- 1
  denv.case$keep[denv.case$hyp==0 & denv.case$serotype=="I3" | denv.case$hyp==1 & denv.case$serotype=="I3" ] <- 0
  denv.case = subset(denv.case, keep==1)
  
  #dat.ts <- ddply(denv.case, .(time, class), summarise, count=sum(count))
  dat.ser <- ddply(denv.case, .(hyp, time, serotype), summarise, count=sum(count))
  
  
  dat.N  <- ddply(denv.case, .(hyp, time), summarise, Ntot=sum(count))
  
  #dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ser <- merge(dat.ser, dat.N, by=c("hyp", "time"))
  
  #dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  dat.ser$proportion <- dat.ser$count/dat.ser$Ntot
  #dat.ts$proportion[is.na(dat.ts$proportion)] <- 0
  dat.ser$proportion[is.na(dat.ser$proportion)] <- 0
  dat.ser$label <- NA
  dat.ser$label[dat.ser$hyp==0] <- "H0: normal\ndemographic\nsimulation"
  dat.ser$label[dat.ser$hyp==1] <- "H1:\nhigh\n2019 FOI"
  dat.ser$label[dat.ser$hyp==2] <- "H2: 2019 strain\nintro + replacement\n+immune waning"
  dat.ser$label[dat.ser$hyp==3] <- "H3: 2019 serotype\nintro + maintenance +\ntertiary case detection"
  dat.ser$label[dat.ser$hyp==4] <- "H4: 3 endemic serotypes\n+increasing tertiary case\ndetection with time"
  
  
  dat.ser$serotype_label <- NA
  dat.ser$serotype_label[dat.ser$serotype=="I1"] <- "serotype-1/strain-type 1"
  dat.ser$serotype_label[dat.ser$serotype=="I2"] <- "serotype-2"
  dat.ser$serotype_label[dat.ser$serotype=="I3" & dat.ser$hyp==2] <- "serotype-1/strain-type 2"
  dat.ser$serotype_label[dat.ser$serotype=="I3" & dat.ser$hyp==3 |dat.ser$serotype=="I3" & dat.ser$hyp==4] <- "serotype-3"
  dat.ser$label <- factor(dat.ser$label, levels = c(unique(dat.ser$label)))
  
  #and plot
  dat.ser$plot_type <- "proportion of cases\nby serotype/strain-type"
  
  p1 <- ggplot(dat.ser) + theme_bw() +
    geom_line(aes(x=time, y=proportion, color=serotype_label),  size=1, position = position_jitter(height = .03)) +
    facet_grid(label~plot_type, switch = "y") + 
    theme(panel.grid = element_blank(), axis.title = element_blank(), 
          legend.position = c(.35,.95), legend.title = element_blank(),
          legend.key.height =  unit(c(.2), "cm"), strip.placement = "outside",#axis.title.y = element_text(size=12), 
          axis.text = element_text(size=12),
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          strip.background = element_rect(fill="white"), strip.text =element_text(size=12)) +
    #coord_cartesian(ylim=c(0,1), xlim=c(2015,2020)) + scale_y_continuous(breaks = c(0,.5,1)) +
    #scale_x_continuous(breaks=c(2016,2018, 2020)) 
    coord_cartesian(ylim=c(0,1), xlim=c(year.start,2020)) + scale_y_continuous(breaks = c(0,.5,1)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) 
  
  
  return(p1)
}
#column 2 - total cases
column.2 <- function(dat, dat.lci, dat.uci, year.start, perc.obs){
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  
  denv.case <- rbind(denv.case, denv.case.tert, denv.case.tert.4)
  
  
  
  #and lci
  dat1.lci = subset(dat.lci, year >= year.start) 
  denv.case.lci = subset(dat1.lci, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert.lci = subset(dat1.lci, hyp==2|hyp==3)
  denv.case.tert.4.lci =subset(dat1.lci,hyp==4)
  
  denv.case.tert.lci = subset(denv.case.tert.lci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert.lci$count <- denv.case.tert.lci$count*perc.obs
  
  denv.case.tert.4.lci = subset(denv.case.tert.4.lci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc.lci <- cbind.data.frame(year=unique(denv.case.tert.4.lci$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4.lci <- merge(denv.case.tert.4.lci, df.perc.lci, by="year")
  denv.case.tert.4.lci$count <- denv.case.tert.4.lci$count*denv.case.tert.4.lci$perc_obs
  denv.case.tert.4.lci <- dplyr::select(denv.case.tert.4.lci, -(perc_obs))
  
  denv.case.lci <- rbind(denv.case.lci, denv.case.tert.lci, denv.case.tert.4.lci)
  
  
  
  
  #and uci
  dat1.uci = subset(dat.uci, year >= year.start) 
  denv.case.uci = subset(dat1.uci, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert.uci = subset(dat1.uci, hyp==2|hyp==3)
  denv.case.tert.4.uci =subset(dat1.uci,hyp==4)
  
  denv.case.tert.uci = subset(denv.case.tert.uci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert.uci$count <- denv.case.tert.uci$count*perc.obs
  
  denv.case.tert.4.uci = subset(denv.case.tert.4.uci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc.uci <- cbind.data.frame(year=unique(denv.case.tert.4.uci$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4.uci <- merge(denv.case.tert.4.uci, df.perc.uci, by="year")
  denv.case.tert.4.uci$count <- denv.case.tert.4.uci$count*denv.case.tert.4.uci$perc_obs
  denv.case.tert.4.uci <- dplyr::select(denv.case.tert.4.uci, -(perc_obs))
  
  denv.case.uci <- rbind(denv.case.uci, denv.case.tert.uci, denv.case.tert.4.uci)
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year), summarise, count=sum(count))
  dat.ts.lci <- ddply(denv.case.lci, .(hyp, year), summarise, count=sum(count))
  dat.ts.uci <- ddply(denv.case.uci, .(hyp, year), summarise, count=sum(count))
  dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total reported cases\nfrom deterministic simulation"
  dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  
  dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  names(dat.ts.lci.merge)[3] <- "lci"
  names(dat.ts.uci.merge)[3] <- "uci"
  
  dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  dat.ts$lci_new <- NA
  dat.ts$uci_new <- NA
  dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  
  
  
  p1 <- ggplot(dat.ts) + theme_bw() + facet_grid(hyp~plot_type) +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    geom_line(aes(x=year, y=count), size=.3) + 
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.text.x.top  = element_text(size=12),
          strip.text.y = element_blank()) +
    #      coord_cartesian(xlim=c(2015,2020), ylim=c(0,13000)) + #scale_y_log10() +
    scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) 
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
          
  return(p1)
}
#column 3 age distribution of cases
column.3 <- function (dat, year.start,perc.obs){
  
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  
  
  denv.case$round <- "secondary"
  denv.case.tert$round <- "tertiary"
  denv.case.tert.4$round <- "tertiary"
  
  
  denv.case <- rbind(denv.case, denv.case.tert, denv.case.tert.4)
  
  denv.case$round <- factor(denv.case$round, levels=c("secondary", "tertiary"))
  #subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(hyp,round, year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  #split by a year
  df.year <- dlply(df.sum,.(hyp, year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year
  df.age <- dlply(df.sum,.(hyp,round, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data, slim.quant=1))
  
  dat.age$plot_type <- "age distribution of\nreported cases by year"
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  colz = c('secondary'="black", "tertiary"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + facet_grid(hyp~plot_type, scales = "free_y") + theme_bw()+
    geom_jitter(aes(x=year, y=age, color=round), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    scale_color_manual(values=colz) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato", size=.8) +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.x.top = element_text(size=12),
          strip.text.y =element_blank()) + 
    scale_x_continuous(breaks=c(c(2000,2005,2010,2015,2020)))
    
  
  return(p1)
  
   }
#column 4 cumulative proportion of cases
column.4 <- function(dat, year.start,perc.obs){
  
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  
  denv.case <- rbind(denv.case, denv.case.tert, denv.case.tert.4)
  
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(hyp, year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  denv.dat$plot_type <- "cumulative proportion\nof cases by age"
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + facet_grid(hyp~plot_type) + theme_bw()+
        geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
        scale_color_viridis_d(direction=-1,option="turbo") +
        theme(panel.grid = element_blank(), 
              strip.background.y  = element_blank(),
              strip.background.x  = element_rect(fill="white"),
          legend.title = element_blank(),
          axis.title = element_blank(), 
          legend.key.size = unit(c(.4), "cm"), 
          legend.position = c(.6,.91),
          legend.text = element_text(size=10),
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.y = element_blank(),
          strip.text.x.top = element_text(size=12)) + 
        guides(color=guide_legend(ncol=3)) 
  
  return(p1)
  
}

unique(cam.sim$hyp)


cam.2019 = subset(cam.sim, hyp!="H3-2007" & hyp!="H1-2007" & hyp!="H2-2007")
cam.lci.2019 = subset(cam.sim.lci, hyp!="H3-2007" & hyp!="H1-2007"& hyp!="H2-2007")
cam.uci.2019 = subset(cam.sim.uci, hyp!="H3-2007" & hyp!="H1-2007"& hyp!="H2-2007")

cam.2007 = subset(cam.sim,hyp=="H1-2007" | hyp=="H2-2007"|  hyp=="H3-2007")
cam.lci.2007 = subset(cam.sim.lci, hyp=="H1-2007" | hyp=="H2-2007"|  hyp=="H3-2007")
cam.uci.2007 = subset(cam.sim.uci, hyp=="H1-2007" | hyp=="H2-2007"|  hyp=="H3-2007")

rm(cam.sim, cam.sim.lci, cam.sim.uci)
cam.2019$hyp <- as.numeric(cam.2019$hyp)
cam.lci.2019$hyp <- as.numeric(cam.lci.2019$hyp)
cam.uci.2019$hyp <- as.numeric(cam.uci.2019$hyp)

cam.2007$hyp[cam.2007$hyp=="H1-2007"] <- "1"
cam.2007$hyp[cam.2007$hyp=="H2-2007"] <- "2"
cam.2007$hyp[cam.2007$hyp=="H3-2007"] <- "3"
cam.2007$hyp<- as.numeric(cam.2007$hyp)

cam.lci.2007$hyp[cam.lci.2007$hyp=="H1-2007"] <- "1"
cam.lci.2007$hyp[cam.lci.2007$hyp=="H2-2007"] <- "2"
cam.lci.2007$hyp[cam.lci.2007$hyp=="H3-2007"] <- "3"
cam.lci.2007$hyp<- as.numeric(cam.lci.2007$hyp)

cam.uci.2007$hyp[cam.uci.2007$hyp=="H1-2007"] <- "1"
cam.uci.2007$hyp[cam.uci.2007$hyp=="H2-2007"] <- "2"
cam.uci.2007$hyp[cam.uci.2007$hyp=="H3-2007"] <- "3"
cam.uci.2007$hyp<- as.numeric(cam.uci.2007$hyp)

col1 <- column.1(dat=cam.2019, year.start = 2000, perc.obs=.1)
col2 <- column.2(dat=cam.2019, dat.lci=cam.lci.2019, dat.uci=cam.uci.2019, year.start = 2000, perc.obs = .1)
col3 <- column.3(dat=cam.2019, year.start = 2000, perc.obs = .1)
col4 <- column.4(dat=cam.2019, year.start = 2000, perc.obs = .1)


#and put them together

Fig5 <- cowplot::plot_grid(col1, col2, col3, col4, nrow = 1, ncol=4, labels=c("A", "B", "C", "D"),rel_widths = c(1.15,1.05,1,1),  label_size = 22, align = "hv", label_x = c(0,0,-.03,0))

ggsave(filename = paste0(homewd, "/final-figures/Fig5.png"),
       plot = Fig5,
       units="mm",  
       width=125, 
       height=90, 
       scale=3, 
       dpi=300)

ggsave(filename = paste0(homewd, "/final-figures/Fig5.pdf"),
       plot = Fig5,
       units="mm",  
       width=125, 
       height=90, 
       scale=3, 
       dpi=300)

#and the supplementary figure (S21)
##this is for the 2007 simulations - same plot as above
#just redo the first column to edit the hypothesis names
#column 1- distribution of serotypes
column.1.2007 <- function(dat, year.start,perc.obs){
  #and get total by time
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  
  denv.case <- rbind(denv.case, denv.case.tert, denv.case.tert.4)
  
  denv.case$serotype <- NA
  denv.case$serotype[denv.case$class=="I12" | denv.case$class=="I32"| denv.case$class=="I132" | denv.case$class=="I312"] <- "I2"
  denv.case$serotype[denv.case$class=="I13" | denv.case$class=="I23" | denv.case$class=="I123" | denv.case$class=="I213"] <- "I3"
  denv.case$serotype[denv.case$class=="I31" | denv.case$class=="I21"| denv.case$class=="I231" | denv.case$class=="I321"] <- "I1"
  
  denv.case$keep <- 1
  denv.case$keep[denv.case$hyp==0 & denv.case$serotype=="I3" | denv.case$hyp==1 & denv.case$serotype=="I3" ] <- 0
  denv.case = subset(denv.case, keep==1)
  
  #dat.ts <- ddply(denv.case, .(time, class), summarise, count=sum(count))
  dat.ser <- ddply(denv.case, .(hyp, time, serotype), summarise, count=sum(count))
  
  
  dat.N  <- ddply(denv.case, .(hyp, time), summarise, Ntot=sum(count))
  
  #dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ser <- merge(dat.ser, dat.N, by=c("hyp", "time"))
  
  #dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  dat.ser$proportion <- dat.ser$count/dat.ser$Ntot
  #dat.ts$proportion[is.na(dat.ts$proportion)] <- 0
  dat.ser$proportion[is.na(dat.ser$proportion)] <- 0
  dat.ser$label <- NA
  #dat.ser$label[dat.ser$hyp==0] <- "H0: normal\ndemographic\nsimulation"
  dat.ser$label[dat.ser$hyp==1] <- "H1: high\n2007 FOI"
  dat.ser$label[dat.ser$hyp==2] <- "H2: 2007 strain\nintro + replacement\n+immune waning"
  dat.ser$label[dat.ser$hyp==3] <- "H3: 2007 serotype\nintro + maintenance +\ntertiary case detection"
  #dat.ser$label[dat.ser$hyp==4] <- "H4: 3 endemic serotypes\n+increasing tertiary case\ndetection with time"
  
  
  dat.ser$serotype_label <- NA
  dat.ser$serotype_label[dat.ser$serotype=="I1"] <- "serotype-1/strain-type 1"
  dat.ser$serotype_label[dat.ser$serotype=="I2"] <- "serotype-2"
  dat.ser$serotype_label[dat.ser$serotype=="I3" & dat.ser$hyp==2] <- "serotype-1/strain-type 2"
  dat.ser$serotype_label[dat.ser$serotype=="I3" & dat.ser$hyp==3 |dat.ser$serotype=="I3" & dat.ser$hyp==4] <- "serotype-3"
  dat.ser$label <- factor(dat.ser$label, levels = c(unique(dat.ser$label)))
  
  #and plot
  dat.ser$plot_type <- "proportion of cases\nby serotype/strain-type"
  
  p1 <- ggplot(dat.ser) + theme_bw() +
    geom_line(aes(x=time, y=proportion, color=serotype_label),  size=1, position = position_jitter(height = .03)) +
    facet_grid(label~plot_type, switch = "y") + 
    theme(panel.grid = element_blank(), axis.title = element_blank(), 
          legend.position = c(.35,.93), legend.title = element_blank(),
          legend.key.height =  unit(c(.2), "cm"), strip.placement = "outside",#axis.title.y = element_text(size=12), 
          axis.text = element_text(size=12),
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          strip.background = element_rect(fill="white"), strip.text =element_text(size=12)) +
    #coord_cartesian(ylim=c(0,1), xlim=c(2015,2020)) + scale_y_continuous(breaks = c(0,.5,1)) +
    #scale_x_continuous(breaks=c(2016,2018, 2020)) 
    coord_cartesian(ylim=c(0,1), xlim=c(year.start,2020)) + scale_y_continuous(breaks = c(0,.5,1)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) 
  
  
  return(p1)
}
column.2.2007 <- function(dat, dat.lci, dat.uci, year.start, perc.obs){
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  #denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  # denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  # #and take increasing proportion of those through time
  # df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  # denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  # denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  # denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  # 
  denv.case <- rbind(denv.case, denv.case.tert)#, denv.case.tert.4)
  
  
  
  #and lci
  dat1.lci = subset(dat.lci, year >= year.start) 
  denv.case.lci = subset(dat1.lci, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert.lci = subset(dat1.lci, hyp==2|hyp==3)
  #denv.case.tert.4.lci =subset(dat1.lci,hyp==4)
  
  denv.case.tert.lci = subset(denv.case.tert.lci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert.lci$count <- denv.case.tert.lci$count*perc.obs
  
  #denv.case.tert.4.lci = subset(denv.case.tert.4.lci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  #df.perc.lci <- cbind.data.frame(year=unique(denv.case.tert.4.lci$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  #denv.case.tert.4.lci <- merge(denv.case.tert.4.lci, df.perc.lci, by="year")
 # denv.case.tert.4.lci$count <- denv.case.tert.4.lci$count*denv.case.tert.4.lci$perc_obs
 # denv.case.tert.4.lci <- dplyr::select(denv.case.tert.4.lci, -(perc_obs))
  
  denv.case.lci <- rbind(denv.case.lci, denv.case.tert.lci)#, denv.case.tert.4.lci)
  
  
  
  
  #and uci
  dat1.uci = subset(dat.uci, year >= year.start) 
  denv.case.uci = subset(dat1.uci, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert.uci = subset(dat1.uci, hyp==2|hyp==3)
  #denv.case.tert.4.uci =subset(dat1.uci,hyp==4)
  
  denv.case.tert.uci = subset(denv.case.tert.uci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert.uci$count <- denv.case.tert.uci$count*perc.obs
  
  #denv.case.tert.4.uci = subset(denv.case.tert.4.uci, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take increasing proportion of those through time
  # df.perc.uci <- cbind.data.frame(year=unique(denv.case.tert.4.uci$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  # denv.case.tert.4.uci <- merge(denv.case.tert.4.uci, df.perc.uci, by="year")
  # denv.case.tert.4.uci$count <- denv.case.tert.4.uci$count*denv.case.tert.4.uci$perc_obs
  # denv.case.tert.4.uci <- dplyr::select(denv.case.tert.4.uci, -(perc_obs))
  # 
  denv.case.uci <- rbind(denv.case.uci, denv.case.tert.uci)#, denv.case.tert.4.uci)
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year), summarise, count=sum(count))
  dat.ts.lci <- ddply(denv.case.lci, .(hyp, year), summarise, count=sum(count))
  dat.ts.uci <- ddply(denv.case.uci, .(hyp, year), summarise, count=sum(count))
  dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total reported cases\nfrom deterministic simulation"
  dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  
  dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  names(dat.ts.lci.merge)[3] <- "lci"
  names(dat.ts.uci.merge)[3] <- "uci"
  
  dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  dat.ts$lci_new <- NA
  dat.ts$uci_new <- NA
  dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  
  
  
  p1 <- ggplot(dat.ts) + theme_bw() + facet_grid(hyp~plot_type) +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    geom_line(aes(x=year, y=count), size=.3) + 
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.text.x.top  = element_text(size=12),
          strip.text.y = element_blank()) +
    #      coord_cartesian(xlim=c(2015,2020), ylim=c(0,13000)) + #scale_y_log10() +
    scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) 
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
  
  return(p1)
}
column.3.2007 <- function (dat, year.start,perc.obs){
  
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  #denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  # denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  # #and take increasing proportion of those through time
  # df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  # denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  # denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  # denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  # 
  # 
  denv.case$round <- "secondary"
  denv.case.tert$round <- "tertiary"
  #denv.case.tert.4$round <- "tertiary"
  
  
  denv.case <- rbind(denv.case, denv.case.tert)#, denv.case.tert.4)
  
  denv.case$round <- factor(denv.case$round, levels=c("secondary", "tertiary"))
  #subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(hyp,round, year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  #split by a year
  df.year <- dlply(df.sum,.(hyp, year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year
  df.age <- dlply(df.sum,.(hyp,round, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data, slim.quant=1))
  
  dat.age$plot_type <- "age distribution of\nreported cases by year"
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  colz = c('secondary'="black", "tertiary"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + facet_grid(hyp~plot_type, scales = "free_y") + theme_bw()+
    geom_jitter(aes(x=year, y=age, color=round), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    scale_color_manual(values=colz) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato", size=.8) +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.x.top = element_text(size=12),
          strip.text.y =element_blank()) + 
    scale_x_continuous(breaks=c(c(2000,2005,2010,2015,2020)))
  
  
  return(p1)
  
}
column.4.2007 <- function(dat, year.start,perc.obs){
  
  dat1 = subset(dat, year >= year.start) 
  denv.case = subset(dat1, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  denv.case.tert = subset(dat1, hyp==2|hyp==3)
  #denv.case.tert.4 =subset(dat1,hyp==4)
  
  denv.case.tert = subset(denv.case.tert, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  #and take 10% of these
  denv.case.tert$count <- denv.case.tert$count*perc.obs
  
  # denv.case.tert.4 = subset(denv.case.tert.4, class=="I123" | class=="I132" | class=="I213" | class == "I231" | class=="I312" | class =="I321")
  # #and take increasing proportion of those through time
  # df.perc <- cbind.data.frame(year=unique(denv.case.tert.4$year), perc_obs=seq(0,1, length.out = length(unique(denv.case.tert.4$year))))
  # denv.case.tert.4 <- merge(denv.case.tert.4, df.perc, by="year")
  # denv.case.tert.4$count <- denv.case.tert.4$count*denv.case.tert.4$perc_obs
  # denv.case.tert.4 <- dplyr::select(denv.case.tert.4, -(perc_obs))
  # 
  denv.case <- rbind(denv.case, denv.case.tert)#, denv.case.tert.4)
  
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(hyp, year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  denv.dat$plot_type <- "cumulative proportion\nof cases by age"
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + facet_grid(hyp~plot_type) + theme_bw()+
    geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
    scale_color_viridis_d(direction=-1,option="turbo") +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          legend.title = element_blank(),
          axis.title = element_blank(), 
          legend.key.size = unit(c(.4), "cm"), 
          legend.position = c(.6,.82),
          legend.text = element_text(size=10),
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.y = element_blank(),
          strip.text.x.top = element_text(size=12)) + 
    guides(color=guide_legend(ncol=3)) 
  
  return(p1)
  
}


col1Supp <- column.1.2007(dat=cam.2007, year.start = 2000, perc.obs=.1)
col2Supp <- column.2.2007(dat=cam.2007, dat.lci=cam.lci.2007, dat.uci=cam.uci.2007, year.start = 2000, perc.obs = .1)
col3Supp <- column.3.2007(dat=cam.2007, year.start = 2000, perc.obs = .1)
col4Supp <- column.4.2007(dat=cam.2007, year.start = 2000, perc.obs = .1)




FigS21 <- cowplot::plot_grid(col1Supp, col2Supp, col3Supp, col4Supp, nrow = 1, ncol=4, labels=c("A", "B", "C", "D"),rel_widths = c(1.15,1.05,1,1),  label_size = 22, align = "hv", label_x = c(0,0,-.03,0))

ggsave(filename = paste0(homewd, "/final-figures/FigS21.png"),
       plot = FigS21,
       units="mm",  
       width=125, 
       height=50, 
       scale=3, 
       dpi=300)

  