rm(list = ls())


library(ggplot2)
#library(bobfunctions2)
library(plyr)
library(dplyr)
library(stringr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/"))

#load the output from the previous trials
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane-uci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007-uci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019-uci.Rdata"))


summarise.age.dist.wane <- function(dat, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #summarise into primary, secondary, tertiary, quaternary cases
  
  dat1$case_type <- (str_extract(dat1$class, "[aA-zZ]+"))
  dat1$case_sum <- nchar(dat1$class)
  
  dat1$serotype <- NA  
  dat1$serotype[dat1$case_type=="I"] <- str_sub(dat1$class[dat1$case_type=="I"], start=-1)
  
  dat1$state <- NA  
  dat1$state[dat1$case_type=="S"] <- "S"
  dat1$state[dat1$case_type=="P" | dat1$case_type=="Pms"] <- "Temporary-Heterotypic-Immunity"
  dat1$state[dat1$case_type=="Pm" ] <- "Pm"
  
  dat1$state[dat1$case_type=="I" & dat1$case_sum==2] <- "Primary-Infection"
  
  # you want to track how often you get "reinfections" within a serotype - 
  # these have to be serotype 4. and then, they need to have previously also experienced
  # infection with serotype 1
  
  dat1$state[dat1$class=="I14"] <- "Primary-Re-Infection"
  dat1$state[dat1$serotype=="4" & dat1$case_sum==4 & dat1$case_type=="I" & dat1$class!="I234" & dat1$class!="I324"] <- "Secondary-Re-Infection"
  dat1$state[dat1$serotype=="4" & dat1$case_sum==5 & dat1$case_type=="I"] <- "Tertiary-Re-Infection"
  
  #and the rest of the cases
  dat1$state[dat1$case_type=="I" & dat1$case_sum==3 & is.na(dat1$state)] <- "Secondary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==4 & is.na(dat1$state)] <- "Tertiary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==5 & is.na(dat1$state) & dat1$serotype!=1] <- "Tertiary-Infection-After-Reinfection"
  
  #and there are a few that still need editing - can't have 1 after 4
  dat1$state[dat1$case_type=="I" & dat1$case_sum==5 & dat1$serotype==1] <- "Not-Possible"
  dat1$state[dat1$class=="I241" | dat1$class=="I341"| dat1$class=="I412"|dat1$class=="I413"|dat1$class=="I421" | dat1$class=="I431"] <- "Not-Possible"
  dat1$state[dat1$class=="I142" | dat1$class=="I143"] <- "Secondary-Infection-After-Reinfection"
  
  
  dat2 = subset(dat1, state!="Not-Possible")
  
  #tracking reinfections with serotype 2
  
  #and sum by year
  df.sum <- ddply(dat2,.(year, age, state), summarise, count = sum(count))
  
  
  
  df.sum$count<- round(df.sum$count,0)
  
  
  df.sum <- df.sum[complete.cases(df.sum),]
  
  #and just focus on infections
  #df.sum.1 <- df.sum
  
  #and return these
  df.sum.I = subset(df.sum, state!="Pm" & state!="Temporary-Heterotypic-Immunity" & state!="S")
  return(df.sum.I)
}

age.out.2007 = summarise.age.dist.wane(dat=out.cam.geno.rep.2007, year.start = min(out.cam.geno.rep.2007$year))
age.out.2007.lci = summarise.age.dist.wane(dat=out.cam.geno.rep.2007.lci, year.start = min(out.cam.geno.rep.2007$year))
age.out.2007.uci = summarise.age.dist.wane(dat=out.cam.geno.rep.2007.uci, year.start = min(out.cam.geno.rep.2007$year))

age.out.2019 = summarise.age.dist.wane(dat=out.cam.geno.rep.2019, year.start = min(out.cam.geno.rep.2019$year))
age.out.2019.lci = summarise.age.dist.wane(dat=out.cam.geno.rep.2019.lci, year.start = min(out.cam.geno.rep.2019$year))
age.out.2019.uci = summarise.age.dist.wane(dat=out.cam.geno.rep.2019.uci, year.start = min(out.cam.geno.rep.2019$year))

age.out.nointro = summarise.age.dist.wane(dat=out.cam.geno, year.start = min(out.cam.geno$year))
age.out.nointro.lci = summarise.age.dist.wane(dat=out.cam.geno.lci, year.start = min(out.cam.geno$year))
age.out.nointro.uci = summarise.age.dist.wane(dat=out.cam.geno.uci, year.start = min(out.cam.geno$year))


#and select what gets counted as symptomatic
select.symptom <- function(df, criteria){
  
   
  
  
  #first, select what you want.
  if(criteria=="Secondary+Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria=="Secondary+Reinfection+ReSecondary"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection" | state == "Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria == "Secondary-Only"){
    df1 = subset(df, state =="Secondary-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if (criteria=="Secondary+Secondary-Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
    
  }else if (criteria=="Secondary-Extension"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection" | state=="Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection" & df1$state!="Secondary-Infection-After-Reinfection"] <- "Repeat-Infection"
    
  }else if (criteria=="Increasing-Tertiary"){
      df1.sec = subset(df, state =="Secondary-Infection" )
      df1.tert = subset(df, state=="Tertiary-Infection" )
      df.perc <- cbind.data.frame(year=unique(df1.tert$year), perc_obs=seq(0,.3, length.out = length(unique(df1.tert$year))))
      df1.tert <- merge(df1.tert, df.perc, by="year", all.x=T)
      df1.tert$count <- df1.tert$count*df1.tert$perc_obs
      df1.tert <- dplyr::select(df1.tert, -(perc_obs))
      
      df1 <- rbind(df1.sec, df1.tert)
      
    
    
  }
  
  
  
  df1$case_type = "symptomatic"
  
  return(df1)
  
  
}

age.out.2007 = subset(age.out.2007, year <2021)
age.sub.2007 = select.symptom(df=age.out.2007,criteria = "Secondary-Extension")
age.sub.2007$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2007.lci = subset(age.out.2007.lci, year <2021)
age.sub.2007.lci = select.symptom(df=age.out.2007.lci,criteria = "Secondary-Extension")
age.sub.2007.lci$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2007.uci = subset(age.out.2007.uci, year <2021)
age.sub.2007.uci = select.symptom(df=age.out.2007.uci,criteria = "Secondary-Extension")
age.sub.2007.uci$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2019 = subset(age.out.2019, year <2021)
age.sub.2019 = select.symptom(df=age.out.2019,criteria = "Secondary-Extension")
age.sub.2019$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2019)"

age.out.2019.lci = subset(age.out.2019.lci, year <2021)
age.sub.2019.lci = select.symptom(df=age.out.2019.lci,criteria = "Secondary-Extension")
age.sub.2019.lci$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2019)"

age.out.2019.uci = subset(age.out.2019.uci, year <2021)
age.sub.2019.uci = select.symptom(df=age.out.2019.uci,criteria = "Secondary-Extension")
age.sub.2019.uci$hyp = "H2: Genotype Replacement\n+ Waning Immunity (2019)"

age.out.nointro = subset(age.out.nointro, year<2021)
age.sub.tert = select.symptom(df=age.out.nointro,criteria = "Increasing-Tertiary")
age.sub.tert$hyp = "H1: Increasing Tertiary\nCase Detection"

age.out.nointro.lci = subset(age.out.nointro.lci, year<2021)
age.sub.tert.lci = select.symptom(df=age.out.nointro.lci,criteria = "Increasing-Tertiary")
age.sub.tert.lci$hyp = "H1: Increasing Tertiary\nCase Detection"

age.out.nointro.uci = subset(age.out.nointro.uci, year<2021)
age.sub.tert.uci = select.symptom(df=age.out.nointro.uci,criteria = "Increasing-Tertiary")
age.sub.tert.uci$hyp = "H1: Increasing Tertiary\nCase Detection"



age.sub.H0 = select.symptom(df=age.out.nointro,criteria = "Secondary-Only")
age.sub.H0$hyp = "H0: Normal Demographic\nSimulation"



age.sub.H0.lci = select.symptom(df=age.out.nointro.lci,criteria = "Secondary-Only")
age.sub.H0.lci$hyp = "H0: Normal Demographic\nSimulation"



age.sub.H0.uci = select.symptom(df=age.out.nointro.uci,criteria = "Secondary-Only")
age.sub.H0.uci$hyp = "H0: Normal Demographic\nSimulation"


#put all the data together
comp.dat <- rbind(age.sub.H0, age.sub.tert, age.sub.2019, age.sub.2007)
comp.dat$hyp <- factor(comp.dat$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))

#and save for fitting
save(comp.dat, file = "comp-dat-sim.Rdata") 

comp.dat.lci <- rbind(age.sub.H0.lci, age.sub.tert.lci, age.sub.2019.lci, age.sub.2007.lci)
comp.dat.lci$hyp <- factor(comp.dat.lci$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
save(comp.dat.lci, file = "comp-dat-sim-lci.Rdata") 
comp.dat.uci <- rbind(age.sub.H0.uci, age.sub.tert.uci, age.sub.2019.uci, age.sub.2007.uci)
comp.dat.uci$hyp <- factor(comp.dat.uci$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
save(comp.dat.uci, file = "comp-dat-sim-lci.Rdata") 


#then, feed into plotting

#helper functions
mean.age <- function(df){
  df$mult <- df$age*df$count
  mean.age <- sum(df$mult)/sum(df$count)
  df2 <- cbind.data.frame(year=unique(df$year), mean_age=mean.age, hyp = unique(df$hyp))
  return(df2)
}
replicate.data.type <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    new.dat$state <- unique(df$state)
    new.dat$hyp <- unique(df$hyp)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    
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

column.1 <- function(dat, dat.lci, dat.uci, year.start){
  dat1 = subset(dat, year >= year.start) 
  dat1.lci = subset(dat.lci, year >= year.start) 
  dat1.uci = subset(dat.uci, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  denv.case.lci = subset(dat1.lci, case_type == "symptomatic")
  denv.case.uci = subset(dat1.uci, case_type == "symptomatic")
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year), summarise, count=sum(count))
  dat.ts.lci <- ddply(denv.case.lci, .(hyp,year), summarise, count=sum(count))
  dat.ts.uci <- ddply(denv.case.uci, .(hyp, year), summarise, count=sum(count))
  #dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # 
  dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  names(dat.ts.lci.merge)[3] <- "lci"
  names(dat.ts.uci.merge)[3] <- "uci"
  # 
  dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  dat.ts$lci_new <- NA
  dat.ts$uci_new <- NA
  dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  # 
  
  #dat.ts$count_dss <- dat.ts$count*perc_dss
  #dat.ts$count_mort <- dat.ts$count_dss*perc_mort #some fraction of dss cases is the total mortality
  
  p1 <- ggplot(dat.ts) + theme_bw() + facet_grid(hyp~plot_type, switch = "y") +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    geom_line(aes(x=year, y=count), size=.3) + 
    #geom_line(aes(x=year, y=count_dss), size=.3, color="navy") + 
    #geom_line(aes(x=year, y=count_mort), size=.3, color="green") + 
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_rect(fill="white"),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.placement = "outside",
          strip.text.x.top  = element_text(size=12),
          strip.text.y.left  = element_text(size=12)) +
    #      coord_cartesian(xlim=c(2015,2020), ylim=c(0,13000)) + #scale_y_log10() +
    #scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) +
    geom_vline(aes(xintercept=2007), linetype=2, color="red") +
    geom_vline(aes(xintercept=2012), linetype=2, color="red") +
    geom_vline(aes(xintercept=2019), linetype=2, color="red")
  
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
          
  return(p1)
}
col1 <- column.1(dat =comp.dat, dat.lci = comp.dat.lci, dat.uci = comp.dat.uci, year.start = 2002)
#column 2 age distribution of cases
column.2 <- function (dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  #denv.case$round <- factor(denv.case$round, levels=c("secondary", "tertiary"))
  #subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
   if(length(unique(denv.case$age))>1){
     denv.case = subset(denv.case, age<max(denv.case$age))  
   }
  # #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
    df.sum = ddply(denv.case,.(hyp, year,age, state),summarise, count=sum(count))  
  
  
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  
  df.year.only <- dlply(df.sum,.(hyp, year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year.only, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and state
  df.age.list <- dlply(df.sum,.(hyp, year, age, state))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age.list, replicate.data.type, slim.quant=1))
  dat.age$state[dat.age$state!="Secondary-Infection"] <- "alternative"
  dat.age$state[dat.age$state=="Secondary-Infection"] <- "secondary"
  dat.age$plot_type <- "age distribution of\nreported cases by year"
  #dat.age$hyp = unique(dat.age$hyp)
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  
  colz = c('secondary'="black", "alternative"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + facet_grid(hyp~plot_type) + theme_bw()+
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
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
col2 <- column.2(dat =comp.dat, year.start = 2002)
#column 3 cumulative proportion of cases
column.3 <- function(dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  
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
col3 <- column.3(dat =comp.dat,  year.start = 2002)


#and put them together

Fig5 <- cowplot::plot_grid(col1, col2, col3, nrow = 1, ncol=3, labels=c("A", "B", "C"),rel_widths = c(1.15,1,1),  label_size = 22, align = "hv", label_x = c(0,-.03,0))

ggsave(filename = paste0(homewd, "/final-figures/Fig5.png"),
       plot = Fig5,
       units="mm",  
       width=110, 
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



plot.age.dist(df = age.out, criteria = "Secondary+Reinfection", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out, criteria = "Secondary+Reinfection+ReSecondary", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out, criteria = "Secondary+Secondary-Reinfection", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out, criteria = "Secondary-Only", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)

plot.age.dist(df = age.out.2019, criteria = "Secondary+Reinfection", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out.2019, criteria = "Secondary+Reinfection+ReSecondary", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out.2019, criteria = "Secondary+Secondary-Reinfection", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.age.dist(df = age.out.2019, criteria = "Secondary-Only", slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)


#Secondary + Secondary Reinfection the best

#and what if mortality is less in repeats within a serotype?
plot.mort.age <- function(df, criteria, cfr_secondary, cfr_repeat, slim.quant, view.plot, save.plot, filename){
  
  #first, select what you want.
  if(criteria=="Secondary+Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria=="Secondary+Reinfection+ReSecondary"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection" | state == "Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria == "Secondary-Only"){
    df1 = subset(df, state =="Secondary-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if (criteria=="Secondary+Secondary-Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
    
  }
  
  
  #then, sum by year
  df1 <- ddply(df1,.(year, age, state), summarise, count = sum(count))
  df1$cfr <- NA
  df1$cfr[df1$state=="Secondary-Infection"] <- cfr_secondary
  df1$cfr[df1$state=="Repeat-Infection"] <- cfr_repeat
  
  #and just get annual total mortality by age
  
  df1$mort <- ceiling(df1$count*df1$cfr)
  
  df.sum <- ddply(df1,.(year, age), summarise, mort = sum(mort))
  df.sum$state <- "deaths"
  
  #split by a year
  df.year <- dlply(df.sum,.(year, state))
  
  # #get mean age
  # mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  # mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  # 
  
  #split by age and year and type
  df.age <- dlply(df.sum,.(state, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type.mort, slim.quant=slim.quant))
  
  #dat.age$state <- factor(dat.age$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #mean.df$state <- factor(mean.df$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  #colz = c("secondary"="black", "tertiary"="royalblue3")
  
  p1 <- ggplot(dat.age) + facet_grid(~state) +
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    #geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + #scale_color_manual(values=colz) +
    coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=100, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  #return(mean.df)
}
plot.mort.total <- function(df, criteria, cfr_secondary, cfr_repeat, slim.quant, view.plot, save.plot, filename){
  
  #first, select what you want.
  if(criteria=="Secondary+Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria=="Secondary+Reinfection+ReSecondary"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection" | state == "Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria == "Secondary-Only"){
    df1 = subset(df, state =="Secondary-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if (criteria=="Secondary+Secondary-Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
    
  }
  
  
  #then, sum by year
  df1 <- ddply(df1,.(year, age, state), summarise, count = sum(count))
  df1$cfr <- NA
  df1$cfr[df1$state=="Secondary-Infection"] <- cfr_secondary
  df1$cfr[df1$state=="Repeat-Infection"] <- cfr_repeat
  
  #and just get annual total mortality by age
  
  df1$mort <- ceiling(df1$count*df1$cfr)
  
  df.sum <- ddply(df1,.(year), summarise, mort = sum(mort))
  df.sum$state <- "deaths"
  # 
  # #split by a year
  # df.year <- dlply(df.sum,.(year, state))
  # 
  # # #get mean age
  # # mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  # # mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  # # 
  # 
  # #split by age and year and type
  # df.age <- dlply(df.sum,.(state, year, age))
  # 
  # 
  # 
  # 
  # dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type.mort, slim.quant=slim.quant))
  # 
  #dat.age$state <- factor(dat.age$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #mean.df$state <- factor(mean.df$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  #colz = c("secondary"="black", "tertiary"="royalblue3")
  
  p1 <- ggplot(df.sum) + facet_grid(~state) +
    geom_line(aes(x=year, y=mort)) #+
  #geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
  #geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + #scale_color_manual(values=colz) +
  #coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=100, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  #return(mean.df)
}

plot.mort.age(df = age.out.2019, criteria = "Secondary+Reinfection", cfr_secondary = 0.05,cfr_repeat = 0.01, slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.mort.age(df = age.out, criteria = "Secondary+Reinfection", cfr_secondary = 0.05,cfr_repeat = 0.00001, slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.mort.total(df = age.out.2019, criteria = "Secondary+Reinfection", cfr_secondary = 0.05,cfr_repeat = 0.00001, slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)
plot.mort.total(df = age.out, criteria = "Secondary+Secondary-Reinfection", cfr_secondary = 0.1,cfr_repeat = 0, slim.quant = 0.05, view.plot = T, save.plot = F, filename = F)

#and plot cases
plot.age.dist <- function(df, criteria, slim.quant, view.plot, save.plot, filename){
  
  #first, select what you want.
  if(criteria=="Secondary+Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection")
    
  }else if(criteria=="Secondary+Reinfection+ReSecondary"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection" | state == "Secondary-Infection-After-Reinfection")
    
  }else if(criteria == "Secondary-Only"){
    df1 = subset(df, state =="Secondary-Infection")
  }else if (criteria=="Secondary+Secondary-Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection")
    
  }
  
  #then, sum by year
  df1 <- ddply(df1,.(year, age), summarise, count = sum(count))
  df1$state = "Apparent-Infection"
  
  #split by a year
  df.year <- dlply(df1,.(year, state))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and type
  df.age <- dlply(df1,.(state, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type, slim.quant=slim.quant))
  
  #dat.age$state <- factor(dat.age$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #mean.df$state <- factor(mean.df$state, levels = c("Primary-Infection", "Secondary-Infection", "Tertiary-Infection", "Quaternary-Infection"))
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  #colz = c("secondary"="black", "tertiary"="royalblue3")
  
  p1 <- ggplot(dat.age) + facet_grid(~state) +
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + #scale_color_manual(values=colz) +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=100, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}

#this is the mortality one: we are going to leave it out and just discuss in the discussion
column.2 <- function(dat, dat.lci, dat.uci, year.start, perc_dss, perc_mort){
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year, state), summarise, count=sum(count))
  #dat.ts.lci <- ddply(denv.case.lci, .( year), summarise, count=sum(count))
  #dat.ts.uci <- ddply(denv.case.uci, .(year), summarise, count=sum(count))
  #dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total dss cases and deaths\nfrom deterministic simulation"
  # dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  # dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  # dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # 
  # dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  # dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  # names(dat.ts.lci.merge)[3] <- "lci"
  # names(dat.ts.uci.merge)[3] <- "uci"
  # 
  # dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  # dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  # dat.ts$lci_new <- NA
  # dat.ts$uci_new <- NA
  # dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  # dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  # dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  # dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  
  
  #dss cases should primarily be in the cases of secondary infection
  #or in cases of tertiary if we are assuming they are easier to detect due 
  #to heightened pathology in older patients
  dat.ts$count_dss <- 0
  dat.ts$count_dss[dat.ts$state=="Secondary-Infection"] <- dat.ts$count[dat.ts$state=="Secondary-Infection"]*perc_dss
  dat.ts$count_dss[dat.ts$state=="Tertiary-Infection"] <- dat.ts$count[dat.ts$state=="Tertiary-Infection"]*perc_dss_tert
  dat.ts$count_mort <- dat.ts$count_dss*perc_mort #some fraction of dss cases is the total mortality
  dat.dss <- ddply(dat.ts, .(plot_type, hyp, year), summarise, count = sum(count_dss))
  dat.mort <- ddply(dat.ts, .(plot_type, hyp, year), summarise, count = sum(count_mort))
  
  dat.dss$count_type <- "dss"
  dat.mort$count_type <- "deaths"
  
  dat.plot = rbind(dat.dss, dat.mort)
  dat.plot$count_type <- factor(dat.plot$count_type, levels = c("dss", "deaths"))
  
  colz = c("dss"="navy", "deaths"="darkred")
  
  p1 <- ggplot(dat.plot) + theme_bw() + facet_grid(hyp~plot_type) +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    #geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    #geom_line(aes(x=year, y=count), size=.3) + 
    geom_line(aes(x=year, y=count, color=count_type, group=count_type), size=.3) + 
    scale_color_manual(values=colz) +
    theme(panel.grid = element_blank(), 
          strip.background.y  =  element_blank(), #element_rect(fill="white"),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          legend.position = c(.15,.9),
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.placement = "outside",
          strip.text.x.top  = element_text(size=12),
          strip.text.y  = element_blank()) +# element_text(size=12)) +
    #      coord_cartesian(xlim=c(2015,2020), ylim=c(0,13000)) + #scale_y_log10() +
    #scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) +
    geom_vline(aes(xintercept=2007), linetype=2, color="red") +
    geom_vline(aes(xintercept=2012), linetype=2, color="red") +
    geom_vline(aes(xintercept=2019), linetype=2, color="red")
  
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
  
  return(p1)
}
