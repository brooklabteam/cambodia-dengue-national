rm(list = ls())


library(ggplot2)
#library(bobfunctions2)
library(plyr)
library(dplyr)
library(stringr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/"))

#load the output from the previous trials
#load(paste0(homewd,"/figure-development/Fig5/sim-new/test.no.intro.small.Rdata"))
#load(paste0(homewd,"/figure-development/Fig5/sim-new/test.2007.intro.small.Rdata"))

load(paste0(homewd,"/figure-development/Fig5/sim-new/test.no.intro.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/test.2007.intro.Rdata"))


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
  dat1$state[dat1$class=="I41"] <- "Not-Possible"
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
  df.sum.I = subset(  df.sum.I, count>0)
  return(df.sum.I)
}

age.out.no.intro = summarise.age.dist.wane(dat=test.no.intro, year.start = min(test.no.intro$year))
age.out.intro = summarise.age.dist.wane(dat=test.2007.intro, year.start = min(test.2007.intro$year))

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
    df.perc <- cbind.data.frame(year=unique(df1.tert$year), perc_obs=seq(0,.15, length.out = length(unique(df1.tert$year))))
    df1.tert <- merge(df1.tert, df.perc, by="year", all.x=T)
    df1.tert$count <- df1.tert$count*df1.tert$perc_obs
    df1.tert <- dplyr::select(df1.tert, -(perc_obs))
    
    df1 <- rbind(df1.sec, df1.tert)
    
    
    
  }
  
  
  
  df1$case_type = "symptomatic"
  
  return(df1)
  
  
}

age.out.no.intro.sub = subset(age.out.no.intro, year <2021)
age.out.no.intro.sub = select.symptom(df=age.out.no.intro.sub,criteria = "Secondary-Extension")
age.out.no.intro.sub$hyp = "no-intro"

age.out.2007.intro.sub = subset(age.out.intro, year <2021)
age.out.2007.intro.sub  = select.symptom(df=age.out.2007.intro.sub ,criteria = "Secondary-Extension")
age.out.2007.intro.sub$hyp = "2007-intro"



age.out.no.intro.sub.tert = subset(age.out.no.intro, year <2021)
age.out.no.intro.sub.tert = select.symptom(df=age.out.no.intro.sub.tert,criteria = "Increasing-Tertiary")
age.out.no.intro.sub.tert$hyp = "no-intro-tert"

age.out.2007.intro.sub.tert = subset(age.out.intro, year <2021)
age.out.2007.intro.sub.tert  = select.symptom(df=age.out.2007.intro.sub.tert,criteria = "Increasing-Tertiary")
age.out.2007.intro.sub.tert$hyp = "2007-intro-tert"





#put all the data together
comp.dat <- rbind(age.out.no.intro.sub, age.out.2007.intro.sub, age.out.no.intro.sub.tert, age.out.2007.intro.sub.tert)
#comp.dat$hyp <- factor(comp.dat$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
comp.dat$hyp <- factor(comp.dat$hyp, levels = c("no-intro",
                                                "2007-intro",
                                                "no-intro-tert", 
                                                "2007-intro-tert"))
comp.dat.test <- comp.dat
#and save for fitting
save(comp.dat.test, file = "comp-dat-sim-test.Rdata") 

# comp.dat.lci <- rbind(age.sub.H0.lci, age.sub.tert.lci, age.sub.2019.lci.b, age.sub.2019.lci.c, age.sub.2019.lci, age.sub.2007.lci)
# comp.dat.lci$hyp <- factor(comp.dat.lci$hyp, levels = c("H0: Normal Demographic\nSimulation", 
#                                                         "H1: Increasing Tertiary\nCase Detection",
#                                                         "H2: Genotype Intro\n+ Normal Immunity (2019)",
#                                                         "H3: Genotype Intro + Increasing\nTertiary Case Detection",
#                                                         "H4: Genotype Replacement\n+ Waning Immunity (2019)", 
#                                                         "H4: Genotype Replacement\n+ Waning Immunity (2007)"))
# 
# save(comp.dat.lci, file = "comp-dat-sim-lci.Rdata") 
# comp.dat.uci <- rbind(age.sub.H0.uci, age.sub.tert.uci, age.sub.2019.uci.b, age.sub.2019.uci.c, age.sub.2019.uci,  age.sub.2007.uci)
# comp.dat.uci$hyp <- factor(comp.dat.uci$hyp, levels = c("H0: Normal Demographic\nSimulation", 
#                                                         "H1: Increasing Tertiary\nCase Detection",
#                                                         "H2: Genotype Intro\n+ Normal Immunity (2019)",
#                                                         "H3: Genotype Intro + Increasing\nTertiary Case Detection",
#                                                         "H4: Genotype Replacement\n+ Waning Immunity (2019)", 
#                                                         "H4: Genotype Replacement\n+ Waning Immunity (2007)"))
# 
# save(comp.dat.uci, file = "comp-dat-sim-uci.Rdata") 
# 
# 
# 
# 
# 
# 
