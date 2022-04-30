

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)


homewd= "/home/rstudio"
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
  df.out$n <- sum(df.out$Nage)
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

#dat.tot <- ddply(dat, .(epiwk, sex), summarise, cases = length(year))





#first, plot total cases as a time series:
dat.tot <- ddply(dat, .(epiwk), summarise, cases = length(year))
dat.tot$epiwk <- as.Date(dat.tot$epiwk)
dat.tot # seems like single week not two weeks
#pA




library(scales)
dat.tot$week = week(dat.tot$epiwk)
dat.tot$year = year(dat.tot$epiwk)
dat.tot$year = factor(dat.tot$year)



colz = hue_pal()(length(unique(dat.tot$year)))
names(colz) <- unique(as.character(dat.tot$year))
colz["2007"] <- "orange"
colz["2012"] <- "red"
colz["2019"] <- "darkred"


colz.long <- rep("black", length(colz))
names(colz.long) <- unique(as.character(dat.tot$year))
colz.long["2007"] <- "orange"
colz.long["2012"] <- "red"
colz.long["2019"] <- "darkred"

pA <- ggplot(data=dat.tot) + geom_line(aes(x=epiwk, cases, color=year), show.legend = F, size=.7) +
  scale_color_manual(values=colz.long) +
  theme_bw() + ylab("reported dengue cases,\nCambodia (2002-2020)")+
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.x = element_blank())


dat.tot$year = factor(dat.tot$year)


colz = hue_pal()(length(unique(dat.tot$year)))
names(colz) <- unique(as.character(dat.tot$year))
colz["2007"] <- "orange"
colz["2012"] <- "tomato2"
colz["2019"] <- "red3"

pS1 <- ggplot(data=dat.tot) + geom_point(aes(x=week, y=cases, color=year)) +
  geom_line(aes(x=week, y=cases, color=year)) +
  theme_bw() + ylab("reported dengue cases,\nCambodia (2002-2020)")+
  xlab("week of year") +
  theme(panel.grid = element_blank(), 
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16)) +
  scale_color_manual(values=colz) 

# pS1




ggsave(file = paste0(homewd, "/FigS1.png"),
       plot=pS1,
       units=c("mm"),
       width=65,
       height=45,
       scale=3,
       dpi=300)




unique(dat$age) #round to years
dat$age <- ceiling(dat$age)
year.split <- dlply(dat, .(year))

dat.sum <- lapply(year.split, sum.yr, age_vect = 1:max(dat$age))
dat.sum <- data.table::rbindlist(dat.sum)
dat.sum$year <- rep(unique(dat$year), each=(max(dat$age)+1))
dat.sum$year <- as.factor(dat.sum$year)

dat.sum$prop_prev <- dat.sum$Nage/dat.sum$n




#and plot as age structured incidence
pB <- ggplot(data=dat.sum) + 
  #geom_point(aes(x=age, y=Nage, color=year), show.legend = F) +
  #geom_line(aes(x=age, y=Nage, color=year, group=year), show.legend = F) + 
  geom_point(aes(x=age, y=prop_prev, color=year), show.legend = F) +
  geom_line(aes(x=age, y=prop_prev, color=year, group=year), show.legend = F) + 
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=prop_prev, color=year), show.legend = F) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=prop_prev, color=year, group=year), show.legend = F) + 
  #ylab("annual case indicence") +
  ylab("annual proportion\nof case indicence") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        legend.text = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(size=14),
        axis.ticks.x = element_blank()) + guides(color=guide_legend(ncol=2))
# pB


# plot pB_blk_lines



colz = hue_pal()(length(unique(dat.tot$year)))
names(colz) <- unique(as.character(dat.tot$year))
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
colz["2019"] <- "red3"


pB_blk_lines <- ggplot(data=dat.sum) + 
  geom_point(aes(x=age, y=prop_prev, color=year), show.legend = F) +
  geom_line(aes(x=age, y=prop_prev, color=year, group=year), show.legend = F) + 
  # geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=prop_prev, group=year), color="black", size=2, show.legend=F) + 
  
  geom_point(data=subset(dat.sum, year==2019), aes(x=age, y=prop_prev), color= "red", show.legend = F) +
  geom_point(data=subset(dat.sum, year==2012), aes(x=age, y=prop_prev), color="tomato", show.legend = F) +
  geom_point(data=subset(dat.sum, year==2007), aes(x=age, y=prop_prev), color="orange", show.legend = F) +
  geom_line(data=subset(dat.sum, year==2019), aes(x=age, y=prop_prev,group=year),color="red",  show.legend = F) + 
  geom_line(data=subset(dat.sum, year==2012), aes(x=age, y=prop_prev, group=year), color="tomato", show.legend = F) + 
  geom_line(data=subset(dat.sum, year==2007), aes(x=age, y=prop_prev, group=year), color="orange", show.legend = F) + 
  ylab("annual proportion\nof case indicence") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  
  # scale_colour_grey(start=0.8,  end = 0)+
  scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), #legend.title = element_blank(),
        # legend.text = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(size=14),
        axis.ticks.x = element_blank()) + guides(color=guide_legend(ncol=2))
pB_blk_lines




pD <- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  
  scale_color_manual(values=colz) +
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))

# pD


# plot pD_blk_lines

pD_blk_lines<- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) +
  
  
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  geom_line(data=subset(dat.sum, year==2007), aes(x=age, y=cum_prop_cases, group=year), color="grey45", size=1) +
  geom_line(data=subset(dat.sum, year==2012), aes(x=age, y=cum_prop_cases, group=year), color="grey29", size=1) +
  geom_line(data=subset(dat.sum, year==2019), aes(x=age, y=cum_prop_cases, group=year), color="black", size=1) +
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  scale_color_manual(values=colz) +
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))# + geom_segment(aes(x = 21.9, y = 0.03, xend = 24.6, yend = 0.03), size=2)
# pD_blk_lines



pD_blk_lines<- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) +
  
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  geom_line(data=subset(dat.sum, year==2007), aes(x=age, y=cum_prop_cases, group=year), color="tan4", size=1) +
  geom_line(data=subset(dat.sum, year==2012), aes(x=age, y=cum_prop_cases, group=year), color="saddlebrown", size=1) +
  geom_line(data=subset(dat.sum, year==2019), aes(x=age, y=cum_prop_cases, group=year), color="coral4", size=1) +
  
  scale_color_manual(values=colz) +
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))# + geom_segment(aes(x = 21.9, y = 0.03, xend = 24.6, yend = 0.03), size=2)
# pD_blk_lines


colz = hue_pal()(length(unique(dat.tot$year)))
names(colz) <- unique(as.character(dat.tot$year))
colz["2007"] <- "orange"
colz["2012"] <- "tomato2"
colz["2019"] <- "red3"

pD_blk_lines<- ggplot(data=subset(dat.sum, year!=2007 & year!=2012 & year!=2019 ))+ geom_point(data=subset(dat.sum, year!=2007 & year!=2012 & year!=2019 ), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year!=2007 & year!=2012 & year!=2019 ), aes(x=age, y=cum_prop_cases, color=year, group=year), show.legend = T) +
  
  geom_point(data=subset(dat.sum, year==2007 ), aes(x=age, y=cum_prop_cases), color="orange") +
  geom_point(data=subset(dat.sum, year==2012 ), aes(x=age, y=cum_prop_cases), color="tomato") +
  geom_point(data=subset(dat.sum, year==2019 ), aes(x=age, y=cum_prop_cases), color="red") +
  
  geom_line(data=subset(dat.sum, year==2007), aes(x=age, y=cum_prop_cases, group=year), color="orange", size=1) +
  geom_line(data=subset(dat.sum, year==2012), aes(x=age, y=cum_prop_cases, group=year), color="tomato", size=1) +
  geom_line(data=subset(dat.sum, year==2019), aes(x=age, y=cum_prop_cases, group=year), color="red", size=1) +
  # scale_color_gradientn(colours = rainbow(5))+
  scale_colour_grey(start=0.8,  end = 0)+
  #  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))# + geom_segment(aes(x = 21.9, y = 0.03, xend = 24.6, yend = 0.03), size=2)
# pD_blk_lines

# gradiant color generator https://colordesigner.io/gradient-generator 

colz = hue_pal()(length(unique(dat.tot$year)))
names(colz) <- unique(as.character(dat.tot$year))
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



pD_blk_lines<- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) +
  
  # geom_point(data=subset(dat.sum, year==2007 ), aes(x=age, y=cum_prop_cases), color="orange") +
  # geom_point(data=subset(dat.sum, year==2012 ), aes(x=age, y=cum_prop_cases), color="tomato") +
  # geom_point(data=subset(dat.sum, year==2019 ), aes(x=age, y=cum_prop_cases), color="red") +
  # 
  # geom_line(data=subset(dat.sum, year==2007), aes(x=age, y=cum_prop_cases, group=year), color="orange", size=1) +
  # geom_line(data=subset(dat.sum, year==2012), aes(x=age, y=cum_prop_cases, group=year), color="tomato", size=1) +
  # geom_line(data=subset(dat.sum, year==2019), aes(x=age, y=cum_prop_cases, group=year), color="red", size=1) +
  # 
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  
  scale_color_manual(values=colz) +
  
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))
# + geom_segment(aes(x = 21.9, y = 0.03, xend = 24.6, yend = 0.03), size=2)
# pD_blk_lines


pD <- ggplot(data=dat.sum)+ geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  geom_point(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(data=subset(dat.sum, year==2007 | year==2012| year==2019), aes(x=age, y=cum_prop_cases, color=year, group=year)) + 
  
  scale_color_manual(values=colz) +
  ylab("annual cumulative\nproportion of cases") +
  theme_bw() + coord_cartesian(xlim = c(0,30)) + 
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size=16),
        legend.position = c(.8,.45)) + guides(color=guide_legend(ncol=2))

# pD



#now calculate annual FOI, using the Cummings method
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
  return(out.all)
  
  
}

#try for one year
comp.all.NS.2002 <- comp.all.NS.multi(lambda.guess = c(.01),
                                      age_vect = 1:35,
                                      dat=year.split[[1]])
#and plot


#and another year
out.2003.likprop <- comp.all.NS.multi(lambda.guess = c(.2),
                                      age_vect = 1:max(dat$age),
                                      dat=year.split[[2]])

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

pFig1 <- cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1,1.1))
pFig1



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

pFig1_blk_lines <- cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1,1.1))# +cowplot::draw_text('epi years of\n\t2007/2012/2019', x = 0.945, y = 0.1, size = 10, hjust = 0.5, vjust = 0.5)




ggsave(file = paste0(homewd, "/Fig1.png"),
       plot=pFig1_blk_lines,
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
         filename = paste0(homewd, "/FigS2.png"))















