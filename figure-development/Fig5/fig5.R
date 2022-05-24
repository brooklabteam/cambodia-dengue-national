rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

#load the transmission tree data
#here for both serotypes


#and load the metadata
homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#load the transmission tree data
#here for both serotypes
all.denv <- read.csv(file=paste0(homewd, "/data/AllDENVtransTreeDat.csv"), header = T, stringsAsFactors = F)
all.denv$distance <- all.denv$distance/1000 #convert to km

#are they in the same season?
all.denv$season <- "yes"
all.denv$season[all.denv$pairtime2-all.denv$pairtime1>.5] <- "no"

#now, only look at those within a season
all.denv <- subset(all.denv, season=="yes")

#dat=subset(all.denv,DENV.serotype=="DENV-1")
mrca_thresh=.5
#and write over to get the transmission trees
geothresh <- list(.2,.4,.6,.8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
geothresh <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
#here's the threshold list of distances over which to compute

#here, decides whether each pair is in a transmission chain, based on the mrca_thresh (in years)
is.trans.chain <- function(dat, mrca_thresh){
  dat$trans_chain <- 0
  dat$trans_chain[dat$tMRCA <= mrca_thresh] <- 1
  return(dat)
  
}
get.prop.chain<- function(thresh, dat){
  tot.chain <- sum(subset(dat, distance<=thresh)$trans_chain)
  N.pairs <- length(subset(dat, distance<=thresh)$trans_chain)
  prop.chain <- tot.chain/N.pairs
  
  CI <- binom.test(x=tot.chain, n=N.pairs, alternative = "two.sided", conf.level = .95)
  
  
  #and get CIs
  
  dat.new <- cbind.data.frame(distance = thresh, prop=prop.chain, prop_lci = CI$conf.int[1], prop_uci= CI$conf.int[2], tot_pairs_in_chain=tot.chain, tot_pairs=N.pairs)
  return(dat.new)
}
make.trans.chains <- function(dat, geothresh, mrca_thresh, character){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #and proportion of transmission chains
  dat.prop = lapply(geothresh, get.prop.chain, dat=dat1)
  prop.df = data.table::rbindlist( dat.prop)
  prop.df$character = character
  prop.df$DENV.serotype=unique(dat$DENV.serotype)
  prop.df$transchain_threshold = mrca_thresh
  prop.df$character <- character
  return(prop.df)
  
}

combine.chain.prop<-function(mrca.thresh){
  out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #head(out.prop1)
  out.prop1$DENV.subtype <- "DENV-1"
  out.prop3$DENV.subtype <- "DENV-2-Cosmopolitan"
  out.prop7$DENV.subtype <- "DENV-2-All"
  out.all <-rbind(out.prop1, out.prop3, out.prop7)
  return(out.all)
  
}
make.chain.diff <- function(mrca.thresh){
out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
head(out.prop1)

out.prop <- out.prop3 #Cosmopolitan is dominant and we subtract from it
out.prop$prop_diff <- out.prop3$prop - out.prop1$prop
out.prop$prop_diff_lci <-  out.prop$prop_lci-out.prop1$prop
out.prop$prop_diff_uci <-  out.prop$prop_uci-out.prop1$prop
out.prop$comp <- "DENV-1"


out.prop2 <- out.prop3 #Cosmopolitan is dominant and we subtract from it
out.prop2$prop_diff <- out.prop2$prop - out.prop7$prop
out.prop2$prop_diff_lci <-  out.prop2$prop_lci-out.prop7$prop
out.prop2$prop_diff_uci <-  out.prop2$prop_uci-out.prop7$prop
#out.prop2$comp <- "DENV-2-Asian-1"#
out.prop2$comp <- "All-DENV-2"

out.prop <- rbind(out.prop, out.prop2)

out.prop <- dplyr::select(out.prop, distance, prop_diff, prop_diff_lci, prop_diff_uci, transchain_threshold, comp)
return (out.prop)
}
mean.trans.chains <- function(dat, geothresh, mrca_thresh){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #get the probability that two cases are from the same chain
  prob <- sum(dat1$trans_chain)/length(dat1$trans_chain)
  prob.uci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[1]
  prob.lci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[2]
  
  eff.chains <- 1/prob #22.264 chains circulating in the region (Denv1)
  eff.chains.lci <- 1/prob.lci #17.019 chains circulating in the region (Denv1)
  eff.chains.uci <- 1/prob.uci #29.298 chains circulating in the region (Denv1)
  
  #and return
  chain.df <- cbind.data.frame(DENV.serotype=unique(dat1$DENV.serotype), DENV.subtype=unique(dat1$DENV.subtype), N_chains=eff.chains, N_chains_lci=eff.chains.lci, N_chains_uci = eff.chains.uci)
  
  return(chain.df)
  #pop size KP: 877,523... this is a little higher than the prediction for rural Thailand and lower
  #than the prediction for Bangkok
  #and get the mean number of transmission chains for this region that are in circulation
  
  
}


out.pt5 = make.chain.diff(.5)
out.pt5 = combine.chain.prop(mrca.thresh=.5)
head(out.pt5)
tail(out.pt5)
#out.prop = rbind(out.pt3, out.pt5, out.1)#, out.3)
#head(out.prop)
out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Cosmopolitan"))


colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")

pB <- ggplot(data=out.pt5) + theme_bw()+
  #facet_grid(dummy_label~.) +
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=DENV.subtype, group=DENV.subtype), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=DENV.subtype, group=DENV.subtype)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        #axis.title.x = element_blank(), axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        legend.title = element_blank(),
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 18),
        axis.text = element_text(size=14), legend.text = element_text(size=12),
        legend.position = c(.73,.84)) + coord_cartesian(ylim=c(0,.8), xlim=c(1,5), expand = F)+
  scale_color_manual(values=colzB) + 
  scale_fill_manual(values=colzB) + 
  #scale_color_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  #scale_fill_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  ylab("Proportion same transmission chain")  +
  xlab("Distance between cases (km)")
  


 
#and age 
#and get these ones that you can divide
out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1" & age<5), geothresh = geothresh, mrca_thresh = .5, character = "age < 5")
out.prop2 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1" & age>=5), geothresh = geothresh, mrca_thresh = .5, character = "age = 5+")

DENV1.prop <- rbind(out.prop1, out.prop2)
DENV1.prop$DENV.subtype <- "DENV-1"

out.prop4 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & age<5  & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = .5, character = "age < 5")
out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & age>=5 & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = .5, character = "age = 5+")

DENV2.cos.prop <- rbind(out.prop4, out.prop3)
DENV2.cos.prop$DENV.subtype <- "DENV-2-Cosmopolitan"

out.prop6 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & age<5  ), geothresh = geothresh, mrca_thresh = .5, character = "age < 5")
out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & age>=5), geothresh = geothresh, mrca_thresh = .5, character = "age = 5+")

DENV2.prop <- rbind(out.prop6, out.prop5)
DENV2.prop$DENV.subtype <- "DENV-2-All"

DENV2.prop <- rbind(DENV2.cos.prop, DENV2.prop)
DENV2.prop$DENV.subtype <- factor(DENV2.prop$DENV.subtype, levels=c("DENV-2-All", "DENV-2-Cosmopolitan"))

out.prop.all <- rbind(DENV1.prop, DENV2.prop)
head(out.prop.all)
out.prop.all$character <- factor(out.prop.all$character, levels=c("age < 5", "age = 5+"))
# 

colz=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")
colz_age <- c("age < 5" = "purple3", "age = 5+" = "turquoise4")
#pC <- ggplot(data=subset(out.prop.all, DENV.subtype=="DENV-1")) + theme_bw()+
pC <- #ggplot(data=out.prop.all) + theme_bw()+
  ggplot(data=subset(out.prop.all, DENV.subtype=="DENV-1")) + theme_bw()+
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  #facet_grid(DENV.subtype~.) +
  facet_grid(~DENV.subtype) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, 
                  fill=character), alpha=.3, show.legend = F) +
  geom_line(aes(x=distance, y=prop, color=character), show.legend = F) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        strip.background = element_rect(fill="white"), 
        strip.text = element_text(size = 18),
        axis.text = element_text(size=14), 
        legend.title = element_blank(),
        #axis.text.x = element_blank(), axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #panel.spacing = unit(4, "lines"),
        plot.margin = unit(c(.1,.3,.1,.9), "cm"),
        legend.text = element_blank())+
        #plot.margin = unit(c(.1,.1,.2,1), "cm"),
        coord_cartesian(ylim=c(0,.8), xlim=c(1,5), expand = F)+
  scale_color_manual(values=colz_age) +
  scale_fill_manual(values=colz_age) +
  ylab("Proportion same transmission chain") +
  xlab("Distance between cases (km)") 
#pC

pD <- #ggplot(data=out.prop.all) + theme_bw()+
  ggplot(data=subset(out.prop.all, DENV.subtype=="DENV-2-Cosmopolitan")) + theme_bw()+
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  #facet_grid(DENV.subtype~.) +
  facet_grid(~DENV.subtype) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=character), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=character)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 18),
        axis.text = element_text(size=14), legend.title = element_blank(),
        #axis.text.x = element_blank(), axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #panel.spacing = unit(4, "lines"),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        legend.text = element_text(size=12),
        #plot.margin = unit(c(.1,.1,.2,1), "cm"),
        legend.position = c(.85,.86)) + coord_cartesian(ylim=c(0,.8), xlim=c(1,5), expand = F)+
  scale_color_manual(values=colz_age) +
  scale_fill_manual(values=colz_age) +
  ylab("Proportion same transmission chain") +
  xlab("Distance between cases (km)") 
#pD





#and then do age for the

#and plot with the salje data
salje.dat <- read.csv(file=paste0(homewd, "/data/salje_chains.csv"), header = T, stringsAsFactors = F)
head(salje.dat)


#and get mean chains
denv.1 = subset(all.denv, DENV.serotype=="DENV-1")
denv.1$DENV.subtype <- "DENV-1"
denv.1.mean = mean.trans.chains(dat=denv.1, geothresh, mrca_thresh)

denv.2 = subset(all.denv, DENV.serotype=="DENV-2")
denv.2$DENV.subtype <- "DENV-2"
denv.2.mean = mean.trans.chains(dat=denv.2, geothresh, mrca_thresh)


denv.cosmo.2 = subset(all.denv, DENV.serotype=="DENV-2" & DENV.subtype=="DENV-2-Cosmopolitan")
denv.cosmo.2$DENV.subtype <- "DENV-2-Cosmopolitan"
denv.cosmo.2.mean = mean.trans.chains(dat=denv.cosmo.2, geothresh, mrca_thresh)

#all.denv.mean <- rbind(denv.1.mean, denv.2.mean, denv.cosmo.2.mean)
all.denv.mean <- rbind(denv.1.mean, denv.2.mean)
#fewer circulating chains for denv-2 vs. 1 and even fewer for co


#salje.dat$locale[salje.dat$locale=="Bangkok"] <- "Salje et al. Bangkok"
salje.dat$locale[salje.dat$locale=="Thai_Countryside"] <- "Rural Thailand"# "Salje et al.\nThai rural"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-1"] <- "Cambodia DENV-1"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2"] <- "Cambodia All-DENV-2"
#all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2-Cosmopolitan"] <- "Cambodia\nDENV-2-Cosmopolitan"

salje.dat$study <- "Salje et al. 2017"
all.denv.mean$study <- "Kampong Speu 2019"
#colznew <- c('Salje et al. Bangkok' = "darkorchid", 'Salje et al.\nThai rural' = "deeppink4", 'Cambodia DENV-1' = "forestgreen", 'Cambodia All-DENV-2' = "navy", 'Cambodia\nDENV-2-Cosmopolitan' = "dodgerblue")
colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia All-DENV-2' = "navy")
shapeznew <- c('Salje et al. 2017' = 21, "Kampong Speu 2019" = 24)
colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia All-DENV-2' = "navy", 'Salje et al. 2017' = "black", "Kampong Speu 2019" = "red")

#and plot
pA <- ggplot(data=salje.dat) + 
  geom_point(aes(x=pop_size, y=eff_chains, fill=locale, 
    shape=study, color=locale), size=4, color="black", show.legend = F) +
  geom_errorbar(aes(x=pop_size, ymin=lci, ymax=uci, color=locale)) +
  scale_y_log10() + 
  scale_x_log10(breaks=c(1e+03, 1e+04, 1e+05, 1e+06, 1e+07), 
                labels=c(1,10,100,1000, 10000)) + theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        plot.margin = unit(c(.3,.2,.1,.5), "cm"),
        #legend.background = element_rect(color="black"),
        legend.position = c(.3,.88)) +
  ylab("Effective # Chains") +
  xlab("Population Size (x1000)") +
  geom_point(data=all.denv.mean,aes(x=877523., y=N_chains, fill=DENV.subtype, shape=study, color=study), size=5,  stroke=2) +
  scale_color_manual(values=colznew, guide=NULL) + scale_fill_manual(values=colznew, guide=NULL) +
  scale_shape_manual(values = shapeznew) 

#pA


#and by sex
#and try by sex

geothresh <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 

out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1" & sex=="male"), geothresh = geothresh, mrca_thresh = .5, character = "male")
out.prop2 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1" & sex=="female"), geothresh = geothresh, mrca_thresh = .5, character = "female")

out.prop1$DENV.subtype<- "DENV-1"
out.prop2$DENV.subtype<- "DENV-1"

out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"  & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan" & sex=="male"), geothresh = geothresh, mrca_thresh = .5, character = "male")
out.prop4 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"  & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan" & sex=="female"), geothresh = geothresh, mrca_thresh = .5, character = "female")
out.prop3$DENV.subtype<- "DENV-2-Cosmopolitan"
out.prop4$DENV.subtype<- "DENV-2-Cosmopolitan"



out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"  & sex=="male"), geothresh = geothresh, mrca_thresh = .5, character = "male")
out.prop6 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"   & sex=="female"), geothresh = geothresh, mrca_thresh = .5, character = "female")
out.prop5$DENV.subtype<- "DENV-2-All"
out.prop6$DENV.subtype<- "DENV-2-All"

out.prop <- rbind(out.prop1, out.prop2, out.prop3, out.prop4, out.prop5, out.prop6)

out.prop$DENV.subtype <- factor(out.prop$DENV.subtype, levels=c("DENV-1","DENV-2-All", "DENV-2-Cosmopolitan" ))
out.prop$character <- factor(out.prop$character, levels=c("female", "male"))

colz=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")
colz <- c("female" = "palevioletred1", "male" = "cyan2")

#denv 1
pE <- #ggplot(data=out.prop) + theme_bw()+
  ggplot(data=subset(out.prop, DENV.subtype=="DENV-1")) + theme_bw()+
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  facet_grid(~DENV.subtype) +
  #facet_grid(DENV.subtype~.) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=character), alpha=.3,show.legend = F) +
  geom_line(aes(x=distance, y=prop, color=character), show.legend = F) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 16),
        axis.text = element_text(size=14), 
        legend.title = element_blank(),
        legend.text = element_blank(),
        plot.margin = unit(c(.1,.3,.1,.9), "cm"))+
        #panel.spacing = unit(4, "lines"),
        #axis.text.x = element_blank(),
        #legend.position = c(.06,.85))+
        #axis.title.x = element_blank(),
        #axis.ticks.x = element_blank()) +
  coord_cartesian(ylim=c(0,.8), xlim=c(1,5), expand = F)+
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) + 
  ylab("Proportion same transmission chain") +
  xlab("Distance between cases (km)") 

#print(pE)
pF <- #ggplot(data=out.prop) + theme_bw()+
  ggplot(data=subset(out.prop, DENV.subtype=="DENV-2-Cosmopolitan")) + theme_bw()+
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  facet_grid(~DENV.subtype) +
  #facet_grid(DENV.subtype~.) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=character), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=character)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 16),
        axis.text = element_text(size=14), legend.title = element_blank(),
        legend.text = element_text(size=12),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        #panel.spacing = unit(4, "lines"),
        #axis.text.x = element_blank(),
        legend.position = c(.86,.85))+
  #axis.title.x = element_blank(),
  #axis.ticks.x = element_blank()) +
  coord_cartesian(ylim=c(0,.8), xlim=c(1,5), expand = F)+
  scale_color_manual(values=colz) + 
  scale_fill_manual(values=colz) + 
  ylab("Proportion same transmission chain") +
  xlab("Distance between cases (km)") 




ptop <- cowplot::plot_grid(pA, pB, nrow = 1, labels = c("a", "b"), label_size = 22, hjust=c(-.5, 0))
pmid <- cowplot::plot_grid(pC, pD, nrow = 1, labels = c("c", "d"), label_size = 22, hjust=c(-.5, 0))
pbottom <- cowplot::plot_grid(pE, pF, nrow = 1, labels = c("e", "f"), label_size = 22, hjust=c(-.5, 0))
pall <- cowplot::plot_grid(ptop, pmid, pbottom, nrow = 3)

ggsave(file = paste0(homewd,"/final-figures/fig5.png"),
        plot=pall,
        units="mm",  
        width=80, 
        height=90, 
        scale=3.5, 
        dpi=300)
