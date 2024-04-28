rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/fit-sim/"))



#lambda fit from the sim model

#and lambda
load("hyp0-fit-lambda.Rdata")
load("hyp1-fit-lambda.Rdata")
load("hyp2-fit-lambda-2007.Rdata")
load("hyp2-fit-lambda-2019.Rdata")
#load("hyp4-fit-lambda.Rdata")

# load("hyp1-refit-lambda.Rdata")
# head(hyp1.refit.lambda)
# load("hyp2-refit-lambda.Rdata")
# head(hyp2.refit.lambda)
# load("hyp3-refit-lambda.Rdata")
# head(hyp3.refit.lambda)
# load("hyp4-refit-lambda.Rdata")
# head(hyp4.refit.lambda)
# load("hyp5-fit-lambda.Rdata")
# hyp5.fit.lambda$sim_type <- "hyp5"
# load("hyp6-fit-lambda.Rdata")
# load("hyp7-fit-lambda.Rdata")



lambda.fit <- rbind(hyp0.fit.lambda, hyp1.fit.lambda, hyp2.fit.lambda.2019, hyp2.fit.lambda.2007)
unique(lambda.fit$convergence)
lambda.fit$label <- NA
lambda.fit$label[lambda.fit$sim_type=="hyp0"] <- "H0: Normal Simulation"
lambda.fit$label[lambda.fit$sim_type=="hyp1"] <- "H1: Increasing Tertiary"
lambda.fit$label[lambda.fit$sim_type=="hyp2-2019"] <- "H2: 2019 Strain Intro + Replacement + Waning Immunity"
lambda.fit$label[lambda.fit$sim_type=="hyp2-2007"] <- "H2: 2007 Strain Intro + Replacement + Waning Immunity"



lambda.fit$label <-factor(lambda.fit$label, levels=unique(lambda.fit$label))
#load the foi fits for cambodia
fit.dat <- read.csv(file = paste0(homewd,"/data/prov-fits-FOI.csv"), stringsAsFactors = F, header = T)
nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 60 before this

#now fit and recover age -FOI fixed
par.dat <- cbind.data.frame(year = 1921:2020,lambda=c(rep(.9,60),fit.dat$lambda[fit.dat$provname=="National"]))
par.dat$lambda[par.dat$year<1999] <- 0.9
par.dat$lci <- c(rep(NA,60), fit.dat$lci[fit.dat$provname=="National"])
par.dat$uci <- c(rep(NA,60), fit.dat$uci[fit.dat$provname=="National"])
par.dat$lci[par.dat$year<1999] <- 0.9
par.dat$uci[par.dat$year<1999] <- 0.9

p1 <- ggplot(lambda.fit) + theme_bw() +
  geom_ribbon(data=par.dat, aes(x=year, ymin=lci, ymax=uci), alpha=.3) +
  geom_line(data=par.dat, aes(x=year, y=lambda)) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=label), alpha=.3) +
  geom_line(aes(x=year, y=lambda, color=label), size=1) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size=14),
        legend.title = element_blank(),
        legend.position = c(.4,.75),
        axis.text = element_text(size=12), axis.title.x = element_blank()) + 
  ylab(bquote("annual force of infection,"~lambda))+ 
  coord_cartesian(ylim = c(0,.7), xlim=c(2000,2020))

print(p1)
setwd(paste0(homewd, "/figure-development/Fig5/fit-sim/fit-sim-wane/out-wane/"))

#and waning immunity
load("hyp0-fit-wane.Rdata")
load("hyp1-fit-wane.Rdata")
load("hyp2-fit-wane-2019.Rdata")
load("hyp2-fit-wane-2007.Rdata")

#fit.dat <- read.csv(file="prov-fits-FOI.csv", header = T, stringsAsFactors = F)
hyp0.fit.wane$sigma <- exp(hyp0.fit.wane$sigma)
hyp1.fit.wane$sigma <- exp(hyp1.fit.wane$sigma)
hyp.fit <- rbind(hyp0.fit.wane, hyp1.fit.wane, hyp2.fit.wane.2019,  hyp2.fit.wane.2007)
hyp.fit <- dplyr::select(hyp.fit, year, sigma, convergence,lci_sigma, uci_sigma, sim_type)

hyp.fit$sim_type <- factor(hyp.fit$sim_type, levels= c("hyp0_wane", "hyp1_wane","hyp2_wane_2019", "hyp2_wane_2007"))

#hyp3.fit.wane <- dplyr::select(hyp3.fit.wane, year, sigma, convergence,lci_sigma, uci_sigma, sim_type)
#hyp.fit<-rbind(hyp.fit, hyp3.fit.wane)

unique(hyp.fit$convergence)

#hyp.fit$sigma <- exp(hyp.fit$sigma)



p2 <- ggplot(hyp.fit) + theme_bw() +#coord_cartesian(ylim=c(0,.01)) +
  geom_point(aes(x=year, y=sigma, color=sim_type), size=2, show.legend = F) +
  geom_line(aes(x=year, y=sigma, color=sim_type), size=1, show.legend = F) +
  #geom_linerange(aes(x=year, ymin=lci_sigma, ymax=uci_sigma, color=sim_type),  show.legend = F) +
  ylab(bquote("annual rate of waning multitypic immunity,"~sigma))+ 
  #geom_point(aes(x=year, y=sigma, color=sim_type), size=1, show.legend = F, size=3) +
  #geom_linerange(aes(x=year, ymin=lci_sigma, ymax=uci_sigma, color=sim_type), linetype=3, size=1,show.legend = F) +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size=12), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank()) 
  
print(p2)

FigS22 <- cowplot::plot_grid(p1, p2, ncol=2, nrow=1, labels = c("A","B"), label_size = 22, align = "hv", label_x = c(0,-.01))


ggsave(filename = paste0(homewd, "/final-figures/FigS22.png"),
       plot = FigS22,
       units="mm",  
       width=110, 
       height=40, 
       scale=3, 
       dpi=300)
