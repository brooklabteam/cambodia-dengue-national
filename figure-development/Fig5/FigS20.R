rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/test-sim/"))


#lambda fit from the sim model

#and lambda
load("hyp1-refit-lambda.Rdata")
head(hyp1.refit.lambda)
load("hyp2-refit-lambda.Rdata")
head(hyp2.refit.lambda)
load("hyp3-refit-lambda.Rdata")
head(hyp3.refit.lambda)
load("hyp4-refit-lambda.Rdata")
head(hyp4.refit.lambda)
load("hyp5-fit-lambda.Rdata")
hyp5.fit.lambda$sim_type <- "hyp5"
load("hyp6-fit-lambda.Rdata")
load("hyp7-fit-lambda.Rdata")



lambda.fit <- rbind(hyp1.refit.lambda, hyp2.refit.lambda, hyp3.refit.lambda, hyp4.refit.lambda, hyp5.fit.lambda, hyp6.fit.lambda, hyp7.fit.lambda)
unique(lambda.fit$convergence)
lambda.fit$label <- NA
lambda.fit$label[lambda.fit$sim_type=="hyp1"] <- "H0: Normal Simulation"
lambda.fit$label[lambda.fit$sim_type=="hyp2"] <- "H1: Elevated FOI in 2019"
lambda.fit$label[lambda.fit$sim_type=="hyp3"] <- "H2: 2019 Strain Intro + Maintenance"
lambda.fit$label[lambda.fit$sim_type=="hyp4"] <- "H3: 2019 Strain Intro + Replacement"
lambda.fit$label[lambda.fit$sim_type=="hyp5"] <- "H4: 2019 Strain Intro + Replacement + Waning Immunity"
lambda.fit$label[lambda.fit$sim_type=="hyp6"] <- "H5: 2019 Serotype Intro + Replacement + Tertiary Case Detection"
lambda.fit$label[lambda.fit$sim_type=="hyp7"] <- "H6: 2019 Sterotype Intro + Maintenance + Tertiary Case Detection"

lambda.fit$label <-factor(lambda.fit$label, levels=unique(lambda.fit$label))
#load the foi fits for cambodia
fit.dat <- read.csv(file = "prov-fits-FOI.csv", stringsAsFactors = F, header = T)
nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 60 before this

#now fit and recover age -FOI fixed
par.dat <- cbind.data.frame(year = 1921:2020,lambda=c(rep(.2,60),fit.dat$lambda[fit.dat$provname=="National"]))
par.dat$lambda[par.dat$lambda<0.05] <- 0.05
par.dat$lci <- c(rep(NA,60), fit.dat$lci[fit.dat$provname=="National"])
par.dat$uci <- c(rep(NA,60), fit.dat$uci[fit.dat$provname=="National"])

p1 <- ggplot(lambda.fit) + theme_bw() +
  geom_ribbon(data=par.dat, aes(x=year, ymin=lci, ymax=uci), alpha=.3) +
  geom_line(data=par.dat, aes(x=year, y=lambda)) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=label), alpha=.3) +
  geom_line(aes(x=year, y=lambda, color=label), size=1) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_text(size=16),
        legend.title = element_blank(),
        legend.position = c(.4,.75),
        axis.text = element_text(size=12), axis.title.x = element_blank()) + 
  ylab(bquote(lambda)) + coord_cartesian(ylim = c(0,1), xlim=c(1980,2020))

print(p1)

#and waning immunity

load("hyp1-fit-wane.Rdata")
load("hyp2-fit-wane.Rdata")
load("hyp3-fit-wane.Rdata")
load("hyp4-fit-wane.Rdata")
load("hyp5-fit-wane.Rdata")
load("hyp6-fit-wane.Rdata")
load("hyp7-fit-wane.Rdata")
#fit.dat <- read.csv(file="prov-fits-FOI.csv", header = T, stringsAsFactors = F)

hyp.fit <- rbind(hyp1.fit.wane, hyp2.fit.wane, hyp3.fit.wane, hyp4.fit.wane, hyp5.fit.wane, hyp6.fit.wane, hyp7.fit.wane)
unique(hyp.fit$convergence)


p2 <- ggplot(hyp.fit) + theme_bw() +
  geom_line(aes(x=year, y=sigma, color=sim_type), size=1, show.legend = F) +
  geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma, fill=sim_type), alpha=.3, show.legend = F) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=13),
        axis.text = element_text(size=12), 
        axis.title.x = element_blank()) + ylab(bquote(sigma)) 
  
print(p2)

FigSX <- cowplot::plot_grid(p1, p2, ncol=2, nrow=1, labels = c("A","B"), label_size = 22, align = "hv", label_x = c(0,-.01))


ggsave(filename = paste0(homewd, "/final-figures/FigS20.png"),
       plot = FigSX,
       units="mm",  
       width=110, 
       height=40, 
       scale=3, 
       dpi=300)
