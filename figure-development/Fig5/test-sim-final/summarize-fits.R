rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/test-sim/"))


#lambda fit from the sim model

#and lambda
load("hyp1-fit-lambda.Rdata")
head(hyp1.fit.lambda)
load("hyp2-fit-lambda.Rdata")
head(hyp2.fit.lambda)
load("hyp3-fit-lambda.Rdata")
head(hyp3.fit.lambda)
hyp3.fit.lambda$sim_type <- "hyp3"
load("hyp4-fit-lambda.Rdata")
head(hyp4.fit.lambda)
hyp4.fit.lambda$sim_type <- "hyp4"
load("hyp5-fit-lambda.Rdata")
hyp5.fit.lambda$sim_type <- "hyp5"

lambda.fit <- rbind(hyp1.fit.lambda, hyp2.fit.lambda, hyp3.fit.lambda, hyp4.fit.lambda, hyp5.fit.lambda)
unique(lambda.fit$convergence)
unique(lambda.fit$sim_type[lambda.fit$convergence==1]) # all except for 5!


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
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=sim_type), alpha=.3) +
  geom_line(aes(x=year, y=lambda, color=sim_type), size=1) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=13),
        axis.text = element_text(size=10), axis.title.x = element_blank()) + 
  ylab(bquote(lambda)) + coord_cartesian(ylim = c(0,1))

print(p1)

#and waning immunity

load("hyp1-fit-wane.Rdata")
load("hyp2-fit-wane.Rdata")
load("hyp3-fit-wane.Rdata")
load("hyp4-fit-wane.Rdata")
load("hyp5-fit-wane.Rdata")
#fit.dat <- read.csv(file="prov-fits-FOI.csv", header = T, stringsAsFactors = F)

hyp.fit <- rbind(hyp1.fit.wane, hyp2.fit.wane, hyp3.fit.wane, hyp4.fit.wane, hyp5.fit.wane)
unique(hyp.fit$convergence)


p2 <- ggplot(hyp.fit) + theme_bw() +
  geom_line(aes(x=year, y=sigma, color=sim_type), size=1) +
  geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma, fill=sim_type), alpha=.3) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size=13),
        axis.text = element_text(size=10), axis.title.x = element_blank()) + ylab(bquote(sigma)) 
  
print(p2)

#save this plot

ggsave(filename = "sigma-fitted-foi.png",
       plot = p1,
       units="mm",  
       width=55, 
       height=40, 
       scale=3, 
       dpi=300)


