rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(sjPlot)

homewd="/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

# here, we ask... is there a relationship between the saliva elisa 
# titer (for mosquitoes) and testing postivie or negative for 
# dengue infection?

# and we ask...are any other variables in our dataset associated with dengue positivity?

# we use IDseq and PAGODAs data from the start of both studies through to the end of 2020

#load the combined metadata
all.dat <- read.csv(file = paste0(homewd, "/data/IDseq_PAGODAS_ALL_metadata_through_2020.csv"), header = T, stringsAsFactors = F)
head(all.dat)
all.dat$date <- as.Date(all.dat$date)

length(unique(all.dat$NIH.ID)) #697 unique in 760 entries, so 63 are duplicate entries

all.dat$dengue_result[all.dat$dengue_result=="neg"] <- 0
all.dat$dengue_result[all.dat$dengue_result=="pos"] <- 1
all.dat$dengue_result <- as.numeric(all.dat$dengue_result)
sum(all.dat$dengue_result) #106 positives

length(unique(all.dat$NIH.ID[all.dat$dengue_result==1])) #106 - no one tests positive twice
dup.id <- cbind.data.frame(NIH.ID=unique(all.dat$NIH.ID[duplicated(all.dat$NIH.ID)]))  #46 repeats

#and get the duplicated in case we want to examine -- they all have different dates
dup.df <-merge(dup.id, all.dat, by="NIH.ID")

all.dat$sex <- as.factor(all.dat$sex)
all.dat$DENV.serotype[all.dat$DENV.serotype=="neg"] <- NA
all.dat$DENV.serotype <- as.factor(all.dat$DENV.serotype)

#don't need path data here
all.dat<- dplyr::select(all.dat, -(WBC), -(platelets))

#and we query impact on dengue results

#then, merge with elisa data from the time of infection
elisa.pag <- read.csv(file = paste0(homewd, "/data/pagodas-elisa-dat.csv"), header = T, stringsAsFactors = F)
elisa.idseq <- read.csv(file = paste0(homewd, "/data/IDseq_saliva_elisa.csv"), header = T, stringsAsFactors = F)

#and slim to just what is needed
elisa.pag <- dplyr::select(elisa.pag,seq.ID, date, elisa)
elisa.idseq <- dplyr::select(elisa.idseq,seq.ID, date, elisa)

elisa.idseq$date <- as.Date(elisa.idseq$date, format = "%m/%d/%y")
elisa.pag$date <- as.Date(elisa.pag$date, format = "%m/%d/%y")
elisa.dat <- rbind(elisa.idseq, elisa.pag)

#and merge
head(elisa.dat)
names(elisa.dat)[1] <- "NIH.ID"
head(all.dat)
all.dat <- merge(all.dat, elisa.dat, by =c("NIH.ID", "date"), all.x = T)
head(all.dat)
tail(all.dat)
length(unique(all.dat$NIH.ID)) #697 - some people (5) reported fevers multiple times 
all.dat$elisa <- as.numeric(all.dat$elisa) #lots of NAs currently

#first, we ask, is higher mosquito biting (in elisas) associated with positive dengue infection?
all.dat$NIH.ID <- as.factor(all.dat$NIH.ID)

elisa.dat = subset(all.dat, !is.na(elisa))

psaliva <- glmer(dengue_result~elisa + (1|NIH.ID), data=elisa.dat, family="binomial")
summary(psaliva)  #no effect - positive trend with dengue result of 1 but not significant


plot_model(psaliva, type="est")


#next, we ask, are any other variables in our dataset associated with dengue positivity?

#first, simplify some of the predictors

length(all.dat$socio[all.dat$socio=="upper"]) #just one
socio.sum <- ddply(all.dat, .(socio), summarise, N=length(socio))
socio.sum

net.sum <- ddply(all.dat, .(net_use), summarise, N=length(socio))
net.sum

housing.sum <- ddply(all.dat, .(housing), summarise, N=length(socio))
housing.sum

insecticide.sum <- ddply(all.dat, .(insecticide_use), summarise, N=length(socio))
insecticide.sum

larvicide.sum <- ddply(all.dat, .(larvicide_use), summarise, N=length(socio))
larvicide.sum

mosquito_coil_use.sum <- ddply(all.dat, .(mosquito_coil_use), summarise, N=length(socio))
mosquito_coil_use.sum

all.dat$socio[all.dat$socio=="upper"] <- "middle"
all.dat$socio[all.dat$socio=="very poor"] <- "lower"

all.dat$net_use[all.dat$net_use=="all the time"] <- "regularly"
all.dat$net_use[all.dat$net_use=="never"] <- "rarely"

all.dat$housing[all.dat$housing=="temporary"] <- "other"
all.dat$housing <- as.factor(all.dat$housing)
all.dat$net_use <- as.factor(all.dat$net_use)
all.dat$mosquito_coil_use <- as.factor(all.dat$mosquito_coil_use)
all.dat$larvicide_use <- as.factor(all.dat$larvicide_use)
all.dat$insecticide_use <- as.factor(all.dat$insecticide_use)
all.dat$socio <- as.factor(all.dat$socio)
all.dat$age <- as.numeric(all.dat$age)


#remove elisa
stat.dat <- dplyr::select(all.dat, dengue_result, age, sex, housing, socio, net_use, insecticide_use, mosquito_coil_use, NIH.ID)
stat.dat <- stat.dat[complete.cases(stat.dat),]
pdengue_global <- lmerTest::lmer(dengue_result~age + sex + housing + socio +
                                        net_use + insecticide_use + mosquito_coil_use + 
                                        (1|NIH.ID),
                                        data=stat.dat,
                                        na.action = na.fail)

summary(pdengue_global) 
# slight sig effect of age (lower ages--so kids--associated with positive infection. might reverse if not active surveillance
# middle/upper socioeconomic status is negatively associated with dengue infection
# living in a house is slightly negatively associated with dengue infection

#try dredge for model selection
library(MuMIn)
library(relaimpo)

out.comp <- dredge(global.model = pdengue_global) #net use and soc and females slightly more likely to encounter
head(out.comp)

top_mod <- lmerTest::lmer(dengue_result~(1|NIH.ID),data=stat.dat,
                                             na.action = na.fail)
summary(top_mod) #random intercept

#very few 'risk factors' of infection were identified. 
#No effect of mosquito saliva and a small effect only of age, housing, 
# and higher socioeconomic status, all negatively associated with 
#dengue infection,

#for socio, we binned "very poor" (N=12) with "lower" (N=202)
#and "upper" (N=1) with "middle" (N=322)