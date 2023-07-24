rm(list = ls())

# Load packages
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(mgcv)
library(RColorBrewer)
library(lubridate)

# Set wd
homewd = "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(paste0(homewd, "ncdp"))


# Load and clean NDCP dataset
NDCP=read_excel("Dengue all years_2002-2020_CompletedV2.xlsx")
NDCP <- dplyr::rename(NDCP,diagnostic=dianostic)
NDCP <- NDCP %>% mutate(diagnosis = recode(NDCP$diagnostic, "DF"="df", "DSS"="dss")) %>% mutate(sex = recode(NDCP$sex, "m"="M"))

# Load mortality dataset
dta0 <- read_excel("Dengue_2002-20020_mortality.xlsx")
dta0 %>% pivot_wider(names_from=Outcome,values_from="No") -> dta0
dta0$Alive <- dta0$Total - dta0$Deaths

# Expand case fatality rates to a new dataframe with a binomial output of 0 for non-fatal and 1 for fatal cases
outdta1 <- NULL

for(i in 1:19){
  # Year 2002 through 2020
  for(j in 1:12){
    # Month 1 through 12
    Outcome <- c(rep(1,dta0$Deaths[12*(i-1)+j]),rep(0,dta0$Alive[12*(i-1)+j]))
    outdta<- as.data.frame(Outcome)
    outdta$Month <- dta0$Month[12*(i-1)+j]
    outdta$Year <- dta0$Year[12*(i-1)+j]
    
    outdta1 <- rbind(outdta1, outdta)
  }
}


# Load climate data
ppt <- read.csv("monthly_ppt.csv") #sumppt is total monthly precipitation averaged across provinces, meanppt is mean daily precipitation. For these analyses we keep sumppt

temp <- read.csv("monthly_temp.csv")  #meant, mint, maxt are average, minimum, and maximum daily temperatures for each month, respectively. For these analyses we keep meant
climate <- full_join(ppt,temp)

# Calculate lag for ppt and temp and create a lagged version of the crude case number dataset

#vThis is the climate dataset
clim.sum <- climate %>% mutate(epiyear=year-2002) %>% mutate(epimonth=12*epiyear+month)

# This is the crude case number dataset
dat.sum <- NDCP %>% group_by(year,month) %>% dplyr::summarise(cases=n()) %>% mutate(epiyear=year-2002) %>% mutate(epimonth=12*epiyear+month)

# Merge the two
merge.dat <- merge(dat.sum,clim.sum,by="epimonth")
merge.melt <- melt(merge.dat,id.vars=c("epimonth"))
case.dat1 = subset(merge.melt,variable=="cases")
case.dat1$variable <- "meant"
case.dat2 = subset(merge.melt,variable=="cases")
case.dat2$variable <- "sumppt"
case.dat <- rbind(case.dat1,case.dat2)


# Plot cases and climate variables
ggplot(data=subset(merge.melt, variable=="meant" |variable=="sumppt")) + 
  geom_line(data=case.dat, aes(x=epimonth, y=value),  size=1, alpha=.2) +
  geom_point(data=case.dat, aes(x=epimonth, y=value), size=3, alpha=.2) +
  geom_point(aes(x=epimonth, y=value, color=variable),  size=3, show.legend = F) +
  geom_line(aes(x=epimonth, y=value, color=variable),  size=1, show.legend = F) +
  facet_grid(variable~., scales = "free") + ylim(c(0, NA)) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     axis.title = element_blank(), 
                     strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"),
                     legend.text = element_text(size=12),
                     plot.margin = unit(c(.2,.1,1.3,1.1), "lines"),
                     axis.text = element_text(size=14))


# Look for cross-correlations
dat.lag <- cbind.data.frame(lag = print(ccf(merge.dat$sumppt, merge.dat$cases))$lag, acf=print(ccf(merge.dat$sumppt, merge.dat$cases))$acf)
dat.lag$variable <- "sumppt"
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)] #Sum ppt precedes cases by 11 months
dat2 = cbind.data.frame(lag = print(ccf(merge.dat$meant, merge.dat$cases))$lag, acf=print(ccf(merge.dat$meant, merge.dat$cases))$acf)
dat2$variable <- "meant"
dat2$lag[dat2$acf==max(dat2$acf)] #Mean temp precedes cases by 3 months
dat.lag <- rbind(dat.lag, dat2)

# Plot acf
# include the optimal lag on plot
max.lag <- dlply(dat.lag, .(variable))
get.lag <- function(df){
  lag = df$lag[df$acf==max(df$acf)]
  df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
  return(df.out)
}
max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))
max.lag$label = paste0("lag=", max.lag$lag, "epimonth")

ggplot(dat.lag) + geom_label(data=max.lag, aes(x=18,y=.4, label=label), label.size = 0) +
  geom_bar(aes(x=lag, y=acf), stat = "identity") + ylim(c(NA,.45)) +
  geom_hline(aes(yintercept=0.09), color="blue", linetype=2) +
  geom_hline(aes(yintercept=-0.09), color="blue", linetype=2) +
  facet_grid(variable~.) + theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                                              legend.title = element_blank(),
                                              axis.title = element_text(size=16),
                                              strip.background = element_rect(fill="white"),
                                              strip.text = element_text(size=14),
                                              legend.text = element_text(size=12),
                                              plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                                              axis.text = element_text(size=14))


# Make a dataframe with lagged climate variables
##Temp lagged 3 months, Precipitation lagged 11 months
merge.shift <- merge.dat[12: length(merge.dat$epimonth),]
merge.shift$precip_lag <- merge.dat$sumppt[1:(length(merge.dat$sumppt)-11)]
merge.shift$meantempLag <- merge.dat$meant[9:(length(merge.dat$meant)-3)]


############################
############################
#### Cara code starts here!!!

# Run glm with fixed effects
m1 <- glm(cases~ year.x + meantempLag + precip_lag, data = merge.shift, family = "poisson")
summary(m1)

# cases are positively associated with high temp, high precip, and increasing year
# we subsequently show that temp and precip are NOT increasing with time (see below). 
# therefore, the case increase is independent of climate.

# I don't think we should put this in the paper but I went ahead and tested the interacting effects
# of the climate data and the dates here
m2 <- glm(cases~ year.x*meantempLag*precip_lag, data = merge.shift, family = "poisson")
summary(m2)


# you can plot model outputs with this great package:
library(sjPlot)

# here for the first model:
plot_model(m1, type = "est") 
# this is the raw version of the plot I had you make down below
# all are above 0, meaning positive associations

plot_model(m1, type = "pred") 
# here are the predicted cases with all the model info, compared with the data

# and here for model 2:
plot_model(m2, type = "est", vline.color = "green")
# a lot more variables here - the red are below 0, 
# so we actually have this interacting effect of high precip 
# and advancing years being associated with fewer cases. 
# same with high temp combined with high precip.

plot_model(m2, type = "pred")
# predictions here, by each variable

plot_model(m2, type = "int")
# pattern here is closest in the last plot with the 3-way interaction.
# temperature trumps precip, so you can get cases where year advances,
# but temp is low and you have few cases despite high vs. low precip.
# at high temp, you always have more cases

# You can compare how many variables are needed in the model using the 'dredge'
# function in MuMIn package

library(MuMIn)

# set your global options - no missing data allowed
options(na.action ="na.fail")

# the 'global model' is the most complex. dredge compares AIC for all subsets
# of this model
out = dredge(global.model = m2)
head(out)

# look at top 10 AIC values
out[1:10,] #AIC values within 4 points of each other are largely equivalent

# The best model is the most complex (all predictors and all interactions)
# but I don't think it is really telling us much more than m1 (model 9 here),
# which is more interpretable. So I would recommend just reporting this as is.
# I mostly ran this to make sure we weren't missing anything.

# Now, to get an annual prediction by year from m1, we have to make a dummy dataset
# that has only one datapoint per month and includes the same columns we had in the 
# data the model was fit to originally. We'll use mean values for temp and precip 
# so these don't vary

# here, I set up the dummy dataset
# first, just one date value per month
#predict.df <- cbind.data.frame(epiyear.x = as.Date(paste0(2002:2020, "-07-01")))
predict.df <- cbind.data.frame(year.x = 2002:2020)

head(predict.df)


# now add mean values for the climate data
predict.df$meantempLag = mean(merge.shift$meantempLag)
predict.df$precip_lag = mean(merge.shift$precip_lag)


head(predict.df)


# Now, an annual prediction by year from m1, using above dummy dataset
# as the input term "newdata" into predict.glm
# Remember that the link function in a poisson model is always "log", 
# so we need to raise these predictions to the power e^x
predict.df$predicted_cases <- exp(predict.glm(m1, newdata = predict.df))

# We can also do the same for upper and lower confidence limits
predict.df$predicted_cases_lci <- exp(predict.glm(m1, newdata = predict.df, type = "link", se.fit = T)$fit -1.96*predict.glm(m1, newdata = predict.df, type = "link", se.fit = T)$se.fit)
predict.df$predicted_cases_uci <- exp(predict.glm(m1, newdata = predict.df, type = "link", se.fit = T)$fit +1.96*predict.glm(m1, newdata = predict.df, type = "link", se.fit = T)$se.fit)
  
# Alternatively--rather than raising to e^x-- we can specify 'type="response"' in the prediction function,
# which I just learned. R will do it automatically. Still, good to understand what is going on!

# Add a "date" column to ease with plotting - I picked the middle of each year
predict.df$date <-  as.Date(paste0(predict.df$year.x, "-07-01"))


# Now plot your predictions on top of cases. I added some styling to make the plot look nicer

p1 <- ggplot(merge.shift) + geom_line(aes(x=date, y=cases)) + 
      geom_ribbon(data=predict.df, aes(x=date, ymin=predicted_cases_lci, ymax=predicted_cases_uci), alpha=.6, fill="red") +
      geom_line(data=predict.df, aes(x=date, y=predicted_cases), color="red") +
      ylab("monthly dengue cases") +
      theme_bw() + 
      theme(axis.title.x = element_blank(), axis.text = element_text(size=16), axis.title.y = element_text(size=18), panel.grid = element_blank())

ggsave(file ="casePlot.png",
       plot=p1,
       units="mm",  
       width=60, 
       height=50, 
       scale=3, 
       dpi=300)
# Note that the CIs are quite narrow. You could make them a different color and get rid of the alph (makes them translucent)
# if you wanted. But I think it is fine as is...

# From m1 (above), we show that cases increase through time, independent of climate predictors. The line is a
# case prediction holding climate variables constant and changing year only

# To put the icing on the cake, in the supplement, you can add some work to show that temp and precip do not
# increase through time. I think the best way to do this is with GAMs. If they think they are too complicated,
# hopefully it won't matter in the supplementary material

# Here is that analysis:

# look at trends in temp and precip through time

head(ppt)
ppt$month[ppt$month<10] <- paste0("0",ppt$month[ppt$month<10])
ppt$date <- as.Date(paste0(ppt$year, "-", ppt$month, "-01"))
pA <- ggplot(data=ppt) + geom_line(aes(x=date, y=sumppt)) +theme_bw() + 
      ylab("total monthly precipitation") +  # add units to this if you have them
        theme(axis.title.x = element_blank(), axis.text = element_text(size=16), 
          axis.title.y = element_text(size=18), panel.grid = element_blank())
print(pA)


head(temp)
temp$month[temp$month<10] <- paste0("0",temp$month[temp$month<10])
temp$date <- as.Date(paste0(temp$year, "-", temp$month, "-01"))
pB <- ggplot(data=temp) + geom_line(aes(x=date, y=meant)) +
      ylab(bquote('mean monthly temperature ('^0*"C)")) + #""
      theme_bw() + 
      theme(axis.title.x = element_blank(), axis.text = element_text(size=16), 
            axis.title.y = element_text(size=18), panel.grid = element_blank())


# now run gam on the two
library(mgcv)
ppt$month <- as.numeric(ppt$month)
mA <- gam(sumppt ~ year + s(month, k=7, bs="cc"), data = ppt)
summary(mA) # no pattern through time. only seasonal variation intra-annually is significant

temp$month <- as.numeric(temp$month)
mB <- gam(meant ~ year + s(month, k=7, bs="cc"), data = temp)
summary(mB) # no pattern through time. only seasonal variation. only seasonal variation intra-annually is significant

# Include these plots together as a supplementary figure and report the GAM results

SuppPlot <- cowplot::plot_grid(pA, pB, ncol=2, labels = c("A", "B"))
print(SuppPlot)

ggsave(file ="suppPlot.png",
       plot=SuppPlot,
       units="mm",  
       width=90, 
       height=40, 
       scale=3, 
       dpi=300)

# You can orient the two plots vertically if you prefer.