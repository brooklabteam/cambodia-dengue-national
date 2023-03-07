rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

dat <- read.csv(file = paste0(homewd, "/data/ageadjust_month.csv"), header = T,stringsAsFactors = F)
head(dat)

#add year, month, and date column
#get year
dat$year <- sapply(strsplit(dat$year_month, "_"), "[", 1)
#get month
dat$month <- sapply(strsplit(dat$year_month, "_"), "[", 2)
#add 0 to the single digit months:
dat$month[dat$month!="10" & dat$month!="11" & dat$month!="12"] <- paste0("0", dat$month[dat$month!="10" & dat$month!="11" & dat$month!="12"])
#make a date column
dat$date <- as.Date(paste0(dat$year, "-", dat$month, "-01"))


#plot 

p1 <- ggplot(data=dat) + theme_bw() +
  geom_line(aes(x=date, y= crude.case.numbers))

print(p1)

p2 <- ggplot(data=dat) + theme_bw() +
  geom_line(aes(x=date, y= age.adjusted.incidence))

print(p2)


# Research Questions:
# (1) Do crude cases numbers change with time?

# make month and year numeric
dat$month <- as.numeric(dat$month)
dat$year <- as.numeric(dat$year)

#add the epi year
dat$epi_year <- 0
dat$epi_year[dat$year==2007 | dat$year==2012 | dat$year==2019] <- 1
dat$epi_year <- as.factor(dat$epi_year)

#with no epi year effect
# spline bs= 'cc' means controls such that january and december are contiguous
m1 <- gam(crude.case.numbers ~ year + s(month, k=7, bs="cc"),
                               family="poisson",data = dat)

summary(m1)
# you do have a significant (slightly) positive slope  of year

# #save this output
# sink(file="case_gam.txt")
# summary(m1)
# sink()


#and here is the internal within-year curve
#you can also save this as an object and extract the data to make a prettier plot
plot(m1) 

#you can plot the trend by year through the "predict" function, dropping
#the effects of month

dat$predict_gam <- predict.gam(m1, type="response", exclude = "s(month)")
dat$predict_gam_lci <- dat$predict_gam - 1.96*predict.gam(m1, type="response", exclude = "s(month)", se.fit = T)$se
dat$predict_gam_uci <- dat$predict_gam + 1.96*predict.gam(m1, type="response", exclude = "s(month)", se.fit = T)$se

# and plot with your plot from before:

p1 <- ggplot(data=dat) + theme_bw() + ylab("crude case numbers") +
  geom_line(aes(x=date, y= crude.case.numbers))

print(p1)

# you are plotting by date on the x-axis, and  this predicts the same value for all 
# dates in a year (because we excluding month), so it looks like a stairstep
# (it's faint because of the scale)

p1b <- ggplot(data=dat) + theme_bw() +
  geom_line(aes(x=date, y= crude.case.numbers)) + ylab("crude case numbers") +
  geom_line(aes(x=date, y=predict_gam), color="red") +
  geom_ribbon(aes(x=date, ymin=predict_gam_lci, ymax=predict_gam_uci), alpha=.3, fill="red")
  
  
print(p1b)

# You could also just make a separate dataframe for the annual predictions

dat.annual <- dplyr::select(dat, year, month, date, predict_gam, predict_gam_lci, predict_gam_uci)
dat.annual <- subset(dat.annual, month ==1)

p1b <- ggplot(data=dat) + theme_bw() +
  geom_line(aes(x=date, y= crude.case.numbers)) + ylab("crude case numbers") +
  geom_line(data=dat.annual, aes(x=date, y=predict_gam), color="red") +
  geom_ribbon(data=dat.annual, aes(x=date, ymin=predict_gam_lci, ymax=predict_gam_uci), 
              alpha=.3, fill="red")

print(p1b)

#or plot just the trend by year on its own:

p1c <- ggplot(data=dat.annual) + theme_bw() + ylab("predicted case\ntrend by year")+
  geom_line(aes(x=date, y=predict_gam), color="red") +
  geom_ribbon( aes(x=date, ymin=predict_gam_lci, ymax=predict_gam_uci), 
              alpha=.3, fill="red")

print(p1c)


#or plot side-by-side

CaseCompiled <- cowplot::plot_grid(p1,p1b,p1c, ncol=1, nrow=3, align = "hv")

#and you can save with this:

ggsave(file ="crude_case_plot.pdf",
      plot=CaseCompiled,
      units="mm",  
      width=60, 
      height=100, 
      scale=3, 
      dpi=300)
        


# I believe Yimei included a random effect of epidemic year
# to allow for a different slope within that year. I don't think
# I would do that, though I tried it here

m1re <- gam(crude.case.numbers ~ year + s(month, k=7, bs="cc") + s(epi_year, bs="re"),
          family="poisson",data = dat)
summary(m1re) 

#random effect is significant but the positive slope is still significant
#this means that the positive trend by year persists despite "removing" the 
#epi years from the data. I would not report this personally. It's weird.



# The script for the age incidence will be largely similar:

# (2) Does age incidence change with time?

m2 <- gam(age.adjusted.incidence ~ year + s(month, k=7, bs="cc"),
          family="gaussian",data = dat)

# You'll need to report all these details in an article supplement, 
# similar to what I did in my Biological Conservation paper.

summary(m2) 
#still a slight positive slope, but less significant than above
# 
# #save
# sink(file="incidence_gam.txt")
# summary(m2)
# sink()


#largely the same internal curve though with wider CIs
#you can also save this as an object and extract the data to make a prettier plot
plot(m2) 


# and you can do the same as above to extract the predictions by year
# in fact, I'll just add to the original dataset
head(dat)

#rename gam prediction columns
names(dat)[8:10]  <- c("predict_gam_cases", "predict_gam_cases_lci", "predict_gam_cases_uci")


#and add in the new predictions
dat$predict_gam_ageincidence <- predict.gam(m2, type="response", exclude = "s(month)")
dat$predict_gam_ageincidence_lci <- dat$predict_gam_ageincidence - 1.96*predict.gam(m2, type="response", exclude = "s(month)", se.fit = T)$se
dat$predict_gam_ageincidence_uci <- dat$predict_gam_ageincidence + 1.96*predict.gam(m2, type="response", exclude = "s(month)", se.fit = T)$se


#and redo the annual df
# You could also just make a separate dataframe for the annual predictions
dat.annual <- dplyr::select(dat, year, month, date, predict_gam_cases, predict_gam_cases_lci, predict_gam_cases_uci, predict_gam_ageincidence, predict_gam_ageincidence_lci, predict_gam_ageincidence_uci)
dat.annual <- subset(dat.annual, month ==1)


#and plot
p2 <- ggplot(data=dat) + theme_bw() + ylab("age adjusted incidence") +
  geom_line(aes(x=date, y= age.adjusted.incidence))

print(p2)


p2b <- ggplot(data=dat) + theme_bw() + ylab("age adjusted incidence") +
  geom_line(aes(x=date, y= age.adjusted.incidence)) +
  geom_line(data=dat.annual, aes(x=date, y=predict_gam_ageincidence), color="red") +
  geom_ribbon(data=dat.annual, aes(x=date, ymin=predict_gam_ageincidence_lci, ymax=predict_gam_ageincidence_uci), 
              alpha=.3, fill="red")

print(p2b)



p2c <- ggplot(data=dat.annual) + theme_bw() + ylab("predicted age adjusted\nincidence by year")+
  geom_line(aes(x=date, y=predict_gam_ageincidence), color="red") +
  geom_ribbon( aes(x=date, ymin=predict_gam_ageincidence_lci, ymax=predict_gam_ageincidence_uci), 
               alpha=.3, fill="red")

print(p2c)


IncidenceCompiled <- cowplot::plot_grid(p2,p2b,p2c, ncol=1, nrow=3, align = "hv")

#and you can save with this:

ggsave(file ="incidence_plot.pdf",
       plot=IncidenceCompiled,
       units="mm",  
       width=60, 
       height=100, 
       scale=3, 
       dpi=300)
