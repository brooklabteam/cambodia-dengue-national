rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(WaveletComp)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)


dat <- read.csv(file = paste0(homewd, "/data/all_case_climate.csv"), header=T, stringsAsFactors = F)
  
head(dat) 
dat$month_date <- as.Date(dat$month_date)

ggplot(dat) + geom_line(aes(x=month_date, y=cases, group=provname)) + 
  facet_wrap(~provname, ncol=5) +
  geom_line(aes(x=month_date, y=temp_C*100, group=provname), color="red")

#subset to national data only
dat.nat <- dplyr::select(dat, month_date, national_cases, temp_C_national, precip_mm_national, oni, mei)
names(dat.nat) <- c("month_date", "cases", "temp_C", "precip_mm", "oni", "mei")
head(dat.nat)
dat.nat <- ddply(dat.nat, .(month_date), summarise, cases=unique(cases), temp_C=unique(temp_C), precip_mm=unique(precip_mm), oni=unique(oni), mei=unique(mei))

ggplot(dat.nat) + geom_point(aes(month_date, cases)) + geom_line(aes(month_date, cases)) +
  geom_line(aes(x=month_date, y=temp_C*100), color="red")


corr.temp <- analyze.coherency(dat.nat, my.pair = c("temp_C", "cases"),
                           loess.span = 0,
                           dt = 1/12, dj = 1/100,
                           lowerPeriod = 2, #shortes possible period in years
                           make.pval = TRUE, n.sim = 100)

wc.image(corr.temp, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)", plot.ridge = T,
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#synchrony intensifies from 2008 onwards (again that 3-4 year period!)
#arrows to the right signify that these are in phase. 
#arrows up signify that x is leading
#arrows down signify that y is leading...
#temp leads cases until 219 when they largely converge
corr.temp$Angle #give lead (in years of x over y per location)

wc.avg(corr.temp, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)")

wc.sel.phases(corr.temp,    spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))

corr.precip <- analyze.coherency(dat.nat, my.pair = c("precip_mm","cases"),
                           loess.span = 0,
                           dt = 1/12, dj = 1/100,
                           lowerPeriod = 2, #shortes possible period in years
                           make.pval = TRUE, n.sim = 100)

wc.avg(corr.precip, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #very sig, around 3

wc.image(corr.precip, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#arrows to the left suggest that these are out of phase... cases follow precip
wc.sel.phases(corr.precip, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                              labels = 2002:2019))



#what about temp and precip
#3-4 years just later on
precip.temp <- analyze.coherency(dat.nat, my.pair = c("temp_C","precip_mm"),
                                   loess.span = 0,
                                   dt = 1/12, dj = 1/100,
                                   lowerPeriod = 2, #shortes possible period in years
                                   make.pval = TRUE, n.sim = 100)


wc.image(precip.temp, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))
#only coherency in precip and temp later in the time series
#out of phase

wc.avg(precip.temp, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #not significant

wc.sel.phases(precip.temp, spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                                                   labels = 2002:2019))


#now check national with oni and mei
corr.mei <- analyze.coherency(dat.nat, my.pair = c("mei","cases"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 lowerPeriod = 2, #shortes possible period in years
                                 make.pval = TRUE, n.sim = 100)

wc.image(corr.mei, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#strong synchrony throughout time series: 2006-2019, in phase, with mei leading cases
wc.avg(corr.mei, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #not significant

wc.sel.phases(corr.mei, spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))



#now check national with  oni
corr.oni <- analyze.coherency(dat.nat, my.pair = c("oni","cases"),
                              loess.span = 0,
                              dt = 1/12, dj = 1/100,
                              lowerPeriod = 2, #shortes possible period in years
                              make.pval = TRUE, n.sim = 100)

wc.image(corr.oni, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
         periodlab= "period (in years)",
         spec.time.axis = list(at = seq(1,12*18, 12), #converting weeks to years
                               labels = 2002:2019))

#strong synchrony throughout time series: 2006-2019, in phase, with oni leading cases
#by how much?
wc.avg(corr.oni, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (years)") #not significant

wc.sel.phases(corr.oni, spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
#basically the same with oni as with mei

ggplot(dat.nat)+ geom_hline(aes(yintercept=0)) +
        geom_vline(aes(xintercept=as.Date("2007-01-01"))) + geom_vline(aes(xintercept=as.Date("2012-01-01"))) + geom_vline(aes(xintercept=as.Date("2019-01-01"))) +
        geom_point(data=subset(dat.nat, oni>0), aes(x=month_date, y=oni), color="red") + 
        geom_line(data=subset(dat.nat, oni>0), aes(x=month_date, y=oni), color="red") +
        geom_point(data=subset(dat.nat, oni<0), aes(x=month_date, y=oni), color="blue") + 
        geom_line(data=subset(dat.nat, oni<0), aes(x=month_date, y=oni), color="blue") 
        

# is temp increasing with time?

head(dat)

ggplot(data=dat) + geom_point(aes(x=month_date, y=temp_C, color=provname)) +
              geom_line(aes(x=month_date, y=temp_C, color=provname))

dat$year <- year(dat$month_date)
dat$month <- month(dat$month_date)
dat$week <- week(dat$month_date)
dat$provname <- as.factor(dat$provname)

#and gam
m1 <- gam(temp_C ~ s(year, by=provname, bs="tp", k=3) + 
                    s(month, k=7, bs="cc"),
          dat = dat)
summary(m1)

#look at slope by province
gam.dat <- plot(m1)#first 25 elements are the info for each province. #26 is the seasonal smoother
plot(gam.dat[[1]]$x, gam.dat[[1]]$fit)
plot(gam.dat[[2]]$x, gam.dat[[2]]$fit)
plot(gam.dat[[3]]$x, gam.dat[[3]]$fit)
plot(gam.dat[[4]]$x, gam.dat[[4]]$fit)
plot(gam.dat[[25]]$x, gam.dat[[25]]$fit)


# #now plot
# 
# #load streicker script and look if there is predictable variation by province
# source("/Users/carabrook/Developer/spillover-virulence/source/mollentze-streicker-2020-functions.R")
# 
# 
# # look at partial effects
# prov.df <- get_partial_effects(fit=m1, var="provname")
# 
# head(prov.df)
# 
# plot.partial(df = prov.df, var="provname", response_var = "temp_C") #huge variation by province
# year.df <- get_partial_effects_continuous(gamFit=m1, var="year")
# plot.partial.cont(df = year.df, var="year", response_var = "temp (C)", log = F, alt_var = "year") #similar to all cases
# 
# prov.temp<- prov.df$effects
# 
# inc.temp = subset(prov.df[[1]], y>0) 
# #these are the capital and the central corridor with the highest pop density
# #this area also has the highest temp and the highest flood risk (surrounds Tonle Sap and Mekong Delta)
# dec.temp = subset(prov.df[[1]], y<0) #these are the exterior of the country


#and map these by slope through time
#and get submap of kampong speu witht the points of cambodia sequences
library(sf)
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)
unique(prov.temp$provname)
prov.temp$provname <- as.character(prov.temp$provname)
setdiff(unique(prov.temp$provname), unique(cam$ADM1_EN))

cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"

#and merge with the slopes
dat.merge <- dplyr::select(prov.temp, provname, y, IsSignificant, provname)
names(dat.merge)[names(dat.merge)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam, dat.merge, by = "ADM1_EN", all.x=T, sort=F)

colz = c("No" = "gray", "Yes" = "black")
# get map
pCam <- ggplot(cam_merge) + geom_sf(aes(fill=y, color=IsSignificant), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', name="partial effect\non slope\nthrough time") +
  scale_color_manual(values = colz) +
  guides(color="none") +ggtitle("Rising temperatures across Cambodia, 2002-2020")

#and plot the multinational time series for temp

dat$proj_temp <- predict.gam(m1, exclude = "s(month)")
#dat$proj_temp <- predict.gam(m1)


ggplot(data=dat) + geom_point(aes(x=month_date, y=temp_C, color=provname)) +
  geom_line(aes(x=month_date, y=temp_C, color=provname)) + 
  geom_line(aes(x=month_date, y=proj_temp), color="red") + 
  facet_wrap(~provname, ncol=5)

#and plot all the slopes only

ggplot(data=dat) +
  geom_line(aes(x=month_date, y=proj_temp, color=provname)) 
