rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
  
head(dat) 

#plot time series of each type by province by year

unique(dat$diagnostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$epiwk <- cut.Date(dat$date, breaks="weeks", start.on.monday = T)
dat$epiwk <- as.Date(as.character(dat$epiwk))
dat$epiwk[dat$epiwk< "2002-01-01"] <- "2002-01-01"
head(dat)

#now attach biwks (2 biweeks per week)
wks = 1:52
biwks <- rep(1:26, each=2)
wk.biwk <- cbind.data.frame(week_report=wks, biwk=biwks)

sort(unique(dat$week_report))
unique(dat$provname)

dat <- merge(dat, wk.biwk, by="week_report")
head(dat) #now plot as the lowest epiweek date per biweek
dat$case=1


#sum by epidemic biweek, type, province
dat.sum <- ddply(dat, .(year, biwk, provname, diagnostic), summarise, cases = sum(case), epiwk = min(epiwk))

head(dat.sum)

#plot infection type by province by year
p1 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=diagnostic)) +
      facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2018-01-01"), as.Date("2021-01-01")))
p1


p2 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=diagnostic)) +
  facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2011-01-01"), as.Date("2014-01-01")))
p2

p3 <- ggplot(data=dat.sum) + geom_line(aes(x=epiwk, y=cases, color=diagnostic)) +
  facet_wrap(~provname, ncol = 5, scales = "free_y") + coord_cartesian(xlim=c(as.Date("2006-01-01"), as.Date("2009-01-01")))
p3

#total cases by year by province
dat.prov <- ddply(dat, .(year, biwk, provname), summarise, cases = sum(case), epiwk=min(epiwk))
head(dat.prov)

#now fill in 0s for the missin biweeks
all.biwk.dat <- cbind.data.frame(year=rep(unique(dat.prov$year), each = (length(unique(dat.prov$biwk))*length(unique(dat$provname)))), biwk=rep(unique(dat.prov$biwk), (length(unique(dat.prov$year))*length(unique(dat$provname)))), provname = rep(unique(dat$provname), each=(length(unique(dat$year))*length(unique(dat$biwk)))))
head(all.biwk.dat)
all.biwk.dat <- all.biwk.dat[!duplicated(all.biwk.dat),]

dat.prov.cases <- dplyr::select(dat.prov, provname, year, biwk, cases)
epiweeks <- dplyr::select(dat.prov,  year, biwk, epiwk)
epiweeks <- epiweeks[!duplicated(epiweeks ),]
epiweeks <- ddply(epiweeks, .(year, biwk), summarise, epiweek = min(epiwk))
dat.prov.final <- merge(all.biwk.dat, dat.prov.cases, by = c("provname", "year", "biwk"), all.x = T, sort = F)
head(dat.prov.final)
dat.prov.final$cases[is.na(dat.prov.final$cases)] <- 0

#and attach epiweeks
dat.prov.final <- merge(dat.prov.final, epiweeks, by = c("year", "biwk"), all.x = T, sort = F)
head(dat.prov.final)
dat.prov.final <- arrange(dat.prov.final, epiweek, provname)

p1 <- ggplot(data=dat.prov) + geom_line(aes(x=epiwk, y=cases, color=provname))
p2 <- ggplot(data=dat.prov) + geom_line(aes(x=epiwk, y=cases, color=provname)) + 
      facet_wrap(provname~., ncol=5, scales = "free_y") + 
      geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
      geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
      geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)


head(dat.prov)
library(mgcv)

unique(dat.prov$biwk)
dat.prov$provname <- as.factor(dat.prov$provname)

#remove the one time series with very limited data
dat.sum <- ddply(dat.prov, .(provname), summarise, min_year = min(year), max_year=max(year))

dat.prov = subset(dat.prov, provname!="Tboung Khmum")

m1 <- gam(cases ~ year + s(year, by=provname, bs="tp", k=3) + 
            s(biwk, k=7, bs="cc") + s(provname, bs="re"),
            dat = dat.prov, family="poisson")
summary(m1)

out <- plot(m1)

unique(dat.prov$provname)
prov.list <- list()
for (i in 1:24){
  prov.list[[i]] <- out[[i]]
}
extract.slope.dat <- function(df){
  dat.out <- cbind.data.frame(year= df$x, slope= df$fit)
  dat.out$slope_uci = dat.out$slope + 1.96*(df$se)
  dat.out$slope_lci = dat.out$slope + 1.96*(df$se)
  dat.out$province <- sapply(strsplit(df$ylab, "provname"), '[',2)
return(dat.out)  
}

prov.list <- lapply(prov.list, extract.slope.dat)
prov.df <- data.table::rbindlist(prov.list)
head(prov.df)
tail(prov.df)


#supplemental plot of the time trends for each province
pS1 <- ggplot(prov.df) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  facet_wrap(~province, ncol=4) + ylab("partial effect on slope cases through time") + geom_hline(aes(yintercept=0))+
  geom_line(aes(x=year, y=slope, color=province)) + 
  geom_ribbon(aes(x=year, ymin=slope_lci, ymax=slope_uci, fill=province), alpha=.3)

prov.sum <- ddply(prov.df, .(province), summarise, median_slope=median(slope), upper_slope=max(slope), lower_slope=min(slope), mean_slope=mean(slope))
ggplot(prov.sum) + geom_point(aes(province, median_slope, color=province)) + geom_linerange(aes(x=province, ymin=lower_slope, ymax=upper_slope, color=province)) + geom_hline(yintercept=0)

pS1 <- ggplot(subset(prov.df, province=="Mondul Kiri")) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  facet_wrap(~province, ncol=4) + ylab("slope cases through time") + geom_hline(aes(yintercept=0))+
  geom_line(aes(x=year, y=slope, color=province)) + 
  geom_ribbon(aes(x=year, ymin=slope_lci, ymax=slope_uci, fill=province), alpha=.3)


dat.prov$projected_cases <- predict.gam(m1, exclude=c("s(biwk)"), type = "response")

pS2 <- ggplot(dat.prov) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  facet_wrap(~provname, ncol=4) + ylab("slope cases through time") +
  geom_line(aes(x=epiwk, y=cases, group=provname)) +
  geom_line(aes(x=epiwk, y=projected_cases, color=provname)) 
  


pS2 <- ggplot(subset(dat.prov, provname=="Mondul Kiri")) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  facet_wrap(~provname, ncol=4) + ylab("slope cases through time") +
  geom_line(aes(x=epiwk, y=cases, group=provname)) +
  geom_line(aes(x=epiwk, y=projected_cases, color=provname)) 



#get slope
diff(prov.df$slope)/diff(prov.df$year)

#then, collect the derivative, combine and take the average
library(gratia)

newDF=cbind.data.frame(year=prov.df$year)
newDF$province = "Phonm Penh"

out = derivatives(m1,
                  newdata = newDF,
                  order=1,
                  type="central",
                  eps = 1e-7,
                  interval = "confidence")

# out_df <- cbind.data.frame(bat_species=bat_spp, 
#                            measurement=measurement, 
#                            day_of_year=days_to_calc, 
#                            slope= out$derivative,
#                            lci = out$lower,
#                            uci=out$upper)


m1b <- gam(cases ~ s(year, by=provname, bs="tp", k=3) + 
            s(biwk, k=7, bs="cc"),
          dat = dat.prov, family="poisson")
summary(m1b)

AIC(m1, m1b)
#plot.gam(m1) #positive trend by year but some deviations
dat.prov$month = month(dat.prov$epiwk)

m1alt <- gam(cases ~ s(year, by=provname, k=3, bs="tp") + s(month, k=7, bs="cc"), dat = dat.prov, family="poisson")
summary(m1alt)
plot(m1alt)
#dat.prov$proj <- predict.gam(m1alt, exclude="s(week)")
dat.prov$proj <- predict.gam(m1alt, exclude="s(month)")

ggplot(dat.prov) + geom_line(aes(x=epiwk, y=proj, color=provname))

#load streicker script and look if there is predictable variation by province
source("/Users/carabrook/Developer/spillover-virulence/source/mollentze-streicker-2020-functions.R")


# look at partial effects
prov.df <- get_partial_effects(fit=m1, var="provname")

head(prov.df)

plot.partial(df = prov.df, var="provname", response_var = "cases") #huge variation by province
year.df <- get_partial_effects_continuous(gamFit=m1, var="year")
plot.partial.cont(df = year.df, var="year", response_var = "all df cases", log = F, alt_var = "year") #similar to all cases


inc.cases = subset(prov.df[[1]], y>0) 
#these are the capital and the central corridor with the highest pop density
#this area also has the highest temp and the highest flood risk (surrounds Tonle Sap and Mekong Delta)
dec.cases = subset(prov.df[[1]], y<0) #these are the exterior of the country

#and check the same for dhf and dss
head(dat)
dat.prov.dhf <- ddply(subset(dat, diagnostic=="dhf"), .(year, biwk, provname), summarise, cases = sum(case), epiwk=min(epiwk))
head(dat.prov.dhf)
#plot

p3 <- ggplot(data=dat.prov.dhf) + geom_line(aes(x=epiwk, y=cases, color=provname)) + 
  facet_wrap(provname~., ncol=5, scales = "free_y") + 
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

dat.prov.dhf$provname <- as.factor(dat.prov.dhf$provname)
dat.prov.dhf[!is.na(dat.prov$epiwk),]$month <- month(dat.prov$epiwk[!is.na(dat.prov$epiwk)])
m2 <- gam(cases ~ s(year, bs="tp", k=3) + s(biwk, k=7, bs="cc") + s(provname, bs="re"), dat = dat.prov.dhf, family="poisson")
summary(m2)
#plot.gam(m2) #positive trend by year but some deviations

#and 

prov.df.dhf <- get_partial_effects(fit=m2, var="provname")
year.df.dhf <- get_partial_effects_continuous(gamFit=m2, var="year")
plot.partial(df = prov.df.dhf, var="provname", response_var = "dhf cases") #similar to all cases
plot.partial.cont(df = year.df.dhf, var="year", response_var = "dhf cases", log = F, alt_var = "year") #similar to all cases

inc.cases.dhf = subset(prov.df.dhf[[1]], y>0) 
dec.cases.dhf = subset(prov.df.dhf[[1]], y<0) #these are the exterior of the country

#above - similar

#and dss
dat.prov.dss <- ddply(subset(dat, diagnostic=="dss"), .(year, biwk, provname), summarise, cases = sum(case), epiwk=min(epiwk))
head(dat.prov.dss)
#plot

#cases much rarer, outbreaks less clear
p4 <- ggplot(data=dat.prov.dss) + geom_line(aes(x=epiwk, y=cases, color=provname)) + 
  facet_wrap(provname~., ncol=5, scales = "free_y") + 
  geom_vline(xintercept = c(as.Date("2007-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2012-01-01")), color="red", linetype=2)+
  geom_vline(xintercept = c(as.Date("2019-01-01")), color="red", linetype=2)

dat.prov.dss$provname <- as.factor(dat.prov.dss$provname)
m3 <- gam(cases ~ s(year, k=3, bs="tp") + s(biwk, k=7, bs="cc") + s(provname, bs="re"), dat = dat.prov.dss, family="poisson")
summary(m3) #negative slope! 
#plot.gam(m3)

prov.df.dss <- get_partial_effects(fit=m3, var="provname")
year.df.dss <- get_partial_effects_continuous(gamFit=m3, var="year")

plot.partial(df = prov.df.dss, var="provname", response_var = "dss cases") 
plot.partial.cont(year.df.dss, var="year", response_var = "dss cases", log=FALSE, alt_var = "year")



#similar to all cases but fewer increasing. but still increasing in the most metropolitan areas

inc.cases.dss = subset(prov.df.dss[[1]], y>0) 
dec.cases.dss = subset(prov.df.dss[[1]], y<0) #these are the exterior of the country

dss.dat <- rbind(inc.cases.dss, dec.cases.dss)
dss.dat$diagnostic = "dss"

dhf.dat <- rbind(inc.cases.dhf, dec.cases.dhf)
dhf.dat$diagnostic = "dhf"

df.dat <- rbind(inc.cases, dec.cases)
df.dat$diagnostic = "df"

#now join 
case.dat <- rbind(df.dat, dhf.dat, dss.dat)

#and map these by slope through time
#and get submap of kampong speu witht the points of cambodia sequences
library(sf)
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)
unique(prov.sum$province)
prov.sum$province <- as.character(prov.sum$province)
setdiff(unique(prov.sum$province), unique(cam$ADM1_EN))

cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"
cam$ADM1_EN[cam$ADM1_EN=="Oddar Meanchey"] <- "Otdar Meanchey"

#and merge with the slopes
#dat.merge <- dplyr::select(case.dat, y, IsSignificant, provname, diagnostic)
dat.merge <- dplyr::select(prov.sum, province, median_slope, upper_slope, lower_slope, mean_slope)
names(dat.merge)[names(dat.merge)=="province"] <- "ADM1_EN"
cam_merge <- merge(cam, dat.merge, by = "ADM1_EN", all.x=T, sort=F)

#colz = c("No" = "gray", "Yes" = "black")
# get map
pCam <- ggplot(cam_merge) + geom_sf(aes(fill=mean_slope), size =.4) + 
   theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  #facet_grid(~diagnostic) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', name="mean partial\neffect on slope of cases\nthrough time") +
  #scale_color_manual(values = colz) +
  guides(color="none")


ggsave(file = paste0(homewd, "/test-figs/CaseTimeTrends.png"),
       plot= pCam,
       units="mm",  
       width=90, 
       height=55, 
       scale=3, 
       dpi=300)





#now do the same for the mean age of infection though time by province
dat.age <- ddply(dat, .(year, provname), summarise, age = mean(age))
head(dat.age)

dat.age$provname <- as.factor(dat.age$provname)
m4 <- gam(age ~ s(year, k=3, bs="tp") + s(provname, bs="re"), dat = dat.age, family="gaussian")
summary(m4) #positive slope! 

prov.age <- get_partial_effects(fit=m4, var="provname")
year.age <- get_partial_effects_continuous(gamFit=m4, var="year")
plot.partial(df = prov.age, var="provname", response_var = "mean age of infection") #similar to all cases
plot.partial.cont(df = year.age, var="year", response_var = "mean age of infection", log = F, alt_var = "year") #similar to all cases

age.dat = subset(prov.age[[1]]) 

#and plot
dat.merge.age <- dplyr::select(age.dat, y, IsSignificant, provname)
names(dat.merge.age)[names(dat.merge.age)=="provname"] <- "ADM1_EN"
cam_merge_age <- merge(cam, dat.merge.age, by = "ADM1_EN", all.x=T, sort=F)

colz = c("No" = "gray", "Yes" = "black")
# get map
pCam_age <- ggplot(cam_merge_age) + geom_sf(aes(fill=y, color=IsSignificant), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', name="slope\nthrough\ntime") +
  scale_color_manual(values = colz) +
  guides(color="none")


ggsave(file = paste0(homewd, "/test-figs/CamAgeTrends.png"),
       plot= pCam_age,
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)


#and plot the mean age of infection across the whole dataset
mean.age.df <- ddply(dat, .(provname), summarise, age = mean(age))
names(mean.age.df) <- c("ADM1_EN", "avg_age")

#save
save(mean.age.df, file = paste0(homewd, "/data/meanagebyprov.Rdata"))
cam_merge_age <- merge(cam_merge_age, mean.age.df, by = "ADM1_EN", all.x=T, sort=F)

pCam_avg_age <- ggplot(cam_merge_age) + geom_sf(aes(fill=avg_age), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_steps(low = 'white', high = 'firebrick', name="mean\nage of\ninfection") 


#and R0 as the per capita birth rate/mean age of infection

#get births per 1000 people. multiple to get per total pop size in cambodia
#then divide by province to get per capita

#353829 = mean number of births per year for cambodia. now get this per capita by province


pop.dat <- read.csv(file=paste0(homewd, "/data/world_bank_cambodia.csv"), header = T, stringsAsFactors = F)
head(pop.dat)

#get population vector
pop.vec <- pop.dat[2,5:ncol(pop.dat)]
names(pop.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
#so one timestep before gives you the population at the beginning of the year
pop.vec <- c(unlist(pop.vec[which(names(pop.vec)=="2001"):which(names(pop.vec)=="2020")]))
#pop.vec[length(pop.vec)] <- pop.vec[length(pop.vec)-1]

#do the same for births - these are births per 1000 people
#get total births
birth.vec <- pop.dat[1,5:ncol(pop.dat)]
names(birth.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="2002"):which(names(birth.vec)=="2020")]
birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year

#this is births per 1000. now take births per capita and scale by the country proportion

#load the proportion by province
load(paste0(homewd, "/data/cambodia_province_proportions.Rdata"))
prop.prov$tot_pop <- mean(pop.vec)*prop.prov$pop_prop
prop.prov$births_per_province <- mean(c(unlist(birth.vec/1000)))*prop.prov$tot_pop
prop.prov$births_per_capita <- prop.prov$births_per_province/prop.prov$tot_pop

#merge mewan age
setdiff(prop.prov$provname, mean.age.df$ADM1_EN)
setdiff(mean.age.df$ADM1_EN,prop.prov$provname)
names(prop.prov)[names(prop.prov)=="provname"] <- "ADM1_EN"
birth.merge <- dplyr::select(prop.prov, ADM1_EN, births_per_capita, tot_pop)

mean.age.df <- merge(mean.age.df, birth.merge, by = "ADM1_EN")
#now get R0
mean.age.df$R0 <- (1/mean.age.df$births_per_capita)/mean.age.df$avg_age

setdiff(mean.age.df$ADM1_EN, cam_merge_age$ADM1_EN)
mean.age.df

cam_merge_R0 <- merge(cam_merge_age, mean.age.df, by = "ADM1_EN", all.x=T, sort=F)

#THIS IS A KEEPER
pCam_R0 <- ggplot(cam_merge_R0) + geom_sf(aes(fill=R0), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_steps(low = 'white', high = 'firebrick', name="R0", limits=c(2.1,6)) 



#and total opop

pCam_pop <- ggplot(cam_merge_R0) + geom_sf(aes(fill=tot_pop), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_steps(low = 'white', high = 'firebrick', name="pop\nsize", limits=c(30000,1500000)) 








#and correlate with temp or precip
climate.dat <- read.csv(file = paste0(homewd, "/data/climate_mean_dat.csv"), header = T, stringsAsFactors = F)
head(climate.dat)
setdiff(climate.dat$adm1_name, mean.age.df$ADM1_EN)
climate.dat$adm1_name[climate.dat$adm1_name=="Oddar Meanchey"] <- "Otdar Meanchey"
climate.dat$adm1_name[climate.dat$adm1_name=="Siemreap"] <- "Siem Reap"

names(climate.dat)[1] <- "ADM1_EN"

climate.dat <- merge(climate.dat, mean.age.df, by="ADM1_EN", all.x = T, sort = F)
head(climate.dat)

pTempR0 <- ggplot(data=climate.dat) + scale_color_continuous(low="white", high="navy") +
           geom_point(aes(x=mean_temp_C, y=R0, color=mean_precip_mm), size=10) + theme_bw() +
            geom_label(aes(x=mean_temp_C, y=R0, label=ADM1_EN), size=2)
pPrecipR0 <- ggplot(data=climate.dat) + theme_bw() +
            geom_point(aes(x=mean_precip_mm, y=R0, color=mean_temp_C), size=10) +
            geom_label(aes(x=mean_precip_mm, y=R0, label=ADM1_EN), size=2) +
            scale_color_continuous(low="skyblue", high="firebrick") 
            
ClimRegress <- pTempR0|pPrecipR0


ggsave(file = paste0(homewd, "/test-figs/ClimateRegressions.png"),
       plot= ClimRegress,
       units="mm",  
       width=100, 
       height=55, 
       scale=3, 
       dpi=300)

#and map mean temp and mean precip

cam_merge_climate <- merge(cam_merge_age, climate.dat, by = "ADM1_EN", all.x=T, sort=F)

pCam_temp <- ggplot(cam_merge_climate) + geom_sf(aes(fill=mean_temp_C), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_steps(low = 'yellow', high = 'firebrick', name="mean\nannual\ntemp (C)", limits=c(25,29))

#and precip
pCam_precip <- ggplot(cam_merge_climate) + geom_sf(aes(fill=mean_precip_mm), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_steps(low = 'white', high = 'navy', name="mean\nannual\nprecip (mm)", limits=c(3,8.5))

#and plot all together

pClimMaps <- pCam_precip | pCam_temp | pCam_R0



ggsave(file = paste0(homewd, "/test-figs/ClimateMaps.png"),
       plot= pClimMaps,
       units="mm",  
       width=88, 
       height=55, 
       scale=3, 
       dpi=300)




#add in population
head(cam_merge_climate)
head(prop.prov)
prop.merge <- dplyr::select(prop.prov, ADM1_EN, tot_pop)
head(prop.merge)

cam_merge_all <- merge(cam_merge_climate, prop.merge, by="ADM1_EN")
head(cam_merge_all)


p1 <- ggplot(data=cam_merge_all) + geom_point(aes(x=tot_pop, y=R0, color=mean_temp_C), size=3) +
      scale_color_steps(low = 'yellow', high = 'firebrick', name="mean\nannual\ntemp (C)", limits=c(25,29))
