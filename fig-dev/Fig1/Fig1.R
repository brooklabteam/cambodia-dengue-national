rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#plot time series and time trend by province

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
head(dat)
hist(dat$age)
sum.dat <- ddply(dat, .(year, age), summarise, N=length(age))
subset(sum.dat, age>50)

p1 <- ggplot(data=dat) + geom_violin(aes(x=year, y=age, group=year), 
            draw_quantiles = c(0,.25,.5,.75)) +
  geom_jitter(aes(x=year, y=age), width=.1, size=.1, alpha=.2) 
  
p1

#and the max age of dengue infection
age.max <- ddply(dat, .(year), summarise, max_age = max(age))
ggplot(age.max) + geom_point(aes(x=year, y = max_age)) + geom_line(aes(x=year, y=max_age))

popdat <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dat.csv"), header=T, stringsAsFactors = F)

#sum by epidemic biweek, type, province
dat.prov <- ddply(dat, .(year, epimonth, provname, diagnostic), summarise, cases = sum(case))
dat.nat <- ddply(dat, .(year, epimonth, diagnostic), summarise, cases = sum(case))
dat.nat$epimonth <- as.Date(dat.nat$epimonth, format = "%m/%d/%y")
dat.prov$epimonth <- as.Date(dat.prov$epimonth, format = "%m/%d/%y")


#add pop
popmerge <- dplyr::select(popdat, year, pop, provname)
dat.prov <- merge(dat.prov, popdat, by=c("year", "provname"), all.x = T)
head(dat.prov)
popnat <- ddply(popmerge, .(year), summarise, pop=sum(pop))

dat.nat <- merge(dat.nat, popnat, by=c("year"), all.x = T)
head(dat.nat)

dat.nat$cases_per_1000 <- (dat.nat$cases/dat.nat$pop)*1000

#ggplot(data=dat.nat) +geom_line(aes(x=epimonth, y=pop))
#ggplot(data=dat.prov) +geom_line(aes(x=epimonth, y=pop, color=provname))

#figure 1 A is the time-trend of this at the national level
colz <- c("df" = "navy", dhf="maroon", dss="darkseagreen4")
Fig1A <- ggplot(data=dat.nat) + theme_bw()+
        ylab("national cases\nper 1000 people") + scale_color_manual(values=colz) +
        theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
              legend.title = element_blank(), legend.position = c(.15,.73),
              axis.text = element_text(size=14), axis.title.y = element_text(size=16),
              legend.text = element_text(size=12),
              plot.margin = unit(c(4.5,1,4.5,1.5),"lines")) +
        geom_line(aes(x=epimonth, y=cases_per_1000, color=diagnostic), linewidth=.8)

#and fit GAM to time trend by province
library(lubridate)

dat.prov$month <- month(dat.prov$epimonth)


#need to eliminate one province with short time series
dat.prov = subset(dat.prov, provname!="Tboung Khmum")
dat.prov$provname <- as.factor(dat.prov$provname)
dat.prov$cases_per_1000 <- (dat.prov$cases/dat.prov$pop)*1000

m1 <- gam(cases_per_1000 ~ year:provname + 
            #s(year, by=provname, bs="tp", k=3) + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
           dat = subset(dat.prov, diagnostic=="df"))
summary(m1) #all increasing


out=summary(m1)

prov.df <- cbind.data.frame(provname = unique(dat.prov$provname), 
                            slope= out$p.coeff[2:length(out$p.coeff)],
                            se = out$se[2:length(out$p.coeff)],
                            p_val = out$p.pv[2:length(out$p.coeff)])

prov.df$slope_uci <- prov.df$slope + 1.96*prov.df$se
prov.df$slope_lci <- prov.df$slope - 1.96*prov.df$se
prov.df$sig <- "sig"
prov.df$sig[prov.df$p_val>0.05] <- "not_sig"
rownames(prov.df) <- c()
prov.df$diagnostic = "df"

#and predict
predict.dat <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat$month <- 1

predict.dat$predict_cases <- predict.gam(m1,newdata = predict.dat,  exclude = "s(month)", type = "response")
predict.dat$predict_cases_lci  <- predict.dat$predict_cases - 1.96*predict.gam(m1, newdata = predict.dat,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat$predict_cases_uci  <- predict.dat$predict_cases + 1.96*predict.gam(m1, newdata = predict.dat,  exclude = "s(month)", type = "response", se.fit = T)$se


predict.dat$diagnostic = "df"



#and repeat for dhf and dss

m2 <- gam(cases_per_1000 ~ year:provname + 
            #s(year, by=provname, bs="tp", k=3) + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = subset(dat.prov, diagnostic=="dhf"))
summary(m2) #more sig. most pos


out=summary(m2)

prov.df.2 <- cbind.data.frame(provname = unique(dat.prov$provname), 
                            slope= out$p.coeff[2:length(out$p.coeff)],
                            se = out$se[2:length(out$p.coeff)],
                            p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.2$slope_uci <- prov.df.2$slope + 1.96*prov.df.2$se
prov.df.2$slope_lci <- prov.df.2$slope - 1.96*prov.df.2$se
prov.df.2$sig <- "sig"
prov.df.2$sig[prov.df.2$p_val>0.05] <- "not_sig"
rownames(prov.df.2) <- c()

prov.df.2$diagnostic <- "dhf"

#and predict
predict.dat.2 <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat.2$month <- 1

predict.dat.2$predict_cases <- predict.gam(m2,newdata = predict.dat.2,  exclude = "s(month)", type = "response")
predict.dat.2$predict_cases_lci  <- predict.dat.2$predict_cases - 1.96*predict.gam(m2, newdata = predict.dat.2,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat.2$predict_cases_uci  <- predict.dat.2$predict_cases + 1.96*predict.gam(m2, newdata = predict.dat.2,  exclude = "s(month)", type = "response", se.fit = T)$se


predict.dat.2$diagnostic = "dhf"


#and dss

m3 <- gam(cases_per_1000 ~ year:provname + 
            #s(year, by=provname, bs="tp", k=3) + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = subset(dat.prov, diagnostic=="dss"))
summary(m3) #more sig. most pos


out=summary(m3)

prov.df.3 <- cbind.data.frame(provname = unique(dat.prov$provname), 
                              slope= out$p.coeff[2:length(out$p.coeff)],
                              se = out$se[2:length(out$p.coeff)],
                              p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.3$slope_uci <- prov.df.3$slope + 1.96*prov.df.3$se
prov.df.3$slope_lci <- prov.df.3$slope - 1.96*prov.df.3$se
prov.df.3$sig <- "sig"
prov.df.3$sig[prov.df.3$p_val>0.05] <- "not_sig"
rownames(prov.df.3) <- c()

prov.df.3$diagnostic <- "dss"

#and predict
predict.dat.3 <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat.3$month <- 1

predict.dat.3$predict_cases <- predict.gam(m3,newdata = predict.dat.3,  exclude = "s(month)", type = "response")
predict.dat.3$predict_cases_lci  <- predict.dat.3$predict_cases - 1.96*predict.gam(m3, newdata = predict.dat.3,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat.3$predict_cases_uci  <- predict.dat.3$predict_cases + 1.96*predict.gam(m3, newdata = predict.dat.3,  exclude = "s(month)", type = "response", se.fit = T)$se

predict.dat.3$diagnostic = "dss"


#and combine
prov.df <- rbind(prov.df,prov.df.2, prov.df.3)
predict.dat <- rbind(predict.dat, predict.dat.2, predict.dat.3)

predict.dat$epiyear <- as.Date(paste0(predict.dat$year, "-01-01"))


#and supplemental figure
FigS1 <- ggplot(predict.dat) + scale_color_manual(values=colz) +scale_fill_manual(values=colz) +
          geom_line(data= subset(dat.prov, diagnostic=="df"), aes(x=epimonth, y=cases_per_1000), color="gray70") +
          geom_line(aes(x=epiyear, y=predict_cases, color=diagnostic)) + theme_bw() +
          geom_ribbon(aes(x=epiyear, ymin=predict_cases_lci, ymax=predict_cases_uci, fill=diagnostic), alpha=.3) + 
          facet_wrap(~provname, ncol=4) + ylab("predicted annual cases per 1000 ppl") +
          theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
                strip.background = element_rect(fill="white"),
          legend.title = element_blank(), legend.position = "bottom",
          axis.text = element_text(size=14), axis.title.y = element_text(size=16),
          legend.text = element_text(size=12)) + coord_cartesian(ylim=c(0,.5))
      

ggsave(paste0(homewd, "/fig-new/FigS1.png"),
       plot = FigS1,
       bg='white',
       width=90, 
       height=80, 
       scale=3.2, 
       dpi=300,
       units = "mm")



### and then Fig1B is the map with the slopes by province
library(sf)
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)

prov.df$provname <- as.character(prov.df$provname)

setdiff(unique(prov.df$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(prov.df$provname)) #"Tboung Khmum"

cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"
cam$ADM1_EN[cam$ADM1_EN=="Oddar Meanchey"] <- "Otdar Meanchey"

#and add the blanks
head(prov.df)
prov.add <- subset(prov.df, provname=="Kandal")
prov.add$provname <- "Tboung Khmum"
prov.add$sig <- "not_sig"
prov.add$slope <- prov.add$se <- prov.add$p_val <- prov.add$slope_lci <- prov.add$slope_uci <- NA

prov.df <- rbind(prov.df, prov.add)
#and merge with the slopes
#dat.merge <- dplyr::select(case.dat, y, IsSignificant, provname, diagnostic)

names(prov.df)[names(prov.df)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam, prov.df, by = "ADM1_EN", all.x=T, sort=F)

cam_sub = subset(cam_merge, is.na(slope))
head(cam_merge)

unique(prov.df$sig)
colorz <- c('sig' = "black", 'not_sig'="grey80")

#library(patchwork)

Fig1B <- ggplot(cam_merge) + geom_sf(aes(fill=slope, color=sig), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,.5,0,1.5),"lines")) +
  facet_grid(~diagnostic) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', 
                       na.value = "grey50",
                       name="slope of cases per 1000\npeople through time") +
  scale_color_manual(values = colorz) +
  guides(color="none")



Fig1 <- cowplot::plot_grid(Fig1A, Fig1B, ncol=2, rel_widths = c(1,1.6), labels = c("A", "B"), label_size = 22, align = "hv", axis="tb", label_x = -.01, label_y = .85)


ggsave(paste0(homewd, "/fig-new/Fig1.png"),
       plot = Fig1,
       bg='white',
       width=100, 
       height=40, 
       scale=3.2, 
       dpi=300,
       units = "mm")


#also try mean age of infection

#now do the same for the mean age of infection though time by province

dat.age <- ddply(subset(dat, provname!="Tboung Khmum"), .(year, provname), summarise, age = mean(age))
head(dat.age)

dat.age$provname <- as.factor(dat.age$provname)

m4 <- gam(age ~ year:provname + s(provname, bs="re"), dat = dat.age, family="gaussian")
summary(m4) #positive slope, with varied steepness


out=summary(m4)

prov.df.4 <- cbind.data.frame(provname = unique(dat.prov$provname), 
                              slope= out$p.coeff[2:length(out$p.coeff)],
                              se = out$se[2:length(out$p.coeff)],
                              p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.4$slope_uci <- prov.df.4$slope + 1.96*prov.df.4$se
prov.df.4$slope_lci <- prov.df.4$slope - 1.96*prov.df.4$se
prov.df.4$sig <- "sig"
prov.df.4$sig[prov.df.4$p_val>0.05] <- "not_sig"
rownames(prov.df.4) <- c()


#and predict
predict.dat.4 <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))


predict.dat.4$predict_cases <- predict.gam(m4,newdata = predict.dat.4,  type = "response")
predict.dat.4$predict_cases_lci  <- predict.dat.4$predict_cases - 1.96*predict.gam(m4, newdata = predict.dat.4,   type = "response", se.fit = T)$se
predict.dat.4$predict_cases_uci  <- predict.dat.4$predict_cases + 1.96*predict.gam(m4, newdata = predict.dat.4,   type = "response", se.fit = T)$se



#and supplemental figure
FigS2 <- ggplot(predict.dat.4) + #scale_color_manual(values=colz) +
  #scale_fill_manual(values=colz) +
  geom_line(data= dat.age, aes(x=year, y=age), color="gray70") +
  geom_line(aes(x=year, y=predict_cases, color=provname)) + theme_bw() +
  geom_ribbon(aes(x=year, ymin=predict_cases_lci, ymax=predict_cases_uci, fill=provname), alpha=.3) + 
  facet_wrap(~provname, ncol=4) + ylab("predicted mean age of infection") +
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
        strip.background = element_rect(fill="white"),
        legend.title = element_blank(), legend.position = "bottom",
        axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        legend.text = element_text(size=12)) 




ggsave(paste0(homewd, "/fig-new/FigS2.png"),
       plot = FigS2,
       bg='white',
       width=90, 
       height=80, 
       scale=3.2, 
       dpi=300,
       units = "mm")



#and space


prov.df.4$provname <- as.character(prov.df.4$provname)

setdiff(unique(prov.df.4$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(prov.df.4$provname)) #"Tboung Khmum"


#and add the blanks
head(prov.df.4)
prov.add <- subset(prov.df.4, provname=="Kandal")
prov.add$provname <- "Tboung Khmum"
prov.add$sig <- "not_sig"
prov.add$slope <- prov.add$se <- prov.add$p_val <- prov.add$slope_lci <- prov.add$slope_uci <- NA

prov.df.4 <- rbind(prov.df.4, prov.add)
#and merge with the slopes
#dat.merge <- dplyr::select(case.dat, y, IsSignificant, provname, diagnostic)

names(prov.df.4)[names(prov.df.4)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam, prov.df.4, by = "ADM1_EN", all.x=T, sort=F)

cam_sub = subset(cam_merge, is.na(slope))
head(cam_merge)

unique(prov.df$sig)
colorz <- c('sig' = "black", 'not_sig'="grey80")

#library(patchwork)

Fig1C <- ggplot(cam_merge) + geom_sf(aes(fill=slope, color=sig), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,.5,0,1.5),"lines")) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', 
                       na.value = "grey50",
                       name="slope of trend of increasing\nmean age of infection") +
  scale_color_manual(values = colorz) +
  guides(color="none")

head(prov.df)
head(prov.df.4)

prov.df.slim <- dplyr::select(subset(prov.df, diagnostic=="df"), ADM1_EN, slope)
prov.df.4.slim <- dplyr::select(prov.df.4, ADM1_EN, slope)

names(prov.df.slim) <- c("provname", "slope_df")
names(prov.df.4.slim) <- c("provname", "slope_age")
prov.merge <- merge(prov.df.slim, prov.df.4.slim, by="provname", all.x = T)       

p1 <- ggplot(data=prov.merge) + geom_point(aes(x=slope_age, y=slope_df, color=provname), size=3) +
      theme_bw() + xlab("slope age increase") + ylab("slope incidence increase") +
      theme(panel.grid = element_blank(),
            legend.title = element_blank(),
            axis.text = element_text(size=14), axis.title = element_text(size=16),
            plot.margin = unit(c(.5,.5,.5,.5),"lines"))

colorMap <- ggplot(cam_merge) + geom_sf(aes(fill=ADM1_EN), size =.4, show.legend = F, color="black") + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,.5,0,1.5),"lines")) +
  guides(color="none")


#what about climate patterns by province?


climdat <- read.csv(file = paste0(homewd, "/data/all_case_climate.csv"), header=T, stringsAsFactors = F)

head(climdat) 
climdat$month_date <- as.Date(climdat$month_date)
climdat$year <- year(climdat$month_date)
climdat$month <- month(climdat$month_date)
climdat$provname <- as.factor(climdat$provname)

m5 <- gam(temp_C ~ year:provname + 
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = climdat)
summary(m5) # all pos but varying

out=summary(m5)

prov.df.5 <- cbind.data.frame(provname = unique(climdat$provname), 
                              slope= out$p.coeff[2:length(out$p.coeff)],
                              se = out$se[2:length(out$p.coeff)],
                              p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.5$slope_uci <- prov.df.5$slope + 1.96*prov.df.5$se
prov.df.5$slope_lci <- prov.df.5$slope - 1.96*prov.df.5$se
prov.df.5$sig <- "sig"
prov.df.5$sig[prov.df.5$p_val>0.05] <- "not_sig"
rownames(prov.df.5) <- c()


#and predict
predict.dat.5 <- cbind.data.frame(provname=rep(unique(climdat$provname), each = length(unique(climdat$year))), year=rep(unique(climdat$year), length(unique(climdat$provname))))

prov.df.5$provname <- as.character(prov.df.5$provname)

setdiff(unique(prov.df.5$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(prov.df.5$provname)) 

cam$ADM1_EN[cam$ADM1_EN=="Otdar Meanchey"] <- "Oddar Meanchey"

#and add the blanks
head(prov.df.5)

names(prov.df.5)[names(prov.df.5)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam, prov.df.5, by = "ADM1_EN", all.x=T, sort=F)

cam_sub = subset(cam_merge, is.na(slope))
head(cam_merge)


colorz <- c('sig' = "black", 'not_sig'="grey80")

#library(patchwork)

FigTemp <- ggplot(cam_merge) + geom_sf(aes(fill=slope, color=sig), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,.5,0,1.5),"lines")) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', 
                       na.value = "grey50",
                       name="slope of trend of increasing\ntemp C") +
  scale_color_manual(values = colorz) +
  guides(color="none")



#and precip

m6 <- gam(precip_mm ~ year:provname + 
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = climdat)
summary(m6) # all pos but varying

out=summary(m6)

prov.df.6 <- cbind.data.frame(provname = unique(climdat$provname), 
                              slope= out$p.coeff[2:length(out$p.coeff)],
                              se = out$se[2:length(out$p.coeff)],
                              p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.6$slope_uci <- prov.df.6$slope + 1.96*prov.df.6$se
prov.df.6$slope_lci <- prov.df.6$slope - 1.96*prov.df.6$se
prov.df.6$sig <- "sig"
prov.df.6$sig[prov.df.6$p_val>0.05] <- "not_sig"
rownames(prov.df.6) <- c()


#and predict
predict.dat.6 <- cbind.data.frame(provname=rep(unique(climdat$provname), each = length(unique(climdat$year))), year=rep(unique(climdat$year), length(unique(climdat$provname))))

prov.df.6$provname <- as.character(prov.df.6$provname)

setdiff(unique(prov.df.6$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(prov.df.6$provname)) 

#cam$ADM1_EN[cam$ADM1_EN=="Otdar Meanchey"] <- "Oddar Meanchey"

#and add the blanks
head(prov.df.6)

names(prov.df.6)[names(prov.df.6)=="provname"] <- "ADM1_EN"
cam_merge2 <- merge(cam, prov.df.6, by = "ADM1_EN", all.x=T, sort=F)


colorz <- c('sig' = "black", 'not_sig'="grey80")

#library(patchwork)

FigPrecip <- ggplot(cam_merge2) + geom_sf(aes(fill=slope, color=sig), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,.5,0,1.5),"lines")) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', 
                       na.value = "grey50",
                       name="slope of trend of increasing\nprecip (mm)") +
  scale_color_manual(values = colorz) +
  guides(color="none")



#and plot together
prov.df.5.slim <- dplyr::select(prov.df.5, ADM1_EN, slope)
prov.df.6.slim <- dplyr::select(prov.df.6, ADM1_EN, slope)
names(prov.df.5.slim ) <- c("provname", "slope_temp")
names(prov.df.6.slim ) <- c("provname", "slope_precip")

prov.merge <- merge(prov.merge, prov.df.5.slim, by="provname", all.x = T)
prov.merge <- merge(prov.merge, prov.df.6.slim, by="provname", all.x = T)




p2 <- ggplot(data=prov.merge) + geom_point(aes(x=slope_temp, y=slope_df, color=provname), size=3) +
  theme_bw() + ylab("slope incidence increase") + xlab("slope temperature change") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.margin = unit(c(.5,.5,.5,.5),"lines"))


p3 <- ggplot(data=prov.merge) + geom_point(aes(x=slope_age, y=slope_precip, color=provname), size=3) +
  theme_bw() + xlab("slope precipitation change") + ylab("slope incidence increase") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.margin = unit(c(.5,.5,.5,.5),"lines"))

