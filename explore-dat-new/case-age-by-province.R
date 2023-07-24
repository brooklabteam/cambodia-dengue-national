rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
  
head(dat) 
min(dat$year)

#plot time series of each type by province by year

unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$epimonth <- cut.Date(dat$date, breaks="months", start.on.monday = T)
dat$epimonth <- as.Date(as.character(dat$epimonth))

head(dat)

class(dat$epimonth)

dat$case=1

#total cases by week by year by province
dat.prov <- ddply(dat, .(epimonth, provname), summarise, cases = sum(case))
head(dat.prov)


#plot cases through time by province and pull out the annual direction
p1 <- ggplot(dat.prov) + geom_line(aes(x=epimonth, y=cases, color=provname))

dat.prov$year <- year(dat.prov$epimonth)
dat.prov$month <- month(dat.prov$epimonth)

dat.list <- dlply(dat.prov, .(provname))

get.annual.direction <- function(df, k){
  
  m1 <- gam(cases ~ year + s(month, k=k, bs="cc"),
            family="poisson",data = df)
  out = summary(m1)
  
  df$predict_gam <- predict.gam(m1, type="response", exclude = "s(month)")
  df$predict_gam_lci <- df$predict_gam - 1.96*predict.gam(m1, type="response", exclude = "s(month)", se.fit = T)$se
  df$predict_gam_uci <- df$predict_gam + 1.96*predict.gam(m1, type="response", exclude = "s(month)", se.fit = T)$se
  
  p1b <- ggplot(data=df) + theme_bw() +
    geom_line(aes(x=epimonth, y= cases)) + ylab("crude case numbers") +
    geom_line(aes(x=epimonth, y=predict_gam), color="red") +
    geom_ribbon(aes(x=epimonth, ymin=predict_gam_lci, ymax=predict_gam_uci), alpha=.3, fill="red")
  
  
  #print(p1b)
  
  #add pval and slop
  df$annual_slope <- out$p.table[2,1]
  df$annual_slope_lci <-df$annual_slope - (1.96*out$p.table[2,2])
  df$annual_slope_uci <-df$annual_slope + (1.96*out$p.table[2,2])
  df$annual_pval <- out$p.table[2,4]
  return(df)
}     

dat.list <- lapply(dat.list, get.annual.direction, k=7)
 
dat.new <- data.table::rbindlist(dat.list)
head(dat.new)
dat.new$direction <- "pos"
dat.new$sig <- "sig"
dat.new$sig[dat.new$annual_pval>0.01] <- "not-sig"
unique(dat.new$sig )
dat.new$direction[dat.new$annual_slope<0] <- "neg"


#now plot
colz = c('pos'="tomato", 'neg'="cornflowerblue", 'not-sig' = "gray30")

p1 <- ggplot(dat.new) + geom_line(aes(x=epimonth, y=cases)) +
      geom_line(aes(x=epimonth, y=predict_gam, color=direction)) +
      geom_ribbon(aes(x=epimonth, ymin=predict_gam_lci, 
                      ymax=predict_gam_uci, fill=direction), alpha=.3) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
      scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
      facet_wrap(~provname, nrow=5, ncol=5)



#do the same for dss and dhf
#total cases by week by year by province
dat.prov.diag <- ddply(dat, .(epimonth, provname, diagnostic), summarise, cases = sum(case))
head(dat.prov.diag)
dat.prov.diag$year <- year(dat.prov.diag$epimonth)
dat.prov.diag$month <- month(dat.prov.diag$epimonth)

dat.prov.diag.list <- dlply(dat.prov.diag, .(provname, diagnostic))

dat.prov.diag <-data.table::rbindlist(lapply(dat.prov.diag.list, get.annual.direction, k=5))

dat.prov.diag$direction <- "pos"
dat.prov.diag$sig <- "sig"
dat.prov.diag$sig[dat.prov.diag$annual_pval>0.01] <- "not-sig"
unique(dat.prov.diag$sig )
dat.prov.diag$direction[dat.prov.diag$annual_slope<0] <- "neg"
dat.prov.diag$direction[dat.prov.diag$sig=="not-sig"] <- "not-sig"


p2 <- ggplot(subset(dat.prov.diag, diagnostic=="df")) + geom_line(aes(x=epimonth, y=cases)) +
  geom_line(aes(x=epimonth, y=predict_gam, color=direction)) +
  geom_ribbon(aes(x=epimonth, ymin=predict_gam_lci, 
                  ymax=predict_gam_uci, fill=direction), alpha=.3) +
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  facet_wrap(~provname, nrow=5, ncol=5) + ggtitle(label="dengue fever")


p3 <- ggplot(subset(dat.prov.diag, diagnostic=="dhf")) + geom_line(aes(x=epimonth, y=cases)) +
  geom_line(aes(x=epimonth, y=predict_gam, color=direction)) +
  geom_ribbon(aes(x=epimonth, ymin=predict_gam_lci, 
                  ymax=predict_gam_uci, fill=direction), alpha=.3) +
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  facet_wrap(~provname, nrow=5, ncol=5) + ggtitle(label="dengue hemorrhagic fever")


p4 <- ggplot(subset(dat.prov.diag, diagnostic=="dss")) + geom_line(aes(x=epimonth, y=cases)) +
  geom_line(aes(x=epimonth, y=predict_gam, color=direction)) +
  geom_ribbon(aes(x=epimonth, ymin=predict_gam_lci, 
                  ymax=predict_gam_uci, fill=direction), alpha=.3) +
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  facet_wrap(~provname, nrow=5, ncol=5) + ggtitle(label="dengue shock syndrome")


#and also join mean age of infection

head(dat)

dat.age <- ddply(dat, .( year, provname), summarise, mean_age=mean(age))
head(dat.age)
dat.age.list <- dlply(dat.age, .(provname))


get.annual.year <- function(df){
  
  m1 <- lm(mean_age ~ year,data = df)
  out = summary(m1)
  
  df$predict <- predict(m1, type="response")
  df$predict_lci <- df$predict - 1.96*predict(m1, type="response", se.fit = T)$se
  df$predict_uci <- df$predict + 1.96*predict(m1, type="response",  se.fit = T)$se
  
  p1b <- ggplot(data=df) + theme_bw() +
    geom_line(aes(x=year, y= mean_age)) + ylab("mean age infection") +
    geom_line(aes(x=year, y=predict), color="red") +
    geom_ribbon(aes(x=year, ymin=predict_lci, ymax=predict_uci), alpha=.3, fill="red")
  
  
  #print(p1b)
  
  #add pval and slop
  df$annual_slope <- out$coefficients[2,1]
  df$annual_slope_lci <-df$annual_slope - (1.96*out$coefficients[2,2])
  df$annual_slope_uci <-df$annual_slope + (1.96*out$coefficients[2,2])
  df$annual_pval <- out$coefficients[2,4]
  return(df)
}  




dat.age <-data.table::rbindlist(lapply(dat.age.list, get.annual.year))

dat.age$direction <- "pos"
dat.age$sig <- "sig"
dat.age$sig[dat.age$annual_pval>0.01] <- "not-sig"
unique(dat.age$sig )
dat.age$direction[dat.age$annual_slope<0] <- "neg"
dat.age$direction[dat.age$sig=="not-sig"] <- "not-sig"


p5 <- ggplot(dat.age) + geom_line(aes(x=year, y=mean_age)) +
  geom_line(aes(x=year, y=predict, color=direction)) +
  geom_ribbon(aes(x=year, ymin=predict_lci, 
                  ymax=predict_uci, fill=direction), alpha=.3) +
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  facet_wrap(~provname, nrow=5, ncol=5) + ggtitle(label="mean age of infection")


#now take each of these directions by province and summarise slopes





