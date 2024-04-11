setwd(homewd)

# plot age distribution by province,
# mean age by province 
# and FOI by province, as multipanels

#then, decide how to arrange all together

dat = read.csv(file = paste0(homewd, "/data/mort-dat.csv"), header=T, stringsAsFactors = F)
head(dat)
library(reshape2)
library(ggplot2)
dat.melt <- melt(dat, id.vars = c("province", "year", "type"))
names(dat.melt)[4:5] <- c("month", "count")
unique(dat.melt$month)
dat.melt$month <- as.numeric(factor(dat.melt$month))
dat.melt$month<- as.character(dat.melt$month)
dat.melt$month[dat.melt$month!="12" & dat.melt$month!="11" & dat.melt$month!="10"] <- paste0("0", dat.melt$month[dat.melt$month!="12" & dat.melt$month!="11" & dat.melt$month!="10"])
dat.melt$year <- gsub(pattern=",", replacement = "", x= dat.melt$year, fixed = T)
dat.melt$year <- as.numeric(dat.melt$year)
dat.melt$date <- as.Date(paste0(dat.melt$year, "-", dat.melt$month, "-01"))
dat.melt$count <- as.numeric(dat.melt$count)
dat.melt$count[is.na(dat.melt$count)] <- 0

#total deaths
ggplot(subset(dat.melt, type=="deaths")) + geom_line(aes(x=date, y=count, color=province)) + facet_wrap(~province)

#and cfr
dat.cases = subset(dat.melt, type=="cases")
dat.deaths = subset(dat.melt, type=="deaths")
dat.deaths$cfr <- dat.deaths$count/dat.cases$count

ggplot(dat.deaths) + geom_line(aes(x=date, y=cfr, color=province)) + facet_wrap(~province) + coord_cartesian(ylim=c(0,.2))


#and also look at age distribution of clinical cases
dat = read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
head(dat)


p1 <- ggplot(subset(dat, diagnostic=="df")) + geom_jitter(aes(x=year, y=age, color=diagnostic), width=.1, size=.1, alpha=.3) + 
  facet_wrap(~provname) + theme_bw() + 
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year, color = diagnostic),  
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p1

p2 <- ggplot(subset(dat, diagnostic=="dhf")) + geom_jitter(aes(x=year, y=age), color="cornflowerblue",  width=.1, size=.1, alpha=.3) + 
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  facet_wrap(~provname) + theme_bw() + 
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year),  color="cornflowerblue", 
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p2


p3 <- ggplot(subset(dat, diagnostic=="dss")) + geom_jitter(aes(x=year, y=age), color="mediumseagreen",  width=.1, size=.1, alpha=.3) + 
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  facet_wrap(~provname) + theme_bw() + 
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year),  color="mediumseagreen",
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p3






p4 <- ggplot(subset(dat, provname=="Battambang")) + geom_jitter(aes(x=year, y=age, color=diagnostic), width=.1, size=.1, alpha=.3) + 
  facet_wrap(~provname) + theme_bw() + facet_grid(~diagnostic) +
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year, color = diagnostic),  position = position_dodge(),
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p4
?geom_violin

library(plyr)
library(dplyr)
rate.df <- dlply(dat,.(year, provname))

tot.sum <- function(df){
  tot.df = sum(df$case[df$diagnostic=="df"])
  tot.df[is.na(tot.df)] <- 0
  tot.dhf = sum(df$case[df$diagnostic=="dhf"])
  tot.dhf[is.na(tot.dhf)] <- 0
  tot.dss = sum(df$case[df$diagnostic=="dss"])
  tot.dss[is.na(tot.dss)] <- 0
  dhf.prop = tot.dhf/(tot.df + tot.dss +  tot.dhf)
  dss.prop = tot.dss/(tot.df + tot.dss +  tot.dhf)
  df.prop = tot.df/(tot.df + tot.dss +  tot.dhf)
  
  out.df = cbind.data.frame(provname=unique(df$provname), year=unique(df$year), prop_dhf=dhf.prop, prop_dss=dss.prop, prop_df =  df.prop)
  return(out.df)
}

out.df <- data.table::rbindlist(lapply(rate.df, tot.sum))

head(out.df)

ggplot(out.df) + geom_line(aes(x=year, y=prop_dhf, color=provname)) + facet_wrap(~provname)

ggplot(out.df) + geom_line(aes(x=year, y=prop_dss, color=provname)) + facet_wrap(~provname)

ggplot(out.df) + geom_line(aes(x=year, y=prop_df, color=provname)) + facet_wrap(~provname)


out.melt <- melt(out.df, id.vars = c("provname", "year"))
head(out.melt)
names(out.melt)[3] <- "diagnostic"

ggplot(out.melt) + geom_line(aes(x=year, y=value, color=diagnostic)) + facet_wrap(~provname)

