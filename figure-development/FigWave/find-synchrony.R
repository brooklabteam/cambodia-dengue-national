rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(WaveletComp)
library(mgcv)
library(reshape2)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)


climdat <- read.csv(file = paste0(homewd, "/data/all_case_climate.csv"), header=T, stringsAsFactors = F)

climdat$month_date <- as.Date(climdat$month_date)
climdat$year <- year(climdat$month_date)

head(climdat) 

#join with pop data
popdat <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dat.csv"), header=T, stringsAsFactors = F)
head(popdat)
popdat <- dplyr::select(popdat, year, pop,provname)

climdat <- merge(climdat, popdat, by=c("provname", "year"), all.x = T)

climdat$cases_per_1000 <- (climdat$cases/climdat$pop)*1000

#ggplot(climdat) +geom_line(aes(x=month_date, y=pop, color=provname))

#first, for each time series per province, run a function to collect:
#(a) the reconstructed period for each timestep, both annual
#(b) and multi-annual
#(c) the average wavelet power per timestep, for annual
#(d) and multi-annual
#(e) the average wavelet coherency with ONI
#(f) the average wavelet coherency with temperature for that province
#(g) the average wavelet cohernecy with precipitation for that province
#(h) the proportion of other provinces with which it shares a statistically significant coherency

prov.split <- dlply(climdat, .(provname))


prov.rank.multi <- function(df2, df1){
  df1a <-  df1
  df2a <- df2
  df1 <- dplyr::select(df1, month_date, cases_per_1000)
  df2 <- dplyr::select(df2, month_date, cases_per_1000)
  
  names(df1)[names(df1)=="cases_per_1000"] <- "cases_this_prov"
  names(df2)[names(df2)=="cases_per_1000"] <- "cases_other_prov"
  
  df.merge <- merge(df1, df2, by="month_date", all.x = T)
  head(df.merge)
  
  df.merge <- df.merge[complete.cases(df.merge),]
  
  corr.prov <- analyze.coherency(df.merge, my.pair = c("cases_this_prov","cases_other_prov"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 12, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 2, #shortest possible period in years (multi)
                                 upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                 make.pval = TRUE, n.sim = 100)
  
  
  # 
  # wc.image(corr.prov, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.prov,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  
  #is there a statistically significant coherence through time?
  #coherence.mat <- corr.precip$Coherence
  pval.mat <- corr.prov$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  provname <- sub(" ", "_", unique(df2a$provname), fixed = T)
  
  df.merge$coherence_multi <- colSums(pval.mat)
  
  df2add <- dplyr::select(df.merge, month_date, coherence_multi)
  
  df1a <- merge(df1a, df2add, by="month_date", all.x = T)
  
  df1a$coherence_prov <- provname
  
  df.out <- dplyr::select(df1a, month_date, coherence_multi, coherence_prov)
  
  return(df.out)
  
  
}
prov.rank.annual <- function(df2, df1){
  
  df1a <-  df1
  df2a <- df2
  df1 <- dplyr::select(df1, month_date, cases_per_1000)
  df2 <- dplyr::select(df2, month_date, cases_per_1000)
  
  names(df1)[names(df1)=="cases_per_1000"] <- "cases_this_prov"
  names(df2)[names(df2)=="cases_per_1000"] <- "cases_other_prov"
  
  df.merge <- merge(df1, df2, by="month_date", all.x = T)
  head(df.merge)
  
  df.merge <- df.merge[complete.cases(df.merge),]
  
  corr.prov <- analyze.coherency(df.merge, my.pair = c("cases_this_prov","cases_other_prov"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 12, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 1/12, #shortest possible period in years (annual)
                                 upperPeriod = 2, #largest possible period (2 years)
                                 make.pval = TRUE, n.sim = 100)
  
  
  # 
  # wc.image(corr.prov, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.prov,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  
  #is there a statistically significant coherence through time?
  #coherence.mat <- corr.precip$Coherence
  pval.mat <- corr.prov$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  provname <- sub(" ", "_", unique(df2a$provname), fixed = T)
  
  df.merge$coherence_annual <- colSums(pval.mat)
  
  df2add <- dplyr::select(df.merge, month_date, coherence_annual)
  
  df1a <- merge(df1a, df2add, by="month_date", all.x = T)
  
  df1a$coherence_prov <- provname
  
  df.out <- dplyr::select(df1a, month_date, coherence_annual, coherence_prov)
  
  return(df.out)
  
  
}
get.mean.period <- function(col.df){
  
  #bind
  dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  mean_period = weighted.mean(x=col.df$period_length, w=col.df$value)
  
  return(mean_period)
}
get.wavelet.dat <- function(dat, dat.all){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual
  
  #(a) annual
  anal.dat.annual <-  analyze.wavelet(dat,
                                     my.series = ncol(dat), #cases per 1000
                                     #loess.span = 0,
                                     dt = 1/12,#this allows for annual timestep
                                     dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                     lowerPeriod = 1/12,#shortest possible period (2)
                                     upperPeriod = 2, #largest possible period (in weeks; here, 10 years)
                                     make.pval = TRUE, n.sim = 100)
  
   
  anal.dat = reconstruct(anal.dat.annual)
  
  
  dat$reconstructed_annual_period <- anal.dat$series$cases_per_1000.r
  
  # (b) multi
  anal.dat.multi <-  analyze.wavelet(dat,
                                    my.series = "cases_per_1000", #cases per 1000
                                    #loess.span = 0,
                                    dt = 1/12,#this allows for annual timestep
                                    dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                    lowerPeriod = 2,#shortest possible period (2)
                                    upperPeriod = 20, #largest possible period (in weeks; here, 10 years)
                                    make.pval = TRUE, n.sim = 100)
  
  multi.dat = reconstruct(anal.dat.multi)
  
  
  #reconstruct(anal.dat.multi, only.ridge = T)
  dat$reconstructed_multi_period <- multi.dat$series$cases_per_1000.r
  
  #and, as bonus, get the average length of the multi-annual dengue cycles
  #col.split = as.list(as.data.frame(anal.dat.multi$Power))
  
  power.df = as.data.frame(anal.dat.multi$Power)
  names(power.df) <- 1:length(dat$reconstructed_multi_period)
  power.df$period <- as.numeric(rownames(power.df))
  
  power.df.melt <- melt(power.df, id.vars = "period")
  names(power.df.melt)[names(power.df.melt)=="variable"] <- "time"
  power.df.melt$time <- as.numeric(as.character(power.df.melt$time))
  
  period.df = cbind.data.frame(period=1:length(unique(power.df.melt$period)), period_length=anal.dat.multi$Period)
  
  power.df.melt <- merge(power.df.melt, period.df, by="period", all.x=T)
  
  col.split <- dlply(power.df.melt, .(time))
  
  #ggplot(power.df.melt) + geom_tile(aes(x=time, y=period, fill=value))
  
  
  mean.period.list <- lapply(col.split, get.mean.period)
  
  dat$multi_period <-  c(unlist(mean.period.list))
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=multi_period))
  
  #(c) then, average wavelet power for annual
  dat$avg_wave_power_annual <- (colSums(anal.dat.annual$Power)/length(anal.dat.annual$Power.avg))

  #(d) average wavelet power for multiannual
  dat$avg_wave_power_multi <- (colSums(anal.dat.multi$Power)/length(anal.dat.multi$Power.avg))
  
  #(e) average wavelet coherency with ONI (multi)
  corr.oni <- analyze.coherency(dat, my.pair = c("oni","cases_per_1000"),
                                loess.span = 0,
                                dt = 1/12, dj = 1/100,
                                window.type.t = 1, window.type.s = 1,
                                window.size.t = 12, #examine coherence year-by-year
                                window.size.s = (1/4), #periods on the order of 
                                lowerPeriod = 2, #shortest possible period in years
                                upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                make.pval = TRUE, n.sim = 100)
  
  # wc.image(corr.oni, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.oni,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.oni$Coherence
  pval.mat <- corr.oni$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  out.mat <- coherence.mat*pval.mat
  #plot(colMeans(out.mat))
  #wc.avg(corr.oni)
  
  
  dat$avg_wave_coherency_oni <- colMeans(out.mat)
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.oni$Power.xy
  power.table.pval <- corr.oni$Power.xy.pval
  power.table.pval <- power.table.pval<0.05
  
  power.mat <- power.table*power.table.pval
  
  #dat$avg_cross_power_oni <-  (colSums(corr.oni$Power.xy)/length(corr.oni$Power.xy.avg))
  dat$avg_cross_power_oni <-  (colSums( power.mat)/length(corr.oni$Power.xy.avg))
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_wave_coherency_oni)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_cross_power_oni)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #and also take the average power level across the time series
  
  #(f) the average wavelet coherency with temperature for that province
  
  corr.temp <- analyze.coherency(dat, my.pair = c("temp_C","cases_per_1000"),
                                loess.span = 0,
                                dt = 1/12, dj = 1/100,
                                window.type.t = 1, window.type.s = 1,
                                window.size.t = 12, #examine coherence year-by-year
                                window.size.s = (1/4), #periods on the order of 
                                lowerPeriod = 2, #shortest possible period in years
                                upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                make.pval = TRUE, n.sim = 100)
  
  # wc.image(corr.temp, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.temp,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.temp$Coherence
  pval.mat <- corr.temp$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  out.mat <- coherence.mat*pval.mat
  #plot(colMeans(out.mat))
  #wc.avg(corr.temp)
  
  
  dat$avg_wave_coherency_temp <- colMeans(out.mat)
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.temp$Power.xy
  power.table.pval <- corr.temp$Power.xy.pval
  power.table.pval <- power.table.pval<0.05
  
  power.mat <- power.table*power.table.pval
  
  
  dat$avg_cross_power_temp <-  (colSums(power.mat)/length(corr.temp$Power.xy.avg))
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_wave_coherency_temp)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_cross_power_temp)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  
  
  #(g) the average wavelet coherency with precip for that province
  
  
  corr.precip <- analyze.coherency(dat, my.pair = c("precip_mm","cases_per_1000"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 12, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 2, #shortest possible period in years
                                 upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                 make.pval = TRUE, n.sim = 100)
  
  # wc.image(corr.precip, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.precip,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.precip$Coherence
  pval.mat <- corr.precip$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  out.mat <- coherence.mat*pval.mat
  #plot(colMeans(out.mat))
  #wc.avg(corr.precip)
  
  
  dat$avg_wave_coherency_precip <- colMeans(out.mat)
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.precip$Power.xy
  power.table.pval <- corr.precip$Power.xy.pval
  power.table.pval <- power.table.pval<0.05
  
  power.mat <- power.table*power.table.pval
  
  
  dat$avg_cross_power_precip <-  (colSums(power.mat)/length(corr.precip$Power.xy.avg))
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_wave_coherency_precip)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_cross_power_precip)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  
  
  #(h) the proportion of other provinces with which it shares a statistically significant coherency
  
  #split all the others and try annual coherence
  df.split <- dlply(dat.all, .(provname))
  
  dat.prov.coherence.annual <- lapply(X=df.split, FUN=prov.rank.annual, df1=dat)
  dat.prov.coherence.annual <- data.table::rbindlist(dat.prov.coherence.annual)
  #head(dat.prov.coherence.annual)
  
  dat.prov.coherence.annual$sig_notsig <- 0
  dat.prov.coherence.annual$sig_notsig[dat.prov.coherence.annual$coherence_annual>0] <- 1
  dat.prov.coherence.annual$sig_notsig[is.na(dat.prov.coherence.annual$coherence_annual)] <- NA
  prov.coher.sum.annual <- ddply(dat.prov.coherence.annual, .(month_date), summarise, num_coherence=sum(sig_notsig, na.rm = T), N_tot=length(unique(coherence_prov)))
  prov.coher.sum.annual$proportion_coherence_annual <- prov.coher.sum.annual$num_coherence/prov.coher.sum.annual$N_tot
  
  #now attach to the rest of the dataset
  
  prov.add.annual <- dplyr::select(prov.coher.sum.annual, month_date, proportion_coherence_annual)
  
  dat <- merge(dat, prov.add.annual, by ="month_date", all.x = T)
  
  #and multi
  dat.prov.coherence.multi <- lapply(X=df.split, FUN=prov.rank.multi, df1=dat)
  dat.prov.coherence.multi <- data.table::rbindlist(dat.prov.coherence.multi)
  #head(dat.prov.coherence.multi)
  
  dat.prov.coherence.multi$sig_notsig <- 0
  dat.prov.coherence.multi$sig_notsig[dat.prov.coherence.multi$coherence_multi>0] <- 1
  dat.prov.coherence.multi$sig_notsig[is.na(dat.prov.coherence.multi$coherence_multi)] <- NA
  prov.coher.sum.multi <- ddply(dat.prov.coherence.multi, .(month_date), summarise, num_coherence=sum(sig_notsig, na.rm = T), N_tot=length(unique(coherence_prov)))
  prov.coher.sum.multi$proportion_coherence_multi <- prov.coher.sum.multi$num_coherence/prov.coher.sum.multi$N_tot
  
  #now attach to the rest of the dataset
  
  prov.add.multi <- dplyr::select(prov.coher.sum.multi, month_date, proportion_coherence_multi)
  
  dat <- merge(dat, prov.add.multi, by ="month_date", all.x = T)
  
  
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=proportion_coherence_multi)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  
  #head(dat)
  
  return(dat)
  
}
get.wavelet.sub <- function(dat, dat.all){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual
  
  #(a) annual
  anal.dat.annual <-  analyze.wavelet(dat,
                                      my.series = ncol(dat), #cases per 1000
                                      #loess.span = 0,
                                      dt = 1/12,#this allows for annual timestep
                                      dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                      lowerPeriod = 1/12,#shortest possible period (2)
                                      upperPeriod = 2, #largest possible period (in weeks; here, 10 years)
                                      make.pval = TRUE, n.sim = 100)
  
  
  anal.dat = reconstruct(anal.dat.annual)
  
  
  dat$reconstructed_annual_period <- anal.dat$series$cases_per_1000.r
  
  # (b) multi
  anal.dat.multi <-  analyze.wavelet(dat,
                                     my.series = "cases_per_1000", #cases per 1000
                                     #loess.span = 0,
                                     dt = 1/12,#this allows for annual timestep
                                     dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                     lowerPeriod = 2,#shortest possible period (2)
                                     upperPeriod = 20, #largest possible period (in weeks; here, 10 years)
                                     make.pval = TRUE, n.sim = 100)
  
  multi.dat = reconstruct(anal.dat.multi)
  
  
  #reconstruct(anal.dat.multi, only.ridge = T)
  dat$reconstructed_multi_period <- multi.dat$series$cases_per_1000.r
  
  
  
  return(dat)
  
}
get.wavelet.sub2 <- function(dat, dat.all){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual
  
  
  corr.temp <- analyze.coherency(dat, my.pair = c("temp_C","cases_per_1000"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 12, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 1/12, #shortest possible period in years
                                 upperPeriod = 2, #largest possible period (in weeks; here, 20 years)
                                 make.pval = TRUE, n.sim = 100)
  
  # wc.image(corr.temp, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.temp,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.temp$Coherence
  pval.mat <- corr.temp$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  out.mat <- coherence.mat*pval.mat
  #plot(colMeans(out.mat))
  #wc.avg(corr.temp)
  
  
  dat$avg_wave_coherency_temp_annual <- colMeans(out.mat)
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.temp$Power.xy
  power.table.pval <- corr.temp$Power.xy.pval
  power.table.pval <- power.table.pval<0.05
  
  power.mat <- power.table*power.table.pval
  
  
  dat$avg_cross_power_temp_annual <-  (colSums(power.mat)/length(corr.temp$Power.xy.avg))
  
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_wave_coherency_temp)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_cross_power_temp)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  
  
  #(g) the average wavelet coherency with precip for that province
  
  
  corr.precip <- analyze.coherency(dat, my.pair = c("precip_mm","cases_per_1000"),
                                   loess.span = 0,
                                   dt = 1/12, dj = 1/100,
                                   window.type.t = 1, window.type.s = 1,
                                   window.size.t = 12, #examine coherence year-by-year
                                   window.size.s = (1/4), #periods on the order of 
                                   lowerPeriod = 1/12, #shortest possible period in years
                                   upperPeriod = 2, #largest possible period (in weeks; here, 20 years)
                                   make.pval = TRUE, n.sim = 100)
  
  # wc.image(corr.precip, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  # 
  # wc.image(corr.precip,which.image = "wc",  
  #          color.key = "interval", n.levels = 250,
  #          siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #          legend.params = list(lab = "wavelet coherence levels"),
  #          spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #          periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # 
  
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.precip$Coherence
  pval.mat <- corr.precip$Coherence.pval
  pval.mat <- pval.mat<0.05
  
  out.mat <- coherence.mat*pval.mat
  #plot(colMeans(out.mat))
  #wc.avg(corr.precip)
  
  
  dat$avg_wave_coherency_precip_annual <- colMeans(out.mat)
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.precip$Power.xy
  power.table.pval <- corr.precip$Power.xy.pval
  power.table.pval <- power.table.pval<0.05
  
  power.mat <- power.table*power.table.pval
  
  
  dat$avg_cross_power_precip_annual <-  (colSums(power.mat)/length(corr.precip$Power.xy.avg))
  
  
  return(dat)
  
}

prov.split.out <- lapply(prov.split, get.wavelet.dat, dat.all=climdat)
#prov.split.out <- lapply(prov.split, get.wavelet.sub, dat.all=climdat)
#prov.split.out <- lapply(prov.split, get.wavelet.sub2, dat.all=climdat)
prov.split.df <- data.table::rbindlist(prov.split.out)
head(prov.split.df)
names(prov.split.df)

# prov.split.df <- dplyr::select(prov.split.df, provname, month_date, avg_wave_coherency_temp_annual, avg_cross_power_temp_annual, avg_wave_coherency_precip_annual, avg_cross_power_precip_annual)
# prov.split.df$month_date <- as.Date(prov.split.df$month_date) 
# unique(prov.split.df$provname)
# 
# dat <- read.csv(file = paste0(homewd, "/data/synchrony_data.csv"), header = T, stringsAsFactors = F)
# dat$month_date <- as.Date(dat$month_date)
# unique(dat$provname)
# 
# dat <- merge(dat, prov.split.df, by=c("provname", "month_date"), all.x = T)
# head(dat)
# 
# write.csv(dat, file = paste0(homewd, "/data/synchrony_data.csv"), row.names = F)
write.csv(prov.split.df, file = paste0(homewd, "/data/synchrony_data.csv"), row.names = F)



