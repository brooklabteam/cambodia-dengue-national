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

dat$month_date <- as.Date(dat$month_date)

head(dat) 

#first, for each time series per province, run a function to collect:
#(a) the reconstructed period for each timestep, both annual
#(b) and multi-annual
#(c) the average wavelet power per timestep, for annual
#(d) and multi-annual
#(e) the average wavelet coherency with ONI
#(f) the average wavelet coherency with temperature for that province
#(g) the average wavelet cohernecy with precipitation for that province
#(h) the proportion of other provinces with which it shares a statistically significant coherency


prov.split <- dlply(dat, .(provname))

get.wavelet.dat <- function(dat, dat.all){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual
  
  #(a) annual
  anal.dat.annual <-  analyze.wavelet(dat,
                                     my.series = 4, #cases
                                     #loess.span = 0,
                                     dt = 1/12,#this allows for annual timestep
                                     dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                     lowerPeriod = 1/12,#shortest possible period (2)
                                     upperPeriod = 2, #largest possible period (in weeks; here, 10 years)
                                     make.pval = TRUE, n.sim = 1000)
  
  anal.dat = reconstruct(anal.dat.annual)
  
  
  dat$reconstructed_annual_period <- anal.dat$series$cases
  
  # (b) multi
  anal.dat.multi <-  analyze.wavelet(dat,
                                    my.series = 4, #cases
                                    #loess.span = 0,
                                    dt = 1/12,#this allows for annual timestep
                                    dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                    lowerPeriod = 2,#shortest possible period (2)
                                    upperPeriod = 20, #largest possible period (in weeks; here, 10 years)
                                    make.pval = TRUE, n.sim = 1000)
  
  multi.dat = reconstruct(anal.dat.multi)
  dat$reconstructed_multi_period <- multi.dat$series$cases
  
  
  #(c) then, average wavelet power for annual
  dat$avg_wave_power_annual <- (colSums(anal.dat.annual$Power)/length(anal.dat.annual$Power.avg))

  #(d) average wavelet power for multiannual
  dat$avg_wave_power_multi <- (colSums(anal.dat.multi$Power)/length(anal.dat.multi$Power.avg))
  
  #(e) average wavelet coherency with ONI (multi)
  corr.oni <- analyze.coherency(dat, my.pair = c("oni","cases"),
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
  #ggplot(dat) + geom_line(aes(x=month_date, y=avg_cross_power)) + geom_vline(aes(xintercept=as.Date("2007-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2012-01-01")), color="red") + geom_vline(aes(xintercept=as.Date("2019-01-01")), color="red")
  #and also take the average power level across the time series
  
  #(f) the average wavelet coherency with temperature for that province
  
  corr.temp <- analyze.coherency(dat, my.pair = c("temp_C","cases"),
                                loess.span = 0,
                                dt = 1/12, dj = 1/100,
                                window.type.t = 1, window.type.s = 1,
                                window.size.t = 12, #examine coherence year-by-year
                                window.size.s = (1/4), #periods on the order of 
                                lowerPeriod = 2, #shortest possible period in years
                                upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                make.pval = TRUE, n.sim = 100)
  
  wc.image(corr.temp, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  
  
  wc.image(corr.temp,which.image = "wc",  
           color.key = "interval", n.levels = 250,
           siglvl.contour = 0.1, siglvl.arrow = 0.05,
           legend.params = list(lab = "wavelet coherence levels"),
           spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  
  
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
  
  
  corr.precip <- analyze.coherency(dat, my.pair = c("precip_mm","cases"),
                                 loess.span = 0,
                                 dt = 1/12, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 12, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 2, #shortest possible period in years
                                 upperPeriod = 20, #largest possible period (in weeks; here, 20 years)
                                 make.pval = TRUE, n.sim = 100)
  
  wc.image(corr.precip, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  
  
  wc.image(corr.precip,which.image = "wc",  
           color.key = "interval", n.levels = 250,
           siglvl.contour = 0.1, siglvl.arrow = 0.05,
           legend.params = list(lab = "wavelet coherence levels"),
           spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
           periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  
  
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
  
  #split all the others and try...
  dat.split <- dlply(dat.all, .(provname))
  
  prov.rank <- function(df2, df1){
    df1a <-  df1
    df2a <- df2
    df1 <- dplyr::select(df1, month_date, cases)
    df2 <- dplyr::select(df2, month_date, cases)
    
    names(df1)[names(df1)=="cases"] <- "cases_this_prov"
    names(df2)[names(df2)=="cases"] <- "cases_other_prov"
    
    df.merge <- merge(df1, df2, by="month_date", all.x = T)
    head(df.merge)
    
    df.merge <- df.merge[complete.cases(df.merge),]
    
    corr.prov <- analyze.coherency(df.merge, my.pair = c("cases_this_prov","cases_other_prov"),
                                     loess.span = 0,
                                     dt = 1/12, dj = 1/100,
                                     window.type.t = 1, window.type.s = 1,
                                     window.size.t = 12, #examine coherence year-by-year
                                     window.size.s = (1/4), #periods on the order of 
                                     lowerPeriod = 2, #shortest possible period in years
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
    
    df2a$coherence <- colSums(pval.mat)
    df2add <- dplyr::select(df2a, month_date, coherence)
    
    df1a <- merge(df1a, df2add, by="month_date", all.x = T)
    
    df1a$coherence_prov <- provname
    
    df.out <- dplyr::select(df1a, month_date, coherence, coherence_prov)
    
    return(df.out)
    
    
  }
  
  
  dat.prov.coherence <- lapply(X=dat.split, FUN=prov.rank, df1=dat)
  dat.prov.coherence <- data.table::rbindlist(dat.prov.coherence)
  head(dat.prov.coherence)
  
  dat.prov.coherence$sig_notsig <- 0
  dat.prov.coherence$sig_notsig[dat.prov.coherence$coherence>0] <- 1
  dat.prov.coherence$sig_notsig[is.na(dat.prov.coherence$coherence)] <- NA
  prov.coher.sum <- ddply(dat.prov.coherence, .(month_date), summarise, num_coherence=sum(sig_notsig, na.rm = T), N_tot=length(unique(coherence_prov)))
  prov.coher.sum$proportion_coherence <- prov.coher.sum$num_coherence/prov.coher.sum$N_tot
  
  #now attach to the rest of the dataset
  
  prov.add <- dplyr::select(prov.coher.sum, month_date, proportion_coherence)
  
  dat <- merge(dat, prov.add, by ="month_date", all.x = T)
  #head(dat)
  
  return(dat)
  
}

prov.split.out <- lapply(prov.split, get.wavelet.dat, dat.all=dat)

