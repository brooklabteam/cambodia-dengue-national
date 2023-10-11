# testing our framework using the simplest SIR model

rm(list=ls())

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/"))
#set wdset wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(dplyr)
library(plyr)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(reshape2)


#helper functions
buildTMat_age_SIR_triple_two <- function(c, Npop, age.classes, surv.biwk, age.brk, foi,age.mult.foi, recov, sigma, age.rate){
  
  #this is a transition matrix set up for 3 serotypes which only allows transitions for two
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,14,14) 
  
  Tmat <- matrix(0,14*nage,14*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- (2*foi[j])
    mat1[2,1] <- foi[j]
    mat1[4,1] <- foi[j]
    #mat1[6,1] <- foi[j]
    
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    mat1[3,3] <- 1-(foi[j])
    mat1[8,3] <- foi[j]
    #mat1[9,3] <- foi[j]
    
    mat1[4,4] <- 1-recov
    mat1[5,4] <- recov
    
    mat1[5,5] <- 1-(foi[j])
    mat1[10,5] <- foi[j]
    #mat1[11,5] <- foi[j]
    
    
    #mat1[6,6] <- 1-recov
    #mat1[7,6] <- recov
    
    
    #mat1[7,7] <-  1-(2*foi[j])
    #mat1[12,7] <- foi[j]
    #mat1[13,7] <- foi[j]
    
    
    mat1[8,8] <- 1-recov
    mat1[14,8] <- recov
    
    #mat1[9,9] <- 1-recov
    #mat1[14,9] <- recov
    
    mat1[10,10] <- 1-recov
    mat1[14,10] <- recov
    
    #mat1[11,11] <- 1-recov
    #mat1[14,11] <- recov
    
    #mat1[12,12] <- 1-recov
    #mat1[14,12] <- recov
    
    #mat1[13,13] <- 1-recov
    #mat1[14,13] <- recov
    
    
    mat1[3,14] <- sigma[j]
    mat1[5,14] <- sigma[j]
    #mat1[7,14] <- sigma[j]
    mat1[14,14] <- 1-(2*sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*14+1):(j*14+14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*14+1):((j-1)*14+14),((j-1)*14+1):((j-1)*14+14)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*14+1):(j*14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
   return(Tmat)
  
}
buildTMat_age_SIR_triple_after <- function(c, Npop, age.classes, surv.biwk, age.brk, foi,age.mult.foi, recov, sigma, age.rate){
  
  #this is a transition matrix set up for 3 serotypes which only allows transitions for two
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,14,14) 
  
  Tmat <- matrix(0,14*nage,14*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- (2*foi[j])
    #mat1[2,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[6,1] <- foi[j]
    
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    mat1[3,3] <- 1-(2*foi[j])
    mat1[8,3] <- foi[j]
    mat1[9,3] <- foi[j]
    
    mat1[4,4] <- 1-recov
    mat1[5,4] <- recov
    
    mat1[5,5] <- 1-(foi[j])
    #mat1[10,5] <- foi[j]
    mat1[11,5] <- foi[j]
    
    
    mat1[6,6] <- 1-recov
    mat1[7,6] <- recov
    
    
    mat1[7,7] <-  1-(foi[j])
    #mat1[12,7] <- foi[j] XXX
    mat1[13,7] <- foi[j]
    
    
    mat1[8,8] <- 1-recov
    mat1[14,8] <- recov
    
    mat1[9,9] <- 1-recov
    mat1[14,9] <- recov
    
    mat1[10,10] <- 1-recov
    mat1[14,10] <- recov
    
    mat1[11,11] <- 1-recov
    mat1[14,11] <- recov
    
    mat1[12,12] <- 1-recov
    mat1[14,12] <- recov
    
    mat1[13,13] <- 1-recov
    mat1[14,13] <- recov
    
    
    #mat1[3,14] <- sigma[j]
    mat1[5,14] <- sigma[j]
    mat1[7,14] <- sigma[j]
    mat1[14,14] <- 1-(2*sigma[j])
    
    
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*14+1):(j*14+14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*14+1):((j-1)*14+14),((j-1)*14+1):((j-1)*14+14)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*14+1):(j*14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_age_SIR_triple <- function(c, Npop, age.classes, surv.biwk, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,14,14) 
  
  Tmat <- matrix(0,14*nage,14*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- (3*foi[j])
    mat1[2,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[6,1] <- foi[j]
    
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    mat1[3,3] <- 1-(2*foi[j])
    mat1[8,3] <- foi[j]
    mat1[9,3] <- foi[j]
    
    mat1[4,4] <- 1-recov
    mat1[5,4] <- recov
    
    mat1[5,5] <- 1-(2*foi[j])
    mat1[10,5] <- foi[j]
    mat1[11,5] <- foi[j]
    
    
    mat1[6,6] <- 1-recov
    mat1[7,6] <- recov
    
    
    mat1[7,7] <-  1-(2*foi[j])
    mat1[12,7] <- foi[j]
    mat1[13,7] <- foi[j]
    
    
    mat1[8,8] <- 1-recov
    mat1[14,8] <- recov
    
    mat1[9,9] <- 1-recov
    mat1[14,9] <- recov
    
    mat1[10,10] <- 1-recov
    mat1[14,10] <- recov
    
    mat1[11,11] <- 1-recov
    mat1[14,11] <- recov
    
    mat1[12,12] <- 1-recov
    mat1[14,12] <- recov
    
    mat1[13,13] <- 1-recov
    mat1[14,13] <- recov
    
    
    mat1[3,14] <- sigma[j]
    mat1[5,14] <- sigma[j]
    mat1[7,14] <- sigma[j]
    mat1[14,14] <- 1-(3*sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*14+1):(j*14+14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*14+1):((j-1)*14+14),((j-1)*14+1):((j-1)*14+14)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*14+1):(j*14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_age_SIR_triple_switch <- function(c, Npop, age.classes, surv.biwk, age.brk, foi,age.mult.foi, recov, sigma, age.rate){
  
  #this is a transition matrix set up for 3 serotypes which only allows transitions for two
  s <- nage <- length(age.classes)
  
  
  
  if (length(foi)==1) foi <- rep(foi,nage) # 
  if (length(foi)==s) foi <- foi
  if (length(foi) > 1 & length(foi)<s){
    foi.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(foi)){
      foi.list[[i]] = rep(foi[i], age.dur[i])
    }
    foi = c(unlist(foi.list))
  } 
  
  
  #then, now that foi is appropriate length of all the age classes, 
  #multiply by the age.mult vector
  foi <- foi*age.mult.foi
  
  if (length(sigma)==1) sigma <- rep(sigma,nage) # 
  if (length(sigma)==s) sigma <- sigma
  if (length(sigma) > 1 & length(sigma)<s){
    sigma.list <- list()
    
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(sigma)){
      sigma.list[[i]] = rep(sigma[i], age.dur[i])
    }
    sigma = c(unlist(sigma.list))
  } 
  
  
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  mat1 <- matrix(0,14,14) 
  
  Tmat <- matrix(0,14*nage,14*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- (2*foi[j])
    #mat1[2,1] <- foi[j]
    mat1[4,1] <- foi[j]
    mat1[6,1] <- foi[j]
    
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    mat1[3,3] <- 1-(2*foi[j])
    mat1[8,3] <- foi[j]
    mat1[9,3] <- foi[j]
    
    mat1[4,4] <- 1-recov
    mat1[5,4] <- recov
    
    mat1[5,5] <- 1-(foi[j])
    #mat1[10,5] <- foi[j]
    mat1[11,5] <- foi[j]
    
    
    mat1[6,6] <- 1-recov
    mat1[7,6] <- recov
    
    
    mat1[7,7] <-  1-(foi[j])
    #mat1[12,7] <- foi[j] XXX
    mat1[13,7] <- foi[j]
    
    
    mat1[8,8] <- 1-recov
    mat1[14,8] <- recov
    
    mat1[9,9] <- 1-recov
    mat1[14,9] <- recov
    
    mat1[10,10] <- 1-recov
    mat1[14,10] <- recov
    
    mat1[11,11] <- 1-recov
    mat1[14,11] <- recov
    
    mat1[12,12] <- 1-recov
    mat1[14,12] <- recov
    
    mat1[13,13] <- 1-recov
    mat1[14,13] <- recov
    
    
    #mat1[3,14] <- sigma[j]
    mat1[5,14] <- sigma[j]
    #mat1[7,14] <- sigma[j]
    mat1[14,14] <- 1-(sigma[j])
    
    
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*14+1):(j*14+14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*14+1):((j-1)*14+14),((j-1)*14+1):((j-1)*14+14)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*14+1):(j*14),((j-1)*14+1):(j*14)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
transform.vect <- function(vec, s, c){
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}
find.biweek = function(t, times){
  
  biwks <- sort(unique(round(revtrunc(times),4)))
  biwks <- c(biwks[2:length(biwks)], biwks[1])
  this.wk = round(revtrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}
get.age.struct = function(pop, s){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=s, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  age.mat <- Reduce('+', age.mat.list)
  return(age.mat)
}
get.age.struct.M = function(pop, c){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=c, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  #zage.mat <- Reduce('+', age.mat.list)
  age.mat <- lapply(age.mat.list, sum)
  age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  names(age.mat.dat) <- c("age", "pop")
  return(age.mat)
}
prev.by.age <- function(dat){
  dat$n_age = sum(dat$count)
  dat$prev = dat$count/dat$n_age
  #and seroprev
  dat$seropos = sum(dat$count[dat$class=="M"], dat$count[dat$class=="R"])
  dat$seroprev = dat$seropos/dat$n_age
  return(dat)
}
revtrunc = function(x){
  newx = x - floor(x)
  return(newx)
}
spline.fit = function(data){
  #replace any Inf vals only the okay values
  #data <- data[complete.cases(data),]
  data$seroprevalence[data$seroprevalence==Inf] <- 1
  data$seroprevalence[data$seroprevalence<0] <- 0
  # if(length(data$age)>1){
  spline1 = with(data, smooth.spline(x=age, y=seroprevalence)) 
  #  return(spline1)
  # }else{
  #  return(NA)
  #  }
}
pred.fit <- function(spline1, ages.new){
  prev.comp = predict(spline1, x=ages.new)
  return(prev.comp$y)
}
get.seas.seroprev = function(dat){
  #get doy in bat calendar
  #convert to biweek
  #calc seroprev by biweek
  
  dat$doy <- yday(dat$date)
  #now correct for birthday of bat in questiob
  if(unique(dat$species=="Pteropus rufus")){
    bday <- yday("2015-10-01") #320
  } else if(unique(dat$species=="Eidolon dupreanum")){
    bday = yday("2015-11-01")
  }
  #write over for doy
  dat$new_doy = NA
  for (i in 1:length(dat$doy)){
    if(dat$doy[i] > bday){
      dat$new_doy[i] <- dat$doy[i] - bday
    } else if (dat$doy[i] <= bday){
      dat$new_doy[i] <- dat$doy[i] + (365-bday)
    }
  }
  
  dat$doy <- dat$new_doy
  dat <- dplyr::select(dat, -(new_doy))
  
  #now convert to biweek
  dat$biwk <- NA
  brk <- seq(0,365, by=14)
  for (i in 1:length(dat$biwk)){
    tmp <- brk[dat$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    dat$biwk[i] <- which(brk==tmp2)
  }
  
  #now calc seroprev by biweek
  biwk.sero <- ddply(dat, .(biwk), summarize, seropos = sum(prev), sero_lci = sum(prev_lci), sero_uci = sum(prev_uci), n=length(prev))  
  biwk.sero$seroprev = biwk.sero$seropos/biwk.sero$n
  biwk.sero$seroprev_lci = biwk.sero$sero_lci/biwk.sero$n
  biwk.sero$seroprev_uci = biwk.sero$sero_uci/biwk.sero$n
  biwk.sero$biwk = as.numeric(biwk.sero$biwk)
  
  return(biwk.sero)
}
get.mod.seas = function(mod.out){
  #convert time to doy
  mod.out$doy = mod.out$time*365
  #then to biwk
  brk <- seq(0,365, by=14) #doy breaks by biweek
  brk <- brk[-length(brk)]
  #brk.seq = 1:length(brk)
  mod.out$biwk <- NA
  for (i in 1:length(mod.out$biwk)){
    tmp <- brk[mod.out$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    mod.out$biwk[i] <- which(brk==tmp2)
  }
  #then return
  return(mod.out)
}
cum.cases = function(df){
  df$cumcases = cumsum(df$count)
  return(df)
}
age.sero.plot = function(dat){
  with(dat, plot(age, seroprevalence, type = "b", ylim=c(0,1)))
}
get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1, 2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}
build.pop.mat = function(surv, surv_juv, s, adult_fec_vector){
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1=adult_fec_vector*surv #birth rates* survival rates for the entire age distribution
  #row1 = c(0, rep((adult_fec*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution
mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}#for matrix slicing. give it the number of rows and columns you want in the resulting matrices
stack.age = function(dat, s){
  dat.new= list()
  for (i in 1:s){
    dat.new[[i]] = dat[i,]
  }
  dat.new = c(unlist(dat.new))
  return(dat.new)
}
stack.class = function(mat.dat, c){
  
  new.dat = list()
  for (i in 1:c){
    new.dat[[i]] = mat.dat[i,]
  }
  new.dat = c(unlist(new.dat))
  return(new.dat)
}
sum.yr.all <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(df$age))
  
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
divide <- function(df, ntyr1){
  dfnew <- (df/ntyr1)
  return(dfnew)
  
}
divide.rate <- function(df, ntyr1){
  dfnew = (1-exp(-(df/ntyr1)))
  
  return(dfnew)
  
}
sim.SIR.age.double.num <- function(yrs, ntyr, age.brk, s,foi, recov, mort,  births, pop_vector, sigma){
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  if(length(mort)<yrs){
    mort = rep(list(mort), yrs)
  }else{
    mort = as.list(mort)
  }
  
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  
  
  
  #otherwise, assume no transformation is needed
  
  #
  #number of epidemic classes
  c=8
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  PM_init = rep(0, s)
  I3_init = rep(0, s)
  I4_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, I2_init, P2_init, I3_init, I4_init, PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  #replicate for births
  #Birth_ts <- N_pop_ts
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #make annual mortality rate into probability of mortality per biweek
  mort.biweek = as.list(c(unlist(lapply(lapply(mort, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  if(length(mort.biweek)!=(length(times)-1) | length(sigma.biweek)!=(length(times)-1) | length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    
    #print(foi.biweek[[i]])
    
    #Tmat <- buildTMat_age_SIR(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.biweek[i])^(1/ntyr), surv.juv.biwk =  (1-mort.biweek[i])^(1/ntyr),	foi = foi, sigma=(sigma/ntyr), recov=recov,  age.rate=1)
    Tmat <- buildTMat_age_SIR_double(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.biweek[[i]]),	foi = foi.biweek[[i]], sigma=sigma.biweek[[i]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    #Fmat <- buildFMatrix(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1)
    
    #Tmat is infection and survival.
    transMat <- Tmat #+ Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    
    #print(sum(N_pop_ts[,i]))
    #print(sum(nt1))
    
    #then, add in births into the 0 age class of the susceptibles 
    births_per_1000pop = birth_vector_biwk[[i]] #these are each biweek per 1000 people
    #we are modeling in thousands of people
    
    
    births_add_biweek <- births_per_1000pop*(sum(N_pop_ts[,i])) #fill in susceptible births into class 0
    
    
    nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
    #print(sum(nt1))
    
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  N.sim.dat$class <- rep(c("S", "I1", "P1", "I2", "P2", "I3", "I4", "Pm"), each = s)
  N.sim.dat$age <- rep(1:s, each = length(years))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  return( N.sim.df)
}
replicate.data <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    
    return(new.dat)
  }
  
  
}
mean.age <- function(df){
  df$mult <- df$age*df$count
  mean.age <- sum(df$mult)/sum(df$count)
  df2 <- cbind.data.frame(year=unique(df$year), mean_age=mean.age)
  return(df2)
}
plot.age.dist.three <- function(dat, save.plot, view.plot, filename, slim.quant, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  # 
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  #split by a year
  df.year <- dlply(df.sum,.(year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year
  df.age <- dlply(df.sum,.(year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data, slim.quant=slim.quant))
  
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  #and plot
  p1 <- ggplot(dat.age) + 
    geom_jitter(aes(x=year, y=age), width=.09, height=.09, size=.09, alpha=.8, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,max(dat$age)))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=80, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}
plot.ages.all <- function(dat, save.plot, view.plot, filename, slim.quant, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  
  if(length(unique(dat$age))>1){
    dat2 = subset(dat1, age<max(dat$age))  
  }
  
  # 
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(dat2,.(year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  
  
  #split by a year
  df.year <- dlply(df.sum,.(year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year
  df.age <- dlply(df.sum,.(year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data, slim.quant=slim.quant))
  
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  #and plot
  p1 <- ggplot(dat.age) + 
    geom_jitter(aes(x=year, y=age), width=.09, height=.09, size=.09, alpha=.8, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,max(dat$age)))#,xlim=c(2015,2020))
  
  
  if(save.plot==TRUE){
    ggsave(file = filename,
           plot= p1,
           units="mm",  
           width=80, 
           height=55, 
           scale=3, 
           dpi=300)
    
  }
  if(view.plot==TRUE){
    print(p1)
  }
  
  
  
  return(mean.df)
}
process.plot.age.cum <-function(dat, year.start){
  
  
  #plot age-structure I3 + I4
  denv.case = subset(dat,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  #denv.case = subset(dat, class=="Pm")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year) 
  print(p1)
  
  return(denv.dat)
}
process.plot.age.cum <-function(dat, year.start){
  
  
  #plot age-structure I3 + I4
  denv.case = subset(dat,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  #denv.case = subset(dat, class=="Pm")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
        scale_color_viridis_d(direction=-1,option="turbo")# + facet_wrap(~year) 
  print(p1)
  
  return(denv.dat)
}
check.equil <- function(dat){
  dat.ts <- ddply(dat, .(time, class), summarise, count=sum(count))
  #and get total by time
  dat.N  <- ddply(dat, .(time), summarise, Ntot=sum(count))
  
  dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  print(p1)
  return(dat.ts)
}
check.equil.cases.only <- function(dat){
  #and get total by time
  denv.case = subset(dat, class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  denv.case$serotype <- NA
  denv.case$serotype[denv.case$class=="I12" | denv.case$class=="I32"] <- "I2"
  denv.case$serotype[denv.case$class=="I13" | denv.case$class=="I23"] <- "I3"
  denv.case$serotype[denv.case$class=="I31" | denv.case$class=="I21"] <- "I1"
  
  #dat.ts <- ddply(denv.case, .(time, class), summarise, count=sum(count))
  dat.ser <- ddply(denv.case, .(time, serotype), summarise, count=sum(count))
  
  
  dat.N  <- ddply(denv.case, .(time), summarise, Ntot=sum(count))
  
  #dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ser <- merge(dat.ser, dat.N, by="time")
  
  #dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  dat.ser$proportion <- dat.ser$count/dat.ser$Ntot
  #dat.ts$proportion[is.na(dat.ts$proportion)] <- 0
  dat.ser$proportion[is.na(dat.ser$proportion)] <- 0
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ser) + geom_line(aes(x=time, y=proportion, color=serotype)) +facet_grid(serotype~.) + coord_cartesian(ylim=c(0,1))
  print(p1)
  return(dat.ser)
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  return(df.sum)
}
process.all <- function(df){
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "S", "I1","P1","I2", "P2","I3","I4", "PM")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  #par(mfrow = c(1,1))
  #with(dat.tot, plot(time, N, type="l")) #ylim =c(0,1.2*max(N))))
  #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  
  # prop.tot = dat.tot
  # 
  # prop.tot$S = prop.tot$S/prop.tot$N
  # prop.tot$I = prop.tot$Ii/prop.tot$N
  # prop.tot$R = prop.tot$R/prop.tot$N
  # 
  
  
  
  
  dat.melt = melt(dat.tot,id.vars = "time")
  
  #head(dat.melt)
  
  #take after burnin if we assume this is a virus at equilibrium
  dat.melt <- subset(dat.melt, time>length(burnin_lambda))
  dat.melt <- subset(dat.melt, !is.na(value))
  
  #dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.melt) = c("time", "class", "count")
  dat.N = subset(dat.melt, class =="N")
  dat.melt = subset(dat.melt, class !="N")
  dat.N <- dplyr::select(dat.N, -(class))
  names(dat.N)[names(dat.N)=="count"] <- "Ntot"
  
  dat.melt <- merge(dat.melt, dat.N, by="time")
  #head(dat.melt)
  #dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.melt$class = factor(dat.melt$class, levels=c("S", "I1","P1","I2", "P2","I3","I4", "PM"))
  dat.melt$proportion <- dat.melt$count/dat.melt$Ntot
  
  dat.melt$time <- dat.melt$time - length(burnin_lambda)
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=count, color=class)) #+ scale_color_manual(values=colz)
  
  #dat.melt$proportion = dat.melt$count/sim_pop
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=proportion, color=class)) #+ scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  #now, get the age-structured incidence - incidence is "I3 + I4"
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("S", "I1","P1","I2", "P2","I3","I4", "PM"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  
  #plot age structured cumulative incidence
  age.incidence = subset(age.dat.tot, class=="I3" | class == "I4")
  
  #and collect everyone else
  age.non = subset(age.dat.tot, class!="I3" & class!="I4")
  #head(age.incidence )
  
  age.incidence$year = trunc(age.incidence$times)
  age.non$year = trunc(age.non$times)
  age.incidence <- arrange(age.incidence, year, age)
  
  age.non <- arrange(age.non, year, age)
  
  names(age.non) <- c("times", "count_non_infectious", "class", "age", "year")
  
  age.non <- dplyr::select(age.non, -(class), -(times))
  
  age.non <- ddply(age.non, .(year, age), summarise, count_non_infectious = sum(count_non_infectious))
  
  
  age.sum <- ddply(age.incidence, .(year, age), summarise, count = sum(count))
  
  age.sum <- merge(age.sum, age.non, by=c("year", "age"))
  
  #and just return this counts of years, ages, and those infected and not
  
  
  return(age.sum)#, seas.prev))
  
}
add.year <- function(dat, add.on){
  dat$year <- trunc(dat$time) + add.on
  return(dat)
}
sim.SIR.age.double.num.intro.clim.maintain <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  
 
  
  
  
  
  
  # if(length(age.mult.foi)==1){
  #   age.mult.foi = rep(age.mult.foi, s)
  # }else if(length(age.mult.foi)>1 & length(age.mult.foi) < s){
  #   age.mult.foi = c(age.mult.foi, rep(age.mult.foi[length(age.mult.foi)], (s-length(age.mult.foi))))
  # }else if(length(age.mult.foi)>s){
  #   age.mult.foi = age.mult.foi[1:s]
  # }
  # 
  # if(length(age.mult.foi)!=s){
  #   print("age multiplier proportions are incorrect")
  # }
  # 
  #otherwise, assume no transformation is needed
  
  #
  #number of epidemic classes
  c=14
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s)# would add some in here initially if you wanted this to be a 3-serotype sim
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, I2_init, P2_init, I3_init, P3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init, PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    biwk1 <- find.biweek(t=i, times=times)
    clim.mod <- clim.vect[biwk1]
    
    #and select the appropriate age vector, based on the year
    age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
 
    
    
    if(nrow(age.mult.sub)>10){
      print("error with age multiplier sub-selection")
    }
    
    age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
    age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
    age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
    
    
    #and the mortality rates
    if(times[i]<min(mort$year)){
      mort.df = subset(mort, year==min(year))
    }else{
      mort.df = subset(mort, year==trunc(times[i]))
    }
    #make annual mortality rate by age into probability of mortality per biweek
    #these are deaths per 1000 people.
    #we multiply by the number in the age classes, so need to divide by 1000 here
    mort.df$deaths_per_cap <- mort.df$count/1000
    mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
    
    
    
    #here are the dynamics before the introduction - 2 strain:
    if(times[i]<((yr.intro) + biwk.intro/26)){ 
      
      
      
      Tmat <- buildTMat_age_SIR_triple_two(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi, sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      #we are modeling in 1000s of people, so no need to divide here
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      
      N_pop_ts[,i] <- nt1
      
      
      
      
    }else if (times[i]==((yr.intro) +biwk.intro/26)) { 
      print(paste0("introduction at time =", times[i]))
      #here are the dynamics at the introduction - converting from 2 to 3 strains
      
      #need to add infecteds to I3 disease status -- here at the lowest age class
      #this is row 6!
      N_pop_ts[6,(i-1)] <- 5
      
      #remove them from the Susceptible class
      N_pop_ts[1,(i-1)] <- (N_pop_ts[1,(i-1)]-5)
      
      #then, let her go from there!
      
      Tmat <- buildTMat_age_SIR_triple(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i] > ((yr.intro) +biwk.intro/26)){ 
      #and here are the dynamics after the introduction - 3 strains only
      #no need to introduce anything
      
      Tmat <- buildTMat_age_SIR_triple(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
    }
    
  }
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  N.sim.dat$class <- rep(c("S", "I1", "P1", "I2", "P2", "I3", "P3", "I12", "I13", "I21", "I23", "I31", "I32", "Pm"), each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$time <- N.sim.df$time-1
  N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.double.num.intro.clim.switch <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  
  
  
  
  
  
  
  # if(length(age.mult.foi)==1){
  #   age.mult.foi = rep(age.mult.foi, s)
  # }else if(length(age.mult.foi)>1 & length(age.mult.foi) < s){
  #   age.mult.foi = c(age.mult.foi, rep(age.mult.foi[length(age.mult.foi)], (s-length(age.mult.foi))))
  # }else if(length(age.mult.foi)>s){
  #   age.mult.foi = age.mult.foi[1:s]
  # }
  # 
  # if(length(age.mult.foi)!=s){
  #   print("age multiplier proportions are incorrect")
  # }
  # 
  #otherwise, assume no transformation is needed
  
  #
  #number of epidemic classes
  c=14
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s)# would add some in here initially if you wanted this to be a 3-serotype sim
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, I2_init, P2_init, I3_init, P3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init, PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    biwk1 <- find.biweek(t=i, times=times)
    clim.mod <- clim.vect[biwk1]
    
    #and select the appropriate age vector, based on the year
    age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
    
    
    
    if(nrow(age.mult.sub)>10){
      print("error with age multiplier sub-selection")
    }
    
    age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
    age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
    age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
    
    
    #and the mortality rates
    if(times[i]<min(mort$year)){
      mort.df = subset(mort, year==min(year))
    }else{
      mort.df = subset(mort, year==trunc(times[i]))
    }
    #make annual mortality rate by age into probability of mortality per biweek
    #these are deaths per 1000 people.
    #we multiply by the number in the age classes, so need to divide by 1000 here
    mort.df$deaths_per_cap <- mort.df$count/1000
    mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
    
    
    
    #here are the dynamics before the introduction - 2 strain:
    if(times[i]<((yr.intro) + biwk.intro/26)){ 
      
      
      
      Tmat <- buildTMat_age_SIR_triple_two(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi, sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      #we are modeling in 1000s of people, so no need to divide here
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      
      N_pop_ts[,i] <- nt1
      
      
      
      
    }else if (times[i]==((yr.intro) +biwk.intro/26)) { 
      print(paste0("introduction at time =", times[i]))
      #here are the dynamics at the introduction - converting from 2 to 3 strains
      
      #need to add infecteds to I3 disease status -- here at the lowest age class
      #this is row 6!
      N_pop_ts[6,(i-1)] <- 5
      
      #remove them from the Susceptible class
      N_pop_ts[1,(i-1)] <- (N_pop_ts[1,(i-1)]-5)
      
      #then, let her go from there!
      
      Tmat <- buildTMat_age_SIR_triple(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i] > ((yr.intro) +biwk.intro/26)){ 
      #and here are the dynamics after the introduction - 3 strains only
      #no need to introduce anything
      
      Tmat <- buildTMat_age_SIR_triple_after(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
    }
    
  }
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  N.sim.dat$class <- rep(c("S", "I1", "P1", "I2", "P2", "I3", "P3", "I12", "I13", "I21", "I23", "I31", "I32", "Pm"), each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$time <- N.sim.df$time-1
  N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
sim.SIR.age.double.num.intro.clim.switch.wane <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
  #check on climate vector - should be length of steps within a year (often 26)
  #(this could be modulated to be the length of the time series if you wanted)
  if(length(clim.vect)==1){
    clim.vect = rep(clim.vect, ntyr)
  }
  
  
  #first, if mort, births, sigma, or foi are shorter than the time series,
  #make them match in length here.
  
  #this will also repeat if they come in as a vector (age-structured)
  
  #make as a list
  if(length(foi)<yrs){
    foi <- rep(list(foi), yrs)
  }else{
    foi <- as.list(foi)
  }
  
  if(length(sigma)< yrs){
    sigma = rep(list(sigma), yrs)
  }else{
    sigma <- as.list(sigma)
  }
  
  # if(length(mort)<yrs){
  #   mort = rep(list(mort), yrs)
  # }else{
  #   mort = as.list(mort)
  # }
  # 
  if(length(births)<yrs){
    births = rep(list(births), yrs)
  }else{
    births = as.list(births)
  }
  
  
  
  
  
  
  
  
  # if(length(age.mult.foi)==1){
  #   age.mult.foi = rep(age.mult.foi, s)
  # }else if(length(age.mult.foi)>1 & length(age.mult.foi) < s){
  #   age.mult.foi = c(age.mult.foi, rep(age.mult.foi[length(age.mult.foi)], (s-length(age.mult.foi))))
  # }else if(length(age.mult.foi)>s){
  #   age.mult.foi = age.mult.foi[1:s]
  # }
  # 
  # if(length(age.mult.foi)!=s){
  #   print("age multiplier proportions are incorrect")
  # }
  # 
  #otherwise, assume no transformation is needed
  
  #
  #number of epidemic classes
  c=14
  
  
  #length of time series
  times <-   seq(0, yrs, by =1/ntyr)
  times <-   times + (year.end-yrs)
  
  #split the births up by biweek too
  birth_vector_biwk = as.list(c(unlist(lapply(lapply(births, divide, ntyr1 =ntyr), rep, ntyr))))
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  #import the stable age structure from the literature
  #mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  #stab.struct = Re(eigen(mat1)$vector[,1])
  #stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  #lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  pop.mat = pop_vector #number of people per age class
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  
  #structure the population for 3 strains even though we only have 2 to start. they will get added later
  
  PM_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s)# would add some in here initially if you wanted this to be a 3-serotype sim
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, I2_init, P2_init, I3_init, P3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init, PM_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  
  N_pop_ts[,1] <- M_pop # by age class
  
  #and fill in for births <- every epidemic class can give birth
  #transform birth vector into one distributed by age and epidemic class
  
  
  # stab.struct <- get.age.struct.M(M_pop, c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  
  # foi gets spread across the biweeks as well
  # (it comes in as 1 value per year)
  # take your Muench-estimated rates of annual infection and convert
  # to biweekly probabilities of infection
  foi.biweek = as.list(c(unlist(lapply(lapply(foi, divide.rate, ntyr1 =ntyr), rep, ntyr)))) #here, it could easily be modulated to be seasonal - (including forced by precip/temp)
  
  
  #and make annual waning immunity rate into probability of waning immunity per biweek
  sigma.biweek = as.list(c(unlist(lapply(lapply(sigma, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
    print("input vectors of unequal length")
  }
  #and your aging probability to biweeks (aging probability per year is 100%)
  #age.biwk = rep(1-exp(-(rep(1, yrs)/ntyr)), each=26)
  
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 2:(length(times))){
    
    biwk1 <- find.biweek(t=i, times=times)
    clim.mod <- clim.vect[biwk1]
    
    #and select the appropriate age vector, based on the year
    age.mult.sub = subset(age.mult.df, year_min<=times[i] & year_max>=times[i])
    
    
    
    if(nrow(age.mult.sub)>10){
      print("error with age multiplier sub-selection")
    }
    
    age.mult.sub$age_max[length(age.mult.sub$age_max)] <- s
    age.mult.sub$dur <- (age.mult.sub$age_max-age.mult.sub$age_min)+1
    age.mult.foi <- c(unlist(mapply(rep, x=as.list(age.mult.sub$age_mult),each = as.list(age.mult.sub$dur))))
    
    
    #and the mortality rates
    if(times[i]<min(mort$year)){
      mort.df = subset(mort, year==min(year))
    }else{
      mort.df = subset(mort, year==trunc(times[i]))
    }
    #make annual mortality rate by age into probability of mortality per biweek
    #these are deaths per 1000 people.
    #we multiply by the number in the age classes, so need to divide by 1000 here
    mort.df$deaths_per_cap <- mort.df$count/1000
    mort.df$deaths_per_cap_per_biwk <- divide.rate(df=mort.df$deaths_per_cap, ntyr1 = ntyr)
    
    
    
    #here are the dynamics before the introduction - 2 strain:
    if(times[i]<((yr.intro) + biwk.intro/26)){ 
      
      
      
      Tmat <- buildTMat_age_SIR_triple_two(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi, sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      #we are modeling in 1000s of people, so no need to divide here
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      
      N_pop_ts[,i] <- nt1
      
      
      
      
    }else if (times[i]==((yr.intro) +biwk.intro/26)){#} & times[i]<=((yr.intro)+1)){ 
      print(paste0("introduction at time =", times[i]))
      print(paste0("rate of waning immunity=",sigma.biweek[[(i-1)]]))
      #here are the dynamics at the introduction - converting from 2 to 3 strains
      
      #need to add infecteds to I3 disease status -- here at the lowest age class
      #this is row 6!
      N_pop_ts[6,(i-1)] <- 5
      
      #remove them from the Susceptible class
      N_pop_ts[1,(i-1)] <- (N_pop_ts[1,(i-1)]-5)
      
      #then, let her go from there!
      
      Tmat <- buildTMat_age_SIR_triple_switch(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
    }else if(times[i] > ((yr.intro) +biwk.intro/26)){ 
      #and here are the dynamics after the introduction - 3 strains only
      #no need to introduce anything
      
      Tmat <- buildTMat_age_SIR_triple_after(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <- births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
    }
    
  }
  
  
  #the first column is the initial conditions. everything after that represents transitions within the year
  #should take the ceiling of each biweek to total within the year
  
  
  # stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  # stab.struct <- c(unlist(stab.struct))
  # stab.struct <- stab.struct/(sum(stab.struct))
  # #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  if(s>1){
    N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  }
  
  N.sim.dat <- as.data.frame(N_pop_ts)
  
  names(N.sim.dat) <- times
  N.sim.dat$class <- rep(c("S", "I1", "P1", "I2", "P2", "I3", "P3", "I12", "I13", "I21", "I23", "I31", "I32", "Pm"), each = s)
  N.sim.dat$age <- rep(1:s, length(unique(N.sim.dat$class)))
  
  #and melt
  N.sim.df <- melt(N.sim.dat, id=c("class", "age"), value.name = "count", variable.name = "time")
  
  N.sim.df$time <- as.numeric(as.character(N.sim.df$time))
  N.sim.df$count <- as.numeric(as.character(N.sim.df$count))
  
  #and add year
  N.sim.df$year <- ceiling(N.sim.df$time)
  N.sim.df$time <- N.sim.df$time-1
  N.sim.df$year <- N.sim.df$year-1
  
  
  
  return( N.sim.df)
}
check.equil <- function(dat){
  dat.ts <- ddply(dat, .(time, class), summarise, count=sum(count))
  #and get total by time
  dat.N  <- ddply(dat, .(time), summarise, Ntot=sum(count))
  
  dat.ts <- merge(dat.ts, dat.N, by="time")
  dat.ts$proportion <- dat.ts$count/dat.ts$Ntot
  
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=proportion, color=class)) #+ coord_cartesian(ylim=c(0,.1))
  print(p1)
  return(dat.ts)
}
plot.cases.annual <- function(dat, year.start){
  dat1 = subset(dat, year>=year.start)
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(year), summarise, count=sum(count))
  dat.ts$tot_count <- dat.ts$count#/.001
  #dat.ts$tot_count <- dat.ts$count/.001
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ts) + theme_bw() +
        geom_vline(aes(xintercept=2007), linetype=2)+
        geom_vline(aes(xintercept=1990), linetype=2)+
        geom_vline(aes(xintercept=2012), linetype=2)+
        geom_vline(aes(xintercept=2019), linetype=2) +
        geom_line(aes(x=year, y=tot_count)) + 
          theme(panel.grid = element_blank())
  print(p1)
  return(dat.ts)
}
plot.cases.biweek <- function(dat, year.start){
  dat1 = subset(dat, year>=year.start)
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts$time <- dat.ts$time + min(dat$year)
  #dat.ts <- ddply(denv.case, .(year), summarise, count=sum(count))
  dat.ts$tot_count <- dat.ts$count/.001
  #and plot
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count, color=class))
  #p1 <- ggplot(dat.ts) + geom_line(aes(x=time, y=count)) #+ coord_cartesian(ylim=c(0,.1))
  p1 <- ggplot(dat.ts) + theme_bw() +
    geom_vline(aes(xintercept=2007), linetype=2)+
    geom_vline(aes(xintercept=2012), linetype=2)+
    geom_vline(aes(xintercept=2019), linetype=2) +
    geom_line(aes(x=time, y=tot_count)) + 
    theme(panel.grid = element_blank())
  print(p1)
  return(dat.ts)
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  return(df.sum)
}
process.all <- function(df){
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "S", "I1","P1","I2", "P2","I3","I4", "PM")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  #par(mfrow = c(1,1))
  #with(dat.tot, plot(time, N, type="l")) #ylim =c(0,1.2*max(N))))
  #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  
  # prop.tot = dat.tot
  # 
  # prop.tot$S = prop.tot$S/prop.tot$N
  # prop.tot$I = prop.tot$Ii/prop.tot$N
  # prop.tot$R = prop.tot$R/prop.tot$N
  # 
  
  
  
  
  dat.melt = melt(dat.tot,id.vars = "time")
  
  #head(dat.melt)
  
  #take after burnin if we assume this is a virus at equilibrium
  dat.melt <- subset(dat.melt, time>length(burnin_lambda))
  dat.melt <- subset(dat.melt, !is.na(value))
  
  #dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.melt) = c("time", "class", "count")
  dat.N = subset(dat.melt, class =="N")
  dat.melt = subset(dat.melt, class !="N")
  dat.N <- dplyr::select(dat.N, -(class))
  names(dat.N)[names(dat.N)=="count"] <- "Ntot"
  
  dat.melt <- merge(dat.melt, dat.N, by="time")
  #head(dat.melt)
  #dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.melt$class = factor(dat.melt$class, levels=c("S", "I1","P1","I2", "P2","I3","I4", "PM"))
  dat.melt$proportion <- dat.melt$count/dat.melt$Ntot
  
  dat.melt$time <- dat.melt$time - length(burnin_lambda)
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=count, color=class)) #+ scale_color_manual(values=colz)
  
  #dat.melt$proportion = dat.melt$count/sim_pop
  #ggplot(data=dat.melt) + geom_line(aes(x=time, y=proportion, color=class)) #+ scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  #now, get the age-structured incidence - incidence is "I3 + I4"
  
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("S", "I1","P1","I2", "P2","I3","I4", "PM"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  
  #plot age structured cumulative incidence
  age.incidence = subset(age.dat.tot, class=="I3" | class == "I4")
  
  #and collect everyone else
  age.non = subset(age.dat.tot, class!="I3" & class!="I4")
  #head(age.incidence )
  
  age.incidence$year = trunc(age.incidence$times)
  age.non$year = trunc(age.non$times)
  age.incidence <- arrange(age.incidence, year, age)
  
  age.non <- arrange(age.non, year, age)
  
  names(age.non) <- c("times", "count_non_infectious", "class", "age", "year")
  
  age.non <- dplyr::select(age.non, -(class), -(times))
  
  age.non <- ddply(age.non, .(year, age), summarise, count_non_infectious = sum(count_non_infectious))
  
  
  age.sum <- ddply(age.incidence, .(year, age), summarise, count = sum(count))
  
  age.sum <- merge(age.sum, age.non, by=c("year", "age"))
  
  #and just return this counts of years, ages, and those infected and not
  
  
  return(age.sum)#, seas.prev))
  
}

# first, sim Cambodia 2 strain to endemic equilibrium for 50 yrs prior to 1960,
# then sim with birth and death rates (and starting at the correct population distribution) 
# through to 2020. Assess how well this recapitulates what we find in the data


#load climate from TSIR first:
clim.dat <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
head(clim.dat)

#then, get the median beta across epidemic periods and provinces
#we just want this to be an input of seasonality to our transmission rate
beta.med <- ddply(clim.dat,.(biweek), summarise, beta = median(beta))
beta.med$clim_vect <- scales::rescale(beta.med$beta, to=c(.5,1.5), from=c(range(beta.med$beta)))



#load birth and death rates for Cambodia
pop.dat <- read.csv(file = paste0(homewd, "/data/pop_data_full.csv"), header = T, stringsAsFactors = F)


birth.dat = subset(pop.dat, metric=="births per\n1000 ppl")
death.dat = subset(pop.dat, metric == "deaths per\n1000 ppl")


mort.dat <-  read.csv(file = paste0(homewd, "/data/cambodia_age_specific_mort_through_time.csv"), header = T, stringsAsFactors = F)
names(mort.dat) <- c("year", seq(0,100,1))
mort.melt <- melt(mort.dat, id.vars = "year", variable.name = "age", value.name = "count")
mort.melt = arrange(mort.melt, year, age)
mort.melt$age <- as.numeric(as.character(mort.melt$age))
mort.melt$year <- as.numeric(as.character(mort.melt$year))


pop.dist <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dist_through_time.csv"), header = T, stringsAsFactors = F)
head(pop.dist)
names(pop.dist) <- c("year", seq(0,100,1))
pop.melt <- melt(pop.dist, id.vars = "year", variable.name = "age", value.name = "count")
pop.melt = arrange(pop.melt, year, age)
pop.melt$age <- as.numeric(as.character(pop.melt$age))
pop.melt$year <- as.numeric(as.character(pop.melt$year))

#load the foi fits for cambodia
fit.dat <- read.csv(file = paste0(homewd, "/data/prov-fits-FOI.csv"), stringsAsFactors = F, header = T)
nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 60 before this

#load age structure
load(paste0(homewd, "/figure-development/Fig3/comp-fits/fit-many-age-mult.Rdata"))
age.fit$year_min <-   0 #for simulation
age.fit$year_max <-   100 #for simulation
age.fit$year_max[age.fit$year_range=="<=2010"] <- 2010
age.fit$year_min[age.fit$year_range=="2011-2018"] <- 2010.000001
age.fit$year_max[age.fit$year_range=="2011-2018"] <- 2019
age.fit$year_min[age.fit$year_range=="2019-2020"] <- 2019.0000001
age.fit$year_max[age.fit$year_range=="2019-2020"] <- 2021



#first, just sim normal
#sim here, using foi at the National level, but replacing the too-low values
# replace those 0 fois with .05
#fit.dat$uci[fit.dat$provname=="National" &fit.dat$uci<.05] <- .05
fit.dat$uci[fit.dat$provname=="National" &fit.dat$year<1999] <- .9

out.cam.norm = sim.SIR.age.double.num.intro.clim.maintain(yrs=100,
                                            ntyr=26,
                                            s=101, 
                                            foi=c(rep(.2,60),fit.dat$uci[fit.dat$provname=="National"]), 
                                            births =  c(rep(40,39),birth.dat$value), # these are per 1000
                                            pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                            recov=1,
                                            age.mult.df=age.fit, 
                                            clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                            age.brk=NA,#this is a placeholder for additional age structure
                                            mort=mort.melt, 
                                            sigma=0, #nothing at first
                                            yr.intro = 10000,
                                            biwk.intro = 10000,
                                            year.end = 2021)

# plot.cases.annual(dat=out.cam.norm, year.start=1920)
# check.equil.cases.only(dat=out.cam.norm)
# plot.age.dist.three(dat=out.cam.norm, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.norm, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum(out.cam.norm, year.start = 2002)
out.cam.norm$sim_type <- "normal"
out.cam.norm$hyp <- 1

foi.up <- c(rep(.2,60),fit.dat$uci[fit.dat$provname=="National"])
foi.up[length(foi.up)-1] <- foi.up[length(foi.up)-1]*1.15

out.cam.upFOI = sim.SIR.age.double.num.intro.clim.maintain(yrs=100,
                                                          ntyr=26,
                                                          s=101, 
                                                          foi=foi.up, 
                                                          births =  c(rep(40,39),birth.dat$value), # these are per 1000
                                                          pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                          recov=1,
                                                          age.mult.df=age.fit, 
                                                          clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                          age.brk=NA,#this is a placeholder for additional age structure
                                                          mort=mort.melt, 
                                                          sigma=0, #nothing at first
                                                          yr.intro = 10000,
                                                          biwk.intro = 10000,
                                                          year.end = 2021)

# check.equil.cases.only(dat=out.cam.upFOI)
# plot.cases.annual(dat=out.cam.upFOI, year.start=1920)
# plot.age.dist.three(dat=out.cam.upFOI, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.upFOI, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum(out.cam.upFOI, year.start = 2002)
out.cam.upFOI$sim_type <- "high_FOI"
out.cam.upFOI$hyp <- 2

#now try with introduction maintained
out.cam.intro.maintain = sim.SIR.age.double.num.intro.clim.maintain(yrs=100,
                                                                ntyr=26,
                                                                s=101, 
                                                                foi=c(rep(.2,60),fit.dat$uci[fit.dat$provname=="National"]), 
                                                                births = c(rep(40,39),birth.dat$value), # these are per 1000
                                                                pop_vector = (pop.melt$count[pop.melt$year==1950]),
                                                                recov=1,
                                                                age.mult.df=age.fit, #nothing here for now... can add in if desired. provides age structure on foi
                                                                clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                age.brk=NA,#this is a placeholder for additional age structure
                                                                mort=mort.melt, 
                                                                sigma=0, #nothing at first
                                                                yr.intro = 2019,
                                                                biwk.intro = 1,
                                                                year.end = 2021)


# check.equil.cases.only(dat=out.cam.intro.maintain)
# plot.cases.annual(dat=out.cam.intro.maintain, year.start=1920)
# plot.age.dist.three(dat=out.cam.intro.maintain, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.intro.maintain, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum(out.cam.intro.maintain, year.start = 2002)
out.cam.intro.maintain$sim_type <- "intro_maintain"
out.cam.intro.maintain$hyp <- 3

#now try with introduction and serotype switch (replacement but no waning immunity)
out.cam.intro.switch = sim.SIR.age.double.num.intro.clim.switch(yrs=100,
                                                  ntyr=26,
                                                  s=101, 
                                                  foi=c(rep(.2,60),fit.dat$uci[fit.dat$provname=="National"]), 
                                                  births = c(rep(40,39),birth.dat$value), # these are per 1000
                                                  pop_vector = (pop.melt$count[pop.melt$year==1950]),
                                                  recov=1,
                                                  age.mult.df=age.fit, #nothing here for now... can add in if desired. provides age structure on foi
                                                  clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                  age.brk=NA,#this is a placeholder for additional age structure
                                                  mort=mort.melt, 
                                                  sigma=0, #nothing at first
                                                  yr.intro = 2019,
                                                  biwk.intro = 1,
                                                  year.end = 2021)

# check.equil.cases.only(dat=out.cam.intro.switch)
# plot.cases.annual(dat=out.cam.intro.switch, year.start=1920)
# plot.age.dist.three(dat=out.cam.intro.switch, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.intro.switch, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum(out.cam.intro.switch, year.start = 2002)
out.cam.intro.switch$sim_type <- "intro_switch"
out.cam.intro.switch$hyp <- 4

out.cam.intro.switch.wane = sim.SIR.age.double.num.intro.clim.switch.wane(yrs=100,
                                                                ntyr=26,
                                                                s=101, 
                                                                foi=c(rep(.2,60),fit.dat$uci[fit.dat$provname=="National"]), 
                                                                births = c(rep(40,39),birth.dat$value), # these are per 1000
                                                                pop_vector = (pop.melt$count[pop.melt$year==1950]),
                                                                recov=1,
                                                                age.mult.df=age.fit, #nothing here for now... can add in if desired. provides age structure on foi
                                                                clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                age.brk=NA,#this is a placeholder for additional age structure
                                                                mort=mort.melt, 
                                                                sigma=c(rep(0,98), .01, rep(0,1)), #nothing at first
                                                                yr.intro = 2019,
                                                                biwk.intro = 1,
                                                                year.end = 2021)

out.cam.intro.switch.wane$sim_type <- "intro_switch_wane"
out.cam.intro.switch.wane$hyp <- 5
# check.equil.cases.only(dat=out.cam.intro.switch.wane)
# plot.cases.annual(dat=out.cam.intro.switch.wane, year.start=1920)
# plot.age.dist.three(dat=out.cam.intro.switch.wane, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.intro.switch.wane, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum(out.cam.intro.switch.wane, year.start = 2002)

#bind all together and save, then plot
cam.sim.uci <- rbind(out.cam.norm, out.cam.upFOI, out.cam.intro.maintain, out.cam.intro.switch, out.cam.intro.switch.wane)

save(cam.sim.uci, file = "cam-sim-uci.Rdata")




