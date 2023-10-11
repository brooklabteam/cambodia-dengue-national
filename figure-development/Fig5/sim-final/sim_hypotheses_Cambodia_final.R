# here is the age-structured 3-strain model, including heterotypic immunity

rm(list=ls())

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig5/sim-final/"))
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
#need to modify pre
buildTMat_three_strain_pre <- function(c, Npop, age.classes, surv.biwk, age.brk,wane_hetero,  foi, age.mult.foi, recov, sigma, age.rate){
  #no transmission with the third strain
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

  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
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
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (2*foi[j])
    mat1[2,1] <- foi[j]
    mat1[5,1] <- foi[j]
    #mat1[8,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    #I2s recover
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov
    
    #I3s recover
    #mat1[8,8] <- 1-recov
    #mat1[9,8] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I1
    mat1[3,3] <- 1- wane_hetero[j]
    mat1[4,3] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[6,6] <- 1- wane_hetero[j]
    mat1[7,6] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    #mat1[9,9] <- 1- wane_hetero[j]
    #mat1[10,9] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[4,4] <- 1-(foi[j])
    mat1[11,4] <- foi[j]
    #mat1[12,4] <- foi[j]
    
    #secondary infections - those previously infected with I2
    mat1[7,7] <- 1-(foi[j])
    mat1[13,7] <- foi[j]
    #mat1[14,7] <- foi[j]
    
    #secondary infections - those previously infected with I3
    #mat1[10,10] <-  1-(2*foi[j])
    #mat1[15,10] <- foi[j]
    #mat1[16,10] <- foi[j]
    
    #recovery from secondary infections - I12
    mat1[11,11] <- 1-recov
    mat1[17,11] <- recov
    
    #recovery from secondary infections - I13
   # mat1[12,12] <- 1-recov
    #mat1[18,12] <- recov
    
    #recovery from secondary infections - I21
    mat1[13,13] <- 1-recov
    mat1[19,13] <- recov
    
    #recovery from secondary infections - I23
    #mat1[14,14] <- 1-recov
    #mat1[20,14] <- recov
    
    #recovery from secondary infections - I31
    #mat1[15,15] <- 1-recov
    #mat1[21,15] <- recov
    
    #recovery from secondary infections - I32
    #mat1[16,16] <- 1-recov
    #mat1[22,16] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I12
    mat1[17,17] <- 1- wane_hetero[j]
    mat1[23,17] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I13
    #mat1[18,18] <- 1- wane_hetero[j]
    #mat1[24,18] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    mat1[19,19] <- 1- wane_hetero[j]
    mat1[25,19] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
   # mat1[20,20] <- 1- wane_hetero[j]
    #mat1[26,20] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    #mat1[21,21] <- 1- wane_hetero[j]
    #mat1[27,21] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    #mat1[22,22] <- 1- wane_hetero[j]
    #mat1[28,22] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12
    #here, they move straight to multi-typic immunity 
    mat1[23,23] <- 1
    #mat1[36,23] <- foi[j]
    
    #tertiary infections - those previously infected with I13
    #mat1[24,24] <- 1-foi[j]
   # mat1[30,24] <- foi[j]
    
    
    #tertiary infections - those previously infected with I21
    #here, they just stay here
    mat1[25,25] <- 1
    #mat1[31,25] <- foi[j]
    
    #tertiary infections - those previously infected with I23
   # mat1[26,26] <- 1-foi[j]
    #mat1[32,26] <- foi[j]
    
    
    #tertiary infections - those previously infected with I31
   # mat1[27,27] <- 1-foi[j]
   # mat1[33,27] <- foi[j]
    
    #tertiary infections - those previously infected with I32
    #mat1[28,28] <- 1-foi[j]
    #mat1[34,28] <- foi[j]
    
    
    #recovery from tertiary infections - I123
   # mat1[29,29] <- 1-recov
    #mat1[35,29] <- recov
    
    # #recovery from tertiary infections - I132
    # mat1[30,30] <- 1-recov
    # mat1[35,30] <- recov
    # 
    # #recovery from tertiary infections - I213
    # mat1[31,31] <- 1-recov
    # mat1[35,31] <- recov
    # 
    # #recovery from tertiary infections - I231
    # mat1[32,32] <- 1-recov
    # mat1[35,32] <- recov
    # 
    # #recovery from tertiary infections - I312
    # mat1[33,33] <- 1-recov
    # mat1[35,33] <- recov
    # 
    # #recovery from tertiary infections - I321
    # mat1[34,34] <- 1-recov
    # mat1[35,34] <- recov
    # 
    
    #waning from tertiary immunity back to homotypic immunity (renewed secondary infections)
    #no waning in this sim
    #mat1[3,35] <- sigma[j]
    #mat1[6,35] <- sigma[j]
    #mat1[9,35] <- sigma[j]
    #mat1[35,35] <- 1-(2*sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_three_strain_pre_wane <- function(c, Npop, age.classes, surv.biwk, age.brk,wane_hetero,  foi, age.mult.foi, recov, sigma, age.rate){
  #no transmission with the third strain
  #waning immunity for just a single serotype - here, serotype 1
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
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
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
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (2*foi[j])
    mat1[2,1] <- foi[j]
    mat1[5,1] <- foi[j]
    #mat1[8,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    #I2s recover
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov
    
    #I3s recover
    #mat1[8,8] <- 1-recov
    #mat1[9,8] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I1
    mat1[3,3] <- 1- wane_hetero[j]
    mat1[4,3] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[6,6] <- 1- wane_hetero[j]
    mat1[7,6] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    #mat1[9,9] <- 1- wane_hetero[j]
    #mat1[10,9] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[4,4] <- 1-(foi[j])
    mat1[11,4] <- foi[j]
    #mat1[12,4] <- foi[j]
    
    #secondary infections - those previously infected with I2
    mat1[7,7] <- 1-(foi[j])
    mat1[13,7] <- foi[j]
    #mat1[14,7] <- foi[j]
    
    #secondary infections - those previously infected with I3
    #mat1[10,10] <-  1-(2*foi[j])
    #mat1[15,10] <- foi[j]
    #mat1[16,10] <- foi[j]
    
    #recovery from secondary infections - I12
    mat1[11,11] <- 1-recov
    mat1[17,11] <- recov
    
    #recovery from secondary infections - I13
    # mat1[12,12] <- 1-recov
    #mat1[18,12] <- recov
    
    #recovery from secondary infections - I21
    mat1[13,13] <- 1-recov
    mat1[19,13] <- recov
    
    #recovery from secondary infections - I23
    #mat1[14,14] <- 1-recov
    #mat1[20,14] <- recov
    
    #recovery from secondary infections - I31
    #mat1[15,15] <- 1-recov
    #mat1[21,15] <- recov
    
    #recovery from secondary infections - I32
    #mat1[16,16] <- 1-recov
    #mat1[22,16] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I12
    #here they go straight to multi-typic immunity and can wane back to P1 or P2
    mat1[17,17] <- 1- wane_hetero[j]
    mat1[35,17] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I13
    #mat1[18,18] <- 1- wane_hetero[j]
    #mat1[24,18] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    #here they go straight to multi-typic immunity and can wane back to P1 and P2
    mat1[19,19] <- 1- wane_hetero[j]
    mat1[35,19] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    # mat1[20,20] <- 1- wane_hetero[j]
    #mat1[26,20] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    #mat1[21,21] <- 1- wane_hetero[j]
    #mat1[27,21] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    #mat1[22,22] <- 1- wane_hetero[j]
    #mat1[28,22] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12
    
    #mat1[23,23] <- 1
    #mat1[36,23] <- foi[j]
    
    #tertiary infections - those previously infected with I13
    #mat1[24,24] <- 1-foi[j]
    # mat1[30,24] <- foi[j]
    
    
    #tertiary infections - those previously infected with I21
    
    #mat1[25,25] <- 1
    #mat1[31,25] <- foi[j]
    
    #tertiary infections - those previously infected with I23
    # mat1[26,26] <- 1-foi[j]
    #mat1[32,26] <- foi[j]
    
    
    #tertiary infections - those previously infected with I31
    # mat1[27,27] <- 1-foi[j]
    # mat1[33,27] <- foi[j]
    
    #tertiary infections - those previously infected with I32
    #mat1[28,28] <- 1-foi[j]
    #mat1[34,28] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    # mat1[29,29] <- 1-recov
    #mat1[35,29] <- recov
    
    # #recovery from tertiary infections - I132
    # mat1[30,30] <- 1-recov
    # mat1[35,30] <- recov
    # 
    # #recovery from tertiary infections - I213
    # mat1[31,31] <- 1-recov
    # mat1[35,31] <- recov
    # 
    # #recovery from tertiary infections - I231
    # mat1[32,32] <- 1-recov
    # mat1[35,32] <- recov
    # 
    # #recovery from tertiary infections - I312
    # mat1[33,33] <- 1-recov
    # mat1[35,33] <- recov
    # 
    # #recovery from tertiary infections - I321
    # mat1[34,34] <- 1-recov
    # mat1[35,34] <- recov
    # 
    
    #waning from tertiary immunity back to homotypic immunity (renewed secondary infections)
    # here only one genotype wanes (genotype 1)
    #mat1[3,35] <- sigma[j] #back to P1
    mat1[6,35] <- sigma[j] # back to P2 only - new P1 strain can now invade and cause new infections in P2 or in P3
    #mat1[9,35] <- sigma[j]
    mat1[35,35] <- 1-(sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_three_strain <- function(c, Npop, age.classes, surv.biwk, wane_hetero, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  
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
  
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
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
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (3*foi[j])
    mat1[2,1] <- foi[j]
    mat1[5,1] <- foi[j]
    mat1[8,1] <- foi[j]
    
    #I1s recover
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    #I2s recover
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov
    
    #I3s recover
    mat1[8,8] <- 1-recov
    mat1[9,8] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I1
    mat1[3,3] <- 1- wane_hetero[j]
    mat1[4,3] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[6,6] <- 1- wane_hetero[j]
    mat1[7,6] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    mat1[9,9] <- 1- wane_hetero[j]
    mat1[10,9] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    mat1[4,4] <- 1-(2*foi[j])
    mat1[11,4] <- foi[j]
    mat1[12,4] <- foi[j]
    
    #secondary infections - those previously infected with I2
    mat1[7,7] <- 1-(2*foi[j])
    mat1[13,7] <- foi[j]
    mat1[14,7] <- foi[j]
    
    #secondary infections - those previously infected with I3
    mat1[10,10] <-  1-(2*foi[j])
    mat1[15,10] <- foi[j]
    mat1[16,10] <- foi[j]
    
    #recovery from secondary infections - I12
    mat1[11,11] <- 1-recov
    mat1[17,11] <- recov
    
    #recovery from secondary infections - I13
    mat1[12,12] <- 1-recov
    mat1[18,12] <- recov
    
    #recovery from secondary infections - I21
    mat1[13,13] <- 1-recov
    mat1[19,13] <- recov
    
    #recovery from secondary infections - I23
    mat1[14,14] <- 1-recov
    mat1[20,14] <- recov
    
    #recovery from secondary infections - I31
    mat1[15,15] <- 1-recov
    mat1[21,15] <- recov
    
    #recovery from secondary infections - I32
    mat1[16,16] <- 1-recov
    mat1[22,16] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I12
    mat1[17,17] <- 1- wane_hetero[j]
    mat1[23,17] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I13
    mat1[18,18] <- 1- wane_hetero[j]
    mat1[24,18] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    mat1[19,19] <- 1- wane_hetero[j]
    mat1[25,19] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    mat1[20,20] <- 1- wane_hetero[j]
    mat1[26,20] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    mat1[21,21] <- 1- wane_hetero[j]
    mat1[27,21] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    mat1[22,22] <- 1- wane_hetero[j]
    mat1[28,22] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12
    mat1[23,23] <- 1-foi[j]
    mat1[29,23] <- foi[j]
    
    #tertiary infections - those previously infected with I13
    mat1[24,24] <- 1-foi[j]
    mat1[30,24] <- foi[j]
    
    
    #tertiary infections - those previously infected with I21
    mat1[25,25] <- 1-foi[j]
    mat1[31,25] <- foi[j]
    
    #tertiary infections - those previously infected with I23
    mat1[26,26] <- 1-foi[j]
    mat1[32,26] <- foi[j]
    
    
    #tertiary infections - those previously infected with I31
    mat1[27,27] <- 1-foi[j]
    mat1[33,27] <- foi[j]
    
    #tertiary infections - those previously infected with I32
    mat1[28,28] <- 1-foi[j]
    mat1[34,28] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    mat1[29,29] <- 1-recov
    mat1[35,29] <- recov
    
    #recovery from tertiary infections - I132
    mat1[30,30] <- 1-recov
    mat1[35,30] <- recov
    
    #recovery from tertiary infections - I213
    mat1[31,31] <- 1-recov
    mat1[35,31] <- recov
    
    #recovery from tertiary infections - I231
    mat1[32,32] <- 1-recov
    mat1[35,32] <- recov
    
    #recovery from tertiary infections - I312
    mat1[33,33] <- 1-recov
    mat1[35,33] <- recov
    
    #recovery from tertiary infections - I321
    mat1[34,34] <- 1-recov
    mat1[35,34] <- recov
    
    
    #waning from tertiary immunity back to homotypic immunity (renewed secondary infections)
    mat1[3,35] <- sigma[j]
    mat1[6,35] <- sigma[j]
    mat1[9,35] <- sigma[j]
    mat1[35,35] <- 1-(3*sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
    }
  }
  
  return(Tmat)
  
}
buildTMat_three_strain_after <- function(c, Npop, age.classes, surv.biwk, wane_hetero, age.brk, foi, age.mult.foi, recov, sigma, age.rate){
  #now eliminate transmission from genotype 1 and allow for waning
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
  
  
  if (length(wane_hetero)==1) wane_hetero <- rep(wane_hetero,nage) # 
  if (length(wane_hetero)==s) wane_hetero <- wane_hetero
  if (length(wane_hetero) > 1 & length(wane_hetero)<s){
    wane_hetero.list <- list()
    age.dur = diff(c(age.brk,s))
    if(sum(age.dur)!=s){print("age breaks don't add up")}
    for (i in 1:length(wane_hetero)){
      wane_hetero.list[[i]] = rep(wane_hetero[i], age.dur[i])
    }
    wane_hetero = c(unlist(wane_hetero.list))
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
  
  
  
  mat1 <- matrix(0,c,c) 
  
  Tmat <- matrix(0,c*nage,c*nage) 
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek
    mat1[] <- 0
    
    #infections in naive susceptibles
    mat1[1,1] <- 1- (2*foi[j])
    #mat1[2,1] <- foi[j]
    mat1[5,1] <- foi[j]
    mat1[8,1] <- foi[j]
    
    #I1s recover -still allowed
    mat1[2,2] <- 1-recov
    mat1[3,2] <- recov
    
    
    #I2s recover
    mat1[5,5] <- 1-recov
    mat1[6,5] <- recov
    
    #I3s recover
    mat1[8,8] <- 1-recov
    mat1[9,8] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I1
    #still allowed
    mat1[3,3] <- 1- wane_hetero[j]
    mat1[4,3] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I2
    mat1[6,6] <- 1- wane_hetero[j]
    mat1[7,6] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I3
    mat1[9,9] <- 1- wane_hetero[j]
    mat1[10,9] <- wane_hetero[j]
    
    #secondary infections - those previously infected with I1
    #still allowed
    mat1[4,4] <- 1-(2*foi[j])
    mat1[11,4] <- foi[j]
    mat1[12,4] <- foi[j]
    
    #secondary infections - those previously infected with I2
    #no more I1
    mat1[7,7] <- 1-(foi[j])
    #mat1[13,7] <- foi[j]
    mat1[14,7] <- foi[j]
    
    #secondary infections - those previously infected with I3
    #no more I1
    mat1[10,10] <-  1-(foi[j])
    #mat1[15,10] <- foi[j]
    mat1[16,10] <- foi[j]
    
    #recovery from secondary infections - I12
    #still allowed
    mat1[11,11] <- 1-recov
    mat1[17,11] <- recov
    
    #recovery from secondary infections - I13
    #still allowed
    mat1[12,12] <- 1-recov
    mat1[18,12] <- recov
    
    #recovery from secondary infections - I21
    #still allowed
    mat1[13,13] <- 1-recov
    mat1[19,13] <- recov
    
    #recovery from secondary infections - I23
    #still allowed
    mat1[14,14] <- 1-recov
    mat1[20,14] <- recov
    
    #recovery from secondary infections - I31
    #still allowed
    mat1[15,15] <- 1-recov
    mat1[21,15] <- recov
    
    #recovery from secondary infections - I32
    #still allowed
    mat1[16,16] <- 1-recov
    mat1[22,16] <- recov
    
    
    #waning out of heterotypic immunity - those previously infected with I12
    #still allowed
    mat1[17,17] <- 1- wane_hetero[j]
    mat1[23,17] <- wane_hetero[j]
    
    
    #waning out of heterotypic immunity - those previously infected with I13
    #still allowed
    mat1[18,18] <- 1- wane_hetero[j]
    mat1[24,18] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I21
    #still allowed
    mat1[19,19] <- 1- wane_hetero[j]
    mat1[25,19] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I23
    #still allowed...but they progress straight to multitypic immunity
    mat1[20,20] <- 1- wane_hetero[j]
    mat1[35,20] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I31
    #still allowed
    mat1[21,21] <- 1- wane_hetero[j]
    mat1[27,21] <- wane_hetero[j]
    
    #waning out of heterotypic immunity - those previously infected with I32
    #still allowed...but they progress straight to multitypic immunity
    mat1[22,22] <- 1- wane_hetero[j]
    mat1[35,22] <- wane_hetero[j]
    
    
    #tertiary infections - those previously infected with I12
    #still allowed
    mat1[23,23] <- 1-foi[j]
    mat1[29,23] <- foi[j]
    
    #tertiary infections - those previously infected with I13
    #still allowed
    mat1[24,24] <- 1-foi[j]
    mat1[30,24] <- foi[j]
    
    
    #tertiary infections - those previously infected with I21
    #still allowed
    mat1[25,25] <- 1-foi[j]
    mat1[31,25] <- foi[j]
    
    #tertiary infections - those previously infected with I23
    #no more - these move to Pm
    mat1[26,26] <- 1-recov
    mat1[35,26] <- recov
    
    
    #tertiary infections - those previously infected with I31
    #still allowed
    mat1[27,27] <- 1-foi[j]
    mat1[33,27] <- foi[j]
    
    #tertiary infections - those previously infected with I32
    #no more - these move to Pm
    mat1[28,28] <- 1-foi[j]
    mat1[35,28] <- foi[j]
    
    
    #recovery from tertiary infections - I123
    #still allowed
    mat1[29,29] <- 1-recov
    mat1[35,29] <- recov
    
    #recovery from tertiary infections - I132
    #still allowed
    mat1[30,30] <- 1-recov
    mat1[35,30] <- recov
    
    #recovery from tertiary infections - I213
    #still allowed
    mat1[31,31] <- 1-recov
    mat1[35,31] <- recov
    
    #recovery from tertiary infections - I231
    #still allowed
    mat1[32,32] <- 1-recov
    mat1[35,32] <- recov
    
    #recovery from tertiary infections - I312
    #still allowed
    mat1[33,33] <- 1-recov
    mat1[35,33] <- recov
    
    #recovery from tertiary infections - I321
    #still allowed
    mat1[34,34] <- 1-recov
    mat1[35,34] <- recov
    
    
    #waning from tertiary immunity back to homotypic immunity (renewed secondary infections)
    #all still possible
    mat1[3,35] <- sigma[j]
    mat1[6,35] <- sigma[j]
    mat1[9,35] <- sigma[j]
    mat1[35,35] <- 1-(3*sigma[j])
    
    
    #put in surv. we first say it is equal across all infectious classes
    
    
    
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class, based on the rate of aging
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]*age.rate[j] #this is transitioning from age class j to j+1
      
      #and some stay in the same age class, but still undergo infection transitions and survival
      Tmat[((j-1)*c+1):((j-1)*c+c),((j-1)*c+1):((j-1)*c+c)] <- mat1*surv.biwk[j]*(1-age.rate[j]) #this is staying within age class j
      
      
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv.biwk[j]
      
      
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
replicate.data.type <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    new.dat$type <- unique(df$type)
    
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
plot.age.dist.three.triple.type <- function(dat, count.type, perc.obs, save.plot, view.plot, filename, slim.quant, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    denv.case.tert$type <- "tertiary"  
    denv.case.tert$count <- denv.case.tert$count*perc.obs
    denv.case <- rbind( denv.case,  denv.case.tert)
    
    
    
  }
  
  # 
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(type, year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  df.sum <- df.sum[complete.cases(df.sum),]
  
  #split by a year
  df.year <- dlply(df.sum,.(year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and type
  df.age <- dlply(df.sum,.(type, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type, slim.quant=slim.quant))
  
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  colz = c("secondary"="black", "tertiary"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + #facet_grid(~type) +
    geom_jitter(aes(x=year, y=age, color=type), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + scale_color_manual(values=colz) +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
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
plot.age.dist.three.triple.type.detect <- function(dat, count.type, perc.obs, save.plot, view.plot, filename, slim.quant, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    denv.case.tert$type <- "tertiary"  
    #detectability increases by year
    mod.df <- cbind.data.frame(year=unique(dat1$year), perc_obs =seq(0,perc.obs, length.out= length(unique(dat1$year))))
    denv.case.tert <- merge(denv.case.tert, mod.df, by="year")
    denv.case.tert$count <- denv.case.tert$count*denv.case.tert$perc_obs
    denv.case.tert <- dplyr::select(denv.case.tert, -(perc_obs))
    denv.case <- rbind( denv.case,  denv.case.tert)
    
  }
  
  # 
  #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(type, year,age),summarise, count=sum(count))
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  
  
  df.sum <- df.sum[complete.cases(df.sum),]
  
  #split by a year
  df.year <- dlply(df.sum,.(year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and type
  df.age <- dlply(df.sum,.(type, year, age))
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age, replicate.data.type, slim.quant=slim.quant))
  
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  colz = c("secondary"="black", "tertiary"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + 
    geom_jitter(aes(x=year, y=age, color=type), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    geom_violin(aes(x=year,y=age, group=year), color="gray60",  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato") + scale_color_manual(values=colz) +
    geom_hline(aes(yintercept=1), color="red") + coord_cartesian(ylim=c(0,80))#,xlim=c(2015,2020))
  
  
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
process.plot.age.cum.triple.type <-function(dat, count.type, perc.obs, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    denv.case.tert$type <- "tertiary"  
    denv.case.tert$count <- denv.case.tert$count*perc.obs
    denv.case <- rbind( denv.case,  denv.case.tert)
    
  }
  
  
  
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
  denv.case = subset(dat,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  denv.case$serotype <- NA
  denv.case$serotype[denv.case$class=="I12" | denv.case$class=="I32" | denv.case$class=="I132" | denv.case$class=="I312" ] <- "I2"
  denv.case$serotype[denv.case$class=="I13" | denv.case$class=="I23" | denv.case$class=="I123" | denv.case$class=="I213"] <- "I3"
  denv.case$serotype[denv.case$class=="I31" | denv.case$class=="I21"| denv.case$class=="I321" | denv.case$class=="I231"] <- "I1"
  
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
sim.SIR.age.triple.num.intro.clim.switch <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
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
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=35
  
  
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
  I321_init = rep(0, s)
  I312_init = rep(0, s)
  I231_init = rep(0, s)
  I213_init = rep(0, s)
  I132_init = rep(0, s)
  I123_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s)# would add some in here initially if you wanted this to be a 3-serotype sim
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, S1_init, I2_init, P2_init, S2_init, I3_init, P3_init, S3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init,P12_init, P13_init, P21_init, P23_init, P31_init, P32_init, S12_init, S13_init, S21_init, S23_init, S31_init, S32_init, I123_init, I132_init, I213_init, I231_init, I312_init, I321_init, PM_init)
  
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
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
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
      
      
      
      Tmat <- buildTMat_three_strain_pre(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi, sigma=sigma.biweek[[(i-1)]], recov=recov,  wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
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
      
      Tmat <- buildTMat_three_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
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
      
      #no need to introduce anything
      #here, we drop transmission for strain on2
      
      Tmat <- buildTMat_three_strain_after(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)), wane_hetero=hetero.biweek[[i-1]])
      
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
  
  N.sim.dat$class <- rep(c("S", "I1", "P1", "S1", "I2", "P2", "S2", "I3", "P3", "S3", "I12", "I13", "I21", "I23", "I31", "I32", "P12", "P13", "P21", "P23", "P31", "P32", "S12", "S13", "S21", "S23", "S31", "S32", "I123", "I132", "I213", "I231", "I312", "I321", "Pm"), each = s)
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
sim.SIR.age.triple.num.intro.clim.maintain <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
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
  
  # 
  if(length(rate_hetero)<yrs){
  rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=35
  
  
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
  I321_init = rep(0, s)
  I312_init = rep(0, s)
  I231_init = rep(0, s)
  I213_init = rep(0, s)
  I132_init = rep(0, s)
  I123_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s)# would add some in here initially if you wanted this to be a 3-serotype sim
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, S1_init, I2_init, P2_init, S2_init, I3_init, P3_init, S3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init,P12_init, P13_init, P21_init, P23_init, P31_init, P32_init, S12_init, S13_init, S21_init, S23_init, S31_init, S32_init, I123_init, I132_init, I213_init, I231_init, I312_init, I321_init, PM_init)
  
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
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
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
    age.mult.sub$dur[age.mult.sub$dur<0] <- 0
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
      
      
      
      Tmat <- buildTMat_three_strain_pre(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi, sigma=sigma.biweek[[(i-1)]], recov=recov,  wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
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
      
      Tmat <- buildTMat_three_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
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
      
      #no need to introduce anything
      #here, we drop transmission for strain on2
      
      Tmat <- buildTMat_three_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov,  age.rate=(1-exp(-(1)/ntyr)), wane_hetero=hetero.biweek[[i-1]])
      
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
  
  N.sim.dat$class <- rep(c("S", "I1", "P1", "S1", "I2", "P2", "S2", "I3", "P3", "S3", "I12", "I13", "I21", "I23", "I31", "I32", "P12", "P13", "P21", "P23", "P31", "P32", "S12", "S13", "S21", "S23", "S31", "S32", "I123", "I132", "I213", "I231", "I312", "I321", "Pm"), each = s)
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
sim.SIR.age.triple.num.intro.clim.endemic <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma){
  
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
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=35
  
  
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
  I321_init = rep(0, s)
  I312_init = rep(0, s)
  I231_init = rep(0, s)
  I213_init = rep(0, s)
  I132_init = rep(0, s)
  I123_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s); I3_init[1] = 5 
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, S1_init, I2_init, P2_init, S2_init, I3_init, P3_init, S3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init,P12_init, P13_init, P21_init, P23_init, P31_init, P32_init, S12_init, S13_init, S21_init, S23_init, S31_init, S32_init, I123_init, I132_init, I213_init, I231_init, I312_init, I321_init, PM_init)
  
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
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
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
    age.mult.sub$dur[age.mult.sub$dur<0] <- 0
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
    
    
      
      #then, let her go from there!
      
      Tmat <- buildTMat_three_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
      
   
    
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
  
  N.sim.dat$class <- rep(c("S", "I1", "P1", "S1", "I2", "P2", "S2", "I3", "P3", "S3", "I12", "I13", "I21", "I23", "I31", "I32", "P12", "P13", "P21", "P23", "P31", "P32", "S12", "S13", "S21", "S23", "S31", "S32", "I123", "I132", "I213", "I231", "I312", "I321", "Pm"), each = s)
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
sim.SIR.age.triple.num.intro.clim.two.wane <- function(yrs, ntyr, age.brk, s,foi, age.mult.df, recov, mort, year.end, rate_hetero, clim.vect, births, pop_vector, sigma, yr.intro, biwk.intro){
  
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
  
  # 
  if(length(rate_hetero)<yrs){
    rate_hetero = rep(list(rate_hetero), yrs)
  }else{
    rate_hetero = as.list(rate_hetero)
  }
  
  
  #
  #number of epidemic classes
  c=35
  
  
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
  I321_init = rep(0, s)
  I312_init = rep(0, s)
  I231_init = rep(0, s)
  I213_init = rep(0, s)
  I132_init = rep(0, s)
  I123_init = rep(0, s)
  S32_init = rep(0, s)
  S31_init = rep(0, s)
  S23_init = rep(0, s)
  S21_init = rep(0, s)
  S13_init = rep(0, s)
  S12_init = rep(0, s)
  P32_init = rep(0, s)
  P31_init = rep(0, s)
  P23_init = rep(0, s)
  P21_init = rep(0, s)
  P13_init = rep(0, s)
  P12_init = rep(0, s)
  I32_init = rep(0, s)
  I31_init = rep(0, s)
  I23_init = rep(0, s)
  I21_init = rep(0, s)
  I13_init = rep(0, s)
  I12_init = rep(0, s)
  S3_init = rep(0, s)
  S2_init = rep(0, s)
  S1_init = rep(0, s)
  P3_init = rep(0, s)
  P2_init = rep(0, s)
  P1_init = rep(0, s)
  I3_init = rep(0, s); I3_init[1] = 5 
  I2_init = rep(0, s); I2_init[1] = 5 #comment out if you just want to check demography
  I1_init = rep(0, s); I1_init[1] = 5 #comment out if you just want to check demography
  S_init = pop.mat - I1_init - I2_init 
  
  
  N_tot = cbind(S_init, I1_init, P1_init, S1_init, I2_init, P2_init, S2_init, I3_init, P3_init, S3_init, I12_init, I13_init, I21_init, I23_init, I31_init, I32_init,P12_init, P13_init, P21_init, P23_init, P31_init, P32_init, S12_init, S13_init, S21_init, S23_init, S31_init, S32_init, I123_init, I132_init, I213_init, I231_init, I312_init, I321_init, PM_init)
  
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
  
  #and the probability of waning heterotypic immunity per biweek
  hetero.biweek = as.list(c(unlist(lapply(lapply(rate_hetero, divide.rate, ntyr1 =ntyr), rep, ntyr))))
  
  #length(mort.biweek)!=(length(times)-1) |
  if( length(sigma.biweek)!=(length(times)-1) | length(hetero.biweek)!=(length(times)-1) |length(birth_vector_biwk)!=(length(times)-1)| length(foi.biweek)!=(length(times)-1)){
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
    age.mult.sub$dur[age.mult.sub$dur<0] <- 0
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
    
    
    if(times[i]<((yr.intro) + biwk.intro/26)){ 
   
    #waning should be 1 year prior to introduction
    Tmat <- buildTMat_three_strain_pre_wane(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
    
    #Tmat is infection and survival.
    transMat <- Tmat  
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,(i-1)]
    
    #then, add in births into the 0 age class of the susceptibles 
    births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
    
    
    births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
    
    
    nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
    
    N_pop_ts[,i] <- nt1
    
    }else if (times[i]==((yr.intro) + biwk.intro/26)){
      print(paste0("introduction at time =", times[i]))
      #here are the dynamics at the introduction - converting from 2 to 3 strains
      
      #need to add infecteds to I3 disease status -- here at the lowest age class
      #this is row 6!
      N_pop_ts[6,(i-1)] <- 5
      
      #remove them from the Susceptible class
      N_pop_ts[1,(i-1)] <- (N_pop_ts[1,(i-1)]-5)
      
      #then, let her go!
      
      Tmat <- buildTMat_three_strain(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
      nt1[1,1] <- nt1[1,1] + births_add_biweek #births get added into the lowest age class of susceptibles
      
      N_pop_ts[,i] <- nt1
      
    }else if(times[i]>((yr.intro) + biwk.intro/26)){
      Tmat <- buildTMat_three_strain_after(c=c, Npop= N_pop_ts[,(i-1)], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort.df$deaths_per_cap_per_biwk),	foi = (foi.biweek[[(i-1)]]*clim.mod), age.mult.foi =age.mult.foi,  sigma=sigma.biweek[[(i-1)]], recov=recov, wane_hetero=hetero.biweek[[i-1]], age.rate=(1-exp(-(1)/ntyr)))
      
      #Tmat is infection and survival.
      transMat <- Tmat  
      
      #move forward in time
      nt1<-(transMat) %*% N_pop_ts[,(i-1)]
      
      #then, add in births into the 0 age class of the susceptibles 
      births_per_1000pop = birth_vector_biwk[[(i-1)]] #these are each biweek per 1000 people
      
      
      births_add_biweek <-births_per_1000pop*((sum(N_pop_ts[,(i-1)]))/1000) #fill in susceptible births into class 0
      
      
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
  
  N.sim.dat$class <- rep(c("S", "I1", "P1", "S1", "I2", "P2", "S2", "I3", "P3", "S3", "I12", "I13", "I21", "I23", "I31", "I32", "P12", "P13", "P21", "P23", "P31", "P32", "S12", "S13", "S21", "S23", "S31", "S32", "I123", "I132", "I213", "I231", "I312", "I321", "Pm"), each = s)
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
plot.cases.annual.triple.type <- function(dat, perc.obs, count.type, year.start){
  dat1 = subset(dat, year>=year.start)

  
  #select only those viewed as "cases"
  denv.case = subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32" | class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  denv.case$type = "secondary"
  
  if(count.type=="tertiary"){
    
    denv.case.tert = subset(dat1, class == "I123" | class=="I132" | class=="I213"  | class=="I231"  | class=="I312"  | class=="I321")
    if(length(unique(denv.case.tert$age))>1){
      denv.case.tert = subset(denv.case.tert, age<max(denv.case$age))  
    }
    
    denv.case.tert$type <- "tertiary"  
    denv.case.tert$count <- denv.case.tert$count*perc.obs
    denv.case <- rbind( denv.case,  denv.case.tert)
    
  }
  
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
clim.dat <- read.csv(file = "beta_TSIR_fit_province.csv", header = T, stringsAsFactors = F)
head(clim.dat)

#then, get the median beta across epidemic periods and provinces
#we just want this to be an input of seasonality to our transmission rate
beta.med <- ddply(clim.dat,.(biweek), summarise, beta = median(beta))
beta.med$clim_vect <- scales::rescale(beta.med$beta, to=c(.5,1.5), from=c(range(beta.med$beta)))



#load birth and death rates for Cambodia
pop.dat <- read.csv(file = "pop_data_full.csv", header = T, stringsAsFactors = F)


birth.dat = subset(pop.dat, metric=="births per\n1000 ppl")
death.dat = subset(pop.dat, metric == "deaths per\n1000 ppl")


mort.dat <-  read.csv(file = "cambodia_age_specific_mort_through_time.csv", header = T, stringsAsFactors = F)
names(mort.dat) <- c("year", seq(0,100,1))
mort.melt <- melt(mort.dat, id.vars = "year", variable.name = "age", value.name = "count")
mort.melt = arrange(mort.melt, year, age)
mort.melt$age <- as.numeric(as.character(mort.melt$age))
mort.melt$year <- as.numeric(as.character(mort.melt$year))


pop.dist <- read.csv(file = "cambodia_pop_dist_through_time.csv", header = T, stringsAsFactors = F)
head(pop.dist)
names(pop.dist) <- c("year", seq(0,100,1))
pop.melt <- melt(pop.dist, id.vars = "year", variable.name = "age", value.name = "count")
pop.melt = arrange(pop.melt, year, age)
pop.melt$age <- as.numeric(as.character(pop.melt$age))
pop.melt$year <- as.numeric(as.character(pop.melt$year))



#load the foi fits for cambodia
fit.dat <- read.csv(file = "prov-fits-FOI.csv", stringsAsFactors = F, header = T)
nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 60 before this

#load age structure
load("fit-many-age-mult.Rdata")
age.fit$year_min <-   0 #for simulation
age.fit$year_max <-   100 #for simulation
age.fit$year_max[age.fit$year_range=="<=2010"] <- 2010
age.fit$year_min[age.fit$year_range=="2011-2018"] <- 2010.000001
age.fit$year_max[age.fit$year_range=="2011-2018"] <- 2019
age.fit$year_min[age.fit$year_range=="2019-2020"] <- 2019.0000001
age.fit$year_max[age.fit$year_range=="2019-2020"] <- 2021

#first, just sim normal
#sim here, using foi at the National level, but replacing the too-low values
fit.dat$lambda[fit.dat$provname=="National" &fit.dat$year<1999] <- .9


# H0
out.cam.norm = sim.SIR.age.triple.num.intro.clim.maintain(yrs=80,
                                             ntyr=26,
                                             s=101, 
                                             foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                             births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                             pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                             recov=1,
                                             age.mult.df=age.fit, 
                                             clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                             age.brk=NA,#this is a placeholder for additional age structure
                                             mort=mort.melt, 
                                             rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                             sigma=0, #nothing at first
                                             yr.intro = 10000,
                                             biwk.intro = 10000,
                                             year.end = 2021)

 
#check.equil(out.cam.norm)
#plot.cases.annual.triple.type(dat=out.cam.norm, year.start=1920, perc.obs=0.05, count.type = "secondary")
#check.equil.cases.only(dat=out.cam.norm)
#plot.age.dist.three.triple.type(dat=out.cam.norm, count.type = "secondary", perc.obs = 0.05, save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
out.cam.norm$sim_type <- "normal-demography"
out.cam.norm$hyp <- "0"


#H4 3 endemic strains and tertiary detection (can also plot increasing detectability with time)
out.cam.intro.3strains = sim.SIR.age.triple.num.intro.clim.endemic(yrs=80,
                                                                    ntyr=26,
                                                                    s=101, 
                                                                    foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                                                    births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                                    pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                                    recov=1,
                                                                    age.mult.df=age.fit, 
                                                                    clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                    age.brk=NA,#this is a placeholder for additional age structure
                                                                    mort=mort.melt, 
                                                                    rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                                    sigma=0, #nothing at first
                                                                    year.end = 2021)

# check.equil(dat=out.cam.intro.3strains)
# check.equil.cases.only(dat=out.cam.intro.3strains)
# # plot.cases.annual.triple.type(dat=out.cam.intro.3strains, count.type = "secondary", perc.obs = 0.05, year.start=1920)
# # plot.cases.annual.triple.type(dat=out.cam.intro.3strains, count.type = "tertiary", perc.obs = 0.05, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains, save.plot=F, view.plot=T, filename=NA,count.type = "secondary", perc.obs = NA, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = .5, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = 1, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type.detect(dat=out.cam.intro.3strains, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = 1, slim.quant=1, year.start = 1920)
# process.plot.age.cum.triple.type(dat=out.cam.intro.3strains, year.start = 2002, count.type = "tertiary", perc.obs = 0.05)
#plot.ages.all(dat=out.cam.intro.3strains, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
out.cam.intro.3strains$sim_type <- "3_strains_full_time_series"
out.cam.intro.3strains$hyp <- "4"

# under three strains, we are increasing the overall FOI - 
# so likely need to decrease the FOI under these assumptions
# regardless, if secondary infections only are counted, 
# this is no different than above. If tertiary infections are counted
# it is too high in age earlier in the time series

#you need the increasing detection - but why would that be?
# and how do you get single year spikes?

# #try with lower FOI
# new.foi <- (((c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]))*2)/3)
# 
# # 3 endemic strains and tertiary detection (can also plot increasing detectability with time)
# out.cam.intro.3strains.low = sim.SIR.age.triple.num.intro.clim.endemic(yrs=80,
#                                                                    ntyr=26,
#                                                                    s=101, 
#                                                                    foi=new.foi, 
#                                                                    births =  c(rep(40,19),birth.dat$value), # these are per 1000
#                                                                    pop_vector =(pop.melt$count[pop.melt$year==1950]), 
#                                                                    recov=1,
#                                                                    age.mult.df=age.fit, 
#                                                                    clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
#                                                                    age.brk=NA,#this is a placeholder for additional age structure
#                                                                    mort=mort.melt, 
#                                                                    rate_hetero = 1/2,#duration of heterotypic immunity (in years)
#                                                                    sigma=0, #nothing at first
#                                                                    year.end = 2021)
# 
# 
# check.equil(dat=out.cam.intro.3strains.low)
# check.equil.cases.only(dat=out.cam.intro.3strains.low)
# # plot.cases.annual.triple.type(dat=out.cam.intro.3strains.low, count.type = "secondary", perc.obs = 0.05, year.start=1920)
# # plot.cases.annual.triple.type(dat=out.cam.intro.3strains.low, count.type = "tertiary", perc.obs = 0.05, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains.low, save.plot=F, view.plot=T, filename=NA,count.type = "secondary", perc.obs = NA, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains.low, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = .5, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.3strains.low, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = 1, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type.detect(dat=out.cam.intro.3strains.low, save.plot=F, view.plot=T, filename=NA,count.type = "tertiary", perc.obs = 1, slim.quant=1, year.start = 1920)
# process.plot.age.cum.triple.type(dat=out.cam.intro.3strains.low, year.start = 2002, count.type = "tertiary", perc.obs = 0.05)
# #plot.ages.all(dat=out.cam.intro.3strains.low, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# out.cam.intro.3strains.low$sim_type <- "3_strains_full_time_series_low_foi"
# out.cam.intro.3strains.low$hyp <- 2



# two strains with a novel serotype invasion
out.cam.intro.maintain = sim.SIR.age.triple.num.intro.clim.maintain(yrs=80,
                                                                    ntyr=26,
                                                                    s=101, 
                                                                    foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                                                    births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                                    pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                                    recov=1,
                                                                    age.mult.df=age.fit, 
                                                                    clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                    age.brk=NA,#this is a placeholder for additional age structure
                                                                    mort=mort.melt, 
                                                                    rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                                    sigma=0, #nothing at first
                                                                    yr.intro = 2019,
                                                                    biwk.intro = 1,
                                                                    year.end = 2021)

# check.equil.cases.only.triple(dat=out.cam.intro.maintain)
# plot.cases.annual.triple(dat=out.cam.intro.maintain, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.maintain,count.type = "secondary", perc.obs = 0.05,  save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.maintain,count.type = "tertiary", perc.obs = 0.5,  save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.ages.all(dat=out.cam.intro.maintain, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum.triple(out.cam.intro.maintain, year.start = 2002)
out.cam.intro.maintain$sim_type <- "intro_maintain_tertiary"
out.cam.intro.maintain$hyp <- "3"


#this can get what we want...very sudden leap in age distribution, but
#unclear how we lose the spike in ages
#also no genomic evidence to suggest this happened



# two strains with a novel serotype invasion - but make it earlier to understand the 
# longer term consequences
out.cam.intro.maintain.mid = sim.SIR.age.triple.num.intro.clim.maintain(yrs=80,
                                                                        ntyr=26,
                                                                        s=101, 
                                                                        foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                                                        births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                                        pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                                        recov=1,
                                                                        age.mult.df=age.fit, 
                                                                        clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                        age.brk=NA,#this is a placeholder for additional age structure
                                                                        mort=mort.melt, 
                                                                        rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                                        sigma=0, #nothing at first
                                                                        yr.intro = 2007,
                                                                        biwk.intro = 1,
                                                                        year.end = 2021)

# check.equil.cases.only.triple(dat=out.cam.intro.maintain.mid)
# plot.cases.annual.triple(dat=out.cam.intro.maintain.mid, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.maintain.mid,count.type = "secondary", perc.obs = 0.05,  save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.maintain.mid,count.type = "tertiary", perc.obs = 0.5,  save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
# plot.ages.all(dat=out.cam.intro.maintain.mid, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum.triple(out.cam.intro.maintain.mid, year.start = 2002)
out.cam.intro.maintain.mid$sim_type <- "intro_maintain_tertiary_mid_season"
out.cam.intro.maintain.mid$hyp <- "H3-2007"

# big spike in age, followed by decline in mean age as 3 get maintained in the system
# we don't see that. might be able to explain if old serotype was eliminated
# again, no evidence of this



#finally, with waning in the intro year and just two strains
#e.g. genotype replacement with waning immunity
out.cam.intro.switch.wane = sim.SIR.age.triple.num.intro.clim.two.wane(yrs=80,
                                                                       ntyr=26,
                                                                       s=101, 
                                                                       foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                                                       births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                                       pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                                       recov=1,
                                                                       age.mult.df=age.fit, 
                                                                       clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                       age.brk=NA,#this is a placeholder for additional age structure
                                                                       mort=mort.melt, 
                                                                       rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                                       sigma=c(rep(0,77),.01,0,0), #nothing at first
                                                                       yr.intro = 2018,
                                                                       biwk.intro = 1,
                                                                       year.end = 2021)

# check.equil.cases.only.triple(dat=out.cam.intro.switch.wane)
# plot.cases.annual.triple(dat=out.cam.intro.switch.wane, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.switch.wane,count.type = "secondary", perc.obs = 0.05,  save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.switch.wane,count.type = "tertiary", perc.obs =1,  save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# #plot.ages.all(dat=out.cam.intro.switch.wane, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum.triple.type(out.cam.intro.switch.wane,count.type = "secondary", perc.obs = 0.05, year.start = 2002)
# process.plot.age.cum.triple.type(out.cam.intro.switch.wane,count.type = "tertiary", perc.obs = 1, year.start = 2002)
out.cam.intro.switch.wane$sim_type <- "intro_switch_wane"
out.cam.intro.switch.wane$hyp <- "2"


#also show in the middle


#finally, with waning in the intro year and just two strains
#e.g. genotype replacement with waning immunity
out.cam.intro.switch.wane.mid = sim.SIR.age.triple.num.intro.clim.two.wane( yrs=80,
                                                                            ntyr=26,
                                                                            s=101, 
                                                                            foi=c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"]), 
                                                                            births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                                            pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                                            recov=1,
                                                                            age.mult.df=age.fit, 
                                                                            clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                                            age.brk=NA,#this is a placeholder for additional age structure
                                                                            mort=mort.melt, 
                                                                            rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                                            sigma=c(rep(0,66),.01,rep(0,13)), #nothing at first
                                                                            yr.intro = 2007,
                                                                            biwk.intro = 1,
                                                                            year.end = 2021)

# check.equil.cases.only.triple(dat=out.cam.intro.switch.wane.mid)
# plot.cases.annual.triple(dat=out.cam.intro.switch.wane.mid, year.start=1920)
# plot.age.dist.three.triple.type(dat=out.cam.intro.switch.wane.mid,count.type = "secondary", perc.obs = 0.05,  save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
# #plot.age.dist.three.triple.type(dat=out.cam.intro.switch.wane.mid,count.type = "tertiary", perc.obs =1,  save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
# #plot.ages.all(dat=out.cam.intro.switch.wane.mid, save.plot=F, view.plot=T, filename=NA, slim.quant=.001, year.start = 1920)
# process.plot.age.cum.triple.type(out.cam.intro.switch.wane.mid, count.type = "secondary", perc.obs = 0.05, year.start = 2002)
out.cam.intro.switch.wane.mid$sim_type <- "intro_switch_wane_mid_time_series"
out.cam.intro.switch.wane.mid$hyp <- "H2-2007"

#yes, you get the spike in the middle - best if you only look at secondary infections in this case



#and high FOI in epidemic year

foi.up <- c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"])
foi.up[length(foi.up)-1] <- foi.up[length(foi.up)-1]*1.1

out.cam.FOIup = sim.SIR.age.triple.num.intro.clim.maintain(yrs=80,
                                                          ntyr=26,
                                                          s=101, 
                                                          foi=foi.up, 
                                                          births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                          pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                          recov=1,
                                                          age.mult.df=age.fit, 
                                                          clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                          age.brk=NA,#this is a placeholder for additional age structure
                                                          mort=mort.melt, 
                                                          rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                          sigma=0, #nothing at first
                                                          yr.intro = 10000,
                                                          biwk.intro = 10000,
                                                          year.end = 2021)


#check.equil(out.cam.norm)
#plot.cases.annual.triple.type(dat=out.cam.norm, year.start=1920, perc.obs=0.05, count.type = "secondary")
#check.equil.cases.only(dat=out.cam.norm)
#plot.age.dist.three.triple.type(dat=out.cam.norm, count.type = "secondary", perc.obs = 0.05, save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
out.cam.FOIup$sim_type <- "high-FOI"
out.cam.FOIup$hyp <- "1"


#and in 2007
foi.up.2007 <- c(rep(.9,40),fit.dat$lambda[fit.dat$provname=="National"])
foi.up.2007[length(foi.up.2007)-13] <- foi.up.2007[length(foi.up.2007)-13]*1.1

out.cam.FOIupmid = sim.SIR.age.triple.num.intro.clim.maintain(yrs=80,
                                                           ntyr=26,
                                                           s=101, 
                                                           foi=foi.up.2007, 
                                                           births =  c(rep(40,19),birth.dat$value), # these are per 1000
                                                           pop_vector =(pop.melt$count[pop.melt$year==1950]), 
                                                           recov=1,
                                                           age.mult.df=age.fit, 
                                                           clim.vect = beta.med$clim_vect, #26 entries, one per biwk, scaled from 0.5-1.5
                                                           age.brk=NA,#this is a placeholder for additional age structure
                                                           mort=mort.melt, 
                                                           rate_hetero = 1/2,#duration of heterotypic immunity (in years)
                                                           sigma=0, #nothing at first
                                                           yr.intro = 10000,
                                                           biwk.intro = 10000,
                                                           year.end = 2021)


#check.equil(out.cam.norm)
#plot.cases.annual.triple.type(dat=out.cam.norm, year.start=1920, perc.obs=0.05, count.type = "secondary")
#check.equil.cases.only(dat=out.cam.norm)
#plot.age.dist.three.triple.type(dat=out.cam.norm, count.type = "secondary", perc.obs = 0.05, save.plot=F, view.plot=T, filename=NA, slim.quant=1, year.start = 1920)
out.cam.FOIupmid$sim_type <- "high-FOI"
out.cam.FOIupmid$hyp <- "H1-2007"








#could get more of a spike if waning immunity lasted longer probably - or if it came a few years before?

#no need to test genotype replacement with no waning immunity since it should
#not look any different than a 2-strain model





#bind all together and save, then plot
cam.sim <- rbind(out.cam.norm, out.cam.intro.3strains, out.cam.intro.maintain, out.cam.intro.maintain.mid, out.cam.intro.switch.wane, out.cam.intro.switch.wane.mid, out.cam.FOIup, out.cam.FOIupmid)
save(cam.sim, file = "cam-sim-final-10-10.Rdata")

