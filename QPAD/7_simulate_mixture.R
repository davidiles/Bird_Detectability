# **************************************
# Repeated simulations across different values of tau and phi
# **************************************

library(tidyverse)
library(ggpubr)
library(detect)
library(ggthemes)

rm(list=ls())

setwd("~/1_Work/Bird_Detectability/QPAD") # <- set to wherever scripts are stored
source("model_fns3.R")

result_df <- expand.grid(sim_rep = 1:100,
                         
                         # Simulation parameters
                         tau_A = 1,
                         phi_A = 0.5,
                         tau_B = 0.8,
                         phi_B = 1.5,
                         pA = 1,
                         Density = 1000
                         
)

for (i in 1:nrow(result_df)){
  
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  tau_A <- result_df$tau_A[i]
  phi_A <- result_df$phi_A[i]
  tau_B <- result_df$tau_B[i]
  phi_B <- result_df$phi_B[i]
  pA <- result_df$pA[i]
  
  set.seed(sim_rep)
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  rint <- c(0.5,1,Inf)
  tint <- seq(1,10)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  
  Yarray <- array(NA,dim=c(1,nrint,ntint))
  rarray <- array(NA,dim=c(1,nrint))
  tarray <- array(NA,dim=c(1,ntint))
  
  rarray[1,1:length(rint)] <- rint
  tarray[1,1:length(tint)] <- tint
  
  Y_A <- Y_B <- 0
  
  # *********************************************************************
  # Simulate observations for bird type A
  # *********************************************************************
  if (pA >0){
    dim <- 10
    N <- round(Density*pA*dim^2) # Number of birds to place on landscape
    
    birds <- data.frame(x = runif(N,-dim/2,dim/2),
                        y = runif(N,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    
    # Remove birds outside maximum distance
    birds <- subset(birds, dist <= (dim/2))
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    N = nrow(birds)
    birds$dist <- birds$dist
    
    # ------------------------------------
    # Simulate bird cues, based on phi_true
    # ------------------------------------
    
    cues <- matrix(NA, nrow=N, ncol = 50)
    for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_A))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    # ------------------------------------
    # Determine which cues are detected, based on tau_true
    # ------------------------------------
    
    cues$p<- exp(-(cues$dist/tau_A)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Isolate first detected EVENT for each bird (bird type = A)
    # ------------------------------------
    
    dat <- subset(cues,detected == 1 & cue_number == 1)
    dat <- dat[!duplicated(dat$bird_id),]
    
    # ------------------------------------
    # Transcription: distance and time bins
    # ------------------------------------
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y_A <- table(dat[,c("rint","tint")])
    Y_A # Data to analyze
  }
  
  # *********************************************************************
  # Simulate observations for bird type A
  # *********************************************************************
  
  if (pA < 1){
    dim <- 10
    N <- round(Density*(1-pA)*dim^2) # Number of birds to place on landscape
    
    birds <- data.frame(x = runif(N,-dim/2,dim/2),
                        y = runif(N,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    
    # Remove birds outside maximum distance
    birds <- subset(birds, dist <= (dim/2))
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    N = nrow(birds)
    birds$dist <- birds$dist
    
    # ------------------------------------
    # Simulate bird cues, based on phi_true
    # ------------------------------------
    
    cues <- matrix(NA, nrow=N, ncol = 50)
    for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_B))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    # ------------------------------------
    # Determine which cues are detected, based on tau_true
    # ------------------------------------
    
    cues$p<- exp(-(cues$dist/tau_B)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Isolate first detected cue for each bird
    # ------------------------------------
    
    dat <- subset(cues,detected == 1)
    dat <- dat[!duplicated(dat$bird_id),]
    
    # ------------------------------------
    # Transcription: distance and time bins
    # ------------------------------------
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y_B <- table(dat[,c("rint","tint")])
    Y_B # Data to analyze
  }
  
  # *********************************************************************
  # Combine Y_A and Y_B
  # *********************************************************************
  
  Y <- Y_A + Y_B
  Yarray[1,1:nrint,1:ntint] <- Y
  
  Ysum <- sum(Yarray)
  result_df$Ysum[i] <- Ysum
  
  if (Ysum == 0) next # Species was never detected at any point counts
  
  # *********************************************************************
  # Fit mixture model
  # *********************************************************************
  
  start <- Sys.time()
  fit_PP <- fit_PP_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0
  )
  end <- Sys.time()
  print(end-start)
  
  tau_est <- exp(fit_PP$coefficients[1])
  phi_est <- exp(fit_PP$coefficients[2])
  gamma_est <- fit_PP$coefficients[3]
  
  offsets_PP <- offsets_PP_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              
                              tau = rep(tau_est,dim(Yarray)[1]),
                              phi = rep(phi_est,dim(Yarray)[1]),
                              gamma = gamma_est
  )
  
  # Estimates
  result_df$tau_est[i] <- tau_est
  result_df$phi_est[i] <- phi_est
  result_df$gamma_est[i] <- gamma_est
  result_df$offset[i] <- offsets_PP$offset_vec[1]
  result_df$Density_est[i] <- Ysum/offsets_PP$offset_vec[1]
  
  print(i)
  
}

result_summary <- result_df %>%
  group_by(tau_A,phi_A,tau_B,phi_B,pA,Density) %>%
  summarize_all(list(mean = mean,med = median))
result_summary %>% as.data.frame()
