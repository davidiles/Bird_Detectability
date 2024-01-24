# **************************************
# Evaluate sensitivity of QPAD to mixtures of birds with multiple values of phi and tau
#  e.g., birds with high and low cue rates
# **************************************

library(tidyverse)
library(ggpubr)
library(detect)
library(ggthemes)
library(ggpubr)

rm(list=ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Bird_Detectability/Analysis/script")

result_df <- expand.grid(sim_rep = seq(1,100),
                         
                         # Cue rates of birds in each group
                         phi_A = 1,
                         phi_B = c(0.1,1,10),
                         
                         # EDR of birds in each group
                         tau_A = 1,
                         tau_B = c(0.5,1,2),
                         
                         # Proportion of birds in each group
                         pA = c(0,0.5,1),
                         Density = 500,
                         
                         Density_est = NA
                         
) %>% arrange(sim_rep,phi_B)

for (i in 1:nrow(result_df)){
  print(i)
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  tau_A <- result_df$tau_A[i]
  phi_A <- result_df$phi_A[i]
  tau_B <- result_df$tau_B[i]
  phi_B <- result_df$phi_B[i]
  pA <- result_df$pA[i]
  det_type <- result_df$det_type[i]
  
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
    
    cues$p <- exp(-(cues$dist/tau_A)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Is first cue detected every time?
    # Is independence assumption met
    # ------------------------------------
    
    dat <- subset(cues,detected == 1) #  & cue_number == 1
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
  # Simulate observations for bird type B 
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
    # Is first cue detected every time?
    # ------------------------------------
    
    dat <- subset(cues,detected == 1) #  & cue_number == 1
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
  # Fit QPAD and estimate density
  # *********************************************************************
  
  # Distance model
  Y_distance = matrix(rowSums(Y),1)
  D_distance = matrix(rint,1)
  nbins_distance <- length(Y_distance)
  
  # Removal model
  Y_removal = matrix(colSums(Y),1)
  D_removal = matrix(tint,1)
  nbins_removal <- length(Y_removal)
  
  # ------------------------------
  # Analysis using conditional maximum likelihood
  # ------------------------------
  
  # Trycatch in case of 'singular' error
  fits <- 0
  fits <- tryCatch({
    fit.q <- detect::cmulti.fit(Y_distance,D_distance, type = "dis")
    tau_MLE = exp(fit.q$coefficients)
    
    fit.p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
    phi_MLE = exp(fit.p$coefficients)
    
    1
  },
  error = function(cond){0})
  
  if (fits == 0) next
  
  # Estimate density
  A_hat = pi*tau_MLE^2
  p_hat = 1-exp(-max(tint)*phi_MLE)
  D_hat <- sum(Y)/(A_hat*p_hat)
  
  # Save results
  result_df$tau_est[i] <- tau_MLE
  result_df$phi_est[i] <- phi_MLE
  result_df$Density_est[i] <- D_hat
  
  # *********************************************************************
  # Plot results
  # *********************************************************************

  if ((i %% 100)==0){
  result_df$pA_label <- paste0("pA = ",result_df$pA)
  result_df$tau_B_label <- paste0("tau_B = ",result_df$tau_B)
  
  result_plot1 <- ggplot(result_df)+

    geom_hline(aes(yintercept = Density), linetype = 2)+

    geom_boxplot(aes(x = factor(phi_B), y = Density_est), fill = "gray", size = 1)+

    ylab("Population Density\n\n(Birds/1000km^2)")+
    xlab("phi_B")+
    ggtitle(paste0("Simulation Results\nphi_A = ",phi_A))+

    facet_grid(pA_label~tau_B_label)+
    theme_bw()

  print(result_plot1)
  
  }
  
  
}
