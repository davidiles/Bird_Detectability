# **************************************
# Repeated simulations across different values of tau and phi
# **************************************

# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
# -------------------------------------------------------------

rm(list=ls())

result_df <- expand.grid(sim_rep = 1:100,
                         tau = 1,
                         phi = 0.5,
                         time_bins = c(1,5),
                         Density = 100,
                         tau_est = NA,
                         phi_est = NA) %>%
  arrange(sim_rep,time_bins)

for (i in 1:nrow(result_df)){
  
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  tau <- result_df$tau[i]
  phi <- result_df$phi[i]
  time_bins <- result_df$time_bins[i]
  
  set.seed(sim_rep)
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  rint <- c(0.5,1,Inf)
  tint <- seq(time_bins,10,by = time_bins)
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
  
  Y <- 0
  
  # *********************************************************************
  # Simulate observations
  # *********************************************************************
  
  dim <- 10
  
  N <- round(Density*dim^2) # Number of birds to place on landscape
  
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
  cue_prob <- plogis(scale(cumsum(rnorm(1000))))
  cues <- matrix(NA, nrow=N, ncol = 50)
  for (bird_id in 1:N){
    cues_potential <- cumsum(rexp(500,phi))
    cues_produced <- sort(rbinom(500,1,cue_prob) * cues_potential)
    cues_produced <- cues_produced[cues_produced>0][1:50]
    cues[bird_id,] <- cues_produced
  }
  
  cues <- cues %>% 
    reshape2::melt() %>% 
    rename(bird_id = Var1, cue_number = Var2, time = value) %>%
    arrange(bird_id,cue_number)
  
  cues$dist <- birds$dist[cues$bird_id]
  
  # ------------------------------------
  # Determine which cues are detected, based on tau_true
  # ------------------------------------
  
  cues$p <- exp(-(cues$dist/tau)^2)  # Probability each cue is detected
  cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
  
  # ------------------------------------
  # Isolate first detected EVENT for each bird (bird type = A)
  # ------------------------------------
  
  dat <- subset(cues,detected == 1)# & cue_number == 1)
  dat <- dat[!duplicated(dat$bird_id),]
  
  # ------------------------------------
  # Transcription: distance and time bins
  # ------------------------------------
  
  # Separate into distance and time bins
  dat$rint <- cut(dat$dist,c(0,rint))
  dat$tint <- cut(dat$time,c(0,tint))
  dat <- na.omit(dat)
  
  Y <- table(dat[,c("rint","tint")])
  Y # Data to analyze
  
  Ysum <- sum(Y)
  result_df$Ysum[i] <- Ysum
  
  if (Ysum == 0) next # Species was never detected at any point counts
  
  # ------------------------------------
  # Analysis: removal model
  # ------------------------------------
  
  Y_removal <- matrix(colSums(Y[,]),1)
  D_removal <- matrix(tint,1)
  nbins_removal <- length(D_removal)
  fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_est = exp(fit_p$coefficients)
  
  # ------------------------------------
  # Analysis: distance model
  # ------------------------------------
  
  Y_distance <- matrix(rowSums(Y[,]),1)
  D_distance <- matrix(rint,1)
  nbins_distance <- length(D_distance)
  fit_q <- detect::cmulti.fit(Y_distance,D_distance, type = "dis")
  tau_est = exp(fit_q$coefficients)
  
  result_df$phi_est[i] <- phi_est
  result_df$tau_est[i] <- tau_est
  
  print(i)
  

  if ((sim_rep %% 10) == 0) {
    ylim=c(0,1)
    result_plot <- ggplot(data = result_df)+
      geom_line(aes(x = time_bins, y = phi_est, col = factor(sim_rep)))+
      coord_cartesian(ylim=ylim)+
      scale_color_manual(values=rep("black",length(unique(result_df$sim_rep))), guide = "none")+
      xlab("Time Bin Duration")+
      ylab("Estimated Phi")+
      theme_bw()
    
    print(result_plot)
  }
  
  
}

