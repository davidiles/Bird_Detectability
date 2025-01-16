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
                         time_bins = c(0.5,1,2,5),
                         Density = 1000,
                         
                         tau_est = NA,
                         phi_est = NA,
                         p_est = NA,
                         D_est = NA,
                         D_true = NA) %>%
  arrange(sim_rep,time_bins)

for (i in 1:nrow(result_df)){
  
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  tau <- result_df$tau[i]
  time_bins <- result_df$time_bins[i]
  
  dat <- data.frame()
  
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
  
  dim <- 10
  
  N <- round(Density*dim^2) # Number of birds to place on landscape
  
  # -------------------------------------------------------
  # Birds with fast cue rate
  # -------------------------------------------------------
  phi_1 <- 20
  birds <- data.frame(x = runif(N*0.5,-dim/2,dim/2),
                      y = runif(N*0.5,-dim/2,dim/2))
  
  # Distances to observer
  birds$dist <- sqrt(birds$x^2 + birds$y^2)
  
  # Remove birds outside maximum distance
  birds <- subset(birds, dist <= (dim/2))
  birds <- birds %>% arrange(dist)
  birds$bird_id = 1:nrow(birds)
  n = nrow(birds)
  birds$dist <- birds$dist
  
  cues <- matrix(NA, nrow=n, ncol = 50)
  for (bird_id in 1:n) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_1))
  cues <- cues %>% 
    reshape2::melt() %>% 
    rename(bird_id = Var1, cue_number = Var2, time = value) %>%
    arrange(bird_id,cue_number)
  
  cues$dist <- birds$dist[cues$bird_id]
  
  cues$p <- exp(-(cues$dist/tau)^2)  # Probability each cue is detected
  cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
  
  cues <- subset(cues,detected == 1 & cue_number == 1)
  cues <- cues[!duplicated(cues$bird_id),]
  
  dat <- rbind(dat, cues)
  
  # -------------------------------------------------------
  # Birds with slow cue rate
  # -------------------------------------------------------
  phi_2 <- 0.1
  birds <- data.frame(x = runif(N*0.5,-dim/2,dim/2),
                      y = runif(N*0.5,-dim/2,dim/2))
  
  # Distances to observer
  birds$dist <- sqrt(birds$x^2 + birds$y^2)
  
  # Remove birds outside maximum distance
  birds <- subset(birds, dist <= (dim/2))
  birds <- birds %>% arrange(dist)
  birds$bird_id = 1:nrow(birds)
  n = nrow(birds)
  birds$dist <- birds$dist
  
  cues <- matrix(NA, nrow=n, ncol = 50)
  for (bird_id in 1:n) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_2))
  cues <- cues %>% 
    reshape2::melt() %>% 
    rename(bird_id = Var1, cue_number = Var2, time = value) %>%
    arrange(bird_id,cue_number)
  
  cues$dist <- birds$dist[cues$bird_id]
  
  cues$p <- exp(-(cues$dist/tau)^2)  # Probability each cue is detected
  cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
  
  cues <- subset(cues,detected == 1 & cue_number == 1) # 
  cues <- cues[!duplicated(cues$bird_id),]
  
  dat <- rbind(dat, cues)
  
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
  
  # Estimate density
  A_est <- pi*tau_est^2 # area sampled by unlimited distance survey
  result_df$p_est[i] <- 1-exp(-10*phi_est) # 10 minute survey
  result_df$D_est[i] <- Ysum/(A_est*result_df$p_est[i])
  print(i)
  
  true_p <- 1-0.5*exp(-10*phi_2)
  
  
  result_summary <- result_df %>%
    group_by(time_bins) %>%
    summarize(p_est = mean(p_est,na.rm = TRUE),
              D_est = mean(D_est,na.rm = TRUE))
  
  if ((sim_rep %% 10) == 0) {
    
    
    result_plot_p <- ggplot(data = result_df)+
      geom_line(aes(x = time_bins, y = p_est, col = factor(sim_rep)))+
      scale_color_manual(values = rep("gray80",length(unique(result_df$sim_rep))), guide = "none")+
      geom_line(data = result_summary, aes(x = time_bins, y = p_est), size = 2, col = "dodgerblue")+
      geom_hline(yintercept = true_p, linetype = 2)+
      xlab("Time bin length")+
      ylab("Estimated p")+
      theme_bw()
    
    print(result_plot_p)
    
    # 
    # result_plot_D <- ggplot(data = result_df)+
    #   geom_line(aes(x = time_bins, y = D_est, col = factor(sim_rep)))+
    #   scale_color_manual(values = rep("gray80",length(unique(result_df$sim_rep))), guide = "none")+
    #   geom_line(data = result_summary, aes(x = time_bins, y = D_est), size = 2, col = "dodgerblue")+
    #   geom_hline(yintercept = mean(result_df$Density), linetype = 2)+
    #   xlab("Time bin length")+
    #   ylab("Estimated Density")+
    #   theme_bw()
    # 
    # print(result_plot_D)
    
  }
  
}

