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
                         time_bins = 1, #c(0.5,1,2,5),
                         Density = 1,
                         D_0 = NA,
                         D_1 = NA) %>%
  arrange(sim_rep,time_bins)

for (i in 1:nrow(result_df)){
  
  print(i)
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  time_bins <- result_df$time_bins[i]
  
  phi_1 <- 1
  phi_2 <- phi_1/10
  
  tau_1 <- 1
  tau_2 <- tau_1 
  
  dat <- data.frame()
  
  set.seed(sim_rep)
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  rint <- c(0.5,1,Inf)
  tint <- seq(time_bins,5,by = time_bins)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  
  n_surveys <- 1000
  Yarray <- array(NA,dim=c(n_surveys,nrint,ntint))
  rarray <- array(NA,dim=c(n_surveys,nrint))
  tarray <- array(NA,dim=c(n_surveys,ntint))
  
  rarray <- matrix(rep(rint,1000),byrow = TRUE, nrow=n_surveys)
  tarray <- matrix(rep(tint,1000),byrow = TRUE, nrow=n_surveys)
  
  Y <- 0
  
  dim <- 10
  
  N <- round(Density*dim^2) # Number of birds to simulate at each point count
  
  # -------------------------------------------------------
  # Birds with fast cue rate
  # -------------------------------------------------------
  
  for (s in 1:n_surveys){
    
    dat <- data.frame()
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
    
    cues$p <- exp(-(cues$dist/tau_1)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    cues <- subset(cues,detected == 1 & cue_number == 1) # 
    cues <- cues[!duplicated(cues$bird_id),]
    
    if (nrow(cues)>0){
      cues$phi <- phi_1
      dat <- rbind(dat, cues)
    }
    
    # -------------------------------------------------------
    # Birds with slow cue rate
    # -------------------------------------------------------
    
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
    
    cues$p <- exp(-(cues$dist/(tau_2))^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    cues <- subset(cues,detected == 1 & cue_number == 1) # & cue_number == 1
    cues <- cues[!duplicated(cues$bird_id),]
    
    if (nrow(cues)>0){
      cues$phi <- phi_2
      dat <- rbind(dat, cues)
    }
    
    # ------------------------------------
    # Transcription: distance and time bins
    # ------------------------------------
    
    Yarray[s,,] <- 0
    if (nrow(dat)>0){
      
      # Separate into distance and time bins
      dat$rint <- cut(dat$dist,c(0,rint))
      dat$tint <- cut(dat$time,c(0,tint))
      dat <- na.omit(dat)
      
      Y <- table(dat[,c("rint","tint")])
      
      Y # Data to analyze
      Yarray[s,,] <- Y
    }
  }
  
  # ------------------------------------
  # Analysis: removal model using finite mixture
  # ------------------------------------
  
  Y_removal <- apply(Yarray,c(1,3),sum)
  #fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  #phi_est = exp(fit_p$coefficients)
  
  Mf0 <- detect::cmulti(Y_removal  | tarray ~ 1, type="fmix")
  
  # Proportion of infrequent singers
  prop_infreq <- plogis(Mf0$coefficients[2])
  
  # Singing rate of infrequent singers
  phi_infreq <- exp(Mf0$coefficients[1])
  
  # Proportion of birds that sing at least once during the survey
  p_mix <- 1 - prop_infreq * exp(-max(tint)*phi_infreq)
  
  # ------------------------------------
  # Analysis: removal model using standard removal modeling
  # ------------------------------------
  
  M0 <- detect::cmulti(Y_removal  | tarray ~ 1, type="rem")
  phi_0 <- exp(M0$coefficients[1])
  p_0 <- 1-exp(-max(tint)*phi_0)
  
  # ------------------------------------
  # Analysis: distance model
  # ------------------------------------
  
  Y_distance <- apply(Yarray,c(1,2),sum)
  a = which(rowSums(Y_distance)==0)
  
  if (length(a)>0){
  fit_q <- detect::cmulti(Y_distance[-a,] | rarray[-a,] ~ 1, type = "dis")
  } else {
    fit_q <- detect::cmulti(Y_distance | rarray ~ 1, type = "dis")
  }
  tau_est = exp(fit_q$coefficients)
  
  # **********************************************************************
  # Density estimates
  # **********************************************************************
  
  # Estimate density for each species
  A <- pi*tau_est^2 # area sampled by unlimited distance survey
  D_1 <- sum(Yarray)/n_surveys/(A*p_mix) %>% as.numeric()
  D_0 <- sum(Yarray)/n_surveys/(A*p_0) %>% as.numeric()
  
  result_df$D_1[i] = D_1
  result_df$D_0[i] = D_0
  
  result_plot <- ggplot(result_df) +
    geom_hline(yintercept = Density, linetype = 2) + 
    geom_hline(yintercept = 0, linetype = 2, col = "transparent") + 
    
    geom_jitter(aes(x = "Standard removal", y = D_0), width = 0.1, height = 0)+
    geom_jitter(aes(x = "Mixture removal", y = D_1), width = 0.1, height = 0)+
    ylab("Density Estimate")+
    xlab("")+
    theme_bw()
  
  print(result_plot)
  
  # result_summary <- result_df %>%
  #   group_by(time_bins) %>%
  #   summarize(p_est = mean(p_est,na.rm = TRUE))
  # 
  # if ((sim_rep %% 10) == 0) {
  # 
  # 
  #   result_plot_p <- ggplot(data = result_df)+
  #     geom_line(aes(x = time_bins, y = p_est, col = factor(sim_rep)))+
  #     scale_color_manual(values = rep("gray80",length(unique(result_df$sim_rep))), guide = "none")+
  #     geom_line(data = result_summary, aes(x = time_bins, y = p_est), size = 2, col = "dodgerblue")+
  #     geom_hline(yintercept = true_p, linetype = 2)+
  #     xlab("Time bin length")+
  #     ylab("Estimated p")+
  #     theme_bw()
  # 
  #   print(result_plot_p)
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
