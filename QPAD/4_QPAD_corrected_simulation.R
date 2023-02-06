library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)
library(jagsUI)


setwd("~/1_Work/Bird_Detectability/QPAD")

# ----------------------------------------
# Script to fit model
# ----------------------------------------

sink("QPAD_corrected.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    phi ~ dunif(0,1)
    tau ~ dunif(0,2)
    
    # ------------------------------
    # Simulation
    # ------------------------------
    
    # For each simulated bird, probability it has generated at least one detectable cue by end of time interval
    for (i in 1:Nsim){
    
      # Simulate cueing for each bird in each time interval
      for (j in 1:ntbin){
      
        # Expected number of cues produced during interval
        lambda[i,j] <- phi*tbin_duration[j]*exp(-(sim_dist[i]/tau)^2)
        
        # Probability bird produces at least one detectable cue
        p[i,j] <- 1 - exp(-lambda[i,j])
        dc[i,j] ~ dbern(p[i,j]) # binary (was bird detected?)
        
      }
      
      # Ensure that each bird only gets detected once
      yy[i,1] <- dc[i,1]
      for (j in 1:(ntbin-1)){
        mx[i,j] <- max(dc[i,1:j]) # Max detection up to this time bin
        yy[i,j+1] <- dc[i,j+1]*(1-mx[i,j]) # If bird was detected previously, set to zero (so we can sum yy in each distance/time bin)
      }
    }
    
    # Sum up first detections in each time and distance bin
    for (i in 1:nrbin){
      for (j in 1:ntbin){
        Ysim[i,j] <- sum(yy[r_from[i]:r_to[i],j])
      }
    }
    
    # -----------------------------------------------------
    # Marginals
    # -----------------------------------------------------
    
    for (i in 1:nrbin){
      p_distance[i] <- sum(Ysim[i,1:ntbin])/sum(Ysim[1:nrbin,1:ntbin])
    }
    for (j in 1:ntbin){
      p_removal[j] <- sum(Ysim[1:nrbin,j])/sum(Ysim[1:nrbin,1:ntbin])
    }
    
    # -----------------------------------------------------
    # Use relative proportions for detectability modeling
    # -----------------------------------------------------
    
    Y_distance[1:nrbin] ~ dmulti(p_distance[],N_det)
    Y_removal[1:ntbin] ~ dmulti(p_removal[],N_det)
    
    # -----------------------------------------------------
    # Derived parameters
    # -----------------------------------------------------
    
    # Density
    A_hat = pi*tau^2
    p_hat = 1-exp(-max(tbin)*phi)
    Density <- N_det/(A_hat*p_hat)
    
    
  }
    
",fill = TRUE)
sink()

# ----------------------------------------
# This script just generates initial values of
# dc[i,j] that are internally consistent with values of tau and phi
#
# Otherwise, bad initial values cause JAGS to fail to initialize
# ----------------------------------------

sink("QPAD_init.jags")
cat("
    model {
    
    # ------------------------------
    # Simulation
    # ------------------------------
    
    # For each simulated bird, probability it has generated at least one detectable cue by end of time interval
    for (i in 1:Nsim){
    
      # Simulate cueing for each bird in each time interval
      for (j in 1:ntbin){
      
        # Expected number of cues produced during interval
        lambda[i,j] <- phi*tbin_duration[j]*exp(-(sim_dist[i]/tau)^2)
        
        # Probability bird produces at least one detectable cue
        p[i,j] <- 1 - exp(-lambda[i,j])
        dc[i,j] ~ dbern(p[i,j]) # binary (was bird detected?)
      }
      
      # Ensure that each bird only gets detected once
      yy[i,1] <- dc[i,1]
      for (j in 1:(ntbin-1)){
        mx[i,j] <- max(dc[i,1:j]) # Max detection up to this time bin
        yy[i,j+1] <- dc[i,j+1]*(1-mx[i,j]) # If bird was detected previously, set to zero (so we can sum yy in each distance/time bin)
      }
      
    }
    
    # Sum up first detections in each time and distance bin
    for (i in 1:nrbin){
      for (j in 1:ntbin){
        Ysim[i,j] <- sum(yy[r_from[i]:r_to[i],j])
      }
    }
    
    # -----------------------------------------------------
    # Marginals
    # -----------------------------------------------------
    
    for (i in 1:nrbin){
      p_distance[i] <- sum(Ysim[i,1:ntbin])/sum(Ysim[1:nrbin,1:ntbin])
    }
    
    for (j in 1:ntbin){
      p_removal[j] <- sum(Ysim[1:nrbin,j])/sum(Ysim[1:nrbin,1:ntbin])
    }
    
    dummy ~ dunif(0,1)
  }
    
",fill = TRUE)
sink()


rm(list=ls())

results_df = data.frame()
for (sim_rep in 5:100){
  #  for (tau_true in c(0.5,1.5)){
  #    for (phi_true in c(0.2,0.4)){
  
  N = 10000 # Number of birds to place on landscape
  dim = 10 # landscape size (dim x dim landscape)
  Density_true = N/dim^2 
  
  tau_true <- 1
  phi_true <- 0.2
  
  # ------------------------------------
  # PART 1: GENERATE 'ACTUAL' DATA
  # ------------------------------------
  
  birds <- data.frame(bird_id = 1:N,
                      x = runif(N,-dim/2,dim/2),
                      y = runif(N,-dim/2,dim/2))
  
  
  # Distances to observer
  birds$dist <- sqrt(birds$x^2 + birds$y^2)
  
  # Simulate bird cues, based on phi_true
  cues <- matrix(NA, nrow=N, ncol = 10)
  for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_true))
  cues <- cues %>% 
    reshape2::melt() %>% 
    rename(bird_id = Var1, cue_number = Var2, time = value) %>%
    arrange(bird_id,cue_number)
  
  cues$dist <- birds$dist[cues$bird_id]
  
  # Determine which cues are detected, based on tau_true
  cues$p<- exp(-(cues$dist/tau_true)^2)  # Probability each cue is detected
  cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
  
  # Isolate first detected cue for each bird
  dat <- subset(cues,detected == 1)
  dat <- dat[!duplicated(dat$bird_id),]
  
  # Transcription: when distance is fixed
  tint <- c(3,5,10)
  rint <- c(0.5,1,Inf)
  
  dat$tint <- cut(dat$time,c(0,tint))
  dat$rint <- cut(dat$dist,c(0,rint))
  Y = table(dat$rint,dat$tint)
  
  # this is the data we will analyze
  #Y_removal <- colSums(Y)
  #Y_distance <- rowSums(Y)
  
  # Distance model
  Y_distance = matrix(rowSums(Y),1)
  D_distance = matrix(rint,1)
  
  # Removal model
  Y_removal = matrix(colSums(Y),1)
  D_removal = matrix(tint,1)
  
  # ------------------------------------
  # PART 2: FIT MODELS WITH STANDARD QPAD
  # ------------------------------------
  
  fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
  tau_MLE = exp(fit.q$coefficients)
  
  fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_MLE = exp(fit.p$coefficients)
  
  # Estimate density
  A_hat = pi*tau_MLE^2
  p_hat = 1-exp(-max(tint)*phi_MLE)
  Density_MLE <- sum(Y)/(A_hat*p_hat)
  
  tau_MLE
  phi_MLE
  # ------------------------------------
  # PART 3: GENERATE 'AUGMENTED' DATA
  #
  # - also specify initial values for JAGS
  # ------------------------------------
  
  # Plausible initial values
  phi_init = 0.5
  tau_init = 1
  
  # Data augmentation
  Nsim = 1000 # Number of sim_birds to place on landscape
  dim = 6 # landscape size (dim x dim landscape)
  
  sim_birds <- data.frame(bird_id = 1:Nsim,
                          x = runif(Nsim,-dim/2,dim/2),
                          y = runif(Nsim,-dim/2,dim/2))
  
  # Random distances to observer
  sim_birds$dist <- sqrt(sim_birds$x^2 + sim_birds$y^2)
  sim_birds$rint <- cut(sim_birds$dist,c(0,rint))
  
  
  # Sort simulated sim_birds by distance
  sim_birds <- sim_birds %>% arrange(dist)
  
  # Group simulated birds into distance bins
  sim_birds$rbin <- cut(sim_birds$dist,c(0,rint))
  sim_birds$rbin_num <- sim_birds$rbin %>% factor() %>% as.numeric()
  
  # For each distance bin, store first and last bird in that bin
  r_from <- r_to <- rint*NA
  for (i in 1:length(rint)){
    r_from[i] <- min(which(sim_birds$rbin_num == i))
    r_to[i] <- max(which(sim_birds$rbin_num == i))
  }
  
  # Prepare to generate initial values
  init_data <- list(
    phi = phi_init,
    tau = tau_init,
    
    # Time bin information
    ntbin = length(tint),
    tbin_duration = diff(c(0,tint)),
    
    # Distance bin information
    nrbin = length(rint),
    
    # Data augmentation
    Nsim = nrow(sim_birds),
    sim_dist = sim_birds$dist,
    
    # Distance that each simulated bird belongs to
    r_from = r_from,
    r_to = r_to
    
  )
  
  # Run once to generate initial values that work
  out_inits <- jags(data = init_data,
                    model.file =  "QPAD_init.jags",
                    parameters.to.save = c("dc","dummy"),
                    inits = NULL,
                    n.chains = 3,
                    n.thin = 1,
                    n.iter = 2,
                    n.burnin = 1,
                    codaOnly=c("dc"))
  
  # ******************************************
  # ******************************************
  # PART 4: WITH INITIAL VALUES CHOSEN, FIT MODEL TO ACTUAL DATA
  # ******************************************
  # ******************************************
  
  jags_data <- list(
    pi = pi,
    N_det = sum(Y),
    Y_distance = rowSums(Y),
    Y_removal = colSums(Y),
    
    # Time bin information
    tbin = tint,
    ntbin = length(tint),
    tbin_duration = diff(c(0,tint)),
    
    # Distance bin information
    nrbin = length(rint),
    
    # Data augmentation
    Nsim = nrow(sim_birds),
    sim_dist = sim_birds$dist,
    
    # Distance that each simulated bird belongs to
    r_from = r_from,
    r_to = r_to
  )
  
  # Fit model
  inits <- function()list(dc = out_inits$sims.list$dc[1,,],
                          phi = phi_init,
                          tau = tau_init)
  
  # Using initial values
  out <- jags(data = jags_data,
              model.file =  "QPAD_corrected.jags",
              parameters.to.save = c("phi","tau","Density"),
              inits = inits,
              n.chains = 3,
              n.thin = 5,
              n.iter = 5000,
              n.burnin = 2000,
              parallel = TRUE)
  
  
  out
  
  phi_est = out$mean$phi
  tau_est = out$mean$tau
  A = pi*tau_est^2
  p = 1-exp(-max(tint)*phi_est)
  Density_est <- sum(Y)/(A*p)
  
  # ----------------------------------
  # STORE RESULTS
  # ----------------------------------
  
  results_df = rbind(results_df,
                     data.frame(sim_rep = sim_rep,
                                tau_true = tau_true,
                                phi_true = phi_true,
                                Density_true = Density_true,
                                tau_label = paste0("EDR = ",tau_true),
                                phi_label = paste0("Cue Rate = ",phi_true),
                                
                                tau_uncorrected = tau_MLE,
                                phi_uncorrected = phi_MLE,
                                Density_uncorrected = Density_MLE,
                                
                                tau_est_q50 = out$q50$tau,
                                tau_est_q025 = out$q2.5$tau,
                                tau_est_q975 = out$q97.5$tau,
                                
                                phi_est_q50 = out$q50$phi,
                                phi_est_q025 = out$q2.5$phi,
                                phi_est_q975 = out$q97.5$phi,
                                
                                Density_est_q50 = out$q50$Density,
                                Density_est_q025 = out$q2.5$Density,
                                Density_est_q975 = out$q97.5$Density
                     ))
  
  # ----------------------------------
  # PLOT RESULTS
  # ----------------------------------
  
  tau_plot <- ggplot(data = results_df, aes(x = sim_rep, y = tau_est_q50, ymin = tau_est_q025, ymax = tau_est_q975))+
    geom_hline(data = results_df, aes(yintercept = tau_true), col = "dodgerblue", size = 2, alpha = 0.5)+
    #geom_errorbar(width = 0, aes(col = "Corrected"))+
    geom_point(aes(col = "Corrected"))+
    geom_point(data = results_df, aes(x = sim_rep, y = tau_uncorrected, col = "QPAD"))+
    scale_color_manual(values = c("black","orangered"), name = "Estimate")+
    facet_grid(phi_label~tau_label)+
    ggtitle("Estimates of EDR")+
    xlab("Simulation #")+
    ylab("EDR")+
    theme_bw()
  
  phi_plot <- ggplot(data = results_df, aes(x = sim_rep, y = phi_est_q50, ymin = phi_est_q025, ymax = phi_est_q975))+
    geom_hline(data = results_df, aes(yintercept = phi_true), col = "dodgerblue", size = 2, alpha = 0.5)+
    #geom_errorbar(width = 0, aes(col = "Corrected"))+
    geom_point(aes(col = "Corrected"))+
    geom_point(data = results_df, aes(x = sim_rep, y = phi_uncorrected, col = "QPAD"))+
    scale_color_manual(values = c("black","orangered"), name = "Estimate")+
    facet_grid(phi_label~tau_label)+
    ggtitle("Estimates of Cue Rate")+
    xlab("Simulation #")+
    ylab("Cue rate")+
    theme_bw()
  
  estimate_plot <- ggarrange(tau_plot,phi_plot,nrow=2)
  print(estimate_plot)  
  
}
#}
#}




