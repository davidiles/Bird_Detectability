library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)
library(jagsUI)
library(reshape2)


setwd("~/1_Work/Bird_Detectability/QPAD")

# ----------------------------------------
# Script to fit model
# ----------------------------------------

sink("QPAD_corrected_vectorized.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    phi ~ dunif(0,2)
    tau ~ dunif(0,3)
    
    # ------------------------------
    # Simulation
    # ------------------------------
    
    # For each simulated bird, probability it has generated at least one detectable cue by end of time interval
    for (i in 1:Nsim){
    
      # Simulate cueing for each bird in each time interval
      for (j in 1:ntint){
      
        # Expected number of cues produced during interval
        lambda[i,j] <- phi*tint_duration[j]*exp(-(sim_dist[i]/tau)^2)
        
        # Probability bird produces at least one detectable cue
        p[i,j] <- 1 - exp(-lambda[i,j])
        dc[i,j] ~ dbern(p[i,j]) # binary (was bird detected?)
      }
      
      # Ensure that each bird only gets detected once
      yy[i,1] <- dc[i,1]
      for (j in 1:(ntint-1)){
        mx[i,j] <- max(dc[i,1:j]) # Max detection up to this time bin
        yy[i,j+1] <- dc[i,j+1]*(1-mx[i,j]) # If bird was detected previously, set to zero (so we can sum yy in each distance/time bin)
      }
      
    }
    
    # Sum up first detections in each time and distance bin
    for (i in 1:nrint){
      for (j in 1:ntint){
        Ysim_matrix[i,j] <- sum(yy[r_from[i]:r_to[i],j])
      }
    }
    
    # Place counts in correct vectorized bins
    for (i in 1:n_bins_vectorized){
        Ysim_vectorized[i] <- Ysim_matrix[Y_rint_numeric[i],Y_tint_numeric[i]]
    }
      
    # -----------------------------------------------------
    # Use relative proportions for detectability modeling
    # -----------------------------------------------------
    
    Y_vectorized[1:n_bins_vectorized] ~ dmulti(Ysim_vectorized[],N_det)
    
  }
    
",fill = TRUE)
sink()

# ----------------------------------------
# This script just generates initial values of
# dc[i,j] that are internally consistent with values of tau and phi
#
# Otherwise, bad initial values cause JAGS to fail to initialize
# ----------------------------------------

sink("QPAD_init_vectorized.jags")
cat("
    model {
    
    dummy ~ dunif(0,1)
    
    # ------------------------------
    # Simulation
    # ------------------------------
    
    # For each simulated bird, probability it has generated at least one detectable cue by end of time interval
    for (i in 1:Nsim){
    
      # Simulate cueing for each bird in each time interval
      for (j in 1:ntint){
      
        # Expected number of cues produced during interval
        lambda[i,j] <- phi*tint_duration[j]*exp(-(sim_dist[i]/tau)^2)
        
        # Probability bird produces at least one detectable cue
        p[i,j] <- 1 - exp(-lambda[i,j])
        dc[i,j] ~ dbern(p[i,j]) # binary (was bird detected?)
      }
      
      # Ensure that each bird only gets detected once
      yy[i,1] <- dc[i,1]
      for (j in 1:(ntint-1)){
        mx[i,j] <- max(dc[i,1:j]) # Max detection up to this time bin
        yy[i,j+1] <- dc[i,j+1]*(1-mx[i,j]) # If bird was detected previously, set to zero (so we can sum yy in each distance/time bin)
      }
      
    }
    
    # Sum up first detections in each time and distance bin
    for (i in 1:nrint){
      for (j in 1:ntint){
        Ysim_matrix[i,j] <- sum(yy[r_from[i]:r_to[i],j])
      }
    }
    
    # Place counts in correct vectorized bins
    for (i in 1:n_bins_vectorized){
        Ysim_vectorized[i] <- Ysim_matrix[Y_rint_numeric[i],Y_tint_numeric[i]]
      }
    
  }
    
",fill = TRUE)
sink()


rm(list=ls())

results_df = data.frame()

for (sim_rep in 1:100){
  
  #  for (tau_true in c(0.5,1.5)){
  #    for (phi_true in c(0.2,0.4)){
  
  N = 100000 # Number of birds to place on landscape (determines overall sample size in point count dataset)
  dim = 10 # landscape size (dim x dim landscape)
  Density_true = N/dim^2 
  
  tau_true <- 0.7
  phi_true <- 0.5
  
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
  
  # Change to vector format (so dcat() can work in JAGS)
  Y_vec <- reshape2::melt(Y) %>%
    rename(rint = Var1, tint = Var2, Y = value)
  
  Y_rint_numeric <- factor(Y_vec$rint) %>% as.numeric()
  Y_tint_numeric <- factor(Y_vec$tint) %>% as.numeric()
  Y_vectorized <- Y_vec$Y
  
  # ------------------------------------
  # PART 2: FIT MODELS WITH STANDARD QPAD
  # ------------------------------------
  
  # Distance model
  Y_distance = matrix(rowSums(Y),1)
  D_distance = matrix(rint,1)
  
  # Removal model
  Y_removal = matrix(colSums(Y),1)
  D_removal = matrix(tint,1)
  
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
  # - reasonably high chance of model not initializing - need to try setting initial values multiple times
  # ------------------------------------
  
  # Trycatch to keep trying initializations until one works
  out = NULL
  
  while (is.null(out)){
    
    out <- tryCatch({
      
      # Plausible initial values
      phi_init = 0.2
      tau_init = 1.2
      
      # Data augmentation
      Nsim = 1000 # Number of sim_birds to place on landscape (size of augmented dataset)
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
      sim_birds$rint <- cut(sim_birds$dist,c(0,rint))
      sim_birds$rint_num <- sim_birds$rint %>% factor() %>% as.numeric()
      
      # For each distance bin, store first and last bird in that bin
      r_from <- r_to <- rint*NA
      for (i in 1:length(rint)){
        r_from[i] <- min(which(sim_birds$rint_num == i))
        r_to[i] <- max(which(sim_birds$rint_num == i))
      }
      
      # Prepare to generate initial values
      init_data <- list(
        phi = phi_init,
        tau = tau_init,
        
        # Time bin information
        ntint = length(tint),
        tint_duration = diff(c(0,tint)),
        
        # Distance bin information
        nrint = length(rint),
        
        # Data augmentation
        Nsim = nrow(sim_birds),
        sim_dist = sim_birds$dist,
        
        # Distance that each simulated bird belongs to
        r_from = r_from,
        r_to = r_to,
        
        # To vectorize
        n_bins_vectorized = length(Y_rint_numeric),
        Y_rint_numeric = Y_rint_numeric,
        Y_tint_numeric = Y_tint_numeric
        
      )
      
      # Run once to generate initial values that work
      out_inits <- jags(data = init_data,
                        model.file =  "QPAD_init_vectorized.jags",
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
        Y_vectorized = Y_vectorized,
        
        # Time bin information
        tint = tint,
        ntint = length(tint),
        tint_duration = diff(c(0,tint)),
        
        # Distance bin information
        nrint = length(rint),
        
        # Data augmentation
        Nsim = nrow(sim_birds),
        sim_dist = sim_birds$dist,
        
        # Distance that each simulated bird belongs to
        r_from = r_from,
        r_to = r_to,
        
        # To vectorize
        n_bins_vectorized = length(Y_rint_numeric),
        Y_rint_numeric = Y_rint_numeric,
        Y_tint_numeric = Y_tint_numeric
      )
      
      # Fit model
      inits <- function()list(dc = out_inits$sims.list$dc[1,,],
                              phi = phi_init,
                              tau = tau_init)
      
      # Using initial values
      out <- jags(data = jags_data,
                  model.file =  "QPAD_corrected_vectorized.jags",
                  parameters.to.save = c("phi","tau","Density"),
                  inits = inits,
                  n.chains = 3,
                  n.thin = 5,
                  n.iter = 5000,
                  n.burnin = 2000,
                  parallel = TRUE)
      
      
      out # 7 min to fit
      
    },
    error = function(cond){NULL})
  } # end while loop
  
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
                                phi_est_q975 = out$q97.5$phi
                     ))
  
  # ----------------------------------
  # PLOT RESULTS
  # ----------------------------------
  mean_df <- results_df %>%
    group_by(tau_label, phi_label) %>%
    summarize(tau_true = mean(tau_true),
              phi_true = mean(phi_true),
              mean_tau_corrected = mean(tau_est_q50),
              mean_tau_uncorrected = mean(tau_uncorrected),
              mean_phi_corrected = mean(phi_est_q50),
              mean_phi_uncorrected = mean(phi_uncorrected)
    )
  
  tau_plot <- ggplot(data = results_df, aes(x = sim_rep, y = tau_est_q50, ymin = tau_est_q025, ymax = tau_est_q975))+
    geom_hline(data = results_df, aes(yintercept = tau_true), col = "dodgerblue", size = 2, alpha = 0.5)+
    #geom_errorbar(width = 0, aes(col = "Corrected"))+
    geom_point(aes(col = "Corrected"))+
    geom_point(data = results_df, aes(x = sim_rep, y = tau_uncorrected, col = "QPAD"))+
    
    # Mean corrected estimate
    geom_hline(data = mean_df, aes(yintercept = mean_tau_corrected, col = "Corrected"), linetype = 2)+
    
    # Mean uncorrected estimate
    geom_hline(data = mean_df, aes(yintercept = mean_tau_uncorrected, col = "QPAD"), linetype = 2)+
    
    scale_color_manual(values = c("black","orangered"), name = "Estimate")+
    facet_grid(phi_label~tau_label)+
    ggtitle("Estimates of EDR")+
    xlab("Simulation #")+
    ylab("EDR")+
    theme_bw()
  #tau_plot
  
  phi_plot <- ggplot(data = results_df, aes(x = sim_rep, y = phi_est_q50, ymin = phi_est_q025, ymax = phi_est_q975))+
    geom_hline(data = results_df, aes(yintercept = phi_true), col = "dodgerblue", size = 2, alpha = 0.5)+
    #geom_errorbar(width = 0, aes(col = "Corrected"), alpha = 0.5)+
    geom_point(aes(col = "Corrected"))+
    geom_point(data = results_df, aes(x = sim_rep, y = phi_uncorrected, col = "QPAD"))+
    
    # Mean corrected estimate
    geom_hline(data = mean_df, aes(yintercept = mean_phi_corrected, col = "Corrected"), linetype = 2)+
    
    # Mean uncorrected estimate
    geom_hline(data = mean_df, aes(yintercept = mean_phi_uncorrected, col = "QPAD"), linetype = 2)+
    
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




