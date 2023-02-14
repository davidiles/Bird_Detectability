library(tidyverse)
library(jagsUI)
library(ggpubr)
library(detect)

sink("analysis_integrated.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    # Relatively wide, flat priors
    phi ~ dunif(0,5)
    tau ~ dunif(0,5)
    
    # ------------------------------
    # Numerical integration
    # ------------------------------
    
    # Calculate CDF across distance and time bins
    for (j in 1:ntint){
      for (i in 1:(nbins-1)){
        CDF[i,j] = ((pi * bins[i+1]^2) - (pi*bins[i]^2)) * (1-exp(-phi*tint[j]*exp(-midpoints[i]^2/tau^2)))
      }
      
      # Sum within each distance bin to calculate binned CDF
      for (i in 1:nrint){
        CDF_binned[i,j] <- sum(CDF[1:rint_indices[i+1],j])
      }
    }
    
    # --------------------------------------------------------
    # Difference CDF to calculate multinomial cell probabilities
    # --------------------------------------------------------
    
    tmp1[1,1:ntint] <- CDF_binned[1,1:ntint]
    for (i in 2:nrint){
      tmp1[i,1:ntint] <- CDF_binned[i,1:ntint] - CDF_binned[i-1,1:ntint]
    }
    
    p_matrix[1:nrint,1] <- tmp1[1:nrint,1]
    for (j in 2:ntint){
      p_matrix[1:nrint,j] <- tmp1[1:nrint,j] - tmp1[1:nrint,j-1]
    }
    
    # --------------------------------------------------------
    # Place  in vectorized cells for use with dmulti()
    # --------------------------------------------------------
    
    for (i in 1:n_bins_vectorized){
        p_vectorized[i] <- p_matrix[Y_rint_numeric[i],Y_tint_numeric[i]]
    }
    
    # -----------------------------------------------------
    # Use relative proportions for detectability modeling
    # Note that JAGS automatically normalizes p_vectorized (p_vectorized/sum(p_vectorized)) so that probabilities sum to 1
    # -----------------------------------------------------
    
    # Observed counts in each distance/time bin
    Y_vectorized[1:n_bins_vectorized] ~ dmulti(p_vectorized[],Total_detected)
    
    # -----------------------------------------------------
    # Density
    # -----------------------------------------------------

    C <- sum(p_vectorized[])/(pi*pow(bins[nbins],2))
    Y_corrected <- sum(Y_vectorized[])/C
    Density <- Y_corrected/(pi*pow(bins[nbins],2))
    
  }
    
",fill = TRUE)
sink()

rm(list=ls())

results_df = data.frame()

for (sim_rep in 1:1000){
  
  # ******************************************
  # PART 1: SIMULATE DATA
  # ******************************************
  
  set.seed(sim_rep)
  
  tau_true = 0.6
  phi_true = 1
  
  N = 10000 # Number of birds to place on landscape (select a high number to provide sufficient sample size)
  dim = 10 # landscape size (metres)
  
  Density_true <- N/dim^2
  
  # ------------------------------------
  # Place birds on landscape around observer (centred on landscape)
  # ------------------------------------
  
  birds <- data.frame(x = runif(N,-dim/2,dim/2),
                      y = runif(N,-dim/2,dim/2))
  
  # Distances to observer
  birds$dist <- sqrt(birds$x^2 + birds$y^2)
  
  # Remove birds outside maximum distance
  birds <- subset(birds, dist <= (dim/2))
  birds <- birds %>% arrange(dist)
  birds$bird_id = 1:nrow(birds)
  N = nrow(birds)
  
  # ------------------------------------
  # Simulate bird cues, based on phi_true
  # ------------------------------------
  
  cues <- matrix(NA, nrow=N, ncol = 50)
  for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_true))
  cues <- cues %>% 
    reshape2::melt() %>% 
    rename(bird_id = Var1, cue_number = Var2, time = value) %>%
    arrange(bird_id,cue_number)
  
  cues$dist <- birds$dist[cues$bird_id]
  
  # ------------------------------------
  # Determine which cues are detected, based on tau_true
  # ------------------------------------
  
  cues$p<- exp(-(cues$dist/tau_true)^2)  # Probability each cue is detected
  cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
  
  # ------------------------------------
  # Isolate first detected cue for each bird
  # ------------------------------------
  
  dat <- subset(cues,detected == 1)
  dat <- dat[!duplicated(dat$bird_id),]
  
  
  # ------------------------------------
  # Transcription: distance and time bins
  # ------------------------------------
  
  rint <- c(0.5,1,2)
  tint <- seq(1,10,1) #c(3,5)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # Separate into distance and time bins
  dat$rint <- cut(dat$dist,c(0,rint))
  dat$tint <- cut(dat$time,c(0,tint))
  dat <- na.omit(dat)
  
  Y <- table(dat[,c("rint","tint")])
  
  Y # Data to analyze
  
  # ******************************************
  # PART 2: FIT MODELS WITH STANDARD QPAD
  # ******************************************
  
  # Distance model
  Y_distance = matrix(rowSums(Y),1)
  D_distance = matrix(rint,1)
  
  # Removal model
  Y_removal = matrix(colSums(Y),1)
  D_removal = matrix(tint,1)
  
  #fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
  #tau_MLE = exp(fit.q$coefficients)
  
  #fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
  #phi_MLE = exp(fit.p$coefficients)
  
  # Estimate density
  #A_hat = pi*tau_MLE^2
  #p_hat = 1-exp(-max(tint)*phi_MLE)
  #Density_MLE <- sum(Y)/(A_hat*p_hat)
  
  # ******************************************
  # PART 3: JOINT ANALYSIS IN JAGS
  # ******************************************
  
  # Prepare for numerical integration (at midpoints)
  
   # Maximum distance for integration
  if (max(rint) != Inf) max_dist = max(rint) else max_dist = dim/2 # If fixed radius point count
  
  bins <- seq(0,max_dist,length.out = 100) # Large number of bins for better approximation (but slower)
  midpoints <- bins[-length(bins)] + diff(bins)
  nbins <- length(bins)
  
  # Determine which bins belong in which rinterval
  rint_indices <- c(1)
  for (i in 1:nrint) rint_indices[i+1] = max(which(midpoints < rint[i]))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Code below just confirms that numerical integration equations are approximately correct
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Ysum_expected = matrix(NA, nrow=nrint, ncol = ntint)

  for (j in 1:ntint){
    Total_detected = c()
    N_tot = c()
    A_tot = c()

    for (i in 1:(nbins-1)){

      P_bin = 1-exp(-phi_true*tint[j]*exp(-midpoints[i]^2/tau_true^2)) # Approximate probability each bird in bin is detected by end of interval
      A_bin = (pi * bins[i+1]^2) - (pi*bins[i]^2)
      N_bin = A_bin*Density_true
      Total_detected_bin = N_bin * P_bin
      Total_detected = c(Total_detected,Total_detected_bin)
    }

    for (k in 1:nrint){
      Ysum_expected[k,j] <- sum(Total_detected[rint_indices[k]:rint_indices[k+1]])
    }

  }

  Ysim_matrix <- Y*NA

  for (i in 1:nrint){
    Ysim_matrix[i,1] <- Ysum_expected[i,1]
    for (j in 2:ntint){
      Ysim_matrix[i,j] <- Ysum_expected[i,j] - sum(Ysum_expected[i,(j-1)])
    }
  }
  (Ysim_matrix/sum(Ysim_matrix)) %>% round(3)
  (Y/sum(Y)) %>% round(3)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # --------------------------------------
  # ANALYZE DATA IN JAGS
  # --------------------------------------
  
  # Vectorize Y (actual data) so we can use multinomial distribution in JAGS
  Y_vec <- reshape2::melt(matrix(Y, nrow = nrint, ncol = ntint)) %>% rename(rint = Var1, tint = Var2, Y = value)
  Y_rint_numeric <- factor(Y_vec$rint) %>% as.numeric()
  Y_tint_numeric <- factor(Y_vec$tint) %>% as.numeric()
  Y_vectorized <- Y_vec$Y
  
  # Data in list for JAGS
  jags_data <- list(
    
    pi = pi,
    
    # Vectorized data
    Total_detected = sum(Y_vectorized),
    Y_vectorized = Y_vectorized,
    n_bins_vectorized = length(Y_rint_numeric),
    Y_rint_numeric = Y_rint_numeric,
    Y_tint_numeric = Y_tint_numeric,
    
    # Time bin information
    tint = tint,
    ntint = length(tint),
    
    # Distance bin information
    nrint = length(rint),
    rint = rint,
    
    # For numerical integration
    bins = bins,
    nbins = nbins,
    midpoints = midpoints,
    rint_indices = rint_indices
  )
  
  # Fit model
  out <- jags(data = jags_data,
              model.file =  "analysis_integrated.jags",
              parameters.to.save = c("phi","tau","Density"),
              n.chains = 3,
              n.thin = 1,
              n.iter = 2000,
              n.burnin = 1000,
              parallel = TRUE)
  
  # ----------------------------------
  # STORE RESULTS FOR THIS SIMULATION
  # ----------------------------------
  
  results_df = rbind(results_df,
                     data.frame(sim_rep = sim_rep,
                                
                                tau_true = tau_true,
                                #tau_MLE = tau_MLE,
                                tau_est_q50 = out$q50$tau,
                                tau_est_q05 = quantile(out$sims.list$tau,0.05),
                                tau_est_q95 = quantile(out$sims.list$tau,0.95),
                                tau_cov90 = quantile(out$sims.list$tau,0.05) < tau_true & quantile(out$sims.list$tau,0.95) > tau_true,
                                
                                phi_true = phi_true,
                                #phi_MLE = phi_MLE,
                                phi_est_q50 = out$q50$phi,
                                phi_est_q05 = quantile(out$sims.list$phi,0.05),
                                phi_est_q95 = quantile(out$sims.list$phi,0.95),
                                phi_cov90 = quantile(out$sims.list$phi,0.05) < phi_true & quantile(out$sims.list$phi,0.95) > phi_true,
                                
                                Density_true = Density_true,
                                #Density_MLE = Density_MLE,
                                Density_est_q50 = out$q50$Density,
                                Density_est_q05 = quantile(out$sims.list$Density,0.05),
                                Density_est_q95 = quantile(out$sims.list$Density,0.95),
                                Density_cov90 = quantile(out$sims.list$Density,0.05) < Density_true & quantile(out$sims.list$Density,0.95) > Density_true
                                
                     ))
  
  # ----------------------------------
  # PLOT RESULTS
  # ----------------------------------
  
  # Summarize results across repeated simulations
  summary_df <- results_df %>%
    group_by(tau_true, phi_true) %>%
    summarize(mean_tau_corrected = mean(tau_est_q50),
              mean_phi_corrected = mean(phi_est_q50),
              mean_Density_corrected = mean(Density_est_q50),
              
              # Credible interval coverage
              tau_cov90 = signif(mean(tau_cov90),2),
              phi_cov90 = signif(mean(phi_cov90),2),
              Density_cov90 = signif(mean(Density_cov90),2),
              
              # Percent bias in each estimate
              tau_bias_percent = mean(100*(tau_est_q50 - tau_true)/tau_true) %>% signif(1),
              phi_bias_percent = mean(100*(phi_est_q50 - phi_true)/phi_true) %>% signif(1),
              Density_bias_percent = mean(100*(Density_est_q50 - Density_true)/Density_true) %>% signif(1)
              
    )
  
  # 
  tau_plot <- ggplot(data = results_df, aes(x = sim_rep, y = tau_est_q50, ymin = tau_est_q05, ymax = tau_est_q95))+
    geom_hline(data = results_df, aes(yintercept = tau_true), size = 2, col = "dodgerblue")+
    geom_errorbar(width=0)+
    geom_point()+
    geom_hline(data = summary_df, aes(yintercept = mean_tau_corrected), linetype = 2)+
    ggtitle(paste0("Tau 90% interval coverage = ",summary_df$tau_cov90))+
    xlab("Simulation #")+
    ylab("tau")+
    theme_bw()
  #tau_plot 
  
  phi_plot <- ggplot(data = results_df, aes(x = sim_rep, y = phi_est_q50, ymin = phi_est_q05, ymax = phi_est_q95))+
    geom_hline(data = results_df, aes(yintercept = phi_true), size = 2, col = "dodgerblue")+
    geom_errorbar(width=0)+
    geom_point()+
    geom_hline(data = summary_df, aes(yintercept = mean_phi_corrected), linetype = 2)+
    ggtitle(paste0("Phi 90% interval coverage = ",summary_df$phi_cov90))+
    xlab("Simulation #")+
    ylab("phi")+
    theme_bw()
  #phi_plot 
  
  Density_plot <- ggplot(data = results_df, aes(x = sim_rep, y = Density_est_q50, ymin = Density_est_q05, ymax = Density_est_q95))+
    geom_hline(data = results_df, aes(yintercept = Density_true), size = 2, col = "dodgerblue")+
    geom_errorbar(width=0)+
    geom_point()+
    geom_hline(data = summary_df, aes(yintercept = mean_Density_corrected), linetype = 2)+
    ggtitle(paste0("Density 90% interval coverage = ",summary_df$Density_cov90))+
    xlab("Simulation #")+
    ylab("Density")+
    theme_bw()
  #Density_plot 
  
  estimate_plot <- ggarrange(tau_plot,phi_plot,Density_plot,nrow=3)
  estimate_plot <- annotate_figure(estimate_plot, top = paste0("Tau = ",tau_true, " , Phi = ", phi_true))
  print(estimate_plot)  
}
