# ******************************************************************
# QPAD CODE FIT WITH JAGS
# - estimates are 'uncorrected' and still biased
# ******************************************************************

library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)
library(jagsUI) # For Bayesian analysis

rm(list=ls())

# ------------------------------------
# QPAD code in JAGS
# ------------------------------------

sink("dist.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    tau ~ dunif(0,5)
    phi ~ dunif(0,2)
    
    # ------------------------------
    # Calculate multinomial cell probabilities for distance modeling
    # ------------------------------
    
    for (i in 1:nbins_distance){
      cdf_rbin[i] <- 1-exp(-(rint[i]/tau)^2)
    }
    
    # Cell probabilities
    p_distance[1] <- cdf_rbin[1]- 0
    for (i in 2:nbins_distance){
      p_distance[i] <- cdf_rbin[i] - cdf_rbin[i-1]
    }
    
    # ------------------------------
    # Calculate multinomial cell probabilities for removal modeling
    # ------------------------------

    for (i in 1:nbins_removal){
      cdf_tbin[i] <- 1-exp(-tint[i]*phi)
    }

    # cell probabilities
    p_removal[1] <- cdf_tbin[1]- 0
    for (i in 2:nbins_removal){
      p_removal[i] <- cdf_tbin[i] - cdf_tbin[i-1]
    }

    # ------------------------------
    # Likelihoods
    # ------------------------------
    
    Y_distance[1,1:nbins_distance] ~ dmulti(p_distance[],N_distance)
    Y_removal[1,1:nbins_removal] ~ dmulti(p_removal[],N_removal)

    }
",fill = TRUE)
sink()


# ------------------------------------
# Conduct repeated simulations
# ------------------------------------

tau = 1
phi = 0.5

results_df <- data.frame()

for (sim_rep in 1:1000){
  
  # ------------------------------------
  # Set up population
  # ------------------------------------
  
  l <- bsims_init(extent = 5)
  p <- bsims_populate(l, density = 10)
  
  # ------------------------------------
  # Animate
  # ------------------------------------
  
  a <- bsims_animate(p, vocal_rate = phi, duration = 10) 
  events <- get_events(a)
  
  # ------------------------------------
  # Which vocalizations are actually detected?
  # ------------------------------------
  
  d <- bsims_detect(a, tau=tau) # 5 heard
  dets <- get_detections(d) 
  
  # ------------------------------------
  # Transcribe
  # ------------------------------------
  
  rint <- c(0.5,1,Inf) # Distance bins
  tint <- c(3,5,10)    # Time bins
  
  x <- bsims_transcribe(d, rint = rint, tint = tint, condition = "det1")
  Y <- get_table(x, "removal")
  
  # **** If using condition = "event1" these numbers are not equal
  length(unique(dets$i)) # number of unique individuals detected
  sum(Y) # Should equal number of unique individuals detected 
  
  # ------------------------------------
  # Analysis
  # ------------------------------------
  
  # Distance model
  Y_distance = matrix(rowSums(Y),1)
  D_distance = matrix(rint,1)
  fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
  tau_MLE = exp(fit.q$coefficients)
  
  # Removal model
  Y_removal = matrix(colSums(Y),1)
  D_removal = matrix(tint,1)
  fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_MLE = exp(fit.p$coefficients)
  
  # ------------------------------------
  # Fit models with JAGS
  # ------------------------------------
  
  jags_data <- list(
    Y_distance = Y_distance,
    N_distance = sum(Y_distance),
    Y_removal = Y_removal,
    N_removal = sum(Y_removal),
    rint = rint,
    nbins_distance = length(rint),
    tint = tint,
    nbins_removal = length(tint)
  )
  
  out <- jags(data = jags_data,
              model.file =  "dist.jags",
              parameters.to.save = c("tau","phi"),
              inits = NULL,
              n.chains = 3,
              n.thin = 1,
              n.iter = 10000,
              n.burnin = 5000,
              parallel = TRUE)
  
  # ------------------------------------
  # Summarize results and append to dataframe
  # ------------------------------------
  
  # Estimate density
  A_hat = pi*tau_MLE^2
  p_hat = 1-exp(-10*phi_MLE)
  D_hat <- sum(Y)/(A_hat*p_hat)
  
  results_df <- rbind(results_df,
                  data.frame(sim_rep = sim_rep,
                             N_detected = sum(Y),
                             D_hat = D_hat,
                             
                             tau_MLE = tau_MLE,
                             tau_q50 = out$q50$tau,
                             tau_q025 = out$q2.5$tau,
                             tau_q975 = out$q97.5$tau,
                             
                             phi_MLE = phi_MLE,
                             phi_q50 = out$q50$phi,
                             phi_q025 = out$q2.5$phi,
                             phi_q975 = out$q97.5$phi ))
  
  # ------------------------------------
  # Plots
  # ------------------------------------
  
  plot_tau <- ggplot(results_df,aes(x = sim_rep, y = tau_q50, ymin = tau_q025, ymax = tau_q975))+
    geom_point()+
    geom_point(data = results_df,aes(x = sim_rep, y = tau_MLE), shape = 4)+
    geom_errorbar(width=0.1)+
    geom_hline(yintercept = tau, size=2,col="dodgerblue",alpha=0.5)+
    geom_hline(yintercept = mean(results_df$tau_q50), linetype = 2)+
    theme_bw()
  #plot_tau
  
  plot_phi <- ggplot(results_df,aes(x = sim_rep, y = phi_q50, ymin = phi_q025, ymax = phi_q975))+
    geom_point()+
    geom_point(data = results_df,aes(x = sim_rep, y = phi_MLE), shape = 4)+
    geom_errorbar(width=0.1)+
    geom_hline(yintercept = phi, size=2,col="dodgerblue",alpha=0.5)+
    geom_hline(yintercept = mean(results_df$phi_q50), linetype = 2)+
    theme_bw()
  #plot_phi
  
  estimate_plot <- ggarrange(plot_tau,plot_phi, nrow=2)
  print(estimate_plot)
}









