library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)

set.seed(999)

tau = 1
phi = 0.5

results_df <- data.frame()

# ------------------------------------
# Conduct repeated simulations
# ------------------------------------

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
  # Summarize results and append to dataframe
  # ------------------------------------
  
  results_df <- rbind(results_df, data.frame(sim_rep = sim_rep,
                                             tau_MLE = tau_MLE,
                                             phi_MLE = phi_MLE))
  
  # ------------------------------------
  # Plot results
  # ------------------------------------
  
  plot_tau <- ggplot(results_df,aes(x = sim_rep, y = tau_MLE))+
    geom_point()+
    geom_hline(yintercept = tau, size=2,col="dodgerblue",alpha=0.5)+
    geom_hline(yintercept = mean(results_df$tau_MLE), linetype = 2)+
    theme_bw()
  
  plot_phi <- ggplot(results_df,aes(x = sim_rep, y = phi_MLE))+
    geom_point()+
    geom_hline(yintercept = phi, size=2,col="dodgerblue",alpha=0.5)+
    geom_hline(yintercept = mean(results_df$phi_MLE), linetype = 2)+
    theme_bw()
  
  estimate_plot <- ggarrange(plot_tau,plot_phi, nrow=2)
  print(estimate_plot)
}




