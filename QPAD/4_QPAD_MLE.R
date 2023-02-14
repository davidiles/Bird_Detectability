library(tidyverse)
library(ggpubr)

rm(list=ls())

# ------------------------------------------------------------
# Multinomial loglikelihood
# ------------------------------------------------------------
# x = counts in each cell; size = total; prob = cell probs
logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}

# ------------------------------------------------------------
# Function to calculate negative log likelihood given tau and phi
# ------------------------------------------------------------

nll_fn <- function(tau_est, phi_est){
  
  # Calculate CDF and p
  f_d = function(dmax){
    
    integrand = substitute(2*pi*dmax *(1-exp(-phi*tmax*exp(-dmax^2/tau^2))),
                           list(phi = phi_est,
                                tau = tau_est,
                                tmax = tmax))
    
    eval(integrand)
    
  }
  
  # Calculate CDF
  CDF_binned <- Y*NA
  for (j in 1:length(tint)){
    tmax = tint[j]
    for (i in 1:length(rint)){
      upper_r = rint[i]
      if (upper_r == Inf) upper_r = max_dist
      CDF_binned[i,j] = integrate(f_d,lower=0.0001,
                                  upper = upper_r, 
                                  subdivisions = 1000)$value
    }
  }
  
  # Difference to calculate multinomial cell probabilities
  tmp1 = CDF_binned
  for (i in 2:nrint){
    tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
  }
  
  p_matrix = tmp1
  for (j in 2:ntint){
    p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
  }
  
  # Normalize
  p_matrix = p_matrix/sum(p_matrix)
  
  # Negative log likelihood
  nll <- -1*logdmultinom(Y, Ysum, p_matrix)
  
  # Check for convergence
  if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2]
  else nll
}

# ------------------------------------------------------------
# Function to calculate offset for density estimation
# ------------------------------------------------------------

offset_fn <- function(tau_est, phi_est){
  
  # Calculate CDF and p
  f_d = function(dmax){
    
    integrand = substitute(2*pi*dmax *(1-exp(-phi*tmax*exp(-dmax^2/tau^2))),
                           list(phi = phi_est,
                                tau = tau_est,
                                tmax = tmax))
    
    eval(integrand)
    
  }
  # Calculate CDF
  CDF_binned <- Y*NA
  for (j in 1:length(tint)){
    tmax = tint[j]
    for (i in 1:length(rint)){
      upper_r = rint[i]
      if (upper_r == Inf) upper_r = max_dist
      CDF_binned[i,j] = integrate(f_d,lower=0.0001,
                                  upper = upper_r, 
                                  subdivisions = 1000)$value
    }
  }
  
  # Difference to calculate multinomial cell probabilities
  tmp1 = CDF_binned
  for (i in 2:nrint){
    tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
  }
  
  p_matrix = tmp1
  for (j in 2:ntint){
    p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
  }
  
  C <- sum(p_matrix)
  C
}






results_df = data.frame()
for (sim_rep in 1:100){
  
  # ******************************************
  # PART 1: SIMULATE DATA
  # ******************************************
  
  set.seed(sim_rep)
  
  tau_true = 1.5
  phi_true = 3
  
  N = 10000 # Number of birds to place on landscape (select a high number to provide sufficient sample size)
  dim = 10 # landscape size (metres)
  
  Density_true <- N/dim^2
  Density_true
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
  
  rint <- c(0.5,1,Inf)
  tint <- c(3,5,10)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # Separate into distance and time bins
  dat$rint <- cut(dat$dist,c(0,rint))
  dat$tint <- cut(dat$time,c(0,tint))
  dat <- na.omit(dat)
  
  Y <- table(dat[,c("rint","tint")])
  
  Y # Data to analyze
  Ysum = sum(Y)
  
  # ******************************************
  # PART 3: JOINT ANALYSIS WITH MLE
  # ******************************************
  
  # Maximum distance for integration
  if (max(rint) != Inf) max_dist = max(rint) else max_dist = dim/2
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  # Maximum likelihood estimation
  res <- suppressWarnings(stats4::mle(nll_fn, 
                                      start = list(tau_est = 1, phi_est = 1),
                                      lower = list(tau_est = 0.01, phi_est = 0.01)))
  
  tau_MLE <- coef(res)[1]
  phi_MLE <- coef(res)[2]
  
  # ----------------------------------
  # Calculate correction factor (i.e., offset) to estimate density
  # ----------------------------------
  
  C <- offset_fn(tau_MLE,phi_MLE)/(pi*max_dist^2)
  Y_corrected <- sum(Y)/C
  Density_MLE <- Y_corrected/(pi*max_dist^2)
  Density_MLE
  
  
  # ----------------------------------
  # STORE RESULTS FOR THIS SIMULATION
  # ----------------------------------
  
  results_df = rbind(results_df,
                     data.frame(sim_rep = sim_rep,
                                
                                tau_true = tau_true,
                                tau_MLE = tau_MLE,
                                phi_true = phi_true,
                                phi_MLE = phi_MLE,
                                Density_true = Density_true,
                                Density_MLE = Density_MLE
                     ))
  
  print(sim_rep)
}

# ----------------------------------
# PLOT RESULTS
# ----------------------------------

# Summarize results across repeated simulations
summary_df <- results_df %>%
  group_by(tau_true, phi_true) %>%
  summarize(mean_tau_corrected = mean(tau_MLE),
            mean_phi_corrected = mean(phi_MLE),
            mean_Density_corrected = mean(Density_MLE),
            
            # Percent bias in each estimate
            tau_bias_percent = mean(100*(tau_MLE - tau_true)/tau_true) %>% signif(1),
            phi_bias_percent = mean(100*(phi_MLE - phi_true)/phi_true) %>% signif(1),
            Density_bias_percent = mean(100*(Density_MLE - Density_true)/Density_true) %>% signif(1)
            
  )

tau_plot <- ggplot(data = results_df, aes(x = sim_rep, y = tau_MLE))+
  geom_hline(data = results_df, aes(yintercept = tau_true), size = 2, col = "dodgerblue")+
  
  geom_point()+
  geom_hline(data = summary_df, aes(yintercept = mean_tau_corrected), linetype = 2)+
  ggtitle("Tau")+
  xlab("Simulation #")+
  ylab("tau")+
  theme_bw()
#tau_plot 

phi_plot <- ggplot(data = results_df, aes(x = sim_rep, y = phi_MLE))+
  geom_hline(data = results_df, aes(yintercept = phi_true), size = 2, col = "dodgerblue")+
  
  geom_point()+
  geom_hline(data = summary_df, aes(yintercept = mean_phi_corrected), linetype = 2)+
  ggtitle("Phi")+
  xlab("Simulation #")+
  ylab("phi")+
  theme_bw()
#phi_plot 
# 
Density_plot <- ggplot(data = results_df, aes(x = sim_rep, y = Density_MLE))+
  geom_hline(data = results_df, aes(yintercept = Density_true), size = 2, col = "dodgerblue")+
  
  geom_point()+
  geom_hline(data = summary_df, aes(yintercept = mean_Density_corrected), linetype = 2)+
  ggtitle("Density")+
  xlab("Simulation #")+
  ylab("Density")+
  theme_bw()
#Density_plot 

estimate_plot <- ggarrange(tau_plot,phi_plot,Density_plot, nrow=3)
estimate_plot <- annotate_figure(estimate_plot, top = paste0("DETECTABILITY PARAMETERS:\nTau = ",tau_true, " , Phi = ", phi_true,"\n\nTRANSCRIPTION:\nrint = ",paste(rint, collapse = " "),"\ntint = ",paste(tint, collapse = " ")))
print(estimate_plot)  
