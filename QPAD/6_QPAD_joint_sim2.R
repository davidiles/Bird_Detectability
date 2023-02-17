# **************************************
# Repeated simulations across different values of tau and phi
#
# Assess independence of estimates via their correlation
# **************************************

library(tidyverse)
library(ggpubr)

rm(list=ls())


setwd("~/1_Work/Bird_Detectability/QPAD") # <- set to wherever scripts are stored
source("joint_fns.R")

result_df <- expand.grid(sim_rep = 1:200,
                         tau = c(0.5,1,1.5),
                         phi = c(0.2,1,5),
                         Density = c(0.1,0.5,1),
                         Ysum = NA,
                         tau_est = NA,
                         phi_est = NA,
                         log_offset = NA,
                         Density_est = NA)

for (sim_rep in 1:nrow(result_df)){
  
  set.seed(sim_rep)
  
  Density <- result_df$Density[sim_rep]
  tau <- result_df$tau[sim_rep]
  phi <- result_df$phi[sim_rep]
  
  # Number of point counts / survey locations to simulate
  nsurvey = 1 
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  # A few different protocols for distance binning
  distance_protocols <- list(p1 = c(0.5,1,Inf))
  
  # A few different protocols for time binning
  time_protocols <- list(p1 = seq(1,10))
  
  # Maximum number of distance bins
  mdbin <- sapply(distance_protocols,function(x) length(x)) %>% max()
  
  # Maximum number of time bins
  mtbin <- sapply(time_protocols,function(x) length(x)) %>% max()
  
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  
  Yarray <- array(NA,dim=c(nsurvey,mdbin,mtbin))
  rarray <- array(NA,dim=c(nsurvey,mdbin))
  tarray <- array(NA,dim=c(nsurvey,mtbin))
  
  for (k in 1:nsurvey){
    
    # Parameters for this survey
    tau_true <- tau
    phi_true <- phi
    Density_true <- Density*1000
    
    # ------------------------------------
    # Place birds on landscape around observer (centred on landscape)
    # ------------------------------------
    dim <- 10 # landscape x and y dimensions (100 metre increments)
    N <- round(Density_true*dim^2) # Number of birds to place on landscape
    
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
    
    rint <- sample(distance_protocols,1)[[1]]
    tint <- sample(time_protocols,1)[[1]]
    nrint <- length(rint)
    ntint <- length(tint)
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y <- table(dat[,c("rint","tint")])
    Y # Data to analyze
    
    
    Yarray[k,1:nrint,1:ntint] <- Y
    rarray[k,1:length(rint)] <- rint
    tarray[k,1:length(tint)] <- tint
    
  }
  
  # ******************************************
  # FIT MODEL TO SIMULATED DATA (only first 1000 point counts)
  # ******************************************
  
  start <- Sys.time()
  fit <- cmulti.fit.joint(Yarray,
                          rarray,
                          tarray,
                          X1 = NULL, # Design matrix for tau
                          X2 = NULL  # Design matrix for phi
  )
  end <- Sys.time()
  print(end-start)
  
  # ******************************************
  # Extract/inspect estimates
  # ******************************************
  
  # Calculate survey-level offsets
  log_offset <- calculate.offsets(fit,
                                   rarray = rarray,
                                   tarray = tarray,
                                   X1 = NULL,
                                   X2 = NULL)
  
  # Estimates
  result_df$Ysum[sim_rep] <- sum(Yarray)
  result_df$tau_est[sim_rep] <- exp(fit$coefficients[1])
  result_df$phi_est[sim_rep] <- exp(fit$coefficients[2])
  result_df$log_offset[sim_rep] <- log_offset
  
  result_df$Density_est[sim_rep] <- sum(Yarray)/exp(log_offset)
  print(sim_rep)
}


# -------------------------------------
# Plot results
# -------------------------------------

result_df$tau_label = paste0("tau = ",result_df$tau)
result_df$phi_label = paste0("phi = ",result_df$phi)

# Confirms that estimates of tau and phi are highly correlated when tau is low
ggplot(data = result_df)+
  geom_point(aes(x = tau_est, y = phi_est))+
  facet_wrap(. ~phi_label+tau_label, scales = "free")+
  theme_bw()

# Is density correlated with phi estimate?
ggplot(data = result_df)+
  geom_hline(yintercept = Density)+
  geom_point(aes(x = phi_est, y = Density_est))+
  facet_wrap(. ~phi_label+tau_label, scales = "free")+
  theme_bw()

# Is density correlated with tau estimate?
ggplot(data = result_df)+
  geom_hline(yintercept = Density)+
  geom_point(aes(x = tau_est, y = Density_est))+
  facet_wrap(. ~phi_label+tau_label, scales = "free")+
  theme_bw()
