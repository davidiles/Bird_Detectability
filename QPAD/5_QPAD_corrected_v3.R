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
      
        # ************************************
        # Expected number of cues produced = phi*tint[j]
        # Expected number of detectable cues = lambda = phi*tint[j]*exp(-(sim_dist[i]/tau)^2)
        # Probability of getting zero detectable cues (under Poisson) = exp(-lambda)
        # Probability of getting at least 1 detectable cue = 1-exp(-lambda)
        # ************************************
        
        # Probability bird produces at least one detectable cue by end of interval (product of expected number of cues and probability each cue is detected)
        p[i,j] <- 1 - exp(-phi*tint[j]*exp(-(sim_dist[i]/tau)^2))
        
      }
    }
    
   # Sum up number of detected birds in each time and distance bin (this is sort of equivalent to the CDF)
    for (i in 1:nrint){
      for (j in 1:ntint){
        ET[i,j] <- sum(p[r_from[i]:r_to[i],j])
      }
    }
   
   # Difference within each distance bin to calculate cell probabilities
    for (i in 1:nrint){
      Ysim_matrix[i,1] <- ET[i,1]
      for (j in 2:ntint){
        Ysim_matrix[i,j] <- ET[i,j] - sum(ET[i,(j-1)])
      }
    } 
    
    # Place counts in correct vectorized bins
    for (i in 1:n_bins_vectorized){
        Ysim_vectorized[i] <- Ysim_matrix[Y_rint_numeric[i],Y_tint_numeric[i]]
      }
      
    # -----------------------------------------------------
    # Use relative proportions for detectability modeling
    # Note that JAGS automatically normalizes Ysim_vectorized (Ysim_vectorized/sum(Ysim_vectorized)) so that probabilities sum to 1
    # -----------------------------------------------------
    
    Y_vectorized[1:n_bins_vectorized] ~ dmulti(Ysim_vectorized[],N_det)
    
  }
    
",fill = TRUE)
sink()

rm(list=ls())

N = 100000 # Number of birds to place on landscape (determines overall sample size in point count dataset)
dim = 10 # landscape size (dim x dim landscape)
Density_true = N/dim^2 

tau_true <- 1
phi_true <- 0.5

# ------------------------------------
# PART 1: GENERATE 'ACTUAL' DATA
# ------------------------------------

birds <- data.frame(x = runif(N,-dim/2,dim/2),
                    y = runif(N,-dim/2,dim/2))

# Distances to observer
birds$dist <- sqrt(birds$x^2 + birds$y^2)

birds <- birds %>% arrange(dist)
bird_id = 1:N

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

# Transcription
tint <- c(3,5,10)
rint <- c(0.5,1,Inf)
dat$tint <- cut(dat$time,c(0,tint))
dat$rint <- cut(dat$dist,c(0,rint))
Y = table(dat$rint,dat$tint)

# --------------------------------------
# Data augmentation (used for integrating across distance bands)
# --------------------------------------

Nsim = 5000 # Choose large number so that *expected* cell probabilities are not affected by sampling variation 
dim = 6     # landscape size (dim x dim landscape)

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


# Change to vector format (so dcat() can work in JAGS)
Y_vec <- reshape2::melt(Y) %>%
  rename(rint = Var1, tint = Var2, Y = value)

Y_rint_numeric <- factor(Y_vec$rint) %>% as.numeric()
Y_tint_numeric <- factor(Y_vec$tint) %>% as.numeric()
Y_vectorized <- Y_vec$Y
Y_vectorized[1] <- Y_vectorized[1] - 20


jags_data <- list(
  
  pi = pi,
  
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
  
  # Vectorized data
  N_det = sum(Y_vectorized),
  Y_vectorized = Y_vectorized,
  n_bins_vectorized = length(Y_rint_numeric),
  Y_rint_numeric = Y_rint_numeric,
  Y_tint_numeric = Y_tint_numeric
)

# Fit model
inits <- function()list(phi = phi_true,tau = tau_true)

# Using initial values
out <- jags(data = jags_data,
            model.file =  "QPAD_corrected_vectorized.jags",
            parameters.to.save = c("phi","tau"),
            #inits = inits,
            n.chains = 3,
            n.thin = 1,
            n.iter = 5000,
            n.burnin = 1000,
            parallel = TRUE)
