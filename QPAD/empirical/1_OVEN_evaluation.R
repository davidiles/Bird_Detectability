library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)
library(jagsUI)
library(reshape2)
library(readxl)

setwd("~/1_Work/Bird_Detectability/QPAD/empirical")

rm(list=ls())


# ----------------------------------------
# Read OVEN data
# ----------------------------------------
file <- "BLPWdetectionhistories.xlsx"
#file <- "OVENdetectionhistory.xlsx"
dat <- read_xlsx(file) %>%
  rename(tbin = 10, dbin = 11)

dat <- dat %>%
  group_by(dbin,tbin) %>%
  summarize(count = n()) %>%
  na.omit()

dat$dbin <- factor(dat$dbin, levels = c("0-49 m","50-100 m",">100 m"))
dat$tbin <- factor(dat$tbin, levels = sort(unique(dat$tbin)))

dat <- dat %>% arrange(dbin,tbin)

# Convert to matrix
tint <- seq(1,10)
rint <- c(50,100,Inf)
Y = dat %>% spread(tbin, count) 
Y <- as.matrix(Y[,2:ncol(Y)])
#colnames(Y) <- tint
rownames(Y) <- levels(dat$dbin)
Y[is.na(Y)] <- 0

# ----------------------------------------
# Fit using standard QPAD
# ----------------------------------------

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

tau_MLE # 92
phi_MLE # 0.39

# ----------------------------------------
# Refit distance model using only first time interval
# ----------------------------------------

# Distance model
Y_distance = matrix(Y[,1],1)
D_distance = matrix(rint,1)

fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
tau_MLE_2 = exp(fit.q$coefficients)
tau_MLE_2 # 82

# ----------------------------------------
# Refit time-removal model using only first distance interval
# ----------------------------------------

# Removal model
Y_removal = matrix(Y[1,],1)
D_removal = matrix(tint,1)

fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
phi_MLE_2 = exp(fit.p$coefficients)

phi_MLE_2 # 0.39

# ----------------------------------------
# Comparisons
# ----------------------------------------


phi_MLE   # Using all distance bins
phi_MLE_2 # Using only nearest distance bin

tau_MLE   # Using all time bins
tau_MLE_2 # Using only first time bin

# ******************************************************
# Analysis with JAGS
# ******************************************************

sink("QPAD_corrected_vectorized.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    phi ~ dunif(0,5)
    tau ~ dunif(0,500)
    
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

# Change to vector format (so dcat() can work in JAGS)
Y_vec <- reshape2::melt(Y) %>%
  rename(rint = Var1, tint = Var2, Y = value)

Y_rint_numeric <- factor(Y_vec$rint) %>% as.numeric()
Y_tint_numeric <- factor(Y_vec$tint) %>% as.numeric()
Y_vectorized <- Y_vec$Y

# --------------------------------------
# Data augmentation (used for integrating across distance bands)
# --------------------------------------

Nsim = 1000 # Choose large number so that *expected* cell probabilities are not affected by sampling variation 
dim = 600     # landscape size (dim x dim landscape)

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
out <- jags(data = jags_data,
            model.file =  "QPAD_corrected_vectorized.jags",
            parameters.to.save = c("phi","tau"),
            n.chains = 3,
            n.thin = 1,
            n.iter = 5000,
            n.burnin = 1000,
            parallel = TRUE)
