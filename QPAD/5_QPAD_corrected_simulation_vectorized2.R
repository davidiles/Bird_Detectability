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
    dummy ~ dunif(0,1)
    #phi ~ dunif(0,2)
    #tau ~ dunif(0,3)
    #phi <- 0.5
    #tau <- 1
    
    # ------------------------------
    # Simulation
    # ------------------------------
    
    # For each simulated bird, probability it has generated at least one detectable cue by end of time interval
    for (i in 1:Nsim){
    
      # Simulate cueing for each bird in each time interval
      for (j in 1:ntint){
      
        # Expected number of cues produced up to end of interval
        lambda[i,j] <- phi*tint[j]*exp(-(sim_dist[i]/tau)^2)
        
        # Probability bird produces at least one detectable cue by end of interval
        p[i,j] <- 1 - exp(-lambda[i,j])
      }
    }
    
   # Sum up number of detected birds in each time and distance bin (this is sort of equivalent to the CDF)
    for (i in 1:nrint){
      for (j in 1:ntint){
        ET[i,j] <- sum(p[r_from[i]:r_to[i],j])
      }
    }
   
   # Difference to calculate cell probabilities
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
    # -----------------------------------------------------
    
    #Y_vectorized[1:n_bins_vectorized] ~ dmulti(Ysim_vectorized[],N_det)
    
  }
    
",fill = TRUE)
sink()

rm(list=ls())

results_df = data.frame()

#for (sim_rep in 1:100){

#  for (tau_true in c(0.5,1.5)){
#    for (phi_true in c(0.2,0.4)){

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

# Transcription: when distance is fixed
tint <- c(3,5,10)
rint <- c(0.5,1,Inf)

dat$tint <- cut(dat$time,c(0,tint))
dat$rint <- cut(dat$dist,c(0,rint))
Y = table(dat$rint,dat$tint)

# ---------------------------------------
# Exploratory
# ---------------------------------------

# For each bird, calculate probability it has been detected by end of each time interval
birds$p_t3 <-  1 - exp(-phi_true*3*exp(-(birds$dist/tau_true)^2))
birds$p_t5 <-  1 - exp(-phi_true*5*exp(-(birds$dist/tau_true)^2))
birds$p_t10 <-  1 - exp(-phi_true*10*exp(-(birds$dist/tau_true)^2))

birds$rint <- cut(birds$dist,c(0,rint))

# Expected totals
ET = birds %>%
  group_by(rint) %>%
  summarize(N_t3 = sum(p_t3),
            N_t5 = sum(p_t5),
            N_t10 = sum(p_t10))

ET = ET[2:ncol(ET)] %>% as.matrix()
rownames(ET) <- rint

# These are equivalent; this is working
ET[,3]
rowSums(Y)

# Calculate cell probabilities
p_matrix <- Y*NA
for (i in 1:nrow(ET)){
  p_matrix[i,1] <- ET[i,1]
  for (j in 2:ncol(ET)){
    p_matrix[i,j] <- ET[i,j] - sum(ET[i,(j-1)])
  }
} 

Y
p_matrix
Y/sum(Y)
p_matrix/sum(p_matrix)
# This works


# -----------------------------------------
# So why not working in JAGS???
# -----------------------------------------

Nsim <- N
ntint <- length(tint)
nrint <- length(rint)
sim_dist <- birds$dist
tau <- tau_true
phi <- phi_true

# Group simulated birds into distance bins
birds$rint <- cut(birds$dist,c(0,rint))
birds$rint_num <- birds$rint %>% factor() %>% as.numeric()

# For each distance bin, store first and last bird in that bin
r_from <- r_to <- rint*NA
for (i in 1:length(rint)){
  r_from[i] <- min(which(birds$rint_num == i))
  r_to[i] <- max(which(birds$rint_num == i))
}




lambda <- p <- matrix(NA,nrow = Nsim, ncol = ntint)
# For each simulated bird, probability it has generated at least one detectable cue by end of time interval
for (i in 1:Nsim){
  
  # Simulate cueing for each bird in each time interval
  for (j in 1:ntint){
    
    # Expected number of cues produced up to end of interval
    lambda[i,j] <- phi*tint[j]*exp(-(sim_dist[i]/tau)^2)
    
    # Probability bird produces at least one detectable cue by end of interval
    p[i,j] <- 1 - exp(-lambda[i,j])
  }
}

# Sum up number of detected birds in each time and distance bin (this is sort of equivalent to the CDF)
ET_test<- Y*NA
for (i in 1:nrint){
  for (j in 1:ntint){
    ET_test[i,j] <- sum(p[r_from[i]:r_to[i],j])
  }
}
ET
ET_test

# Difference to calculate cell probabilities
p_matrix_test <- ET*NA
for (i in 1:nrint){
  p_matrix_test[i,1] <- ET_test[i,1]
  for (j in 2:ntint){
    p_matrix_test[i,j] <- ET_test[i,j] - sum(ET_test[i,(j-1)])
  }
} 

# Place counts in correct vectorized bins

for (i in 1:n_bins_vectorized){
  Ysim_vectorized[i] <- Ysim_matrix[Y_rint_numeric[i],Y_tint_numeric[i]]
}








# Determine number of birds detected by end of each time interval
cues$tint <- cut(cues$time,c(0,tint)) 
cues$rint <- cut(cues$dist,c(0,rint))
cues <- na.omit(cues)

n_t3 <- cues %>%
  subset(time <= 3) %>%
  group_by(bird_id,rint) %>%
  summarize(n_cues_detected = sum(detected)) %>%
  
  # Remove birds that were not detected
  subset(n_cues_detected > 0) %>%
  group_by(rint) %>%
  summarize(n = n())
  
  

# Change to vector format (so dcat() can work in JAGS)
Y_vec <- reshape2::melt(Y) %>%
  rename(rint = Var1, tint = Var2, Y = value)

Y_rint_numeric <- factor(Y_vec$rint) %>% as.numeric()
Y_tint_numeric <- factor(Y_vec$tint) %>% as.numeric()
Y_vectorized <- Y_vec$Y
Y_vectorized[1] <- Y_vectorized[1] - 20

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
# ------------------------------------

# Plausible initial values
phi_init = 0.2
tau_init = 1.2

# Data augmentation
Nsim = 50000 # Number of sim_birds to place on landscape (size of augmented dataset)
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
inits <- function()list(phi = phi_true,
                        tau = tau_true)

# Using initial values
out <- jags(data = jags_data,
            model.file =  "QPAD_corrected_vectorized.jags",
            parameters.to.save = c("phi","tau","Ysim_matrix","p","dummy"),
            #inits = inits,
            n.chains = 1,
            n.thin = 1,
            n.iter = 2,
            n.burnin = 1,
            parallel = TRUE,
            codaOnly=c("phi","tau","Ysim_matrix","p"))


