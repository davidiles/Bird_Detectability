library(tidyverse)

# This script confirms that the equations used for integration
# correctly describe the cell probabilities for the multinomial
# distribution (i.e., relative proportions of birds in each distance/time bin)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 1: SIMULATE DATA
#   - birds produce cues at rate phi_true
#   - detectability per cue is determined by half-normal function
#   - observers only record the first cue they detect from each bird
#   - observers place birds into discrete time and distance bins
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())


tau_true = 1.5
phi_true = 0.2

N = 500000 # Number of birds to place on landscape (select a high number to provide sufficient sample size)
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

cues <- matrix(NA, nrow=N, ncol = 100)
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

rint <- c(0.5,1,2,Inf)
tint <- seq(1,5,1)
nrint <- length(rint)
ntint <- length(tint)

# Separate into distance and time bins
dat$rint <- cut(dat$dist,c(0,rint))
dat$tint <- cut(dat$time,c(0,tint))
dat <- na.omit(dat)

Y <- table(dat[,c("rint","tint")])

Y # Data to analyze


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 2: NUMERICAL INTEGRATION
# If the equation below is correct, we should be able to recreate
# the relative proportions of birds in each bin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_d = function(dmax){
  integrand = substitute(2*pi*dmax *(1-exp(-phi_true*tmax*exp(-dmax^2/tau_true^2))),
                 list(phi_true = phi_true,
                      tau_true = tau_true,
                      tmax = tmax))
  
  eval(integrand)
  
}

# Equation 3 (dropped density term)
max_dist <- max(rint)
if (max_dist==Inf) max_dist = 10
tmax <- max(tint)
ptot <- integrate(f_d,lower=0.0001,upper = max_dist, subdivisions = 1000)$value

# Calculate cell probabilities
p_matrix <- Y*NA
rint2 <- c(0,rint)
for (i in 2:length(rint2)){
  tmax <- tint[1]
  p <- integrate(f_d,lower=rint2[i-1]+0.0001,upper = rint2[i], subdivisions = 1000)$value
  p_matrix[i-1,1] <- p/ptot
}

for (j in 2:length(tint)){
  for (i in 2:length(rint2)){
    tmax <- tint[j-1]
    pt1 <- integrate(f_d,lower=rint2[i-1]+0.0001,upper = rint2[i], subdivisions = 1000)$value
    tmax <- tint[j]
    pt2 <- integrate(f_d,lower=rint2[i-1]+0.0001,upper = rint2[i], subdivisions = 1000)$value
    
    p_matrix[i-1,j] = (pt2-pt1)/ptot
  }
}

# Should be approximately equal if equation/integral is correct
p_matrix %>% round(3)
(Y/sum(Y)) %>% round(3)

