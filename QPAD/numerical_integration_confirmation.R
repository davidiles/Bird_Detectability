
# ******************************************
# PART 1: SIMULATE DATA
# ******************************************

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

cues <- matrix(NA, nrow=N, ncol = 20)
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
tint <- c(3,5,10)
nrint <- length(rint)
ntint <- length(tint)

# Separate into distance and time bins
dat$rint <- cut(dat$dist,c(0,rint))
dat$tint <- cut(dat$time,c(0,tint))
dat <- na.omit(dat)

Y <- table(dat[,c("rint","tint")])

Y # Data to analyze


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NUMERICAL INTEGRATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_d = function(dmax){
  integrand = substitute(2*pi*dmax *(1-exp(-phi_true*tmax*exp(-dmax^2/tau_true^2))),
                 list(phi_true = phi_true,
                      tau_true = tau_true,
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
    CDF_binned[i,j] = integrate(f_d,lower=0.0001,upper = upper_r, subdivisions = 1000)$value
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

p_matrix = p_matrix/sum(p_matrix)

# Should be equal (if Y has a large sample size)
p_matrix %>% round(3)
(Y/sum(Y)) %>% round(3)

