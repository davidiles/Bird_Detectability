# **************************************
# Simulate 1000 point counts that were collected under 30 different protocols
#
# Survey-level density has a quadratic relationship with a 'mean annual temperature' covariate
#
# Tau is negatively affected by a 'forest cover' covariate
# Phi has a quadratic relationship with a 'day of year' covariate
# **************************************

library(tidyverse)
library(ggpubr)

rm(list=ls())

setwd("~/1_Work/Bird_Detectability/QPAD") # <- set to wherever scripts are stored
source("joint_fns.R")

# ----------------------------------------------------------
# Specify density relationship
# ----------------------------------------------------------

# Number of point counts / survey locations to simulate
nsurvey = 1000 

# Mean annual temperature at each survey
covariate.MAT <- runif(nsurvey,0,25)
Density <- exp(log(0.2) + 0.15*covariate.MAT -0.008*covariate.MAT^2)

plot(Density~covariate.MAT)

# ----------------------------------------------------------
# Simulate relationship between tau and forest cover
# ----------------------------------------------------------

covariate.FC <- plogis(rnorm(nsurvey,1,2))
tau_betas <- c(log(0.7),-0.5)

# Tau for each survey
tau <- exp(tau_betas[1] + tau_betas[2]*covariate.FC) 
plot(tau~covariate.FC)

# Design matrix for tau
zFC <- scale(covariate.FC)
X1 <- model.matrix(~zFC)
colnames(X1) <- c("tau_int","tau_b1")

# ----------------------------------------------------------
# Simulate relationship between phi and day-of-year covariate
# ----------------------------------------------------------

covariate.DOY <- round(runif(nsurvey,120,160))
phi_betas <- c(-15.5,0.22,-0.0008)

# Phi for each survey
phi <- exp(phi_betas[1] + phi_betas[2]*covariate.DOY + phi_betas[3]*covariate.DOY^2) 
plot(phi~covariate.DOY)

# Design matrix for phi (scaled covariate for better model convergence)
zDOY <- scale(covariate.DOY)
X2 <- model.matrix(~zDOY + I(zDOY^2))
colnames(X2) <- c("phi_int","phi_b1","phi_b2")


# ----------------------------------------------------------
# Simulate data collection at each point count
# ----------------------------------------------------------

# A few different protocols for distance binning
distance_protocols <- list(p1 = c(0.5,1,Inf),
                           p2 = c(0.5,1,2),
                           p3 = c(0.5,1),
                           p4 = c(0.5,Inf),
                           p5 = Inf,
                           p6 = 1)

# A few different protocols for time binning
time_protocols <- list(p1 = c(3,5,10),
                       p2 = 3,
                       p3 = 5,
                       p4 = 10,
                       p5 = seq(1,10,1))

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
  tau_true <- tau[k]
  phi_true <- phi[k]
  Density_true <- Density[k]
  
  
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
  
  # Randomly select upper 
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
  print(k)
  
  rarray[k,1:length(rint)] <- rint
  tarray[k,1:length(tint)] <- tint
  
}


Ysum <- apply(Yarray,1,sum,na.rm = TRUE) # Total counts at each point count location

# ******************************************
# FIT MODEL TO SIMULATED DATA
# ******************************************

fit <- cmulti.fit.joint(Yarray,
                        rarray,
                        tarray,
                        X1 = X1, # Design matrix for tau
                        X2 = X2  # Design matrix for phi
)

# ******************************************
# Extract/inspect estimates
# ******************************************
fit$coefficients

tau_est <- exp(X1 %*% fit$coefficients[1:2])

plot(tau_est ~ covariate.FC, col = "dodgerblue", pch = 19) # Estimate
points(tau~ covariate.FC, pch = 19) # Truth

# -----------------
# Estimates of phi
# -----------------

plot(phi_est ~ covariate.DOY, col = "dodgerblue", pch = 19) # Estimate
points(phi~ covariate.DOY, pch = 19) # Truth


# ******************************************
# Calculate detectability offsets for each survey
# ******************************************

# Calculate survey-level offsets
log_offsets <- calculate.offsets(fit)

# ******************************************
# Estimate density at each survey location, using GLM with effect of mean annual temperature
# ******************************************

zMAT <- scale(covariate.MAT) # z-standardize
glm1 <- glm(Ysum ~ zMAT + I(zMAT^2) + offset(log_offsets), family = poisson(link="log"))

# Predict density
bhat <- coef(glm1)
Dhat <- exp(bhat[1] + bhat[2]*zMAT + bhat[3]*zMAT^2)

# Estimated density at each survey location (after correcting for detectability)
ggplot()+
  geom_point(aes(x = covariate.MAT,y=Dhat, col = "Estimate"))+
  geom_line(aes(x = covariate.MAT,y=Density, col = "Truth"))+
  xlab("Mean Annual Temperature")+
  ylab("Density (birds/ha)")+
  scale_color_manual(values=c("dodgerblue","black"), name = "")+
  ggtitle("Predicted vs True Density")+
  theme_bw()

