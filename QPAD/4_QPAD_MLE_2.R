library(tidyverse)
library(ggpubr)

rm(list=ls())

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

#results_df = data.frame()
#for (sim_rep in 1:100){

# ******************************************
# PART 1: SIMULATE DATA
# ******************************************

nsurvey = 100 # Number of point counts to simulate

# Two covariates affect tau linearly
tau_intercept <- log(1.5)
#tau_beta1 <- -0.05
#tau_beta2 <- 0.1
#x1 <- runif(nsurvey,-3,3)
#x2 <- runif(nsurvey,-3,3)
#tau <- exp(tau_intercept + tau_beta1*x1 + tau_beta2*x2)
#X1 <- model.matrix(~x1+x2)
#colnames(X1) <- c("tau_b0","tau_b1","tau_b2")

# One covariate affects phi, but isn't modeled
phi_intercept <- log(0.2)
#phi_beta1 <- 0.1
#phi_beta2 <- 0
#x1 <- runif(nsurvey,-3,3)
#x2 <- runif(nsurvey,-3,3)
#phi <- exp(phi_intercept + phi_beta1*x1 + phi_beta2*x2)


tau <- rep(exp(tau_intercept),nsurvey)
phi <- rep(exp(phi_intercept),nsurvey)

# -------------------------------------------------
# Covariate design matrices for tau and phi
# -------------------------------------------------
#head(X1)
#head(X2)

# -------------------------------------------------
# Point count data are stored in arrays
# -------------------------------------------------
Dvec <- rep(NA,nsurvey) # True density at each survey location

Yarray <- array(NA,dim=c(nsurvey,mdbin,mtbin))
rarray <- array(NA,dim=c(nsurvey,mdbin))
tarray <- array(NA,dim=c(nsurvey,mtbin))

for (k in 1:nsurvey){
  
  tau_true <- tau[k]
  phi_true <- phi[k]
  
  dim = 10 # landscape x and y dimensions (100 metre increments)
  D <- runif(1,50,100)  # Target density (birds per ha)
  N <- round(D*dim^2)    # Number of birds to place on landscape
  Density_true <- N/dim^2 # birds per ha
  Dvec[k] <- Density_true
  
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

# ******************************************
# PART 3: JOINT ANALYSIS WITH MLE
# ******************************************

setwd("~/1_Work/Bird_Detectability/QPAD")
source("joint_fns.R")
fit <- cmulti.fit.joint(Yarray,
                        rarray,
                        tarray,
                        X1 = NULL, # Implies intercept only model
                        X2 = NULL # Implies intercept only model
                        )

fit$coefficients

# Calculate survey-level offsets
log_offsets <- calculate.offsets(fit)

# Use offsets to estimate true density at each survey location
Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
Dhat <- Ysum*exp(log_offsets)

ggplot()+
  geom_point(aes(x = Dvec,y=Dhat), col = "dodgerblue")+
  geom_abline(slope=1,intercept=0)+
  ylim(range(c(Dvec,Dhat)))+
  xlim(range(c(Dvec,Dhat)))+
  xlab("True density at survey location")+
  ylab("Estimated density at survey location")+
  theme_bw()
  
