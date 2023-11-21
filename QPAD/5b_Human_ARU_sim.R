
# **************************************
# Simulate 10000 point counts that were collected under 30 different protocols
#
# Survey-level density has a quadratic relationship with a 'mean annual temperature' covariate
#
# Tau is negatively affected by a 'forest cover' covariate
# Phi has a quadratic relationship with a 'day of year' covariate
# **************************************

library(tidyverse)
library(ggpubr)

# Install/load David Hope's C++ version of the 'joint model'

##### remotes::install_github("dhope/CmultiJoint.dev")
library(CmultiJoint.dev) 

rm(list=ls())

set.seed(999)

# Number of point counts / survey locations to simulate
nsurvey = 10000 

# ----------------------------------------------------------
# Specify relationship between density and mean annual temperature
# ----------------------------------------------------------

# Mean annual temperature at each survey
covariate.MAT <- runif(nsurvey,0,25)

# Parameters controlling the quadratic
A <- -0.01     # steepness
k <- log(1)  # maximum of quadratic  
h <- 10        # x location of vertex

# Phi for each survey
Density <- exp(A*(covariate.MAT-h)^2 + k)

plot(Density~covariate.MAT)

# ----------------------------------------------------------
# Simulate relationship between tau and forest cover
# ----------------------------------------------------------

covariate.FC <- plogis(rnorm(nsurvey,1,2))
tau_betas <- c(log(1),-1)

# Tau for each survey
tau <- exp(tau_betas[1] + tau_betas[2]*covariate.FC) 
plot(tau~covariate.FC, ylim=c(0,max(tau)))

# ----------------------------------------------------------
# Simulate relationship between phi and day-of-year covariate
# ----------------------------------------------------------

covariate.DOY <- round(runif(nsurvey,120,160))

# Parameters controlling the quadratic
A <- -0.002     # steepness
k <- log(1.5)  # maximum of quadratic  
h <- 130       # x location of vertex

# Phi for each survey
phi <- exp(A*(covariate.DOY-h)^2 + k)
plot(phi~covariate.DOY, ylim=c(0,max(phi)))

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

# Maximum number of distance bins in the dataset
mdbin <- sapply(distance_protocols,function(x) length(x)) %>% max()

# Maximum number of time bins in the dataset
mtbin <- sapply(time_protocols,function(x) length(x)) %>% max()

# -------------------------------------------------
# Data from each point count is stored in a series of arrays
# -------------------------------------------------

Yarray <- array(NA,dim=c(nsurvey,mdbin,mtbin))
rarray <- array(NA,dim=c(nsurvey,mdbin))
tarray <- array(NA,dim=c(nsurvey,mtbin))

# Simulate point counts (takes a little while for large sample sizes and/or high densities)
for (k in 1:nsurvey){
  
  # Parameters for this survey
  tau_true <- tau[k]
  phi_true <- phi[k]
  Density_true <- Density[k]
  
  # Determine if this will be a human or ARU survey
  ARU <- sample(c(TRUE,FALSE),1,prob=c(0.5,0.5))
  if (ARU){
    rint <- Inf
    tint <- seq(1,10)
  } else{
    rint <- c(0.5,1,Inf)
    tint <- seq(1,10)
  }
  
  # Randomly select sampling protocol
  nrint <- length(rint)
  ntint <- length(tint)
  
  rarray[k,1:length(rint)] <- rint
  tarray[k,1:length(tint)] <- tint
  
  # ------------------------------------
  # Place birds on landscape around observer (centred on landscape)
  # ------------------------------------
  dim <- 10 # landscape x and y dimensions (100 metre increments)
  N <- rpois(1,Density_true*dim^2) # Number of birds to place on landscape
  
  if (N>0){
    birds <- data.frame(x = runif(N,-dim/2,dim/2),
                        y = runif(N,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    
    # ------------------------------------
    # Simulate bird cues, based on phi_true
    # ------------------------------------
    
    cues <- matrix(NA, nrow=N, ncol = 100)
    
    # For each bird, simulate time between cues
    for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_true))
    
    # Convert to a dataframe
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    # Distance of each cue from the observer
    cues$dist <- birds$dist[cues$bird_id]
    
    # ------------------------------------
    # Determine which cues are detected, based on tau_true
    # ------------------------------------
    
    cues$p<- exp(-(cues$dist/tau_true)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Isolate first detected cue for each bird
    # ------------------------------------
    
    # Only select cues that were detected
    dat <- subset(cues,detected == 1)
    
    # Only select the first cue detected for each bird
    dat <- dat[!duplicated(dat$bird_id),]
    
    # ------------------------------------
    # Transcribe the observed data
    # ------------------------------------
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y <- table(dat[,c("rint","tint")])
    Yarray[k,1:nrint,1:ntint] <- Y
    
  } else{
    Yarray[k,1:nrint,1:ntint] <- 0
  }
  
}

# Total number of birds detected at each point count
Ysum <- apply(Yarray,1,sum,na.rm=TRUE)

# Number of point counts that detected at least one bird
sum(Ysum>0)

# Take a look at a few simulated point counts

Yarray[1,,] # Data from first point count
rarray[1,]  # Distance bins used for the first point count
tarray[1,]  # Time bins used for the first point count

Yarray[2,,] # Data from 4th point count
rarray[2,]  # Distance bins used for the 4th point count
tarray[2,]  # Time bins used for the 4th point count

# -------------------------------------------------
# Create covariate design matrices
# -------------------------------------------------

# Design matrix for tau (scaled covariate for better model convergence)
zFC <- scale(covariate.FC)
X1 <- model.matrix(~zFC)
colnames(X1) <- c("tau_int","tau_b1")

# Design matrix for phi (scaled covariate for better model convergence)
zDOY <- scale(covariate.DOY)
X2 <- model.matrix(~zDOY + I(zDOY^2))
colnames(X2) <- c("phi_int","phi_b1","phi_b2")


# -------------------------------------------------
# Fit the joint model
# -------------------------------------------------

start <- Sys.time()
fit <- cmulti_fit_joint(Yarray,
                        rarray,
                        tarray,
                        X1 = NULL, # Design matrix for tau
                        X2 = NULL # Design matrix for phi
)
end <- Sys.time()
print(end-start) # 1.7 minutes on my personal laptop

coef_estimates <- fit$coefficients

# Coefficient estimates for tau (on log scale)
tau_coefficients <- coef_estimates[1:ncol(X1)]
names(tau_coefficients) <- colnames(X1)
tau_coefficients

# Coefficient estimates for phi (on log scale)
phi_coefficients <- coef_estimates[(ncol(X1)+1):length(coef_estimates)]
names(phi_coefficients) <- colnames(X2)
phi_coefficients

# -------------------------------------------------
# Generate predictions from the fitted model
# -------------------------------------------------

# Predictions of tau across range of forest cover
tau_hat <- exp(X1 %*% tau_coefficients)
tau_plot <- ggplot()+
  geom_line(aes(x = covariate.FC, y = tau_hat, col = "Estimate"))+
  geom_line(aes(x = covariate.FC, y = tau, col = "Truth"))+
  scale_color_manual(values=c("dodgerblue","black"), name = "")+
  ylab(expression(tau))+
  theme_bw()+
  ggtitle(expression(Estimate~of~tau))
print(tau_plot)

# Predictions of phi across days of the year
phi_hat <- exp(X2 %*% phi_coefficients)
phi_plot <- ggplot()+
  geom_line(aes(x = covariate.DOY, y = phi_hat, col = "Estimate"))+
  geom_line(aes(x = covariate.DOY, y = phi, col = "Truth"))+
  scale_color_manual(values=c("dodgerblue","black"), name = "")+
  ylab(expression(phi))+
  theme_bw()+
  ggtitle(expression(Estimate~of~phi))
print(phi_plot)

# -------------------------------------------------
# Generate detectability offsets from fitted model
# -------------------------------------------------

log_offsets <- calculate_offsets(fit,
                                 rarray = rarray,
                                 tarray = tarray,
                                 X1 = X1,
                                 X2 = X2)

length(log_offsets) # One for each survey

# -------------------------------------------------
# Use detectability offsets in a density model
# -------------------------------------------------

dat <- data.frame(Y = apply(Yarray,1,sum,na.rm = TRUE),
                  MAT = covariate.MAT,
                  log_off = log_offsets)
head(dat)

glm1 <- glm(Y ~ MAT + I(MAT^2) + offset(log_off), family = poisson(link="log"), data = dat)

# Generate predictions of density across values of MAT
pred_df <- data.frame(MAT = sort(covariate.MAT),log_off = 0)
pred_df$D_hat <- predict(glm1, newdata = pred_df, type = "response")

Density_plot <- ggplot()+
  geom_line(data = pred_df, aes(x = MAT, y = D_hat, col = "Estimate"))+
  geom_line(aes(x = covariate.MAT, y = Density, col = "Truth"))+
  scale_color_manual(values=c("dodgerblue","black"), name = "")+
  ylab("Density")+
  theme_bw()+
  ggtitle("Estimate of Density")
print(Density_plot)

# ***********************************************************
# Fit with "independent" QPAD approach
# ***********************************************************



# -------------------------------------------------
# Bootstrapping if uncertainty in estimates is needed
# -------------------------------------------------

# bootreps <- 1000
# Dhat_boot <- tau_boot <- phi_boot <- matrix(NA,nrow=bootreps,ncol=500)
# 
# for (b in 1:bootreps){
#   
#   bootsamps <- sample(1:nsurvey,nsurvey,replace=TRUE)
#   Yboot <- Yarray[bootsamps,,]
#   rboot <- rarray[bootsamps,]
#   tboot <- tarray[bootsamps,]
#   X1boot <- X1[bootsamps,]
#   X2boot <- X2[bootsamps,]
#   
#   start <- Sys.time()
#   fit <- cmulti_fit_joint(Yboot,
#                           rboot,
#                           tboot,
#                           X1 = X1boot, # Design matrix for tau
#                           X2 = X2boot  # Design matrix for phi
#   )
#   end <- Sys.time()
#   print(end-start)
#   log_offsets <- calculate_offsets(fit,
#                                    rarray = rboot,
#                                    tarray = tboot,
#                                    X1 = X1boot,
#                                    X2 = X2boot)
#   
#   # Predictions of tau across range of FC
#   FCmn <- mean(covariate.FC)
#   FCsd <- sd(covariate.FC)
#   predFC <- seq(min(covariate.FC),max(covariate.FC),length.out = ncol(Dhat_boot))
#   predFCz <- (predFC-FCmn)/FCsd
#   predX1 <- model.matrix(~predFCz)
#   tauhat <- exp(predX1 %*% fit$coefficients[1:ncol(X1)])
#   tau_boot[b,] <- tauhat
#   
#   # Predictions of phi across range of DOY
#   DOYmn <- mean(covariate.DOY)
#   DOYsd <- sd(covariate.DOY)
#   predDOY <- seq(min(covariate.DOY),max(covariate.DOY),length.out = ncol(Dhat_boot))
#   predDOYz <- (predDOY-DOYmn)/DOYsd
#   predX2 <- model.matrix(~predDOYz+I(predDOYz^2))
#   phihat <- exp(predX2 %*% fit$coefficients[(ncol(X1)+1):length(fit$coefficients)])
#   phi_boot[b,] <- phihat
#   
#   # Predictions of density across range of values of MAT
#   MATmn <- mean(covariate.MAT)
#   MATsd <- sd(covariate.MAT)
#   MATboot <- covariate.MAT[bootsamps]
#   zMATboot <- (MATboot - MATmn)/MATsd
#   dat <- data.frame(Y = apply(Yboot,1,sum,na.rm = TRUE),
#                     zMAT = zMATboot,
#                     log_off = log_offsets)
#   
#   glm1 <- glm(Y ~ zMAT + I(zMAT^2) + offset(log_off), family = poisson(link="log"), data = dat)
#   
#   predMAT <- seq(0,25,length.out = ncol(Dhat_boot))
#   predMATz <- (predMAT-MATmn)/MATsd
#   pred_df <- data.frame(zMAT = predMATz,log_off = 0)
#   Dhat <- predict(glm1, newdata = pred_df, type = "response")
#   Dhat_boot[b,] <- Dhat
#   
#   # Summarize results
#   
#   Dhat_mean <- apply(Dhat_boot,2,mean,na.rm = TRUE)
#   Dhat_lcl <- apply(Dhat_boot,2,function(x)quantile(x,0.05,na.rm = TRUE))
#   Dhat_ucl <- apply(Dhat_boot,2,function(x)quantile(x,0.95,na.rm = TRUE))
#   
#   Density_plot <- ggplot()+
#     geom_ribbon(aes(x = predMAT,ymin=Dhat_lcl,ymax=Dhat_ucl), 
#                 fill = "dodgerblue",
#                 col="transparent",
#                 alpha = 0.5)+
#     geom_line(aes(x = predMAT,y = Dhat_mean),col="dodgerblue")+
#     geom_line(aes(x = covariate.MAT, y = Density))+
#     theme_bw()
#   print(Density_plot)
#   print(b)
# }


