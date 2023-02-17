# **************************************
# Repeated simulations across different values of tau and phi
#
# Assess independence of estimates via their correlation
# **************************************

library(tidyverse)
library(ggpubr)
library(detect)
library(ggthemes)

rm(list=ls())

setwd("~/1_Work/Bird_Detectability/QPAD") # <- set to wherever scripts are stored
source("joint_fns.R")

result_df <- expand.grid(sim_rep = 1:100,
                         tau = seq(0.2,2,0.4),
                         phi = c(0.1,0.5,2.5),
                         Density = c(0.1,0.5,2.5),
                         
                         tau_est_joint = NA,
                         phi_est_joint = NA,
                         log_offset_joint = NA,
                         Density_est_joint = NA,
                         
                         tau_est_indep = NA,
                         phi_est_indep = NA,
                         log_offset_indep = NA,
                         Density_est_indep = NA
                         
                         
                         
)

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
  
  result_df$Ysum[sim_rep] <- sum(Yarray)
  
  # ******************************************
  # FIT JOINT MODEL AND ESTIMATE DENSITY
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
  
  # Calculate survey-level offsets
  log_offset <- calculate.offsets(fit,
                                  rarray = rarray,
                                  tarray = tarray,
                                  X1 = NULL,
                                  X2 = NULL)
  
  # Estimates
  result_df$tau_est_joint[sim_rep] <- exp(fit$coefficients[1])
  result_df$phi_est_joint[sim_rep] <- exp(fit$coefficients[2])
  result_df$log_offset_joint[sim_rep] <- log_offset
  result_df$Density_est_joint[sim_rep] <- (sum(Yarray)/exp(log_offset))/1000
  
  # ******************************************
  # FIT INDEPENDENT DISTANCE AND REMOVAL MODELS
  # ******************************************
  tryCatch({
    Y_distance <- apply(Yarray,2,sum) %>% matrix(.,1)
    fit.q <- cmulti.fit(Y_distance,rarray, type = "dis")
    tau_indep = exp(fit.q$coefficients)
    
    Y_removal <- apply(Yarray,3,sum) %>% matrix(.,1)
    fit.p <- cmulti.fit(Y_removal,tarray, type = "rem")
    phi_indep = exp(fit.p$coefficients)
    
    # Estimate density
    A_hat = pi*tau_indep^2
    p_hat = 1-exp(-max(tarray)*phi_indep)
    D_hat <- sum(Yarray)/(A_hat*p_hat)
    
    result_df$tau_est_indep[sim_rep] <- tau_indep
    result_df$phi_est_indep[sim_rep] <- phi_indep
    result_df$Density_est_indep[sim_rep] <- D_hat/1000
  },
  error = function(e){
    
  }
  )
  print(sim_rep)
  
}

result_df$tau_label = paste0("Tau = ",result_df$tau)
result_df$phi_label = paste0("Phi = ",result_df$phi)
result_df$Density_label = paste0("True Density = ",result_df$Density)
save(result_df,file="results/QPAD_joint_sim3.R")

median_results <- result_df %>%
  group_by(tau_label,phi_label,Density_label) %>%
  summarize_all(median,na.rm=TRUE)

result_plot <- ggplot(median_results)+
  geom_hline(yintercept=0,col="white")+
  geom_line(aes(x = tau, y = Density_est_joint, 
                col = "Joint", 
                linetype = "Joint"), size = 1)+
  geom_line(aes(x = tau, y = Density_est_indep, 
                col = "Indep", 
                linetype = "Indep"), size = 1)+
  geom_hline(aes(yintercept = Density, 
                 col = "Truth", 
                 linetype = "Truth"), size = 1)+
  scale_color_manual(values=c("orangered","dodgerblue","black"),name="Model")+
  scale_linetype_manual(values=c(1,1,3),name="Model")+
  ylab("Density Estimate")+
  xlab("Tau")+
  theme_bw()+
  
  facet_grid(Density_label~phi_label, scales = "free")+
  labs(title='Simulation Results',
       subtitle='Median across 100 repeated simulations')
result_plot

png("results/QPAD_joint_sim3.png", height=6,width=6,units="in",res=600)
print(result_plot)
dev.off()