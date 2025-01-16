
# 1) All birds have the same exponential distribution that describes the time between cues

# 2) All birds have a constant interval between cues, but each bird's constant interval is drawn from an exponential


#  How do these different assumptions affect goodness of fit and estimates of cue rate?

# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
library(viridis)
# -------------------------------------------------------------

rm(list=ls())

result_df <- expand.grid(scenario = c("A","B"),
                         phi = c(0.2,1,5),
                         tau = 1,
                         reps = 1:1000,
                        Density_Est = NA)


Density <- 100

# Conduct simulations
for (i in 1:nrow(result_df)){
  
  print(i)
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  dat <- data.frame()
  
  time_bins <- 1 # 1 minute time bins
  rint <- c(0.5,1,Inf)
  tint <- seq(time_bins,10,by = time_bins)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  
  n_surveys <- 1 # Number of point counts
  
  Yarray <- array(NA,dim=c(n_surveys,nrint,ntint))
  rarray <- array(NA,dim=c(n_surveys,nrint))
  tarray <- array(NA,dim=c(n_surveys,ntint))
  
  rarray <- matrix(rep(rint,n_surveys),byrow = TRUE, nrow=n_surveys)
  tarray <- matrix(rep(tint,n_surveys),byrow = TRUE, nrow=n_surveys)
  
  Y <- 0
  
  dim <- 10
  
  N <- round(Density*dim^2) # Number of birds to simulate at each point count
  
  dat <- data.frame()
  
  # -------------------------------------------------------
  # Simulate cues
  # -------------------------------------------------------
  
  phi = result_df$phi[i]
  tau = result_df$tau[i]
  scenario = result_df$scenario[i]
    
  for (s in 1:n_surveys){
    
    birds <- data.frame(x = runif(N,-dim/2,dim/2),
                        y = runif(N,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    
    # Remove birds outside maximum distance
    birds <- subset(birds, dist <= (dim/2))
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    n = nrow(birds)
    birds$dist <- birds$dist
    
    cues <- matrix(NA, nrow=n, ncol = 1000)
    
    start_survey <- 0 #runif(1,0,5)
    
    # scenario = "B"
    # Choose one of two ways to simulate bird cues
    if (scenario == "A"){
      
      # All birds have the same (variable) cue rate as each other
      
      # Exponential
      for (bird_id in 1:n) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi)) - start_survey # The interval between cues described by phi
    
      
      }
    
    if (scenario == "B"){
      
      # Each bird has a constant cue rate over time, but different birds have different rates
      intervals_per_bird <- rexp(n,phi)
      
      # Using exponential
      for (bird_id in 1:n) cues[bird_id,] <- cumsum(rep(intervals_per_bird[bird_id],ncol(cues))) - start_survey # The interval between cues described by phi
      
      start_per_bird <- runif(n,-10,0)
      intervals_per_bird <- rnorm(n,0.3,0.05)
      for (bird_id in 1:n) cues[bird_id,] <- start_per_bird[bird_id] + cumsum(rep(intervals_per_bird[bird_id],ncol(cues))) - start_survey # The interval between cues described by phi
      
    }
    
    #cues[,1] %>% hist(breaks = seq(0,50,0.1)) #**** This distribution is indistinguishable between the two scenarios
    
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number) %>%
      subset(time > 0)
    
    
    first_cues <- cues %>% group_by(bird_id) %>% summarize(first_cue = min(cue_number))
    cues <- left_join(cues,first_cues)
    #hist(subset(cues,cue_number == first_cue)$time)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    cues$p <- exp(-(cues$dist/tau)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # Assume observers only record the first cue actually produced by the bird
    cues <- subset(cues,detected == 1 & time <= max(tint)) # & cue_number == first_cue)
    
    # hist(cues$time) #**** This distribution is indistinguishable between the two scenarios
    
    # Assume observers only record a single cue for each bird
    cues <- cues[!duplicated(cues$bird_id),]
    
    if (nrow(cues)>0){
      cues$phi <- phi
      dat <- rbind(dat, cues)
    }
    
    # -------------------------------------------------------
    # Transcription: distance and time bins
    # -------------------------------------------------------
    
    Yarray[s,,] <- 0
    
    if (nrow(dat)>0){
      
      # Separate into distance and time bins
      dat$rint <- cut(dat$dist,c(0,rint))
      dat$tint <- cut(dat$time,c(0,tint))
      dat <- na.omit(dat)
      
      Y <- table(dat[,c("rint","tint")])
      
      Y # Data to analyze
      Yarray[s,,] <- Y
    }
    
  } # s
  
  # ------------------------------------
  # Analysis: removal model using standard removal modeling
  # ------------------------------------
  
  Y_removal <- apply(Yarray,c(1,3),sum)
  M0 <- detect::cmulti.fit(Y_removal, tarray , type="rem")
  phi_est <- exp(M0$coefficients[1])
  p_est <- 1-exp(-max(tint)*phi_est)
  
  # ------------------------------------
  # Analysis: distance model
  # ------------------------------------
  
  Y_distance <- apply(Yarray,c(1,2),sum)
  fit_q <- detect::cmulti.fit(Y_distance , rarray, type = "dis")
  tau_est = exp(fit_q$coefficients)
  
  # **********************************************************************
  # Density estimates
  # **********************************************************************
  
  # Estimate density for each species
  A <- pi*tau_est^2 # area sampled by unlimited distance survey
  D_est <- sum(Yarray)/n_surveys/(A*p_est) %>% as.numeric()
  
  result_df$Density_Est[i] <- D_est
  
  plot1 <- ggplot(result_df)+geom_histogram(aes(x = Density_Est, fill = scenario)) + theme_bw()+ facet_grid(scenario~phi)
  print(plot1)
}

result_df$Density_bias <- 100*(result_df$Density_Est-Density)/Density

ggplot(result_df)+
  geom_histogram(aes(x = log(Density_Est), fill = scenario)) + 
  theme_bw()+ 
  facet_wrap(scenario~phi, scales = "free_x")+
  geom_vline(xintercept = log(Density))
