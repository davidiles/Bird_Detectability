# **************************************
# Repeated simulations across different values of tau and phi
# **************************************

# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
library(viridis)
# -------------------------------------------------------------

rm(list=ls())

result_df <- data.frame(phi_male = runif(1000,0.1,10),
                         phi_adjust_female = sample(c(0.01,0.1,0.5,1),1000,replace = TRUE),
                         tau_male = sample(c(0.5,1,2),1000,replace=TRUE),
                         tau_adjust_female = sample(c(0.25,0.5,1),1000,replace = TRUE),
                         Density_Est = NA,
                         proportion_females = NA) %>%
  mutate(phi_female = phi_male * phi_adjust_female,
         tau_female = tau_male * tau_adjust_female)


Density <- 1000

# Conduct simulations
for (i in 1:nrow(result_df)){
  
  print(i)
  
  # Cue rate of males and females
  phi_male <- result_df$phi_male[i]
  phi_female <- result_df$phi_female[i]
  
  # EDR of males and females
  tau_male <- result_df$tau_male[i]
  tau_female <- result_df$tau_female[i]
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  dat <- data.frame()
  
  time_bins <- 1 # 1 minute time bins
  rint <- c(0.5,1,Inf)
  tint <- seq(time_bins,5,by = time_bins)
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
  # Simulate males
  # -------------------------------------------------------
  
  for (s in 1:n_surveys){
    
    birds <- data.frame(x = runif(N*0.5,-dim/2,dim/2),
                        y = runif(N*0.5,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    
    # Remove birds outside maximum distance
    birds <- subset(birds, dist <= (dim/2))
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    n = nrow(birds)
    birds$dist <- birds$dist
    
    cues <- matrix(NA, nrow=n, ncol = 50)
    for (bird_id in 1:n) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_male))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    cues$p <- exp(-(cues$dist/tau_male)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    cues <- subset(cues,detected == 1 & cue_number == 1) # 
    cues <- cues[!duplicated(cues$bird_id),]
    
    if (nrow(cues)>0){
      cues$phi <- phi_male
      cues$sex <- "male"
      dat <- rbind(dat, cues)
    }
    
    # -------------------------------------------------------
    # Simulate females
    # -------------------------------------------------------
    
    birds <- data.frame(x = runif(N*0.5,-dim/2,dim/2),
                        y = runif(N*0.5,-dim/2,dim/2))
    
    # Distances to observer
    birds$dist <- sqrt(birds$x^2 + birds$y^2)
    
    # Remove birds outside maximum distance
    birds <- subset(birds, dist <= (dim/2))
    birds <- birds %>% arrange(dist)
    birds$bird_id = 1:nrow(birds)
    n = nrow(birds)
    birds$dist <- birds$dist
    
    cues <- matrix(NA, nrow=n, ncol = 50)
    for (bird_id in 1:n) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_female))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    cues$p <- exp(-(cues$dist/(tau_female))^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    cues <- subset(cues,detected == 1 & cue_number == 1) # & cue_number == 1
    cues <- cues[!duplicated(cues$bird_id),]
    
    if (nrow(cues)>0){
      cues$phi <- phi_female
      cues$sex <- "female"
      dat <- rbind(dat, cues)
    }
    
  } # s
  
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
  
  # # ------------------------------------
  # # Analysis: removal model using finite mixture
  # # ------------------------------------
  # 
  # Y_removal <- apply(Yarray,c(1,3),sum)
  # #fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  # #phi_est = exp(fit_p$coefficients)
  # 
  # Mf0 <- detect::cmulti(Y_removal  | tarray ~ 1, type="fmix")
  # 
  # # Proportion of infrequent singers
  # prop_infreq <- plogis(Mf0$coefficients[2])
  # 
  # # Singing rate of infrequent singers
  # phi_infreq <- exp(Mf0$coefficients[1])
  # 
  # # Proportion of birds that sing at least once during the survey
  # p_mix <- 1 - prop_infreq * exp(-max(tint)*phi_infreq)
  # 
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
  result_df$proportion_females[i] <- sum(dat$sex == "female")/nrow(dat)
  
  #   
  # result_plot <- ggplot(result_df) +
  #   geom_hline(yintercept = Density, linetype = 2) + 
  #   geom_hline(yintercept = 0, linetype = 2, col = "transparent") + 
  #   
  #   geom_jitter(aes(x = "Standard removal", y = D_0), width = 0.1, height = 0)+
  #   geom_jitter(aes(x = "Mixture removal", y = D_1), width = 0.1, height = 0)+
  #   ylab("Density Estimate")+
  #   xlab("")+
  #   theme_bw()
  # 
  # print(result_plot)
  
}

result_df$Density_bias <- 100*(result_df$Density_Est-Density)/Density

ggplot(result_df)+
  geom_point(aes(x = proportion_females, y = Density_bias, col = phi_male))+
  ylab("Percent Bias in Density Estimate")+
  xlab("Proportion of females in dataset")+
  scale_color_gradientn(colors = viridis(10), name = "Male cue rate\n\n(cues/min)")+
  theme_bw()


