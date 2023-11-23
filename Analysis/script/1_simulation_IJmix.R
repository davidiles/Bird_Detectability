# **************************************
# Repeated simulations across different values of tau and phi
# **************************************

library(tidyverse)
library(ggpubr)
library(detect)
library(ggthemes)
library(ggpubr)

rm(list=ls())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Bird_Detectability/Analysis/script")
source("0_model_fns.R")

result_df <- expand.grid(sim_rep = 1:100,
                         
                         # Do birds follow Indep or Joint detection process?
                         det_type = c("Indep","Joint"),
                         
                         # Properties of birds in group A
                         tau_A = 1,
                         phi_A = 3,
                         
                         # Properties of birds in group B
                         tau_B = 1,
                         phi_B = 0.5,
                         
                         # Proportion of birds in each group
                         pA = 0.5,
                         Density = 1000,
                         
                         Density_est_IJ = NA,
                         Density_est_PP = NA,
                         Density_est_IJmix = NA,
                         Density_est_PPmix = NA,
                         
                         best_model = NA
                         
) %>%
  arrange(sim_rep,pA)

for (i in 1:nrow(result_df)){
  
  sim_rep = result_df$sim_rep[i]
  
  Density <- result_df$Density[i]
  tau_A <- result_df$tau_A[i]
  phi_A <- result_df$phi_A[i]
  tau_B <- result_df$tau_B[i]
  phi_B <- result_df$phi_B[i]
  pA <- result_df$pA[i]
  det_type <- result_df$det_type[i]
  
  set.seed(sim_rep)
  
  # ----------------------------------------------------------
  # Simulate data collection at each point count
  # ----------------------------------------------------------
  
  rint <- c(0.5,1,Inf)
  tint <- seq(1,10)
  nrint <- length(rint)
  ntint <- length(tint)
  
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  
  Yarray <- array(NA,dim=c(1,nrint,ntint))
  rarray <- array(NA,dim=c(1,nrint))
  tarray <- array(NA,dim=c(1,ntint))
  
  rarray[1,1:length(rint)] <- rint
  tarray[1,1:length(tint)] <- tint
  
  Y_A <- Y_B <- 0
  
  # *********************************************************************
  # Simulate observations for bird type A
  # *********************************************************************
  
  if (pA >0){
    dim <- 10
    N <- round(Density*pA*dim^2) # Number of birds to place on landscape
    
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
    for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_A))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    # ------------------------------------
    # Determine which cues are detected, based on tau_true
    # ------------------------------------
    
    cues$p <- exp(-(cues$dist/tau_A)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Is first cue detected every time?
    # ------------------------------------
    
    if (det_type == "Indep"){
      dat <- subset(cues,detected == 1 & cue_number == 1)
      dat <- dat[!duplicated(dat$bird_id),]
    } else {
      dat <- subset(cues,detected == 1)
      dat <- dat[!duplicated(dat$bird_id),]
      
    }
    
    # ------------------------------------
    # Transcription: distance and time bins
    # ------------------------------------
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y_A <- table(dat[,c("rint","tint")])
    Y_A # Data to analyze
  }
  
  # *********************************************************************
  # Simulate observations for bird type B 
  # *********************************************************************
  
  if (pA < 1){
    dim <- 10
    N <- round(Density*(1-pA)*dim^2) # Number of birds to place on landscape
    
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
    for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi_B))
    cues <- cues %>% 
      reshape2::melt() %>% 
      rename(bird_id = Var1, cue_number = Var2, time = value) %>%
      arrange(bird_id,cue_number)
    
    cues$dist <- birds$dist[cues$bird_id]
    
    # ------------------------------------
    # Determine which cues are detected, based on tau_true
    # ------------------------------------
    
    cues$p<- exp(-(cues$dist/tau_B)^2)  # Probability each cue is detected
    cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
    
    # ------------------------------------
    # Is first cue detected every time?
    # ------------------------------------
    
    if (det_type == "Indep"){
      dat <- subset(cues,detected == 1 & cue_number == 1)
      dat <- dat[!duplicated(dat$bird_id),]
    } else {
      dat <- subset(cues,detected == 1)
      dat <- dat[!duplicated(dat$bird_id),]
      
    }
    
    # ------------------------------------
    # Transcription: distance and time bins
    # ------------------------------------
    
    # Separate into distance and time bins
    dat$rint <- cut(dat$dist,c(0,rint))
    dat$tint <- cut(dat$time,c(0,tint))
    dat <- na.omit(dat)
    
    Y_B <- table(dat[,c("rint","tint")])
    Y_B # Data to analyze
  }
  
  # *********************************************************************
  # Combine Y_A and Y_B
  # *********************************************************************
  
  Y <- Y_A + Y_B
  Yarray[1,1:nrint,1:ntint] <- Y
  
  Ysum <- sum(Yarray)
  result_df$Ysum[i] <- Ysum
  
  if (Ysum == 0) next # Species was never detected at any point counts
  
  # *********************************************************************
  # Fit IJ model
  # *********************************************************************
  
  start <- Sys.time()
  fit_IJ <- fit_IJ_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0
  )
  end <- Sys.time()
  print(end-start)
  
  tau_est_IJ <- exp(fit_IJ$coefficients[1])
  phi_est_IJ <- exp(fit_IJ$coefficients[2])
  
  offsets_IJ <- offsets_IJ_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              
                              tau = rep(tau_est_IJ,dim(Yarray)[1]),
                              phi = rep(phi_est_IJ,dim(Yarray)[1]))
  
  result_df$offset_IJ[i] <- offsets_IJ$offset_vec[1]
  result_df$Density_est_IJ[i] <- Ysum/offsets_IJ$offset_vec[1]
  
  # *********************************************************************
  # Fit PP model (phi and tau estimated jointly)
  # *********************************************************************
  
  start <- Sys.time()
  fit_PP <- fit_PP_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0
  )
  end <- Sys.time()
  print(end-start)
  
  tau_est_PP <- exp(fit_PP$coefficients[1])
  phi_est_PP <- exp(fit_PP$coefficients[2])
  
  offsets_PP <- offsets_PP_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              
                              tau = rep(tau_est_PP,dim(Yarray)[1]),
                              phi = rep(phi_est_PP,dim(Yarray)[1]))
  
  result_df$offset_PP[i] <- offsets_PP$offset_vec[1]
  result_df$Density_est_PP[i] <- Ysum/offsets_PP$offset_vec[1]
  
  # *********************************************************************
  # Fit IJmix model
  # *********************************************************************
  
  start <- Sys.time()
  fit_IJmix <- fit_IJmix_fn(Yarray = Yarray, 
                            rarray = rarray, 
                            tarray = tarray,
                            tau_A_inits = 0,
                            phi_A_inits = 0,
                            phi_B_inits = 0,
                            pA_inits = 0
  )
  end <- Sys.time()
  print(end-start)
  
  tau_A_est_IJmix <- exp(fit_IJmix$coefficients[1])
  phi_A_est_IJmix <- exp(fit_IJmix$coefficients[2])
  phi_B_est_IJmix <- exp(fit_IJmix$coefficients[3])
  pA_est_IJmix <- plogis(fit_IJmix$coefficients[4])
  
  offsets_IJmix <- offsets_IJmix_fn(Yarray = Yarray,
                                    rarray = rarray, 
                                    tarray = tarray,
                                    
                                    tau_A = rep(tau_A_est_IJmix,dim(Yarray)[1]),
                                    phi_A = rep(phi_A_est_IJmix,dim(Yarray)[1]),
                                    phi_B = rep(phi_B_est_IJmix,dim(Yarray)[1]),
                                    pA = rep(pA_est_IJmix,dim(Yarray)[1]))
  
  result_df$offset_IJmix[i] <- offsets_IJmix$offset_vec[1]
  result_df$Density_est_IJmix[i] <- Ysum/offsets_IJmix$offset_vec[1]
  
  # *********************************************************************
  # Fit PPmix model
  # *********************************************************************
  
  start <- Sys.time()
  fit_PPmix <- fit_PPmix_fn(Yarray = Yarray, 
                            rarray = rarray, 
                            tarray = tarray,
                            tau_A_inits = 0,
                            phi_A_inits = 0,
                            phi_B_inits = 0,
                            pA_inits = 0
  )
  end <- Sys.time()
  print(end-start)
  
  tau_A_est_PPmix <- exp(fit_PPmix$coefficients[1])
  phi_A_est_PPmix <- exp(fit_PPmix$coefficients[2])
  phi_B_est_PPmix <- exp(fit_PPmix$coefficients[3])
  pA_est_PPmix <- plogis(fit_PPmix$coefficients[4])
  
  offsets_PPmix <- offsets_PPmix_fn(Yarray = Yarray,
                                    rarray = rarray, 
                                    tarray = tarray,
                                    
                                    tau_A = rep(tau_A_est_PPmix,dim(Yarray)[1]),
                                    phi_A = rep(phi_A_est_PPmix,dim(Yarray)[1]),
                                    phi_B = rep(phi_B_est_PPmix,dim(Yarray)[1]),
                                    pA = rep(pA_est_PPmix,dim(Yarray)[1]))
  
  result_df$offset_PPmix[i] <- offsets_PPmix$offset_vec[1]
  result_df$Density_est_PPmix[i] <- Ysum/offsets_PPmix$offset_vec[1]
  
  # *********************************************************************
  # Calculate BIC of each model
  # *********************************************************************
  
  BIC_IJ <- 2*log(1000)-2*fit_IJ$loglik
  BIC_PP <- 2*log(1000)-2*fit_PP$loglik
  BIC_IJmix <- 4*log(1000)-2*fit_IJmix$loglik
  BIC_PPmix <- 4*log(1000)-2*fit_PPmix$loglik
  BIC_vals <- c("IJ" = BIC_IJ, "PP" = BIC_PP, "IJmix" = BIC_IJmix, "PPmix" = BIC_PPmix)
  best_model <- BIC_vals[which(BIC_vals == min(BIC_vals))] %>% names()
  result_df$best_model[i] <- best_model
  
  # *********************************************************************
  # Plot results
  # *********************************************************************
  
  result_df$scenario_label <- paste0("Prop 'A' birds = ",result_df$pA)
  result_plot1 <- ggplot(result_df)+
    
    geom_hline(aes(yintercept = Density), linetype = 2)+
    
    geom_boxplot(aes(x = "Ind", y = Density_est_IJ, fill = "Ind", col = "Ind"), alpha = 0.5, size = 1)+
    geom_boxplot(aes(x = "Jnt", y = Density_est_PP, fill = "Jnt", col = "Jnt"), alpha = 0.5, size = 1)+
    geom_boxplot(aes(x = "IndMix", y = Density_est_IJmix, fill = "IndMix", col = "IndMix"), alpha = 0.5, size = 1)+
    geom_boxplot(aes(x = "JntMix", y = Density_est_PPmix, fill = "JntMix", col = "JntMix"), alpha = 0.5, size = 1)+
    
    scale_color_manual(values=c("dodgerblue","gray50","black","orangered"), name = "Model", guide = "none")+
    scale_fill_manual(values=c("dodgerblue","gray50","black","orangered"), name = "Model", guide = "none")+
    
    ylab("Population Density\n\n(Birds/1000km^2)")+
    xlab("Model")+
    ggtitle("Simulation Results")+
    
    facet_grid(det_type~scenario_label)+
    theme_bw()
  
  print(result_plot1)
  
  # # Best model selected by BIC
  # best_model_summary <- result_df %>%
  #   group_by(scenario_label,best_model)%>%
  #   summarize(n = n())%>%
  #   na.omit()
  # 
  # best_model_summary$best_model <- factor(best_model_summary$best_model,
  #                                         levels = c("IJ","IJPP","IJPP5","PP"),
  #                                         labels = c("Indep","In+Jt3","In+Jt5","Joint"))
  # result_plot2 <- ggplot(best_model_summary)+
  #   geom_bar(aes(x = best_model, y = n, fill = best_model, col = best_model), 
  #            stat = "identity", alpha = 0.5, size = 1)+
  #   scale_fill_manual(values=c("dodgerblue","gray50","black","orangered"), drop=FALSE,name = "Model", guide = "none")+
  #   scale_color_manual(values=c("dodgerblue","gray50","black","orangered"), drop=FALSE,name = "Model", guide = "none")+
  #   facet_grid(.~scenario_label)+
  #   theme_bw()+
  #   xlab("Model")+
  #   ylab("# times selected as best")+
  #   ggtitle("Frequency model was selected 'best' with BIC")
  # 
  #result_plot <- ggarrange(result_plot1,result_plot2,nrow=2, align = "hv")
  #print(result_plot)
}


# *********************************************************************
# Plot results
# *********************************************************************
result_df$scenario_label <- paste0("Prop 'indep' birds = ",result_df$pA)
result_plot <- ggplot(result_df)+
  
  geom_hline(aes(yintercept = Density), linetype = 2)+
  
  geom_boxplot(aes(x = "Indep", y = Density_est_IJ, fill = "Indep", col = "Indep"), alpha = 0.5)+
  geom_boxplot(aes(x = "Joint", y = Density_est_PP, fill = "Joint", col = "Joint"), alpha = 0.5)+
  geom_boxplot(aes(x = "Indep+Joint (3 param)", y = Density_est_IJPP, fill = "Indep+Joint3", col = "Indep+Joint3"), alpha = 0.5)+
  geom_boxplot(aes(x = "Indep+Joint (5 param)", y = Density_est_IJPP5, fill = "Indep+Joint5", col = "Indep+Joint5"), alpha = 0.5)+
  
  scale_color_manual(values=c("dodgerblue","gray50","black","orangered"), name = "Model", guide = "none")+
  scale_fill_manual(values=c("dodgerblue","gray50","black","orangered"), name = "Model", guide = "none")+
  
  ylab("Population Density\n\n(Birds/1000km^2)")+
  xlab("Model")+
  ggtitle("Simulation Results")+
  
  facet_grid(.~scenario_label)+
  theme_bw()
print(result_plot)

png(filename = "../results/Sim_Results.png", width = 10, height = 4, units = "in", res = 600)
print(result_plot)
dev.off()
