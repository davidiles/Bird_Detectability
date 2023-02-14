# ******************************************************************
# CUSTOM SIMULATION SCRIPT
# ******************************************************************
# - identical to 'vanilla' version of bSims
# - runs a bit faster because it lacks more complex situations
# - estimates q and p using cmulti()
# - illustrates that estimates are biased, unless each bird can only make a single cue during point counts
# ******************************************************************

library(ggpubr)
library(tidyverse)
library(detect)

rm(list=ls())

N = 50000 # Number of birds to place on landscape (select a high number to allow sufficient sample size)
dim = 1000 # landscape size (metres)
Density = N/dim^2 
Density
results_df <- data.frame()

# Conduct simulations across a range of tau and phi
for (tau in seq(50,200,50)){
  for (phi in c(0.1,1)){
    for (sim_rep in 1:10){
      
      print(sim_rep)
      
      # ------------------------------------
      # Place birds on landscape
      # ------------------------------------
      
      birds <- data.frame(bird_id = 1:N,
                          x = runif(N,-dim/2,dim/2),
                          y = runif(N,-dim/2,dim/2))
      
      
      # Distances to observer
      birds$dist <- sqrt(birds$x^2 + birds$y^2)
      
      # ------------------------------------
      # Simulate bird cues, based on phi
      # ------------------------------------
      
      cues <- matrix(NA, nrow=N, ncol = 10)
      for (bird_id in 1:N) cues[bird_id,] <- cumsum(rexp(ncol(cues),phi))
      cues <- cues %>% 
        reshape2::melt() %>% 
        rename(bird_id = Var1, cue_number = Var2, time = value) %>%
        arrange(bird_id,cue_number)
      
      cues$dist <- birds$dist[cues$bird_id]
      
      # ------------------------------------
      # Determine which cues are detected, based on tau
      # ------------------------------------
      
      # *****************************
      # QPAD IS UNBIASED IF THE FOLLOWING SNIPPET IS USED FIRST
      # cues <- subset(cues, cue_number == 1) # Each bird is only allowed to produce a single cue during the point count
      # *****************************
      
      cues$p<- exp(-(cues$dist/tau)^2)  # Probability each cue is detected
      cues$detected <- rbinom(nrow(cues),1,cues$p) # binary variable: was cue actually detected?
      
      # ------------------------------------
      # Isolate first detected cue for each bird
      # ------------------------------------
      
      dat <- subset(cues,detected == 1)
      dat <- dat[!duplicated(dat$bird_id),]
      
      # ------------------------------------
      # Transcription: distance and time bins
      # ------------------------------------
      
      rint <- c(50,100,Inf)
      tint <- c(3,5,10)
      
      # Separate into distance and time bins
      dat$rint <- cut(dat$dist,c(0,rint))
      dat$tint <- cut(dat$time,c(0,tint))
      dat <- na.omit(dat)
      
      Y <- table(dat[,c("rint","tint")])
      
      # Distance model
      Y_distance = matrix(rowSums(Y),1)
      D_distance = matrix(rint,1)
      nbins_distance <- length(Y_distance)
      
      # Removal model
      Y_removal = matrix(colSums(Y),1)
      D_removal = matrix(tint,1)
      nbins_removal <- length(Y_removal)
      
      # ------------------------------
      # Analysis using conditional maximum likelihood
      # ------------------------------
      
      # Trycatch in case of 'singular' error
      fits <- 0
      fits <- tryCatch({
        fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
        tau_MLE = exp(fit.q$coefficients)
        
        fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
        phi_MLE = exp(fit.p$coefficients)
        
        1
      },
      error = function(cond){0})
      
      if (fits == 0) next
      
      # Estimate density
      A_hat = pi*tau_MLE^2
      p_hat = 1-exp(-max(tint)*phi_MLE)
      D_hat <- sum(Y)/(A_hat*p_hat)
      
      results_df <- rbind(results_df,
                          data.frame(sim_rep = sim_rep,
                                     tau = tau,
                                     phi = phi,
                                     N_detected = sum(Y),
                                     D_hat = D_hat,
                                     
                                     tau_MLE = tau_MLE,
                                     
                                     phi_MLE = phi_MLE))
      
    }
  }
}


# --------------------------------------------------------------
# Plot results for a single combination of tau and phi
# --------------------------------------------------------------
tau_to_plot <- 100
phi_to_plot <- 0.1
results_to_plot <- subset(results_df, tau == tau_to_plot & phi == phi_to_plot)

EDR_bias <- median((results_to_plot$tau_MLE - tau_to_plot)/tau_to_plot * 100) %>% round(1)
plot_tau <- ggplot(results_to_plot,aes(x = sim_rep, y = tau_MLE))+
  geom_point(col="orangered")+
  geom_hline(yintercept = tau_to_plot*0.8,col="transparent")+
  geom_hline(yintercept = tau_to_plot,col="black")+
  geom_hline(yintercept = mean(results_to_plot$tau_MLE), linetype = 2, col = "orangered")+
  ylab("Tau (metres)")+
  xlab("Simulation #")+
  ggtitle(paste0("Tau\nMedian bias = +",EDR_bias,"%"))+
  theme_bw()

cuerate_bias <- median((results_to_plot$phi_MLE - phi_to_plot)/phi_to_plot * 100) %>% round(1)
plot_phi <- ggplot(results_to_plot,aes(x = sim_rep, y = phi_MLE))+
  geom_point(col="orangered")+
  geom_hline(yintercept = phi_to_plot/0.8,col="transparent")+
  geom_hline(yintercept = phi_to_plot,col="black")+
  geom_hline(yintercept = mean(results_to_plot$phi_MLE), linetype = 2, col = "orangered")+
  ylab("Phi (cues per minute)")+
  xlab("Simulation #")+
  ggtitle(paste0("Phi\nMedian bias = ",cuerate_bias,"%"))+
  theme_bw()

Density_bias <- median((results_to_plot$D_hat - Density)/Density * 100) %>% round(1)
plot_Density <- ggplot(results_to_plot,aes(x = sim_rep, y = D_hat))+
  geom_point(col="orangered")+
  geom_hline(yintercept = median(results_to_plot$D_hat)*0.8,col="transparent")+
  geom_hline(yintercept = Density,col="black")+
  geom_hline(yintercept = median(results_to_plot$D_hat), linetype = 2, col = "orangered")+
  ylab("Density")+
  xlab("Simulation #")+
  ggtitle("Density")+
  ggtitle(paste0("Density\nMedian bias = +",Density_bias,"%"))+
  theme_bw()
plot_Density

estimate_plot <- ggarrange(plot_tau,plot_phi,plot_Density, nrow=3)
estimate_plot <- annotate_figure(estimate_plot, top = paste0("True tau = ",tau_to_plot, " , True phi = ", phi_to_plot))
print(estimate_plot)



# # --------------------------------------------------------------
# # Plot median results across combinations of tau and phi
# # --------------------------------------------------------------
# 
# results_df$percent_bias_Density <- (results_df$D_hat - Density)/Density * 100
# results_summarized <- results_df %>%
#   group_by(tau,phi) %>%
#   summarize(percent_bias_Density = median(percent_bias_Density),
#             D_hat = median(D_hat),
#             tau_MLE = median(tau_MLE),
#             phi_MLE = median(phi_MLE))
# 
# plot_Density_bias <- ggplot(data = results_summarized, 
#                             aes(x = tau, y = percent_bias_Density, col = factor(phi)))+
#   geom_line()+
#   theme_bw()+
#   scale_color_manual(name = "Phi", values = c("dodgerblue","orangered"))+
#   ylab("Percent bias in density estimate")+
#   xlab("Tau")+
#   ggtitle("Bias in density estimates across\ncombinations of Tau and Phi")+
#   scale_y_continuous(limits = c(0,max(results_summarized$percent_bias_Density)))
# 
# plot_Density_bias
