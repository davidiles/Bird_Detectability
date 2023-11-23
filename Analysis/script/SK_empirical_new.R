# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
# -------------------------------------------------------------

rm(list=ls())

theme_set(theme_bw())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Bird_Detectability/Analysis/script")
source("0_model_fns_new.R")

`%!in%` <- Negate(`%in%`)

# -------------------------------------------------------------
# Read data
# -------------------------------------------------------------

dat <- read.csv("../data/Survey_Detection_histories.csv")

# -------------------------------------------------------------
# Format data
# -------------------------------------------------------------

dat <- dat %>% dplyr::rename(Species = SPECIES)

# Survey type
dat$SurveyType[dat$SurveyType != "ARU"] <- "HUM"
dat$SurveyType <- dat$SurveyType %>% 
  as.factor() %>%
  relevel("HUM")

# ----------------------
# Remove some weird data that is probably mislabeled
# ----------------------

# Remove "unlimited distance" human surveys that have no  (are these actually ARU data?) # 1 survey
dat <- subset(dat, !(DISTANCE == "Unlimited" & SurveyType == "HUM"))

# Remove any distance information in ARU surveys (are these actually human surveys?) # 38 surveys
dat <- subset(dat, !(DISTANCE != "Unlimited" & SurveyType == "ARU"))

# Only include 3 and 10 minute counts
dat <- subset(dat, Duration %in% c(3,10))

# Distance binning
dat <- subset(dat, DISTANCE != "")
dat$DistBin <- 1
dat$DistBin[dat$DISTANCE == "50-100 m"] <- 2
dat$DistBin[dat$DISTANCE == ">100 m"] <- 3
dat$nDistBin <- 1
dat$nDistBin[which(dat$DISTANCE %in% c("0-49 m","50-100 m", ">100 m"))] <- 3

# Time Binning
dat <- subset(dat, Time.interval != "")
dat$TimeBin <- factor(dat$Time.interval, levels = c("0-1 min","1-2 min","2-3 min","3-4 min","4-5 min","5-6 min","6-7 min","7-8 min","8-9 min","9-10 min"))
dat$TimeBin <- as.numeric(dat$TimeBin)

dat <- subset(dat, nDistBin == 3 & Duration == 10)

# Remove species that commonly do flyovers, flock, unknowns, or are otherwise inappropriate for point counts
aa <- table(dat$Species)
SPP <- names(aa)
SPP <- SPP[!(SPP %in% c("AMCR","AMGO","BCFR","BLJA","CANG","COLO","COGO","COME","CORA","EVGR","FRGU","GRAJ",
                        "RESQ","WOSP","WWCR","WISN","SACR","PISI","UNBI","LEYE","GRYE","WOFR","UNKN",
                        "UNBL","UNDU","UNFI","UNGU","UNTH","UNVI","UNWA","UNWO","UNYE",
                        "RECR"))]

dat <- dat[dat$Species %in% SPP,]

# -------------------------------------------------------------
# Only include species with at least 50 detections in human point counts
# -------------------------------------------------------------

species_totals <- dat %>%
  group_by(Species,SurveyType) %>%
  summarize(TotalCount = sum(Count))

species_to_include <- species_totals %>%
  subset(SurveyType == "HUM") %>%
  rename(Species = Species) %>%
  arrange(desc(TotalCount)) %>%
  subset(TotalCount > 500)

dat <- subset(dat, Species %in% species_to_include$Species)

# -------------------------------------------------------------
# Survey information
# -------------------------------------------------------------

dat$SurveyID <- paste0(dat$Unique.ID,dat$SurveyDate,dat$Time,dat$SurveyType,dat$Observer) %>% as.factor() %>% as.numeric()
nsurvey <- length(unique(dat$SurveyID))

survey_info <- dat %>%
  group_by(SurveyID) %>%
  summarize(Duration = mean(Duration),
            SurveyType = SurveyType[1])
# Survey summary
table(survey_info[,c("Duration","SurveyType")])

human_surveys <- survey_info$SurveyID[survey_info$SurveyType == "HUM"]
ARU_surveys <- survey_info$SurveyID[survey_info$SurveyType == "ARU"]

# -------------------------------------------------------------
# Create rarray and tarrays
# -------------------------------------------------------------

tarray <- array(NA,dim = c(nsurvey,length(unique(dat$TimeBin))))
rarray <- array(NA,dim = c(nsurvey,3))

for (k in 1:nsurvey){
  
  dat_survey <- subset(dat, SurveyID == k)
  
  max_dur <- survey_info$Duration[survey_info$SurveyID == k]
  
  tarray[k,1:max_dur] <- seq(1,max_dur)
  rarray[k,1] <- Inf
  if (max(dat_survey$nDistBin) == 3) rarray[k,] <- c(0.5,1,Inf)
  
}

# All surveys are the same format, so simplify the analysis
rarray <- array(rarray[1,],dim=c(1,ncol(rarray)))
tarray <- array(tarray[1,],dim=c(1,ncol(tarray)))

# **********************************************************************
# **********************************************************************
# Loop through species and conduct analyses
# **********************************************************************
# **********************************************************************

species_results_offsets <- data.frame()
species_results_p <- data.frame()
species_estimates <- data.frame()

for (sp in species_to_include$Species){
  
  # Generate Yarray for this species
  Yarray <- array(0,dim=c(1,length(rarray),length(tarray)))
  
  # Fill Yarray with counts
  for (k in 1:nsurvey){
    dat_survey <- subset(dat, SurveyID == k)
    
    max_dur <- survey_info$Duration[survey_info$SurveyID == k]
    
    # Select data for this species
    dat_survey_sp <- subset(dat_survey, Species == sp)
    
    # If this survey did not detect the species, skip
    if (nrow(dat_survey_sp)==0) next
    
    for (i in 1:nrow(dat_survey_sp)){
      Yarray[1,dat_survey_sp$DistBin[i],dat_survey_sp$TimeBin[i]] <- Yarray[1,dat_survey_sp$DistBin[i],dat_survey_sp$TimeBin[i]] + dat_survey_sp$Count[i]
    }
    
  }
  
  # ****************************************************************
  # Fit INDEP model
  # ****************************************************************
  
  start_Ind <- Sys.time()
  fit_Ind <- fit_Ind_fn(Yarray = Yarray, 
                        rarray = rarray, 
                        tarray = tarray,
                        tau_inits = 0,
                        phi_inits = 0,
                        method = "SANN"
  )
  
  tau_est_Ind <- exp(fit_Ind$coefficients[1])
  phi_est_Ind <- exp(fit_Ind$coefficients[2])
  
  offsets_Ind <- offsets_Ind_fn(Yarray = Yarray,
                                rarray = rarray, 
                                tarray = tarray,
                                tau = rep(tau_est_Ind,dim(Yarray)[1]),
                                phi = rep(phi_est_Ind,dim(Yarray)[1]))
  end_Ind <- Sys.time()
  
  # ****************************************************************
  # Fit JOINT model
  # ****************************************************************
  
  start_Jnt <- Sys.time()
  fit_Jnt <- fit_Jnt_fn(Yarray = Yarray, 
                        rarray = rarray, 
                        tarray = tarray,
                        tau_inits = 0,
                        phi_inits = 0,
                        method = "SANN"
  )
  
  tmp <- fit_Jnt_fn(Yarray = Yarray, 
                    rarray = rarray, 
                    tarray = tarray,
                    X1 = NULL,
                    X2 = NULL,
                    fit_Jnt$coefficients
  )
  
  tau_est_Jnt <- exp(fit_Jnt$coefficients[1])
  phi_est_Jnt <- exp(fit_Jnt$coefficients[2])
  
  offsets_Jnt <- offsets_Jnt_fn(Yarray = Yarray,
                                rarray = rarray, 
                                tarray = tarray,
                                tau = rep(tau_est_Jnt,dim(Yarray)[1]),
                                phi = rep(phi_est_Jnt,dim(Yarray)[1]))
  
  
  end_Jnt <- Sys.time()
  
  # ****************************************************************
  # Fit Indep MIXTURE model with 5 parameters
  # ****************************************************************
  
  start_IndMix <- Sys.time()
  fit_IndMix <- fit_IndMix_fn(Yarray = Yarray, 
                              rarray = rarray, 
                              tarray = tarray,
                              tau_A_inits = 0,
                              tau_B_inits = 0,
                              phi_A_inits = log(2),
                              phi_B_inits = log(0.5),
                              pA_inits = 0,
                              method = "SANN"
  )
  
  tmp <- fit_IndMix_fn(Yarray = Yarray, 
                    rarray = rarray, 
                    tarray = tarray,
                    X1 = NULL,
                    X2 = NULL,
                    X3 = NULL,
                    X4 = NULL,
                    X5 = NULL,
                    fit_IndMix$coefficients
  )
  tau_A_est_IndMix <- exp(fit_IndMix$coefficients[1])
  tau_B_est_IndMix <- exp(fit_IndMix$coefficients[2])
  phi_A_est_IndMix <- exp(fit_IndMix$coefficients[3])
  phi_B_est_IndMix <- exp(fit_IndMix$coefficients[4])
  pA_est_IndMix <- plogis(fit_IndMix$coefficients[5])
  
  offsets_IndMix <- offsets_IndMix_fn(Yarray = Yarray,
                                      rarray = rarray, 
                                      tarray = tarray,
                                      
                                      tau_A = rep(tau_A_est_IndMix,dim(Yarray)[1]),
                                      tau_B = rep(tau_B_est_IndMix,dim(Yarray)[1]),
                                      phi_A = rep(phi_A_est_IndMix,dim(Yarray)[1]),
                                      phi_B = rep(phi_B_est_IndMix,dim(Yarray)[1]),
                                      pA = rep(pA_est_IndMix,dim(Yarray)[1]))
  
  end_IndMix <- Sys.time() # 14 seconds
  
  # ****************************************************************
  # Fit Joint MIXTURE model with 5 parameters
  # ****************************************************************
  
  start_JntMix <- Sys.time()
  fit_JntMix <- fit_JntMix_fn(Yarray = Yarray, 
                              rarray = rarray, 
                              tarray = tarray,
                              tau_A_inits = 0,
                              tau_B_inits = 0,
                              phi_A_inits = log(2),
                              phi_B_inits = log(0.5),
                              pA_inits = 0,
                              method = "SANN"
  )
  
  tau_A_est_JntMix <- exp(fit_JntMix$coefficients[1])
  tau_B_est_JntMix <- exp(fit_JntMix$coefficients[2])
  phi_A_est_JntMix <- exp(fit_JntMix$coefficients[3])
  phi_B_est_JntMix <- exp(fit_JntMix$coefficients[4])
  pA_est_JntMix <- plogis(fit_JntMix$coefficients[5])
  
  offsets_JntMix <- offsets_JntMix_fn(Yarray = Yarray,
                                      rarray = rarray, 
                                      tarray = tarray,
                                      
                                      tau_A = rep(tau_A_est_JntMix,dim(Yarray)[1]),
                                      tau_B = rep(tau_B_est_JntMix,dim(Yarray)[1]),
                                      phi_A = rep(phi_A_est_JntMix,dim(Yarray)[1]),
                                      phi_B = rep(phi_B_est_JntMix,dim(Yarray)[1]),
                                      pA = rep(pA_est_JntMix,dim(Yarray)[1]))
  
  end_JntMix <- Sys.time() # 14 seconds
  
  
  # ****************************************************************
  # Save offsets and model fit summary for this species
  # ****************************************************************
  
  sp_offsets <- data.frame(Species = sp,
                           nsurvey = length(unique((subset(dat, Species == sp)$SurveyID))),
                           
                           loglik_Ind = fit_Ind$loglik,
                           loglik_Jnt = fit_Jnt$loglik,
                           loglik_IndMix = fit_IndMix$loglik,
                           loglik_JntMix = fit_JntMix$loglik,
                           
                           off_Ind = offsets_Ind$offset_vec[1],
                           off_Jnt = offsets_Jnt$offset_vec[1],
                           off_IndMix = offsets_IndMix$offset_vec[1],
                           off_JntMix = offsets_JntMix$offset_vec[1])
  
  species_results_offsets <- rbind(species_results_offsets,sp_offsets)
  
  # ****************************************************************
  # Save predictions of proportion of birds detected in each distance/time bin
  # ****************************************************************
  
  Y <- Yarray[1,,] # Observed
  p_Ind <- offsets_Ind$parray[1,,]
  p_Jnt <- offsets_Jnt$parray[1,,]
  p_IndMix <- offsets_IndMix$parray[1,,]
  p_JntMix <- offsets_JntMix$parray[1,,]
  
  colnames(Y) <- colnames(p_Ind) <- colnames(p_Jnt) <- colnames(p_IndMix) <- colnames(p_JntMix) <- c("0-1 min","1-2 min","2-3 min","3-4 min","4-5 min","5-6 min","6-7 min","7-8 min","8-9 min","9-10 min")
  rownames(Y) <- rownames(p_Ind) <- rownames(p_Jnt) <- rownames(p_IndMix) <- rownames(p_JntMix) <- c("0-49 m","50-100 m", ">100 m")
  
  p_obs <- reshape2::melt(Y/sum(Y)) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "Observed",
           error = 0)
  
  p_Ind <- reshape2::melt(p_Ind) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "Ind",
           error = p_obs$p - p)
  
  p_Jnt <- reshape2::melt(p_Jnt) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "Jnt",
           error = p_obs$p - p)
  
  p_IndMix <- reshape2::melt(p_IndMix) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "IndMix",
           error = p_obs$p - p)
  
  p_JntMix <- reshape2::melt(p_JntMix) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "JntMix",
           error = p_obs$p - p)
  
  species_results_p <- rbind(species_results_p,
                             p_obs,p_Ind,p_Jnt,p_JntMix,p_IndMix,p_JntMix)
  
  print(sp)
  
}

# ****************************************************************
# Manually calculate BIC
# ****************************************************************

species_results_offsets$BIC_Ind <- 2*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_Ind
species_results_offsets$BIC_Jnt <- 2*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_Jnt
species_results_offsets$BIC_IndMix <- 5*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_IndMix
species_results_offsets$BIC_JntMix <- 5*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_JntMix

# ****************************************************************
# Plot offsets, and log likelihoods of each model
# ****************************************************************

sp_off_long <- species_results_offsets[,c("Species","off_Ind","off_Jnt","off_IndMix","off_JntMix")] %>%
  pivot_longer(
    cols = starts_with("off_"),
    names_to = "Model",
    names_prefix = "off_",
    values_to = "Offset",
    values_drop_na = TRUE
  )

sp_BIC_long <- species_results_offsets[,c("Species","BIC_Ind","BIC_Jnt","BIC_IndMix","BIC_JntMix")] %>%
  pivot_longer(
    cols = starts_with("BIC_"),
    names_to = "Model",
    names_prefix = "BIC_",
    values_to = "BIC",
    values_drop_na = TRUE
  )

sp_results <- full_join(sp_off_long,sp_BIC_long)
sp_results$Best = NA

# For each species, identify which model has lowest BIC
for (sp in unique(sp_results$Species)){
  sp_dat <- subset(sp_results, Species == sp)
  best <- sp_dat$Model[which(sp_dat$BIC == min(sp_dat$BIC))]
  sp_results$Best[which(sp_results$Species == sp & sp_results$Model == best)] <- "Best"
}

sp_results$Model <- factor(sp_results$Model,levels = c("Ind","Jnt","IndMix","JntMix"), labels = c("Ind","Jnt","IndMix","JntMix"))

# Arrange species by relative commonness
sp_results <- sp_results %>%
  full_join(species_to_include) %>%
  arrange(TotalCount)
sp_results$Species <- factor(sp_results$Species, 
                             levels = unique(sp_results$Species))

colpal <- viridis::inferno(6)
ggplot()+
  geom_bar(data = sp_results, aes(x = Species, y = Offset, fill = Model,
                                  col = Best), 
           stat = "identity", 
           position = position_dodge(),
           size = 2)+
  geom_text(data = sp_results, aes(x = Species, y = -0.75,
                                   size = Model,
                                   alpha = Best,
                                   label = round(BIC,2)), 
            stat = "identity", 
            position = position_dodge(width = 1),
            hjust=0,vjust=0.5,
            angle = 90)+
  
  scale_size_manual(values = c(3,3,3,3,3), guide = FALSE)+
  scale_alpha_manual(values = c(1), na.value = 0.5,guide = FALSE)+
  scale_fill_manual(values = colpal[c(2,3,4,5,6)])+
  scale_color_manual(values="black",na.value = "transparent", guide = FALSE)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Model Comparison")

# ****************************************************************
# Relative density predictions
# ****************************************************************

# Compare Best model to Ind model
alternative <- subset(sp_results, Best == "Best") %>%
  full_join(species_to_include) %>%
  arrange(TotalCount)

base <- subset(sp_results, Model == "Ind") %>%
  full_join(species_to_include) %>%
  arrange(TotalCount)

alternative$Species <- factor(alternative$Species, levels = alternative$Species)
base$Species <- factor(base$Species, levels = base$Species)

ggplot()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = alternative$Offset, x = base$Offset, label = alternative$Species))+
  geom_text_repel(aes(y = alternative$Offset, x = base$Offset, label = alternative$Species))

# Percent difference in relative density estimate
# offset represents the relative density of birds, given an observed count of 1
alternative$Density <- 1/alternative$Offset
base$Density <- 1/base$Offset

percent_diff <- 100*(alternative$Density/base$Density - 1)

plot_diff <- ggplot()+
  
  geom_bar(aes(x = alternative$Species, 
               y = percent_diff),stat = "identity",
           fill = "dodgerblue")+
  xlab("Species")+
  ylab("% Difference")+
  ggtitle("")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(breaks = seq(-100,500,10))+
  ggtitle("Change in density estimate when using Best model compared to using Ind model\n\n(10 minute, unlimited distance point count)")
plot_diff 

# ****************************************************************
# Examine predicted versus observed proportions in each distance/time bin
# ****************************************************************

best_model_p <- data.frame()
for (sp in unique(sp_results$Species)){
  best_model <- subset(sp_results, Species == sp & Best == "Best")
  best_model_p <- rbind(best_model_p,subset(species_results_p,Model == best_model$Model))
}

best_model_p$Distance <- factor(best_model_p$Distance, levels = rev(c("0-49 m","50-100 m", ">100 m")))
species_results_p$Distance <- factor(species_results_p$Distance, levels = rev(c("0-49 m","50-100 m", ">100 m")))

lim <- max(abs(species_results_p$error))

p_Ind_plot <- ggplot(subset(species_results_p, Model == "Ind") ,
                     aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("QPAD Model")

p_best_plot <- ggplot(best_model_p,
                      aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Best Model")
