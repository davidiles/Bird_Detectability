# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
# -------------------------------------------------------------

rm(list=ls())

theme_set(theme_bw())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Bird_Detectability/Analysis/script")
source("0_model_fns.R")

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
  
  start_IJ <- Sys.time()
  fit_IJ <- fit_IJ_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0
  )
  
  tau_est_IJ <- exp(fit_IJ$coefficients[1])
  phi_est_IJ <- exp(fit_IJ$coefficients[2])
  
  offsets_IJ <- offsets_IJ_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              tau = rep(tau_est_IJ,dim(Yarray)[1]),
                              phi = rep(phi_est_IJ,dim(Yarray)[1]))
  end_IJ <- Sys.time()
  
  # ****************************************************************
  # Fit JOINT model
  # ****************************************************************
  
  start_PP <- Sys.time()
  fit_PP <- fit_PP_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0
  )
  
  tau_est_PP <- exp(fit_PP$coefficients[1])
  phi_est_PP <- exp(fit_PP$coefficients[2])
  
  offsets_PP <- offsets_PP_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              tau = rep(tau_est_PP,dim(Yarray)[1]),
                              phi = rep(phi_est_PP,dim(Yarray)[1]))
  end_PP <- Sys.time()
  
  
  # ****************************************************************
  # Fit INDEP-JOINT MIXTURE model with 3 parameters
  # ****************************************************************
  
  start_IJPP <- Sys.time()
  fit_IJPP <- fit_IJPP_fn(Yarray = Yarray, 
                      rarray = rarray, 
                      tarray = tarray,
                      tau_inits = 0,
                      phi_inits = 0,
                      pA_inits = 0
  )
  
  tau_est_IJPP <- exp(fit_IJPP$coefficients[1])
  phi_est_IJPP <- exp(fit_IJPP$coefficients[2])
  pA_est_IJPP <- plogis(fit_IJPP$coefficients[3])
  
  offsets_IJPP <- offsets_IJPP_fn(Yarray = Yarray,
                              rarray = rarray, 
                              tarray = tarray,
                              tau = rep(tau_est_IJPP,dim(Yarray)[1]),
                              phi = rep(phi_est_IJPP,dim(Yarray)[1]),
                              pA = rep(pA_est_IJPP,dim(Yarray)[1]))
  end_IJPP <- Sys.time()
  
  # ****************************************************************
  # Fit INDEP-JOINT MIXTURE model with 5 parameters
  # ****************************************************************
  
  start_IJPP5 <- Sys.time()
  tries <- 0
  loglik_IJPP5 <- -Inf
  while(tries < 50){
  fit_IJPP5 <- fit_IJPP5_fn(Yarray = Yarray, 
                          rarray = rarray, 
                          tarray = tarray,
                          tau_A_inits = rnorm(1,log(tau_est_IJ),0.1),
                          tau_B_inits = rnorm(1,log(tau_est_IJ),0.1),
                          phi_A_inits = rnorm(1,log(phi_est_IJ),0.1),
                          phi_B_inits = rnorm(1,log(phi_est_IJ),0.1),
                          pA_inits = qlogis(runif(1,0.1,0.9))
  )
  if ("try-error" %in% class(fit_IJPP5)){
    tries <- tries+1
    next
  }
  
  if (fit_IJPP5$convergence == 0){
    if (fit_IJPP5$loglik > loglik_IJPP5){
      fit_IJPP5_best <- fit_IJPP5
      loglik_IJPP5 <- fit_IJPP5_best$loglik
      print("Updating model")
    }
  }
  tries <- tries+1
  
  }
  
  tau_A_est_IJPP5 <- exp(fit_IJPP5$coefficients[1])
  tau_B_est_IJPP5 <- exp(fit_IJPP5$coefficients[2])
  phi_A_est_IJPP5 <- exp(fit_IJPP5$coefficients[3])
  phi_B_est_IJPP5 <- exp(fit_IJPP5$coefficients[4])
  pA_est_IJPP5 <- plogis(fit_IJPP5$coefficients[5])
  
  offsets_IJPP5 <- offsets_IJPP5_fn(Yarray = Yarray,
                                  rarray = rarray, 
                                  tarray = tarray,
                                  tau_A = rep(tau_A_est_IJPP5,dim(Yarray)[1]),
                                  tau_B = rep(tau_B_est_IJPP5,dim(Yarray)[1]),
                                  phi_A = rep(phi_A_est_IJPP5,dim(Yarray)[1]),
                                  phi_B = rep(phi_B_est_IJPP5,dim(Yarray)[1]),
                                  pA = rep(pA_est_IJPP5,dim(Yarray)[1]))
  end_IJPP5 <- Sys.time()
  
  # ****************************************************************
  # Save offsets and model fit summary for this species
  # ****************************************************************
  
  sp_offsets <- data.frame(Species = sp,
                           nsurvey = length(unique((subset(dat, Species == sp)$SurveyID))),
                           
                           convergence_IJ = fit_IJ$convergence,
                           convergence_PP = fit_PP$convergence,
                           convergence_IJPP = fit_IJPP$convergence,
                           convergence_IJPP5 = fit_IJPP5$convergence,
                           
                           loglik_IJ = fit_IJ$loglik,
                           loglik_PP = fit_PP$loglik,
                           loglik_IJPP = fit_IJPP$loglik,
                           loglik_IJPP5 = fit_IJPP5$loglik,
                           
                           off_IJ = offsets_IJ$offset_vec[1],
                           off_PP = offsets_PP$offset_vec[1],
                           off_IJPP = offsets_IJPP$offset_vec[1],
                           off_IJPP5 = offsets_IJPP5$offset_vec[1])
  
  species_results_offsets <- rbind(species_results_offsets,sp_offsets)
  
  # ****************************************************************
  # Save predictions of proportion of birds detected in each distance/time bin
  # ****************************************************************
  
  Y <- Yarray[1,,] # Observed
  p_IJ <- offsets_IJ$parray[1,,]
  p_PP <- offsets_PP$parray[1,,]
  p_IJPP <- offsets_IJPP$parray[1,,]
  p_IJPP5 <- offsets_IJPP5$parray[1,,]
  
  colnames(Y) <- colnames(p_IJ) <- colnames(p_PP) <- colnames(p_IJPP) <- colnames(p_IJPP5) <- c("0-1 min","1-2 min","2-3 min","3-4 min","4-5 min","5-6 min","6-7 min","7-8 min","8-9 min","9-10 min")
  rownames(Y) <- rownames(p_IJ) <- rownames(p_PP) <- rownames(p_IJPP) <- rownames(p_IJPP5) <- c("0-49 m","50-100 m", ">100 m")
  
  p_obs <- reshape2::melt(Y/sum(Y)) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "Observed",
           error = 0)
  
  p_IJ <- reshape2::melt(p_IJ) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "IJ",
           error = p_obs$p - p)
  
  p_PP <- reshape2::melt(p_PP) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "PP",
           error = p_obs$p - p)
  
  p_IJPP <- reshape2::melt(p_IJPP) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "IJPP",
           error = p_obs$p - p)
  
  p_IJPP5 <- reshape2::melt(p_IJPP5) %>%
    rename(Distance = Var1, Time = Var2, p = value) %>%
    mutate(Species = sp,
           Model = "IJPP5",
           error = p_obs$p - p)
  
  species_results_p <- rbind(species_results_p,
                             p_obs,p_IJ,p_PP,p_IJPP,p_IJPP5)
  
  print(sp)
  
}

# # ****************************************************************
# # Take a look at actual parameter estimates from IJPP (are they believable?)
# # ****************************************************************
# 
# species_coefficients_IJPP <- species_coefficients_IJPP %>%
#   full_join(species_to_include) %>%
#   arrange(TotalCount)
# species_coefficients_IJPP$Species <- factor(species_coefficients_IJPP$Species, 
#                                              levels = unique(species_coefficients_IJPP$Species))
# 
# ggplot(subset(species_coefficients_IJPP, Parameter %in% c("tau_A","tau_B","phi")), 
#        aes(x = Species, y = Estimate, fill = Parameter))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   scale_fill_manual(values = viridis::viridis(5)[c(2,3,4)])+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# ggplot(subset(species_coefficients_IJPP, Parameter %in% c("pA")), 
#        aes(x = Species, y = Estimate, fill = Parameter))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   scale_fill_manual(values = "dodgerblue")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   coord_cartesian(ylim=c(0,1))
# 

ggplot(species_results_offsets,aes(x = off_IJ,y = off_IJPP5, label = Species))+
  geom_point()+
  geom_text_repel()+
  geom_abline(slope=1,intercept=0)

# ****************************************************************
# Manually calculate BIC
# ****************************************************************

species_results_offsets$BIC_IJ <- 2*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_IJ
species_results_offsets$BIC_PP <- 2*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_PP
species_results_offsets$BIC_IJPP <- 3*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_IJPP
species_results_offsets$BIC_IJPP5 <- 5*log(species_results_offsets$nsurvey)-2*species_results_offsets$loglik_IJPP5

# ****************************************************************
# Plot offsets, and log likelihoods of each model
# ****************************************************************

sp_off_long <- species_results_offsets[,c("Species","off_IJ","off_PP","off_IJPP","off_IJPP5")] %>%
  pivot_longer(
    cols = starts_with("off_"),
    names_to = "Model",
    names_prefix = "off_",
    values_to = "Offset",
    values_drop_na = TRUE
  )

sp_BIC_long <- species_results_offsets[,c("Species","BIC_IJ","BIC_PP","BIC_IJPP","BIC_IJPP5")] %>%
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
sp_results$Model[sp_results$Model == "IJPP"] <- "IJPP3"
sp_results$Model <- factor(sp_results$Model,levels = c("IJ","PP","IJPP3","IJPP5"), labels = c("IJ","PP","IJPP3","IJPP5"))

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
  ggtitle("Comparison of offsets from:\n\n(a) IJ\n(b) PP\n(c) IJPP with 3 parameters\n(d) IJPP with 5 parameters")

# ****************************************************************
# Relative density predictions
# ****************************************************************

# Compare IJPP5 model to IJ model
alternative <- subset(sp_results, Model == "IJPP5") %>%
  full_join(species_to_include) %>%
  arrange(TotalCount)

base <- subset(sp_results, Model == "IJ") %>%
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
  ggtitle("Change in density estimate when using IJ+PP model compared to using IJ model\n\n(10 minute, unlimited distance point count)")
plot_diff 


# ****************************************************************
# Examine predicted versus observed proportions in each distance/time bin
# ****************************************************************

species_results_p$Distance <- factor(species_results_p$Distance, levels = rev(c("0-49 m","50-100 m", ">100 m")))

lim <- max(abs(species_results_p$error))

p_IJ_plot <- ggplot(subset(species_results_p, Model == "IJ") ,
                    aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Model IJ")

p_PP_plot <- ggplot(subset(species_results_p, Model == "PP") ,
                    aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Model PP")

p_IJPP_plot <- ggplot(subset(species_results_p, Model == "IJPP") ,
                       aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Model IJPP")

p_IJPP5_plot <- ggplot(subset(species_results_p, Model == "IJPP5") ,
                      aes(x = Time, y = Distance, fill = error))+
  geom_raster()+
  coord_equal()+
  xlab("Time Interval")+
  ylab("Distance Bin")+
  facet_wrap(Species~.)+
  scale_fill_gradientn(colors = c("red","white","blue"), limits = c(-lim,lim),
                       name = "Error\n(observed - expected)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ggtitle("Model IJPP5")

p_IJ_plot
p_PP_plot
p_IJPP_plot
p_IJPP5_plot
