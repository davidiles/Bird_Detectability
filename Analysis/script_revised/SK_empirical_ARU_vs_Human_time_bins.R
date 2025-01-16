# -------------------------------------------------------------
library(tidyverse) # For general data manipulation
library(ggrepel)
library(detect)
library(ggpubr)
# -------------------------------------------------------------

rm(list=ls())

theme_set(theme_bw())

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Landbirds/Bird_Detectability/Analysis/script")

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
# Fix data labels
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

dat <- subset(dat, Duration == 10)

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
  subset(TotalCount > 100)

dat <- subset(dat, Species %in% species_to_include$Species)

# -------------------------------------------------------------
# Create 2-minute and 5-minute bins
# -------------------------------------------------------------
dat$TimeBin_2 <- ceiling(dat$TimeBin/2)
dat$TimeBin_5 <- ceiling(dat$TimeBin/5)

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

HUM_surveys <- survey_info$SurveyID[survey_info$SurveyType == "HUM"]
ARU_surveys <- survey_info$SurveyID[survey_info$SurveyType == "ARU"]

# **********************************************************************
# **********************************************************************
# Loop through species and conduct analyses
# **********************************************************************
# **********************************************************************

results <- data.frame()
for (sp in species_to_include$Species){
  
  # **********************************************************************
  # Fit removal model to ARU data with 1-minute bins
  # **********************************************************************
  
  # Generate Yarray for this species
  Yarray <- array(0,dim=c(1,1,10))
  
  # Fill Yarray with counts
  for (k in 1:length(ARU_surveys)){
    
    id = ARU_surveys[k]
    dat_survey <- subset(dat, SurveyID == id)
    
    max_dur <- survey_info$Duration[survey_info$SurveyID == id]
    
    # Select data for this species
    dat_survey_sp <- subset(dat_survey, Species == sp)
    
    # If this survey did not detect the species, skip
    if (nrow(dat_survey_sp)==0) next
    
    for (i in 1:nrow(dat_survey_sp)){
      Yarray[1,1,dat_survey_sp$TimeBin[i]] <- Yarray[1,1,dat_survey_sp$TimeBin[i]] + dat_survey_sp$Count[i]
    }
    
  }
  
  Y_removal <- matrix(Yarray[1,1,],1)
  
  # Fit removal model
  D_removal <- matrix(sort(unique(dat$TimeBin)),1)
  nbins_removal <- length(Y_removal)
  fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_ARU_1 = exp(fit_p$coefficients)
  
  # **********************************************************************
  # Fit removal model to ARU data with 2-minute bins
  # **********************************************************************
  
  # Generate Yarray for this species
  Yarray <- array(0,dim=c(1,1,10/2))
  
  # Fill Yarray with counts
  for (k in 1:length(ARU_surveys)){
    
    id = ARU_surveys[k]
    dat_survey <- subset(dat, SurveyID == id)
    
    max_dur <- survey_info$Duration[survey_info$SurveyID == id]
    
    # Select data for this species
    dat_survey_sp <- subset(dat_survey, Species == sp)
    
    # If this survey did not detect the species, skip
    if (nrow(dat_survey_sp)==0) next
    
    for (i in 1:nrow(dat_survey_sp)){
      Yarray[1,1,dat_survey_sp$TimeBin_2[i]] <- Yarray[1,1,dat_survey_sp$TimeBin_2[i]] + dat_survey_sp$Count[i]
    }
    
  }
  
  Y_removal <- matrix(Yarray[1,1,],1)
  
  # Fit removal model
  D_removal <- matrix(seq(2,10,2),1)
  nbins_removal <- length(Y_removal)
  fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_ARU_2 = exp(fit_p$coefficients)
  
  # **********************************************************************
  # Fit removal model to ARU data with 5-minute bins
  # **********************************************************************
  
  # Generate Yarray for this species
  Yarray <- array(0,dim=c(1,1,10/5))
  
  # Fill Yarray with counts
  for (k in 1:length(ARU_surveys)){
    
    id = ARU_surveys[k]
    dat_survey <- subset(dat, SurveyID == id)
    
    max_dur <- survey_info$Duration[survey_info$SurveyID == id]
    
    # Select data for this species
    dat_survey_sp <- subset(dat_survey, Species == sp)
    
    # If this survey did not detect the species, skip
    if (nrow(dat_survey_sp)==0) next
    
    for (i in 1:nrow(dat_survey_sp)){
      Yarray[1,1,dat_survey_sp$TimeBin_5[i]] <- Yarray[1,1,dat_survey_sp$TimeBin_5[i]] + dat_survey_sp$Count[i]
    }
    
  }
  
  Y_removal <- matrix(Yarray[1,1,],1)
  
  # Fit removal model
  D_removal <- matrix(seq(5,10,5),1)
  nbins_removal <- length(Y_removal)
  fit_p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_ARU_5 = exp(fit_p$coefficients)
  
  results <- rbind(results,data.frame(Species = sp,
                                      time_bin_duration = c(1,2,5),
                                      phi_est = c(phi_ARU_1,phi_ARU_2,phi_ARU_5)))
  print(sp)
  
  phi_comparison <- ggplot(results, aes(x = time_bin_duration, y = phi_est, col = Species))+
    geom_line()+
    scale_color_manual(values = rep("black",length(unique(results$Species))),guide = "none")+
    ggtitle("Comparison of cue rate (phi) estimated using different time bin lengths")+
    theme_bw()+
    geom_text(data = subset(results, time_bin_duration == 5), aes(x = time_bin_duration+0.1, y = phi_est, label = Species), size = 2)+
    scale_y_continuous(trans = "log10")
  print(phi_comparison)
}

