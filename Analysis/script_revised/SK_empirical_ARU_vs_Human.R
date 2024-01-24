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

HUM_surveys <- survey_info$SurveyID[survey_info$SurveyType == "HUM"]
ARU_surveys <- survey_info$SurveyID[survey_info$SurveyType == "ARU"]

# -------------------------------------------------------------
# Create rarray and tarrays
# -------------------------------------------------------------

tarray <- array(NA,dim = c(nsurvey,length(unique(dat$TimeBin))),
                    dimnames = list(NULL,sort(unique(dat$TimeBin))))
rarray <- array(NA,dim = c(nsurvey,3))

for (k in 1:nsurvey){
  
  dat_survey <- subset(dat, SurveyID == k)
  
  max_dur <- survey_info$Duration[survey_info$SurveyID == k]
  
  tarray[k,1:max_dur] <- seq(1,max_dur)
  rarray[k,1] <- Inf
  if (max(dat_survey$nDistBin) == 3) rarray[k,] <- c(0.5,1,Inf)
  
}

# **********************************************************************
# **********************************************************************
# Loop through species and conduct analyses
# **********************************************************************
# **********************************************************************

results <- data.frame()
for (sp in species_to_include$Species){
  
  # **********************************************************************
  # Fit removal model to ARU data
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
  fit.p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_ARU = exp(fit.p$coefficients)
  
  # **********************************************************************
  # Fit model to HUMAN data
  # **********************************************************************
  
  # Generate Yarray for this species
  Yarray <- array(0,dim=c(1,3,10))
  
  # Fill Yarray with counts
  for (k in 1:length(HUM_surveys)){
    
    id = HUM_surveys[k]
    dat_survey <- subset(dat, SurveyID == id)
    
    max_dur <- survey_info$Duration[survey_info$SurveyID == id]
    
    # Select data for this species
    dat_survey_sp <- subset(dat_survey, Species == sp)
    
    # If this survey did not detect the species, skip
    if (nrow(dat_survey_sp)==0) next
    
    for (i in 1:nrow(dat_survey_sp)){
      Yarray[1,dat_survey_sp$DistBin[i],dat_survey_sp$TimeBin[i]] <- Yarray[1,dat_survey_sp$DistBin[i],dat_survey_sp$TimeBin[i]] + dat_survey_sp$Count[i]
    }
    
  }
  
  Y_removal <- matrix(colSums(Yarray[1,,]),1)
  
  # Fit removal model
  D_removal <- matrix(sort(unique(dat$TimeBin)),1)
  nbins_removal <- length(Y_removal)
  fit.p <- detect::cmulti.fit(Y_removal,D_removal, type = "rem")
  phi_HUM = exp(fit.p$coefficients)
  
  results <- rbind(results,data.frame(Species = sp,
                                      phi_HUM = phi_HUM,
                                      phi_ARU = phi_ARU))
}

lim = range(results[,c(2,3)])
ggplot(results, aes(x = phi_HUM, y = phi_ARU, label = Species))+
  geom_abline(slope=1,intercept=0)+
  geom_point(col = "gray70")+
  geom_text_repel()+
  coord_cartesian(xlim=lim,ylim=lim)+
  theme_bw()

# Key result: phi is higher when estimated from ARU