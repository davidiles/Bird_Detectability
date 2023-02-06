##### Reformat data from 'Code from Amundson.R' for use in model fitting with STAN code

library(tidyverse)

## Fixed effect matrices
# Abundance: we use only three of the 8 provided habitats to alleviate correlation caused by the 8-category simplex
Xab  <- sitecovs %>% select(shrub, dec, spruce) %>% as.matrix %>% scale()
Xavl <- as.matrix(jdate)                                                # Note: Julian date is already centered and scaled
Xpev <- as.matrix(model.matrix(y ~ factor(observer))[,-1]) %>% scale()  # Observers split into two groups by Amundson et al
Xpst <- as.matrix(sound) %>% scale()                                    # 0 < sound < 1 : measures percent cover of dense veg

# REQUIRED DATA STRUCTURE FOR STAN MODEL CODE
dat = list(
  n_surv  = nsurveys,
  surveys = 1:nsurveys, # Index of survey number... deprecated
  tau     = 5,          # Maximum detection time
  maxdist = max(brks),  # Maximum observable distance
  bins    = 50,         # Number of bins for step approximation
  ab_prior     = 0,     # Prior expected site abundance: median = 1.0, 95% CI = (0.23, 4.35); Actual _COUNTS_ for GCSP: 0.49
  ab_prior_sd  = 0.75,
  avl_prior    = -1.4,  # Prior Pr(t < 5): median = 0.71, 95% CI = (0.25, 0.995) --- hard-to-detect group
  avl_prior_sd = 0.75,
  n_bab  = ncol(Xab),   # Number of fixed effects (habitat types).
  n_bavl = ncol(Xavl),  # jdate (Julian date)
  n_bpev = ncol(Xpev),  # obsnum (observer group ID)
  n_bpst = ncol(Xpst),  # sound
  n_rab  = 1,  # Number of random effects.  This one is for Transect ID
  n_ravl = 0,
  n_rpev = 0,
  n_rpst = 0,
  n_rabs  = array(ntrans, dim=1), # Number of levels for each random effect
  n_ravls = array(0, dim=0),
  n_rpevs = array(0, dim=0),
  n_rpsts = array(0, dim=0),
  vab_id  = rep(1, times=169),    # Index associating level with effect... all 169 are from the first and only effect (transect)
  vavl_id = array(0, dim=0),
  vpev_id = array(0, dim=0),
  vpst_id = array(0, dim=0),
  Xab  = Xab,                     # Covariate matrices by survey.  No intercept terms.
  Xavl = Xavl,
  Xpev = Xpev,
  Xpst = Xpst,
  Zab  = matrix(tran),            # Matrix of random effect levels by survey (row).  Number of columns = n_rab.
  Zavl = matrix(0, nrow=length(y), ncol=0),
  Zpev = matrix(0, nrow=length(y), ncol=0),
  Zpst = matrix(0, nrow=length(y), ncol=0),
  n_bysurv = gcsites$count,       # Counts by survey
  n_bird   = sum(gcsites$count),  # Total count across all surveys
  time_obs = specobsonly$elapsed_secs / 60, # Observed detection times in minutes
  dist_obs = specobsonly$exact_distance,    # Observed detection distances in meters
  surv_obs = specobsonly$PTPnum,  # SurveyID for each detected bird
  groups = 2                      # Number of heterogeneous groups to use in fitting model
)

saveRDS(dat, "gcspdata.rds")
