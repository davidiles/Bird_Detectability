library(bSims)
library(ggpubr)
library(tidyverse)
library(detect)
library(jagsUI)
library(reshape2)
library(readxl)

setwd("~/1_Work/Bird_Detectability/QPAD/empirical")

rm(list=ls())


# ----------------------------------------
# Read OVEN data
# ----------------------------------------
#file <- "BLPWdetectionhistories.xlsx"
file <- "OVENdetectionhistory.xlsx"
dat <- read_xlsx(file) %>%
  rename(tbin = 10, dbin = 11)

dat <- dat %>%
  group_by(dbin,tbin) %>%
  summarize(count = n()) %>%
  na.omit()

dat$dbin <- factor(dat$dbin, levels = c("0-49 m","50-100 m",">100 m"))
dat$tbin <- factor(dat$tbin, levels = sort(unique(dat$tbin)))

dat <- dat %>% arrange(dbin,tbin)

# Convert to matrix
tint <- seq(1,10)
rint <- c(50,100,Inf)
Y = dat %>% spread(tbin, count) 
Y <- as.matrix(Y[,2:ncol(Y)])
#colnames(Y) <- tint
rownames(Y) <- levels(dat$dbin)
Y[is.na(Y)] <- 0

# ----------------------------------------
# Fit using standard QPAD
# ----------------------------------------

# Distance model
Y_distance = matrix(rowSums(Y),1)
D_distance = matrix(rint,1)

# Removal model
Y_removal = matrix(colSums(Y),1)
D_removal = matrix(tint,1)

fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
tau_MLE = exp(fit.q$coefficients)

fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
phi_MLE = exp(fit.p$coefficients)

tau_MLE # 92
phi_MLE # 0.39

# ----------------------------------------
# Refit distance model using only first time interval
# ----------------------------------------

# Distance model
Y_distance = matrix(Y[,1],1)
D_distance = matrix(rint,1)

fit.q <- cmulti.fit(Y_distance,D_distance, type = "dis")
tau_MLE_2 = exp(fit.q$coefficients)
tau_MLE_2 # 82

# ----------------------------------------
# Refit time-removal model using only first distance interval
# ----------------------------------------

# Removal model
Y_removal = matrix(Y[1,],1)
D_removal = matrix(tint,1)

fit.p <- cmulti.fit(Y_removal,D_removal, type = "rem")
phi_MLE_2 = exp(fit.p$coefficients)

phi_MLE_2 # 0.39

# ----------------------------------------
# Comparisons
# ----------------------------------------


phi_MLE   # Using all distance bins
phi_MLE_2 # Using only nearest distance bin

tau_MLE   # Using all time bins
tau_MLE_2 # Using only first time bin


A_hat = pi*tau_MLE^2
p_hat = 1-exp(-max(tint)*phi_MLE)
D_hat <- sum(Y)/(A_hat*p_hat)

A_hat_2 = pi*tau_MLE_2^2
p_hat_2 = 1-exp(-max(tint)*phi_MLE_2)
D_hat_2 <- sum(Y)/(A_hat_2*p_hat_2)


