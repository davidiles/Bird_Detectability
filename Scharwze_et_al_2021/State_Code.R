library(plyr)
library(tidyverse)
library(rstan)

dat <- readRDS("gcspdata.rds")

# Data alterations -- 'dat' has both IJ and PP predictors.  This model is strictly IJ, so 
#                     PP predictors (in the IJ+PP model) must be recoded as IJ for use here.
dat.state <- dat
dat.state$n_bpst <- 2
dat.state$n_bpev <- NULL
dat.state$Xpst   <- cbind(dat$Xpev, dat$Xpst)
dat.state$Xpev   <- NULL

options(mc.cores=parallel::detectCores()) # Enables parallel computation

model = I("Survey")   # Name of model to fit... We now call it 'IJ', but at first it was 'Survey'
# Which parameters to record and save to file
record.list <- c(recordlist(method=model, surv_level=TRUE),                     
                 c("log_lik", "mean_pavl", "mean_pst", "mean_pp"))

# Define stan model object ... may take a couple of minutes to compile
m.state <- stan_model("StanModels/Survey.stan", auto_write=rstan_options("auto_write"=TRUE))

gcsp.state <- sampling(m.state, data=dat.state, chains=4, iter=5000, 
                       thin=4, pars=record.list, init="random")

saveRDS(gcsp.state, "gcsp_state.rds")

