library(plyr)
library(tidyverse)
library(rstan)

dat <- readRDS("gcspdata.rds")

# Data alterations -- 'dat' has both IJ and PP predictors.  This model is strictly PP, so 
#                     IJ predictors (in the IJ+PP model) must be recoded as PP for use here.
dat.event <- dat
dat.event$n_bpev <- 2
dat.event$n_bpst <- NULL
dat.event$Xpev   <- cbind(dat$Xpev, dat$Xpst)
dat.event$Xpst   <- NULL

options(mc.cores=parallel::detectCores()) # Enables parallel computation

model = I("Event")   # Name of model to fit.  "Event" is our old name for the PP model.
# Which parameters to record and save to file
record.list <- c(recordlist(method=model, surv_level=TRUE),                     
                 c("log_lik", "mean_pavl", "mean_pev", "mean_pp"))

# Define stan model object ... may take a couple of minutes to compile
m.event <- stan_model("StanModels/Event.stan", auto_write=rstan_options("auto_write"=TRUE))

# Load initial values for MCMC sampling
# Convergence of this model can be susceptible to initial values -- the defaults used
#   in rstan are not the best for fitting this model.
init <- readRDS("inits.rds")

### Fit model
gcsp.event <- sampling(m.event, data=dat.event, chains=4, iter=5000, 
                    thin=4, pars=record.list, init=init)

saveRDS(gcsp.event, "gcsp_event.rds")
