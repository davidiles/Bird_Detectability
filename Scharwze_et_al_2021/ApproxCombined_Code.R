library(tidyverse)
library(rstan)

# Read in data
dat <- readRDS("gcspdata.rds")

### Fit the IJ+PP Model using the step function approximation to distance functions

options(mc.cores=parallel::detectCores()) # Enables parallel computation

model = I("Combined_Approx")   # Name of model to fit
# Which parameters to record and save to file
record.list <- c(recordlist(method=model, surv_level=TRUE),                     
                 c("log_lik", "mean_pavl", "mean_pst", "mean_pev", "mean_pp"))

# Define stan model object ... may take a couple of minutes to compile
m.combined <- stan_model("StanModels/Combined_Approx.stan", auto_write=rstan_options("auto_write"=TRUE))

### Fit model
# Note: initial values can be troublesome.  See Stan documentation for default initial values.
#       Generally, Stan defaults are not fully appropriate for model intercept terms, especially for the
#       various intcpt_* terms (see Combined_Approx.stan for model variable names).
init <- readRDS("inits_IJPP.rds")
gcspfit <- sampling(m.combined, data=dat, chains=4, iter=5000, 
                    thin=4, pars=record.list, init=init)

saveRDS(gcspfit, "gcspfit.rds")
