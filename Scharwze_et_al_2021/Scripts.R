########## Helper Scripts
# Many of these functions are not be used in the GCSP analysis but are holdovers
#   from other models and analyses.



# Generate list of parameters to retrieve, manipulate, etc.
recordlist <- function(method="Combined", surv_level=F){
  # Base list: all hyperparameters and global_p
  stopifnot(method %in% c("Combined", "Event", "Survey", "Combined_Approx"))
  record.list <- c("intcpt_ab", "bab", "sigma_ab", 
                   "intcpt_avl", "intcpt_g", "int_vec", "bavl", "sigma_avl", 
                   "gamma", "mixing", "global_p")
  if(method != "Survey") record.list <- c(record.list, "intcpt_pev", "bpev", "sigma_pev")
  if(method != "Event")  record.list <- c(record.list, "intcpt_pst", "bpst", "sigma_pst")
  if(surv_level) {
    record.list <- c(record.list, "rab", "ravl", "unobserved", "totN")
    if(method != "Survey") record.list <- c(record.list, "rpev")
    if(method != "Event")  record.list <- c(record.list, "rpst")
  }
  return(record.list)
}

### Calculate p(det), p_avl, p_pcpt given fixed effects and uniform abundance for each model type
# Standard half-normal perceptibility function
g_r <- function(distance, sigma) exp(-(distance/sigma)^2)   # Distance function

# Pr(avl|phi)
p_avail <- function(phi, duration=10, log=T){
  if(log) phi <- exp(phi)
  pexp(duration, phi)
}

# p(det|phi,sig) for event model
# Parameters on the log-scale [log(phi) and log(sigma)], 'sig' is scaled to the maximum distance w
p_event <- function(phi, sig, duration=10, log=T){
  library(gsl) # Loads the exponential integral function
  if(log){
    phi <- exp(phi)
    sig <- exp(sig)
  }
  return(1 + sig^2*(expint_E1(duration*phi) - expint_E1(duration*phi*g_r(1,sig)))) 
}

# p(det|phi,sig) for survey model
# Parameters on the log-scale [log(phi) and log(sigma)], 'sig' is scaled to the maximum distance w
p_survey <- function(phi, sig, duration=10, log=T){
  if(log){
    phi <- exp(phi)
    sig <- exp(sig)
  }
  return(sig^2 * (1-g_r(1,sig)) * pexp(duration,phi))
}

# p(det|phi,sig_ev,sig_st) for event model
# Parameters on the log-scale [log(phi) and log(sigma)], 'sig' is scaled to the maximum distance w
p_combined <- function(phi, sigev, sigst, duration=10, log=T){
  if(log){
    phi <- exp(phi)
    sigev <- exp(sigev)
    sigst <- exp(sigst)
  }
  pdet_r <- function(r, phi, sigst, sigev, duration) 2*r*g_r(r,sigst)*pexp(duration, phi*g_r(r,sigev))
  pdet <- integrate(f=pdet_r, lower=0, upper=1, phi=phi, sigst=sigst, sigev=sigev, duration=duration)$value
  return(pdet) 
}

# p(det|avail,sig) a.k.a p_pcpt for the event model
p_pcpt_e <- function(phi, sig, duration=10, log=T){
  return(p_event(phi, sig, duration=duration, log=log) / p_avail(phi, duration=duration, log=log))
}

# p(det|avail,sig) a.k.a p_pcpt for the state model
p_pcpt_s <- function(sig, log=T){
  if(log) sig <- exp(sig)
  return(sig^2 * (1-g_r(1,sig)))
}

# p(det|avail,sig) a.k.a p_pcpt for the combined model
p_pcpt_c <- function(phi, sigev, sigst, duration=10, log=T){
  return(p_combined(phi, sigev, sigst, duration=duration, log=log) / p_avail(phi, duration=duration, log=log))
}



##### Calculate intercept parameters required to attain desired p_avl, p_pst, or p_pev
# Calculate desired intcpt_avl given target Pr(avail)
intcpt_a_fxn <- function(pavail) log(-0.1*log(1-pavail))
# Calculate desired intcpt_pst given target Pr(det|availability)
survey_intcpt_p <- function(pp_target){
  sigma <- optimize(f = function(sig) abs(sig^2 * (1-g_r(1,sig)) - pp_target), interval = c(0,5))$minimum
  return(log(sigma))
}
# Calculate desired intcpt_pev given the target pdet and Pr(avail)
event_intcpt_p <- function(pdet, pavail){
  phi <- intcpt_a_fxn(pavail)
  intcpt_p <- optimize(f = function(sig) abs(p_event(phi, sig, log=T) - pdet), interval=c(-5,5))$minimum
  return(intcpt_p)
}



########## Analyses Scripts

### This is a function to fit a given dataset with a given .stan model via MCMC sampling
### It will fit a model repeatedly until a target n_eff is achieved or a defined number of attempts have been made
### If a good n_eff is not achieved, the last model is returned (even though it's bad)
### A similar filter could be imposed for Rhat values, but my experience is that the n_eff filter is far more restrictive
# Function Terms:
# dat    = dataset...   see SimData_CT.R for format of inputs
# model  = name of model.stan ("Combined", "Event", or "Survey")
# chains = number of chains to sample
# iter   = number of iterations pers chain (including warmup)
# thin   = degree of thinning of samples
# n_eff  = minumum acceptable number of effective samples for a model fit
# ntries = number of times to attempt fitting a dataset with acceptable n_eff before giving up
# record.plus = additional parameters to record besides the built-in record.list.  Other possibilities include:
# c("log_lambda", "log_phi_vec", "log_phi_avail", "log_psig", "pdet", "p_surv")
# Thin=5 was the product of trial-and-error and reviewing AC plots for a no-effects model

#fitModel <- function(dat, model, chains=4, iter=75, thin=5, record.plus="log_lik"){
fitModel <- function(dat, model, chains=4, iter=6250, thin=5, NEFF=1000, ntries=4, 
                     record.plus="log_lik", initlist="random"){
  library(rstan)
  options(mc.cores=parallel::detectCores()) # Enables parallel computation
  # If the input dataset is from a simulation, use only the data, not the true values, 
  if(identical(names(dat), c("inputs", "truevalues"))) dat <- dat$inputs
  # Define stan model object
  m <- stan_model(paste0("StanModels/", model, ".stan"), auto_write=rstan_options("auto_write"=TRUE))
  # Identify which parameters to record
  record.list <- c(recordlist(method=model, surv_level=TRUE), record.plus)
  tf <- FALSE; fitnumb <- 0
  # while(!tf){    # Allows us to refit models if the number of effective samples is inadequate
  # Quadratically increase iterations after each fail
  fit <- sampling(m, data=dat, chains=chains, iter=ceiling(iter*1.5^fitnumb), 
                  thin=thin, pars=record.list, init=initlist)  # Fit Model
  # # Calculate minimum ESS for parameters of interest
  # parmsupp <- parm_supps(fit)
  # df <- getNeff(fit, parmsupp)
  # n_eff <- subset(df, parameter %in% c("intcpt_ab", "intcpt_avl", "global_p", "p_avl", "p_pcpt"))$ESS
  # tf <- (min(n_eff)>NEFF)            # Stop if the minimum sample size is above 1000
  # fitnumb <- fitnumb+1
  # if(fitnumb > ntries-1) tf <- TRUE  # Stop if we've fit the model too many times
  # }
  return(fit)
}

### Function to create a complete initial values list from a list of some parameter initial values
makeInitial <- function(target=analyses$target[indx], filename="initlist_c.Rdata"){
  initnm <- gsub("_[^_]+$", "", target)  # Drop '_rep##' from the target name
  load(filename)                         # Load initial values list-of-lists.  Sub-list names encode: Model_Data_nobs_Avl_Pcpt
  initvalues <- eval(parse(text = paste0("initlist_c$", initnm)))  # Extract the analysis-appropriate list of initial values
  # Assign zero-dimension arrays to all of the following initial values... it's the price I pay for having a very flexible .stan script
  null.vals <- c("intcpt_g", "bab", "bavl", "bpev", "bpst", "sigma_ab", "sigma_avl", "sigma_pev", "sigma_pst", "rab", "ravl", "rpev", "rpst", "gamma")
  for(k in 1:length(null.vals)){
    for(j in 1:length(initvalues)) initvalues[[j]][[null.vals[k]]] <- array(0, dim=0)
  }
  return(initvalues)
}

### Function to generate effective sample sizes for desired parameters
getNeff <- function(fit, parmsupp, exclude=c("p_dataset", "int_vec[1]", "mixing[1]")){
  # Examining n_eff
  library(rstan)
  library(mcmcse)
  nlist     <- recordlist(fit@model_name)
  postSamp  <- as.data.frame(fit, pars=c(nlist, "lp__"))
  names(postSamp)[names(postSamp)=="global_p"] <- "p_dataset"   # BUG fix
  postSamp  <- cbind(postSamp, parmsupp)                        # Posterior samples and supplementary quantities
  postSamp  <- postSamp[,-which(names(postSamp) %in% exclude)]  # Exclude desired parameters
  nms       <- names(postSamp)
  out.df    <- data.frame(parameter=nms, ESS=NA)
  for(j in 1:length(nms)) out.df$ESS[j] <- ess(postSamp[,j])
  return(out.df)
}



########## Summaries scripts

# A function to gather posterior samples for all parameters into a unified long-format data table
# Survey-level parameters are not by default included (e.g. "unobserved", "totN", "rab", "ravl", "rpev", "rpst", etc.)
posterior_samples <- function(fit, record.plus=NULL){
  library(reshape2)
  library(rstan)
  record.list <- c(recordlist(method=fit@model_name, surv_level=FALSE), record.plus)
  tmp <- as.data.frame(fit, pars=record.list)
  return(melt(tmp, id.vars=NULL))
}

# This is a function to generate posterior samples for global_p, p_dataset, p_globalp_avl, and p_pcpt
# It also corrects the naming inconsistency (from the .stan models vis-a-vis SimData) in global_p and p_dataset
parm_supps <- function(fit){  
  # Extract samples for parameters needed to calculated p_avl, p_pcpt, pdet
  parvec <- c("intcpt", "global_p")
  samples <- as.data.frame(fit, pars=names(fit)[grepl(paste0("^",paste(parvec, collapse="|^")), names(fit))])
  
  # Calculate p_avail (same for all models)
  samples$p_avl <- with(samples, p_avail(intcpt_avl))
  
  ### Calculating pdet two ways and naming them consistently with SimData.R
  # SimData stores the parameter-based estimate as 'global_p', and stores the dataset-based estimate as 'p_dataset'
  # Unfortunately, the model code stores the dataset-based estimate as 'global_p', when it should be 'p_dataset' (BUG)
  # See definitions for both terms in the equations below.
  # 'global_p' (the SimData way) is more correct, but the difference between methods is minimal.
  # Also, 'global_p' is tricky to implement when there are fixed effects and nigh impossible with random effects, whereas
  # 'p_dataset' is easy to do for any model
  
  # Rename the model output from 'global_p' to 'p_dataset'
  samples$p_dataset <- samples$global_p   # pdet estimated from model code: n / N(rep) --- i.e., based on the posterior estimate of N
  # Calculate the appropriate 'global_p' = fxn(intcpt_avl, intcpt_pev, intcpt_pst)
  if(fit@model_name=="Event") {
    samples$global_p <- with(samples, p_event(intcpt_avl, intcpt_pev))
  } else if(fit@model_name=="Survey") {
    samples$global_p <- with(samples, p_survey(intcpt_avl, intcpt_pst))
  } else if(grepl("Combined", fit@model_name)) {
    for(j in 1:nrow(samples)) samples$global_p[j] <- with(samples, p_combined(intcpt_avl[j], intcpt_pev[j], intcpt_pst[j]))
  }
  
  # Calculate p_pcpt
  samples$p_pcpt <- with(samples, global_p / p_avl)
  
  return(samples[,c("global_p", "p_dataset", "p_avl", "p_pcpt")])
}

# A function to calculate posterior p-values each parameter: Pr(theta < true_theta)
# Survey-level parameters are not by default included (e.g. "unobserved", "totN", "rab", "ravl", "rp", etc.)
# The global_p from this function is REALLY p_dataset, because it is based on model output (see explanation
# in the parm_supps function)
posterior_p <- function(fit, dataset, record.plus=NULL){
  library(rstan)
  record.list <- c(recordlist(method=fit@model_name, surv_level=FALSE), record.plus)
  tmp <- as.data.frame(fit, pars=record.list)
  truth <- data.frame(dataset=deparse(substitute(dataset)),
                      variable=names(tmp), trueval=NA, postp=NA)
  for(i in 1:length(names(tmp))){ 
    # Leverage the square-bracket nomenclature from stan summary output to index true values from the original data simulation
    truth$trueval[i] <- eval(parse(text=paste("dataset$truevalues$",names(tmp)[i],sep="")))
    # Calculate posterior p-value
    truth$postp[i] <- sum(tmp[,i] < truth$trueval[i]) / nrow(tmp)
  }
  return(truth)
}

# Function to extract true simulation values and calculate posterior p-values for 
# global_p (based on parameters), p_dataset, p_avl, and p_pcpt
postpsupp <- function(parmsupp, dataset){
  library(rstan)
  true.intcpt.avl <- dataset$truevalues$intcpt_avl
  true.intcpt.pev <- dataset$truevalues$intcpt_pev
  true.intcpt.pst <- dataset$truevalues$intcpt_pst
  truth <- data.frame(dataset=deparse(substitute(dataset)),
                      variable=names(parmsupp), trueval=NA, postp=NA)
  truth$trueval[truth$variable=="global_p"]  <- dataset$truevalues$global_p
  truth$trueval[truth$variable=="p_dataset"] <- dataset$truevalues$p_dataset
  truth$trueval[truth$variable=="p_avl"]     <- p_avail(dataset$truevalues$intcpt_avl)
  truth$trueval[truth$variable=="p_pcpt"]    <- dataset$truevalues$global_p / p_avail(dataset$truevalues$intcpt_avl)
  for(i in 1:nrow(truth)) truth$postp[i] <- sum(parmsupp[,i] < truth$trueval[i]) / nrow(parmsupp)
  return(truth)
}

# Wrapper functions to gather/calculate posterior fit statistics from log-likelhood values.
loo_calc <- function(fit){
  library(loo)
  library(rstan)
  return(loo(extract_log_lik(fit)))
}

waic_calc <- function(fit){
  library(loo)
  library(rstan)
  return(waic(extract_log_lik(fit)))
}

# Calculate DIC using the minimum sampled deviance as an estimate of the deviance at the MAP
dic_calc <- function(fit){
  library(loo)
  library(rstan)
  ll <- rowSums(extract_log_lik(fit))  # overall log-likelihood for each MCMC iteration
  dev <- -2*ll
  return(2*mean(dev)-min(dev))
}

# A function to generate an rstan summary table for posterior samples from all parameters
# Survey-level parameters are not by default included (e.g. "unobserved", "totN", "rab", "ravl", "rp", etc.)
sum_tables <- function(fit, record.plus=NULL){
  library(rstan)
  record.list <- c(recordlist(method=fit@model_name, surv_level=FALSE), record.plus)
  return(data.frame(summary(fit, pars=record.list)$summary))
}

# A function to generate additional summary tables for global_p, p_dataset, p_avl, p_pcpt 
sum_tables_supp <- function(parmsupp){
  library(coda)
  return(cbind(summary(mcmc(parmsupp))[[1]], summary(mcmc(parmsupp))[[2]], NA, NA)[,c(1,4,2,5:11)])
}





########## Output scripts

# A function to write csv's to a directory other than the working directory
# csvwrite <- function(object, targetname, dir) write.csv(object, file=paste0(dir,"/",targetname))
# csvwrite <- function(object, targetname) write.csv(object, file=paste0("Outputs/", targetname))  # Doesn't work
csvwrite <- function(object, targetname) write.csv(object, file=paste0("./", targetname))  # Works (but not to desired destination)
rdswrite <- function(object, targetname) saveRDS(object, file=paste0("./", targetname))

########## Plotting Scripts

# Helper function to trim plot names... lops off the first, second, and (if rep=F) last terms
# (where terms in the input object are separated by underscores)
# 'x' is a target from a remakeGenerator table such as 'datasets' or 'analyses' or 'summaries'
trimfxn <- function(x, rep=FALSE){
  library(plyr)
  x.sep  <- strsplit(x,"_")
  x.glue <- laply(x.sep, function(x,rep) paste(x[-c(1,2,if(!rep){length(x)})],collapse="_"), rep=rep)
  return(x.glue)
}

# Function to plot posterior densities along with true values.
# Requires that the original data (truevals) has the structure: truevals$truevalues$...
post_dens <- function(postdraws, truevals){
  library(ggplot2)
  record.list <- unique(postdraws$variable) # List of posterior quantities
  truth <- data.frame(variable=record.list, trueval=NA)
  # Leverages the square-bracket nomenclature from stan summary output 
  # to index true values from the original data simulation
  for(i in 1:length(record.list)) 
    truth$trueval[i] <- eval(parse(text=paste("truevals$truevalues$",record.list[i],sep="")))
  nm <- deparse(substitute(postdraws))      # Collect name of submitted postdraws argument for titling purposes.
  print(
    ggplot(data=postdraws, aes(value)) +
      geom_density() +
      geom_vline(data=truth, aes(xintercept=trueval), size=1.25, col="magenta") +
      theme_bw() +
      facet_wrap(~variable, scale="free") +
      ggtitle(nm)
  )
}




##### Functions to help generate posterior predictive plots

# Calculate the proportion of counts for each distance-time bin by posterior draw
#   The input 'postp' should be generated from postp_dataset()
#   Output is a single data frame with counts and proportions by distance-time-bin across all 
#     posterior predictive draws.
postp_fxn <- function(postp, db=distbins, tb=timebins){
  bind_rows(postp, .id="Rep") %>%
    mutate(distbin = cut(distances, breaks = db, include.lowest = TRUE),
           timebin = cut(times,     breaks = tb, include.lowest = TRUE)) %>%
    group_by(Rep) %>%
    mutate(Count = n()) %>%
    group_by(Rep, distbin, timebin, Count) %>%
    summarize(N = n()) %>%
    ungroup() %>%
    mutate(prop = N/Count)
}

# Calculate posterior predictive p-values (proportion of posterior predictive simulations < observed)
#   by distance-time-bin
# Includes 0.025, 0.25, 0.5, 0.75, 0.975 quantiles
pval_fxn <- function(binpost, BD=bindat){
  binpost %>%
    left_join(BD %>% 
                rename(obs_prop = "prop", obs_N = "N") %>%
                select(distbin, timebin, obs_N, obs_prop)) %>%
    group_by(distbin, timebin) %>%
    summarize(p_N = 0.5*(sum(N < obs_N) + sum(N <= obs_N))/length(Rep),           # Note: average (p_pred < obs) and (p_pred <= obs)
              p_prop = 0.5*(sum(prop < obs_prop) + sum(prop <= obs_prop))/length(Rep),
              N_L95 = quantile(N, 0.025),          # Calculate specific quantiles
              N_U95 = quantile(prop, 0.975),
              prop_L95 = quantile(N, 0.025),
              prop_U95 = quantile(prop, 0.975),
              N_L50 = quantile(N, 0.25),
              N_U50 = quantile(prop, 0.75),
              prop_L50 = quantile(N, 0.25),
              prop_U50 = quantile(prop, 0.75),
              N_med = quantile(prop, 0.5),
              prop_med = quantile(N, 0.5)
    )
}

### Code to simulate a dataset from a complete set of hyperparameter values.
#   Output is a dataset in the same format as the GCSP data
postp_dataset <- function(method="Event", n_surv=600, 
                          intcpt_ab= 1.25, bab=NULL, n_lvlab=NULL, sigma_ab=NULL,
                          intcpt_avl=-1.5, bavl=NULL, n_lvlavl=NULL, sigma_avl=NULL, 
                          intcpt_pev=-0.05, bpev=NULL, n_lvlpev=NULL, sigma_pev=NULL, 
                          intcpt_pst=-0.05, bpst=NULL, n_lvlpst=NULL, sigma_pst=NULL, 
                          mixing=1, avl_diff=0, maxdist=300, tau=10,
                          Xab = NULL, Xavl = NULL, Xpev = NULL, Xpst = NULL,
                          Zab = NULL, Zavl = NULL, Zpev = NULL, Zpst = NULL)
{

  library(plyr)
  
  # Extract dimensions
  n_bab    <- length(bab)
  n_bavl   <- length(bavl)
  n_bpev   <- length(bpev)
  n_bpst   <- length(bpst)
  n_rdmab  <- length(n_lvlab)
  n_rdmavl <- length(n_lvlavl)
  n_rdmpev <- length(n_lvlpev)
  n_rdmpst <- length(n_lvlpst)
  groups   <- length(avl_diff)
  
  stopifnot(method %in% c("Combined", "Event", "Survey")) 
  stopifnot(sum(mixing)==1) 
  stopifnot(avl_diff[1]==0)
  stopifnot(length(mixing)==length(avl_diff))
  stopifnot(n_rdmab==length(sigma_ab),
            n_rdmavl==length(sigma_avl),
            n_rdmpev==length(sigma_pev),
            n_rdmpst==length(sigma_pst))
  
  if(method=="Event")  intcpt_pst <- Inf  # In a pure event model, all individuals are perceivable
  if(method=="Survey") intcpt_pev <- Inf  # In a pure state model, all perceivable individuals are detected when available
  
  
  ### Helper functions
  # Distance function
  g.r <- function(distances, sigma) exp(-(distances/sigma)^2)
  # Distribution for covariates
  covar_dist <- function(n_surv) scale(runif(n_surv))[,1]    # Randomly set Xab, Xavl, Xpev, Xpst ~ Unif.
  # Values are standardized to mean=0 and variance=1, matching our approach for real data
  
  ### Build abundance effects.  Observational unit = survey
  Xab <- cbind(rep(1, n_surv), Xab)
  rab <- array(0, dim=0)            # Instantiate a zero-length array for random effect sizes (Stan compatible)
  if(n_rdmab > 0) {
    rab <- rnorm(sum(n_lvlab), mean=0, sd=rep(sigma_ab, times=n_lvlab)) # Generate random effect sizes
  }  
  # loglambda = log(expected abundance)
  loglambda <- Xab %*% as.matrix(c(intcpt_ab, bab))
  if(n_rdmab > 0) for(i in 1:n_rdmab) loglambda <- loglambda + rab[Zab[,i]]
  # Actual survey-level abundance
  N <- rpois(n_surv, exp(loglambda))
  
  ### Build detection effects.  Observational unit = bird, but covariates are survey-level.  
  # Create survey-level matrix first following same steps as above
  Xavl <- cbind(rep(1, n_surv), Xavl)
  ravl <- array(0, dim=0)                  # A zero-length array for random effect sizes (like abundance above).
  if(n_rdmavl > 0) {
    ravl <- rnorm(sum(n_lvlavl), mean=0, sd=rep(sigma_avl, times=n_lvlavl)) # Generate random effect sizes
  }
  # logphi = log(availability rate) assuming no group heterogeneity... heterogeneity added below
  logphi <- Xavl %*% as.matrix(c(intcpt_avl, bavl))
  if(n_rdmavl > 0) for(i in 1:n_rdmavl) logphi <- logphi + ravl[Zavl[,i]]
  # Create bird-level matrix
  grp    <- sample(1:groups, sum(N), replace=T, prob=mixing) # Assign each bird to a group
  logphi <- rep(logphi[N>0], times=N[N>0])          # logphi for each bird (ignoring group)
  logphi <- logphi + avl_diff[grp]                  # Incorporate group detection heterogeneity on the intercept term
  survey <- factor(1:n_surv)                        # Index of surveyIDs for survey-level data
  survey.bird <- rep(survey[N>0], times=N[N>0])     # Index of surveyIDs for bird-level data
  
  ### Build event perceptibility effects.
  Xpev <- cbind(rep(1, n_surv), Xpev)
  rpev <- array(0, dim=0)          # Instantiate a zero-length array for random effect sizes (Stan compatible)
  if(n_rdmpev > 0) {
    rpev <- rnorm(sum(n_lvlpev), mean=0, sd=rep(sigma_pev, times=n_lvlpev)) # Generate random effect sizes
  }  
  # logpsig = log(distance parameter standardized to maximum detectable distance)
  log_pev_sig <- Xpev %*% as.matrix(c(intcpt_pev, bpev))
  if(n_rdmpev > 0) for(i in 1:n_rdmpev) log_pev_sig <- log_pev_sig + rpev[Zpev[,i]]
  
  ### Build state perceptibility effects.
  Xpst <- cbind(rep(1, n_surv), Xpst)
  rpst <- array(0, dim=0)          # Instantiate a zero-length array for random effect sizes (Stan compatible)
  if(n_rdmpst > 0) {
    rpst <- rnorm(sum(n_lvlpst), mean=0, sd=rep(sigma_pst, times=n_lvlpst)) # Generate random effect sizes
  }  
  # logpsig = log(distance parameter standardized to maximum detectable distance)
  log_pst_sig <- Xpst %*% as.matrix(c(intcpt_pst, bpst))
  if(n_rdmpst > 0) for(i in 1:n_rdmpst) log_pst_sig <- log_pst_sig + rpst[Zpst[,i]]
  
  ### Data -- observed and unobserved
  distances <- sqrt(runif(sum(N), min=0, max=maxdist^2))                              # RNG for distance from a uniform spatial distribution (P.I.T.)
  
  state <- rbinom(sum(N), 1, g.r(distances, exp(log_pst_sig[survey.bird])*maxdist))   # Whether or not each bird is state-perceivable = g_s(r)
  times <- rexp(sum(N),  
                g.r(distances, exp(log_pev_sig[survey.bird])*maxdist)*exp(logphi))    # Event-level model: rate = g(r)*phi
  times[!state] <- Inf                                                                # Non-perceivable birds are never detected
  
  df <- data.frame(survey = survey.bird, Time = times)
  n  <- plyr::ddply(df[df$Time < tau,], .(survey), summarize, n=length(Time))         # Counts for surveys where birds are present
  n  <- join(data.frame(survey=survey), n, by="survey"); n[is.na(n)] <- 0             # Counts by survey (includes zeros)
  
  dfr = data.frame(survey = survey.bird[times <= tau],
                   distances = distances[times <= tau], 
                   times = times[times <= tau])
  
  return(dfr)
}
