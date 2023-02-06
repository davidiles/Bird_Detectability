##### The following data handling code is minimally altered from the code used in the paper:

# Amundson, C. L., C. M. Handel, D. R. Ruthrauff, T. L. Tibbitts, and R. E. Gill, Jr. 2018. Montane-breeding bird distribution and abundance across national parks of southwestern Alaska. Journal of Fish and Wildlife Management 9(1):180-207. doi:10.3996/062017-JFWM-050

# With data and code available to download at:
# Amundson, C. L., Handel, C. M., Ruthrauff, D. R., Tibbitts, T. L. and Gill, R. E., Jr , 2018, Data for Montane-breeding Bird Distribution and Abundance across National Parks of Southwestern Alaska, 2004-2008: U.S. Geological Survey data release, https://doi.org/10.5066/F7MW2GD3.

# The creation of distance classes and time intervals is not required for this analysis
# Code to fit a JAGS model has been switched off by introducing a toggle variable 'fit.JAGS'



####### C. Amundson USGS AK Science Center, AK montane parks JAGS model
#camundson@usgs.gov
# Note, model will mostly likely not run in winBUGS or openBUGS
# ID variables for jags file
# y = individuals per site
# tran = transect ID for each site
# observer = condense 5 observers to 2 observer groups – see main document methods
# park = national park ID per site
# jdate = fixed covariate
# water = fixed covariate
# dec = fixed covariate
# mixed = fixed covariate
# spruce = fixed covariate
# dshrub = fixed covariate
# baresnow = fixed covariate
# shrub = fixed covariate
# sound = fixed covariate
# tinterval = time interval category ONLY FOR OBSERVATIONS − here coded for equal time intervals
# dclass = distance class ONLY FOR OBSERVATIONS
# surveyid = site ID associated with each observation ONLY FOR OBSERVATIONS


# Change this to 1 if I want to fit the Amundson model in JAGS
fit.JAGS <- 0
maxdist  <- 280
species.to.model <- "GCSP"

# Begin data formatting
# Import data from local directory.  Data can be obtained from https://doi.org/10.5066/F7MW2GD3
obsdata  <- read.csv("montaneBird_distribution_SW_obsdata_amundson_2004_2008.csv")
sitecovs <- read.csv("montaneBird_distribution_SW_sitecovs_amundson_2004_2008.csv")


# activate packages
library(plyr)

# need to make sure each data frame is ordered according to survey ID - 
# need a numeric ID by site starting at 1
sitecovs$PTP <- as.factor(sitecovs$PTP)
obsdata$PTP  <- as.factor(obsdata$PTP)
sitecovs <- arrange(sitecovs,PTP)
obsdata  <- arrange(obsdata,PTP)
sitecovs$PTPnum <- seq(1,nrow(sitecovs))
# Site information
water    <- sitecovs$water + sitecovs$wetland
dec      <- sitecovs$dec
mixed    <- sitecovs$mixed
spruce   <- sitecovs$spruce
dshrub   <- sitecovs$dshrub
wetland  <- sitecovs$wetland
baresnow <- sitecovs$baresnow
shrub    <- sitecovs$shrub
sound    <- sitecovs$sound
# use GCSP as example
spdat <- subset(obsdata, species==species.to.model)
# assign breaks in distance bins, these were determined by examining density plots (see Methods) 
# Cut points used for remaining species can be obtained by contacting C. Amundson, USGS
brks <- c(92,147,212,280)
# Don't use observations with missing distance or time removal information
gc <- subset(spdat,exact_distance !="NA"&!is.na(elapsed_secs))
# subset to males only
# separate out observations of > 1 male into two rows of data.
#   note: code leverages the fact that max(number) = 2.  Needs editing if max(number) > 2 -- Adam
gc.x <- subset(gc, number>1) # here 1 observation
gcs  <- merge(gc.x,gc,all=T)
gcs$number <- 1
spec <- subset(gcs, exact_distance <= maxdist & sex=="M")

############# create distance classes and time interval bins
spec$dclass <-
  with(spec, ifelse(exact_distance <= brks[1], 1, ifelse(
    exact_distance > brks[1] & exact_distance <= brks[2], 2, ifelse(
      exact_distance > brks[2] & exact_distance <= brks[3], 3, 4))))
spec$pass <-
  as.factor(with(spec, ifelse(elapsed_secs < 100, 1, ifelse(
    elapsed_secs >= 100 & elapsed_secs < 200, 2, 3))))
# Need to get zero counts for sites where no GCSP were observed
# subset observation data to necessary columns: PTP, elapsed_secs, number, exact_distance, dclass, pass
specm <- spec[, c("PTP", "elapsed_secs", "number", "exact_distance", "dclass", "pass")]
# merge with site data and keep all values, now sites without birds are NA
specobs <- merge(specm, sitecovs, by="PTP", all=T)
# time removal
specobs$pass <- as.numeric(specobs$pass)
# create numeric site ID to maintain consistent sorting among datasets
specobs <- arrange(specobs,PTPnum)
# subset data to only observations - to get dclass and tinterval values connected with surveyID
specobsonly <- subset(specobs, number>0)
specobsonly <- arrange(specobsonly,PTPnum)
# now tally counts per site for site-level of model
specsite  <- ddply(specobs,~PTP,summarize,count=sum(number),PTPnum=mean(PTPnum))
specsites <- merge(specsite,sitecovs)
specsites$count[is.na(specsites$count)]=0 # make NAs zero counts
# need to create transect variable that is unique for each plot (combine plot and tran)
t <- ddply(specsites, .(plot,transect), summarize, test=mean(count))
t <- arrange(t,plot,transect)
t$tran  <- 1:169
gcsites <- merge(specsites,t,all=T)
gcsites <- arrange(gcsites,PTPnum)



# now organize data for bugs file
# observation only data
tinterval <- specobsonly$pass # time interval for each observation
dclass    <- specobsonly$dclass # distance class for each observation
surveyid = specobsonly$PTPnum # site ID for each observation
J = max(unique(as.numeric(tinterval))) # number of time intervals
# site level data
y <- gcsites$count # site level counts
nsurveys   <- length(y) # number of survey points
park       <- as.numeric(gcsites$park) # park ID
tran       <- gcsites$tran # transect ID
ntrans     <- length(unique(tran)) # number of transects
observer   <- gcsites$obsnum # observer ID
nobservers <- length(unique(observer)) # number of observers
jdate      <- as.numeric(scale(gcsites$jdate)) # season date; standardize to facilitate convergence
ncov       <- 7 # number of covariates in abundance model\
# now for distance component
nobs    <- length(tinterval) # number of observations
brk     <- brks/100 # here scaled breaks by 100 to facilitate convergence
maxd    <- max(brk) # max distance - here 2.8 m because its scaled
# NOTE: sigma parameter estimates will also be scaled to this max distance
nbreaks <- length(brk) # number of distance bins
breaks  <- c(0, brk) # list cutpoints, add zero
delta   <- c(brk[1],brk[2]-brk[1],brk[3]-brk[2],brk[4]-brk[3]) # number of meters (or fractions of meters) in each bin
mdpts   <-c(brk[1]/2,(brk[1]+brk[2])/2,(brk[2]+brk[3])/2,(brk[3]+brk[4])/2) # really midpoints of each bin

if(fit.JAGS==1){
  library(jagsUI)
  
  # Bundle data for JAGS model
  win.data <-
    list(y=y,surveyid=surveyid,dclass=dclass,nsurveys=nsurveys,nobs=nobs,delta=delta,nbreaks=nbreaks, mdpts=mdpts,
         maxd=maxd,J=J,tinterval=tinterval,dshrub=dshrub,dec=dec,shrub=shrub, mixed=mixed,spruce=spruce,water=water,
         baresnow=baresnow,wetland=wetland,tran=tran, ncov=ncov,ntrans=ntrans,jdate=jdate,observer=observer,
         nobservers=nobservers,sound=sound)
  
  # Initial values for counts must be close to predictions
  Nst <- y+1
  # Inits function
  inits <- function(){list(N=Nst,emp.new=Nst,sigma0=runif(2,10,20), beta0=runif(169,-1,1),
                           beta0.a=runif(1,-1,1),beta.a1=runif(1,-1,1),beta.p1=runif(1,-1,1),
                           beta=runif(ncov,-1,1),sd.obs=runif(1,0,20),
                           sd.tran=runif(1,0,2))} #
  # parameters to estimate
  params<-c("bayes.p","dens","beta","sigma0","meanpcap","beta0.a",
            "beta.p1","beta.a1","meanpdet","meansig","beta0","aniaN","katmN","laclN",
            "aniadens","katmdens","lacldens","pdet","pavail","N","mu.tran",
            "meanpavail","sd.tran","sd.obs","mu.obs","totN","fit","fit.new")
  # MCMC settings
  nc <- 3
  ni <- 1000 # increase iterations to reach convergence
  nb <- 250
  nt <- 2
  # Start Gibbs Sampler
  # Code to create JAGS model starts on page 7. habmod.txt must exist in your working directory
  # for the following code to work
  fit.spec <- jags(data=win.data, inits=inits, parameters.to.save=params, model.file="habmod.txt",
                   n.chains=nc, n.burnin=nb, n.iter=ni, codaOnly =
                     c("N","pavail","pdet","fit","fit.new","beta0"),
                   n.thin=nt, modules = c("basemod","bugs","glm","dic"),parallel=TRUE)
  fit.spec
  
  
  # This is the 'Bayesian P-value' that has always confused me.  It tests whether the sum of residuals 
  #   (on the data scale) is larger in the model fit than in the posterior predictive distribution
  # GOF within jagsUI, should be about the same as bayes.p from output
  pp.check(fit.spec,"fit","fit.new")
  
  
  # Summarize output and plot predictions
  specpdet   <- colMeans(fit.spec$sims.list$pdet)
  specpa     <- colMeans(fit.spec$sims.list$pavail)
  specN      <- colMeans(fit.spec$sims.list$N)
  specdens   <- (colMeans(fit.spec$sims.list$N)*2)/(maxd*maxd*100*100*3.14159/10000)
  specCIpdet <- t(apply(fit.spec$sims.list$pdet,2,quantile,probs=c(0.025,0.975)))
  specCIN    <- t(apply(fit.spec$sims.list$N,2,quantile,probs=c(0.025,0.975)))
  specCIpa   <- t(apply(fit.spec$sims.list$pavail,2,quantile,probs=c(0.025,0.975)))
  specCIdens <- (specCIN*2)/(maxd*maxd*100*100*3.14159/10000)
  # some quick and dirty summary plots
  plot(as.factor(observer), specpdet)
  plot(sound,specpdet) # intercept offset for each observer group
  plot(sitecovs$jdate,specpa)
  scatter.smooth(shrub,specdens)
}

