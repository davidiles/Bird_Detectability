library(plyr)
library(tidyverse)
library(rstan)

# Directory: Adam/IAState/DistRemPaper/Code

### Plot of the observed distribution of distances faceted by minute of detection
source("Scripts.R")

# Load GCSP data
dat <- readRDS("gcspdata.rds")

# Transform data for plotting
dfr <- data.frame(Time     = dat$time_obs,
                  TimeF    = cut(dat$time_obs, seq(0,5,1), include.lowest=T, labels=F),
                  # The following is a trick to triplicate the data, allowing the smoothed function
                  #   to not disappear to zero at r = 0 and r = maxdist
                  Dist     = dat$dist_obs,
                  DistNeg  = -dat$dist_obs,
                  DistPlus = 2*dat$maxdist-dat$dist_obs) %>% gather(DistCat, Distance, 3:5)
# Set factor levels for minute-to-detection
dfr$TimeF <- factor(dfr$TimeF, labels=paste0("Minute ", 1:5))
# Labels for average detection distance
textlab <- data.frame(TimeF=levels(dfr$TimeF), label=paste0("Mean = ", c(130, 144, 148, 167, 173)))
# Generate plot
ggplot(dfr, aes(x=Distance, color=TimeF)) + 
  stat_density(geom="line", size=2, adjust=0.25) +
  geom_text(data=textlab, aes(x=140, y=0.00005, label=label), size=4.5, color="black") +
  facet_grid(~TimeF) + 
  xlab("Distance") + ylab("Conditional\nDensity") + 
  coord_cartesian(xlim=c(12,268)) +
  theme_bw() + 
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text.x = element_text(size=12)) + 
  scale_color_manual(values=paste0("gray", seq(25,73,12)))

# Save output
ggsave("DistByTime.png", width=7, height=2.7, units="in")



### Create the summary of GCSP plots
#   Reminder of which predictor goes with which perceptibility
#   State pcpt: proportion dense habitat
#   Event pcpt: observer group

# Parameters to plot
par.keep <- c("global_p", "mean_pavl", "mean_pp", "mean_pst", "mean_pev",
              "bab[1]", "bab[2]", "bab[3]", "bavl[1]", "bpst[1]", "bpst[2]", "bpev[1]", "bpev[2]")

# Load summary tables from model outputs (from 'ExtractFits.R')
summaries <- readRDS("Outputs.rds")$summaries %>%
  filter(model %in% c("gcspfit", "gcsp_state", "gcsp_event") &
           parameter %in% par.keep) %>%
  gather("Quantile", "qpost", contains("X"), mean) %>%
  left_join(data.frame(model = c("gcspfit", "gcsp_state", "gcsp_event"),
                       .id   = c("Comb_orig", "State", "Event"),
                       stringsAsFactors=FALSE)) %>%
  select(.id, parameter, Quantile, qpost) %>%
  rename(model = .id) %>%
  # Because, in State and Event models, both perceptibility parameters were either bpst or bpev
  mutate(parameter = replace(parameter, parameter == "bpev[2]", "bpst[1]")) %>%
  mutate(parameter = replace(parameter, parameter == "bpst[1]" & model == "State", "bpev[1]")) %>%
  mutate(parameter = replace(parameter, parameter == "bpst[2]", "bpst[1]"))
  
# Load estimates of abundance density
density <- readRDS("Outputs.rds")$density %>% .[1:3] %>%
  ldply(function(x) c(quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(x))) %>%
  setNames(c("model", unique(summaries$Quantile))) %>%
  mutate(parameter = "Density") %>%
  gather("Quantile", "qpost", contains("X"), mean)

# Labels for facets in plot  
labvec <- c("Density", "Abund: Shrub", "Abund: Deciduous", "Abund: Spruce",
            "Avail: Julian Date", "Pcpt: Observer", "Pcpt: Dense Cover", "Pr(Detection)",
            "Availability", "PP Pcpt", "Perceptibility", "IJ Pcpt")
# Combine everything into one data frame, carefully numbering the parameter label factors so as
#   to generate the desired output faceting
plotdf <- bind_rows(summaries, density) %>%
  mutate(model = factor(model, labels=c("IJ+PP", "PP", "IJ"))) %>%
  # Associate parameter names with descriptive label
  left_join(data.frame(
    parameter = c("Density", sort(unique(summaries$parameter))),
    Label = labvec,
    stringsAsFactors=F)) %>%
  mutate(Label = factor(Label, levels = labvec[c(1,11,8,12,9,10,2,5,3,7,4,6)]))

# Partition into two datasets for left-hand and right-hand sets of plots (six each)
plotdf1 <- plotdf %>% filter(!grepl("\\[", parameter))
plotdf2 <- plotdf %>% filter(grepl("\\[", parameter))
  
# Vertical dashed lines at zero
vlines <- data.frame(
  Label     = labvec[c(2:7)],
  qpost     = 0,
  stringsAsFactors = FALSE) %>%
  mutate(Label = factor(Label, levels = labvec[c(1,11,8,12,9,10,2,5,3,7,4,6)]))

# Very manually built function for making x-axis ticks
facetbreaks <- function(x){
  m <- min(x)
  M <- max(x)
  if(m>0.2 & m<0.35 & M>0.5) ticks <- c(0.35, 0.5)
  if(m<0.2 & M>0.5 & M<0.6) ticks <- seq(0.2, 0.5, 0.1)
  if(m<0.2 & M>0.6) ticks <- seq(0.2, 0.6, 0.2)
  if(M > 0.8) ticks <- seq(0.6, 1, 0.2)
  if(M < 0.25) ticks <- seq(0.08, 0.2, 0.06)
  return(ticks)
}

p1 <- ggplot() +
  geom_line(data=plotdf1 %>% filter(grepl("\\.5", Quantile)), 
            aes(x=qpost, y=model), size=1.2) +
  geom_line(data=plotdf1 %>% filter(grepl("X.5", Quantile)), 
            aes(x=qpost, y=model), size=2.2, color="gray55") +
  geom_point(data=plotdf1 %>% filter(Quantile=="mean"), 
             aes(x=qpost, y=model), size=2.5) +
  facet_wrap(~Label, nrow=3, scale="free_x") +
  ylab("Model") +
  scale_x_continuous(NULL, breaks=facetbreaks) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12)
  )

p2 <- ggplot() +
  geom_vline(data=vlines, aes(xintercept=qpost), linetype="dashed") +
  geom_line(data=plotdf2 %>% filter(grepl("\\.5", Quantile)), 
            aes(x=qpost, y=model), size=1.2) +
  geom_line(data=plotdf2 %>% filter(grepl("X.5", Quantile)), 
            aes(x=qpost, y=model), size=2.2, color="gray55") +
  geom_point(data=plotdf2 %>% filter(Quantile=="mean"), 
             aes(x=qpost, y=model), size=2.5) +
  facet_wrap(~Label, nrow=3, scale="free_x") +
  ylab(NULL) + xlab(NULL) + xlim(-0.6, 0.8) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank())

library(ggpubr)
ggarrange(p1, p2, ncol=2, widths=c(14,10))

ggsave("Estimates.pdf", width=6.1, height=4.2, units="in")





### Part IV: Posterior predictive p-values on binned data

# TO MAKE THIS WORK, I HAVE TO GENERATE POSTERIOR PREDICTIVE TIME/DISTANCE BINNED COUNTS FROM POSTERIOR SAMPLES!!
# Would have been simpler to code into the original Stan code

dat <- readRDS("gcspdata.rds")

fit.ev <- readRDS("gcsp_event.rds")
fit.st <- readRDS("gcsp_state.rds")
fit.co <- readRDS("gcspfit.rds")

# ### Generate posterior predictive draws from the combined model

fit <- rstan::extract(fit.co)       # Posterior draws
pars <- names(fit)                  # Parameter names vector
for(p in pars) assign(paste0("f_", p), fit[[p]])  # Assign object to posterior draws of each parameter

postp.co <- vector('list', length(f_intcpt_ab))
for(indx in 1:length(f_intcpt_ab)){
  postp.co[[indx]] <- postp_dataset(method = "Combined", n_surv = dat$n_surv,
                                    intcpt_ab=f_intcpt_ab[indx], bab=f_bab[indx,], n_lvlab=dat$n_rabs, sigma_ab=f_sigma_ab[indx],
                                    intcpt_avl=f_intcpt_avl[indx], bavl=f_bavl[indx,], n_lvlavl=NULL, sigma_avl=NULL,
                                    intcpt_pev=f_intcpt_pev[indx], bpev=f_bpev[indx,], n_lvlpev=NULL, sigma_pev=NULL,
                                    intcpt_pst=f_intcpt_pst[indx], bpst=f_bpst[indx,], n_lvlpst=NULL, sigma_pst=NULL,
                                    mixing=f_mixing[indx,], avl_diff=c(0,f_intcpt_g[indx]), maxdist=dat$maxdist, tau=dat$tau,
                                    Xab = dat$Xab, Xavl = dat$Xavl, Xpev = dat$Xpev, Xpst = dat$Xpst,
                                    Zab = dat$Zab, Zavl = NULL, Zpev = NULL, Zpst = NULL)
}

### Generate posterior predictive draws from the event model
dat <- readRDS("gcspdata.rds")

# Data alterations
dat <- dat
dat$n_bpev <- 2
dat$n_bpst <- NULL
dat$Xpev   <- cbind(dat$Xpev, dat$Xpst)
dat$Xpst   <- NULL

fit <- rstan::extract(fit.ev)       # Posterior draws
pars <- names(fit)                  # Parameter names vector
for(p in pars) assign(paste0("f_", p), fit[[p]])  # Assign object to posterior draws of each parameter

postp.ev <- vector('list', length(f_intcpt_ab))
for(indx in 1:length(f_intcpt_ab)){
  postp.ev[[indx]] <- postp_dataset(method = "Event", n_surv = dat$n_surv,
                                    intcpt_ab=f_intcpt_ab[indx], bab=f_bab[indx,], n_lvlab=dat$n_rabs, sigma_ab=f_sigma_ab[indx],
                                    intcpt_avl=f_intcpt_avl[indx], bavl=f_bavl[indx,], n_lvlavl=NULL, sigma_avl=NULL,
                                    intcpt_pev=f_intcpt_pev[indx], bpev=f_bpev[indx,], n_lvlpev=NULL, sigma_pev=NULL,
                                    intcpt_pst=NULL, bpst=NULL, n_lvlpst=NULL, sigma_pst=NULL,
                                    mixing=f_mixing[indx,], avl_diff=c(0,f_intcpt_g[indx]), maxdist=dat$maxdist, tau=dat$tau,
                                    Xab = dat$Xab, Xavl = dat$Xavl, Xpev = dat$Xpev, Xpst = dat$Xpst,
                                    Zab = dat$Zab, Zavl = NULL, Zpev = NULL, Zpst = NULL)
}


### Generate posterior predictive draws from the state model
dat <- readRDS("gcspdata.rds")

# Data alterations
dat$n_bpst <- 2
dat$n_bpev <- NULL
dat$Xpst   <- cbind(dat$Xpev, dat$Xpst)
dat$Xpev   <- NULL

# Load fit and define parameter
fit <- rstan::extract(fit.st)       # Posterior draws
pars <- names(fit)                  # Parameter names vector
for(p in pars) assign(paste0("f_", p), fit[[p]])  # Assign object to posterior draws of each parameter

postp.st <- vector('list', length(f_intcpt_ab))
for(indx in 1:length(f_intcpt_ab)){
  postp.st[[indx]] <- postp_dataset(method = "Survey", n_surv = dat$n_surv,
                                    intcpt_ab=f_intcpt_ab[indx], bab=f_bab[indx,], n_lvlab=dat$n_rabs, sigma_ab=f_sigma_ab[indx],
                                    intcpt_avl=f_intcpt_avl[indx], bavl=f_bavl[indx,], n_lvlavl=NULL, sigma_avl=NULL,
                                    intcpt_pev=NULL, bpev=NULL, n_lvlpev=NULL, sigma_pev=NULL,
                                    intcpt_pst=f_intcpt_pst[indx], bpst=f_bpst[indx,], n_lvlpst=NULL, sigma_pst=NULL,
                                    mixing=f_mixing[indx,], avl_diff=c(0,f_intcpt_g[indx]), maxdist=dat$maxdist, tau=dat$tau,
                                    Xab = dat$Xab, Xavl = dat$Xavl, Xpev = dat$Xpev, Xpst = dat$Xpst,
                                    Zab = dat$Zab, Zavl = NULL, Zpev = NULL, Zpst = NULL)
}

### Read in posterior predictives for each model.
#   For postp.co, it's a list of 2500 posterior predictive data frames with columns for survey, distance, and time
#   A total of about 1.3M rows!
## I recommend saving copies of the above posterior predictive samples to save future computation time
# postp.st <- readRDS("state_postp.rds")
# postp.ev <- readRDS("event_postp.rds")
# postp.co <- readRDS("combo_postp.rds")

# Quintiles of observed distances and times
# distbins <- quantile(dat$dist_obs, seq(0,1,0.2), type=4); distbins[1] <- 0
distbins <- c(0, 79, 115.5, 155, 202, 280)    # A more even split (by distance) is achieved by using 79 instead of 80
timebins <- quantile(dat$time_obs, seq(0,1,0.2), type=4); timebins[1] <- 0

# Bin the data / posterior predictive draws according to distance and time
bindat     <- data.frame(distances=dat$dist_obs, times=dat$time_obs) %>% postp_fxn()
binpost.st <- postp_fxn(postp.st)
binpost.ev <- postp_fxn(postp.ev)
binpost.co <- postp_fxn(postp.co)

# Calculate posterior predictive p-values for 
pval.st <- pval_fxn(binpost.st)
pval.ev <- pval_fxn(binpost.ev)
pval.co <- pval_fxn(binpost.co)

# A hand-coded simplification of labels
levels(binpost.st$timebin) <- levels(binpost.ev$timebin) <- levels(binpost.co$timebin) <- 
  levels(pval.st$timebin) <- levels(pval.ev$timebin) <- levels(pval.co$timebin) <- 
  c("[0,0.19]", "(0.19,0.47]", "(0.47,0.96]", "(0.96,2.33]", "(2.33,5]")

# Plot proportion of observed/modeled distances and times in each bin
#  The diagonal pattern in the data is curiously strong.
#  If we do an equivalent plot for a random posterior draw from the combined model, we find a much noisier pattern.
data_plot <- function(BD, title){
  ggplot(BD, aes(x=distbin, y=timebin, fill=-prop, label=round(prop,2))) + 
    geom_raster() + 
    geom_text(data = BD %>% filter(prop >= 0.04), color=I("white")) + 
    geom_text(data = BD %>% filter(prop < 0.04), color=I("black")) + 
    # geom_text(data=BD %>% filter(prop > 0.06), color=I("white")) + 
    scale_fill_gradient2(midpoint=-0.04, limits=c(-0.096,0), high="white", mid="gray50", low="black") +
    # scale_fill_gradient2(midpoint=-0.04, limits=c(-0.096,0), high="#F7FBFF", mid="#6BAED6", low="#08306B") +
    ggtitle(title) + ylab("Time") + xlab("Distance") +
    theme_bw() + theme(legend.position="none")
}

p.obs <- data_plot(bindat, title="GCSP data")
p.s   <- data_plot(binpost.st %>% group_by(distbin, timebin) %>% summarize(N=mean(N), prop=mean(prop)), title="IJ Model Mean Posterior")
p.e   <- data_plot(binpost.ev %>% group_by(distbin, timebin) %>% summarize(N=mean(N), prop=mean(prop)), title="PP Model Mean Posterior")
p.c   <- data_plot(binpost.co %>% group_by(distbin, timebin) %>% summarize(N=mean(N), prop=mean(prop)), title="IJ+PP Model Mean Posterior")


# Dark squares had relatively more observed values than posterior predictives, while
#   light squares had relatively fewer observed values.
# Note the intervals... in distance, the middle three bins are narrower
#   In time, the bottom three rows are all within the first minute
plot_fxn <- function(dfr, title){
  ggplot(dfr, aes(x=distbin, y=timebin, fill=p_prop, label=round(p_prop,2))) + 
    geom_raster() + 
    geom_text(data = dfr %>% filter(p_prop > 0.025 & p_prop < 0.5), color="white") + 
    geom_text(data = dfr %>% filter(p_prop >= 0.5 & p_prop < 0.975)) + 
    geom_label(data = dfr %>% filter(p_prop < 0.025 | p_prop > 0.975), color=I("black"), fontface="bold", fill=I("white")) + 
    scale_fill_gradient2(low="black", high="white", mid="gray50", midpoint=0.5) +
    # geom_label(data = dfr %>% filter(p_prop < 0.025 | p_prop > 0.975), color=I("red"), fontface="bold", fill=I("white")) + 
    # scale_fill_gradient2(low="blue", high="red", midpoint=0.5) +
    theme_bw() + theme(legend.position="none") + 
    xlab("Distance") + ylab("Time") +
    ggtitle(title)
}

p.st <- plot_fxn(pval.st, title="IJ Model p-values")
p.ev <- plot_fxn(pval.ev, title="Cue-Based Model p-values")
p.co <- plot_fxn(pval.co, title="Animal-Cue Model p-values")

library(ggpubr)
ggarrange(p.obs, NULL, p.s, p.st, p.e, p.ev, p.c, p.co, nrow=4, ncol=2)
ggsave("postp.png", width=8, height=8, units="in")





# Just for a sense of the bin scaling with respect to observations
# Plot observed times vs. observed distances with a grid representing quintiles
ggplot(data.frame(dist_obs = dat$dist_obs, time_obs = dat$time_obs), 
       aes(x=dist_obs, y=time_obs)) + geom_jitter(width=1, height=0.02) +
  geom_hline(yintercept=timebins) + geom_vline(xintercept=distbins) +
  theme_bw() + ylab("Observed Time") + xlab("Observed Distance")
