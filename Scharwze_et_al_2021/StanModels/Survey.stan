// This model fits datasets using an independent joint model (IJ).


functions{
  real g_r(real r, real sigma){ 
    return exp(-(r/sigma)^2);
  }

  // Mean state perceptibility given a half-normal distance function and uniform species density
  //   Note: r and maxdist are irrelevant for this calculation
  real p_state(real sigma){
    return sigma^2 * (1 - g_r(1, sigma));
  }
}



data {
// Data dimensions
int<lower=0> n_surv;    // Number of surveys 
int<lower=0> n_bird;    // Total birds counted
int<lower=1,upper=2> groups; // Number of behavioral groups to be modeled (only two for now)

// Methodologically defined parameters
real<lower=0> tau;      // Duration of surveys
real<lower=0> maxdist;  // Maximum distance used
real ab_prior;          // Prior mean on log-abundance
real ab_prior_sd;       // Std deviation for prior log-abundance
real avl_prior;         // Prior mean on log-availability rate
real avl_prior_sd;      // Std deviation for prior log-availability rate

// Fixed Effects
int<lower=0> n_bab;         // Number of abundance fixed effects
int<lower=0> n_bavl;        // Number of availability fixed effects
int<lower=0> n_bpst;        // Number of state perceptibility fixed effects
matrix[n_surv,n_bab]  Xab;  // Abundance fixed effect covariates
matrix[n_surv,n_bavl] Xavl; // Availability fixed effect covariates
matrix[n_surv,n_bpst] Xpst; // State perceptibility fixed effect covariates

// Random Effects
int<lower=0> n_rab;                 // Number of abundance random effects
int<lower=0> n_ravl;                // Number of availability random effects
int<lower=0> n_rpst;                // Number of state perceptibility random effects
int<lower=0> n_rabs[n_rab];         // Number of levels for each abundance random effect --- an (n_rab)-length vector
int<lower=0> n_ravls[n_ravl];       // Number of levels for each availability random effect --- an (n_ravl)-length vector
int<lower=0> n_rpsts[n_rpst];       // Number of levels for each state perceptibility random effect --- an (n_rpst)-length vector
int<lower=1> vab_id[sum(n_rabs)];   // Effect-category ID for each effect level --- this is a vector of indices.
     // length(vab_id) = total number of abundance random effect levels across all random effects = sum(n_rabs).
int<lower=1> vavl_id[sum(n_ravls)]; // Effect-category ID for each effect level --- this is a vector of indices.
int<lower=1> vpst_id[sum(n_rpsts)]; // Effect-category ID for each effect level --- this is a vector of indices.
int Zab[n_surv,n_rab];              // Matrix of random effect levels associated with each survey
int Zavl[n_surv,n_ravl];       
int Zpst[n_surv,n_rpst];       

// Counts and times
int<lower=0>  n_bysurv[n_surv];               // Counts by survey
real<lower=0> time_obs[n_bird];               // Observed detection times
real<lower=0> dist_obs[n_bird];               // Observed distances
int<lower=1, upper=n_surv> surv_obs[n_bird];  // Survey ID associated with each bird detected
}



transformed data {
}



parameters {
// Fixed effects
real intcpt_ab;         // Abundance intercept
real intcpt_avl;        // Availability intercept
real intcpt_pst;        // State perceptibility intercept
vector<lower=0>[groups-1] intcpt_g; // Difference between hard- and easy-to-detect intercepts
vector[n_bab] bab;      // Abundance fixed effects
vector[n_bavl] bavl;    // Availability fixed effects
vector[n_bpst] bpst;    // State perceptibility fixed effects

// Random effects
vector<lower=0>[n_rab]  sigma_ab;  // Abundance random effect standard deviation(s)
vector<lower=0>[n_ravl] sigma_avl; // Availability random effect standard deviation(s)
vector<lower=0>[n_rpst] sigma_pst; // State perceptibility random effect standard deviation(s)
vector[sum(n_rabs)]  rab;          // Vector of estimated abundance random effects for all levels
vector[sum(n_ravls)] ravl;         // Vector of estimated availability random effects for all levels
vector[sum(n_rpsts)] rpst;         // Vector of estimated state perceptibility random effects for all levels

// Other parameters
vector<lower=0,upper=1>[groups-1] gamma; // Mixing parameter. gamma = hard-to-detect, (1-gamma) = easy-to-detect 
}



transformed parameters {
vector[n_surv] log_lambda;                     // log(expected survey abundance)
vector[n_surv] log_phi_vec;                    // log(fixed and rdm effects on availability... ignores intercepts, which vary by group)
matrix[n_surv,groups] log_phi_avail;           // log(survey availability rate parameter) for each group
vector[n_surv] log_pst_sig;                    // log(distance parameter)
matrix<lower=0>[n_surv,groups] phi_avail;      // Survey availability rate parameter
matrix<lower=0, upper=1>[n_surv,groups] pdet;  // Probability of detection by survey x group
vector<lower=0, upper=1>[n_surv] p_surv;       // Probability of detection by survey across groups
vector<lower=0, upper=1>[groups] mixing;       // Mixing vector, sums to 1
row_vector[groups] int_vec;                    // Vector of detection intercepts

// expected abundance
for (s in 1:n_surv) log_lambda[s] = intcpt_ab;                        // Intercepts
if(n_bab > 0) log_lambda = log_lambda + Xab*bab;                      // Fixed effects
if(n_rab > 0)
  for (s in 1:n_surv)
    for (i in 1:n_rab) log_lambda[s] = log_lambda[s] + rab[Zab[s,i]]; // Random effects

// Availability rate
if(n_bavl > 0)
  log_phi_vec = Xavl*bavl;                                                     // Fixed effects
else 
  for(s in 1:n_surv) log_phi_vec[s] = 0;
if(n_ravl > 0) 
  for (s in 1:n_surv)
    for (i in 1:n_ravl) log_phi_vec[s] = log_phi_vec[s] + ravl[Zavl[s,i]];     // Random effects
int_vec[1] = intcpt_avl;
if (groups==2) int_vec[2] = intcpt_avl + intcpt_g[1];
log_phi_avail = rep_matrix(log_phi_vec, groups) + rep_matrix(int_vec, n_surv); // Intercepts
phi_avail = exp(log_phi_avail);

// State perceptibility parameter
for (s in 1:n_surv) log_pst_sig[s] = intcpt_pst;                           // Intercepts
if(n_bpst > 0) log_pst_sig = log_pst_sig + Xpst*bpst;                      // Fixed effects
if(n_rpst > 0)
  for (s in 1:n_surv)
    for (i in 1:n_rpst) log_pst_sig[s] = log_pst_sig[s] + rpst[Zpst[s,i]]; // Random effects

// Detection probability by survey-group
for(s in 1:n_surv){
  for(g in 1:groups){                                                          
    pdet[s,g] = exp(2*log_pst_sig[s]) * exponential_cdf(tau, phi_avail[s,g]) * (1 - g_r(1, exp(log_pst_sig[s])));
  }
}

// Mixing and detection probability by survey
if(groups==1) mixing[1]=1;
if(groups==2){
  mixing[1] = gamma[1];
  mixing[2] = 1-gamma[1];
}
p_surv = pdet * mixing;
}



model {
// Fixed effect priors
if (n_bab > 0)  bab  ~ normal(0,1); 
if (n_bavl > 0) bavl ~ normal(0,1); 
if (n_bpst > 0) bpst ~ normal(0,1);
intcpt_ab  ~ normal(ab_prior, ab_prior_sd);  
intcpt_avl ~ normal(avl_prior, avl_prior_sd);
if (groups==2) intcpt_g ~ exponential(0.5);
intcpt_pst ~ normal(0.35,1);      // Prior Pr(det|avail): median = 0.79, 95% CI = (0.040, 0.995)

// MC estimates (n=10000) of prior on Pr(detection): median = 0.47, 95% CI = (0.02, 0.98)

// Random effect priors
if(n_rab > 0)  for (i in 1:sum(n_rabs))  rab[i]  ~ normal(0,sigma_ab[vab_id[i]]);
if(n_ravl > 0) for (i in 1:sum(n_ravls)) ravl[i] ~ normal(0,sigma_avl[vavl_id[i]]);
if(n_rpst > 0) for (i in 1:sum(n_rpsts)) rpst[i] ~ normal(0,sigma_pst[vpst_id[i]]);
if(n_rab > 0)  sigma_ab  ~ student_t(3,0,1);
if(n_ravl > 0) sigma_avl ~ student_t(3,0,1);
if(n_rpst > 0) sigma_pst ~ student_t(3,0,1);

// Other parameter priors
if (groups==2) gamma ~ beta(1,1);

// Data models
n_bysurv ~ poisson_log(log(p_surv) + log_lambda);
// Increment the log-likelihood for the joint distribution of (r,t):
for(b in 1:n_bird){            // log(f_{r,t}|det) = log(f_r|det) + log(f_t|det) in this model
  // log(f_r|det) = log(2r/w^2 * g(r)) = constant + log(r) + log(g(r))
  target += log(dist_obs[b]) + 
            log(g_r(dist_obs[b], exp(log_pst_sig[surv_obs[b]])*maxdist));  
}
if (groups==2) for(b in 1:n_bird){
  // f_t|det = [gamma*dexp(t,phi1) + (1-gamma)*dexp(t,phi2)] / pdet
  target += log(mixing[1] * exp(exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],1])) +
                mixing[2] * exp(exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],2]))) -
            log(p_surv[surv_obs[b]]);
}
if (groups==1) for(b in 1:n_bird){
  target += exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],1]) -
            log(p_surv[surv_obs[b]]);
}
}



generated quantities {
real<lower=0> unobserved[n_surv];   // Sampled uncounted birds (sampling in this section refers to posterior predictive)
real<lower=0> totN[n_surv];         // Sampled total birds
vector<upper=0>[n_surv] log_lik;    // LogLikelhood for LOO / WAIC calculations.  Survey is the observational unit.
real<lower=0,upper=1> pavl_surv[n_surv]; // Survey-level availability
real<lower=0,upper=1> pst_surv[n_surv];  // Survey-level state perceptibility
real<lower=0> ppcpt_surv[n_surv];   // Survey-level pr(det|available)... removing upper bound b/c of step-approximation used for p_surv
real<lower=0,upper=1> global_p;     // Sampled overall detection probability
real<lower=0> mean_pavl;            // Average over all birds (totN) of availability
real<lower=0> mean_pst;             // Average over all birds of state perceptibility
real<lower=0> mean_pp;              // Average over all birds of total perceptibility = pdet / availability

// Surprisingly, neither TotN nor global_p can be vectorized:
for (s in 1:n_surv) {
  unobserved[s] = poisson_rng(exp(log_lambda[s] + log(1-p_surv[s]))); // Uncounted ~ Po(\lambda*(1-p))
  totN[s]       = n_bysurv[s] + unobserved[s];
  log_lik[s]    = poisson_log_lpmf(n_bysurv[s] | log(p_surv[s]) + log_lambda[s]);
  if(groups == 1) pavl_surv[s]  = mixing[1]*exponential_cdf(tau, phi_avail[s,1]);
  if(groups == 2) pavl_surv[s]  = mixing[1]*exponential_cdf(tau, phi_avail[s,1]) +
                                  mixing[2]*exponential_cdf(tau, phi_avail[s,2]);
  pst_surv[s]   = p_state(exp(log_pst_sig[s]));
  ppcpt_surv[s] = p_surv[s] / pavl_surv[s];
}
// Calculate log-likelihood for (n, r, t | \lambda, \varphi, \sigma)... same as above
if(groups==2) for(b in 1:n_bird){
  log_lik[surv_obs[b]] = log_lik[surv_obs[b]] + log(2) - 2*log(maxdist) + 
              log(dist_obs[b]) + log(g_r(dist_obs[b], exp(log_pst_sig[surv_obs[b]])*maxdist)) +
              log(mixing[1] * exp(exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],1])) +
                  mixing[2] * exp(exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],2]))) -
              log(p_surv[surv_obs[b]]);
}
if(groups==1) for(b in 1:n_bird){
  log_lik[surv_obs[b]] = log_lik[surv_obs[b]] + log(2) - 2*log(maxdist) + 
              log(dist_obs[b]) + log(g_r(dist_obs[b], exp(log_pst_sig[surv_obs[b]])*maxdist)) +
              exponential_lpdf(time_obs[b] | phi_avail[surv_obs[b],1]) -
              log(p_surv[surv_obs[b]]);
}

global_p  = sum(n_bysurv) / sum(totN);
mean_pavl = dot_product(pavl_surv, totN) / sum(totN);
mean_pst  = dot_product(pst_surv, totN) / sum(totN);
mean_pp   = dot_product(ppcpt_surv, totN) / sum(totN);
}

