
cmulti.fit.joint <- function()(Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                               rarray, # distance intervals for each point count
                               tarray, # time intervals for each point count
                               X1,     # Design matrix for distance model
                               X2,     # Design matrix for removal model
                               maxdistint = 2) 

logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}

# ------------------------------------
# Calculate maximum distance for integration
# ------------------------------------
max_r <- apply(rarray,1,max,na.rm = TRUE)
max_r[max_r == Inf] <- maxdistint

# ------------------------------------
# Drop point counts with no observations
# ------------------------------------
Ysum <- apply(Yarray,1,sum)
nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)

# ----------------------------------------
# Rework this next section as necessary
# ----------------------------------------

pifun <- function(r, tau) 1 - exp(-(r/tau)^2)

nll.fun <- function(param) {
  
  # Calculate cumulative prob
  if (is.null(X)){
    CP <- pifun(D, poisson("log")$linkinv(param))
  } else{
    CP <-pifun(D, poisson("log")$linkinv(drop(X %*% param)))
  }
  
  # Calculate cell probabilities
  P <- CP - cbind(0, CP[, -k, drop = FALSE])
  Psum <- rowSums(P, na.rm = TRUE)
  PPsum <- P/Psum
  
  
  nll <- -sum(sapply(1:n, function(i) {
    logdmultinom(Y[i, Yok[i, ]], Ysum[i], PPsum[i, 
                                                Yok[i, ]])
  }))
  if (nll %in% c(NA, NaN, Inf, -Inf)) 
    nlimit[2]
  else nll
}

res <- optim(inits, nll.fun, method = method, hessian = TRUE, 
             ...)
rval <- list(coefficients = res$par, vcov = try(.solvenear(res$hessian)), 
             loglik = -res$value)


                                                                                                phi, c) 1 - c * exp(-t * phi))