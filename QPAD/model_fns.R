# *************************************************************************************
# IJ + PP model
# *************************************************************************************

fit_IJPP_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                         rarray, # distance intervals for each point count
                         tarray, # time intervals for each point count
                         
                         X1 = NULL,
                         X2 = NULL,
                         X3 = NULL, 
                         
                         tau_A_inits = NULL,
                         tau_B_inits = NULL,
                         phi_inits = NULL,
                         
                         method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X3)[1] <- "log_phi"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_inits) != ncol(X3)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X3))
    names(phi_inits) <- phi_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    tau_B <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_A_params)+1):(length(tau_A_params)+length(tau_B_params))]))
    phi <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params)+length(tau_B_params) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_B[k]
      phi_k <- phi[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * (exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-exp(-(dmax/tau_A)^2))*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2)))),
                               list(phi = phi_k, 
                                    tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
  
}

offsets_IJPP_fn <- function (Yarray = Yarray,
                             rarray, # distance intervals for each point count
                             tarray, # time intervals for each point count
                             
                             # Estimates for each survey
                             tau_A = NULL,
                             tau_B = NULL,
                             phi = NULL) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_k <- phi[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * (exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-exp(-(dmax/tau_A)^2))*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2)))),
                             list(phi = phi_k, 
                                  tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}

# *************************************************************************************
# PP model
# *************************************************************************************

fit_PP_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                       rarray, # distance intervals for each point count
                       tarray, # time intervals for each point count
                       
                       X1 = NULL,
                       X2 = NULL,
                       X3 = NULL, 
                       
                       tau_A_inits = NULL,
                       tau_B_inits = NULL,
                       phi_inits = NULL,
                       
                       method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X3)[1] <- "log_phi"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_inits) != ncol(X3)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X3))
    names(phi_inits) <- phi_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    tau_B <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_A_params)+1):(length(tau_A_params)+length(tau_B_params))]))
    phi <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params)+length(tau_B_params) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_B[k]
      phi_k <- phi[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * (1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))),
                               list(phi = phi_k, 
                                    tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
  
}

offsets_PP_fn <- function (Yarray = Yarray,
                           rarray, # distance intervals for each point count
                           tarray, # time intervals for each point count
                           
                           # Estimates for each survey
                           tau_A = NULL,
                           tau_B = NULL,
                           phi = NULL) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_k <- phi[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * (1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))),
                             list(phi = phi_k, 
                                  tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}

# *************************************************************************************
# IJ model
# *************************************************************************************

fit_IJ_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                       rarray, # distance intervals for each point count
                       tarray, # time intervals for each point count
                       
                       X1 = NULL,
                       X2 = NULL,
                       X3 = NULL, 
                       
                       tau_A_inits = NULL,
                       tau_B_inits = NULL,
                       phi_inits = NULL,
                       
                       method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X3)[1] <- "log_phi"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_inits) != ncol(X3)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X3))
    names(phi_inits) <- phi_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    tau_B <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_A_params)+1):(length(tau_A_params)+length(tau_B_params))]))
    phi <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params)+length(tau_B_params) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_B[k]
      phi_k <- phi[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * (exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax))),
                               list(phi = phi_k, 
                                    tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
}

offsets_IJ_fn <- function (Yarray = Yarray,
                           rarray, # distance intervals for each point count
                           tarray, # time intervals for each point count
                           
                           # Estimates for each survey
                           tau_A = NULL,
                           tau_B = NULL,
                           phi = NULL) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_k <- phi[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * (exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax))),
                             list(phi = phi_k, 
                                  tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}

# *************************************************************************************
# IJ + PP model 2
# *************************************************************************************

fit_IJPP2_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                          rarray, # distance intervals for each point count
                          tarray, # time intervals for each point count
                          
                          X1 = NULL,
                          X2 = NULL,
                          X3 = NULL, 
                          X4 = NULL, 
                          
                          tau_A_inits = NULL,
                          tau_B_inits = NULL,
                          phi_inits = NULL,
                          pA_inits = NULL,
                          
                          method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X3)[1] <- "log_phi"
  }
  
  if (!is.null(X4)){
    X4 <- X4[Ykeep, ]
    pA_params <- colnames(X4)
  } else {
    X4 <- matrix(1,nrow = nsurvey,ncol=1)
    pA_params <- colnames(X4)[1] <- "logit_pA"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_inits) != ncol(X3)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X3))
    names(phi_inits) <- phi_params
  }
  
  if (length(pA_inits) != ncol(X4)) pA_inits <- NULL
  if (is.null(pA_inits)) {
    pA_inits <- rep(0, ncol(X4))
    names(pA_inits) <- pA_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_inits,pA_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    tau_B <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_A_params) + 1):(length(tau_A_params) + length(tau_B_params))]))
    phi <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params) + length(tau_B_params) + 1):((length(tau_A_params) + length(tau_B_params) + length(phi_params)))]))
    pA <-   binomial("logit")$linkinv(drop(X4 %*% params[(length(tau_A_params) + length(tau_B_params) + length(phi_params) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_B[k]
      phi_k <- phi[k]
      pA_k = pA[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-pA)*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))) ),
                               list(tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    phi = phi_k, 
                                    pA = pA_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
  
}

offsets_IJPP2_fn <- function (Yarray = Yarray,
                              rarray, # distance intervals for each point count
                              tarray, # time intervals for each point count
                              
                              # Estimates for each survey
                              tau_A = NULL,
                              tau_B  = NULL,
                              phi  = NULL,
                              pA  = NULL,
                              ...) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_k <- phi[k]
    pA_k <- pA[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-pA)*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))) ),
                             list(tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  phi = phi_k, 
                                  pA = pA_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}

# *************************************************************************************
# IJ + PP model 3
# Similar to IJPP2, but place constraint on tau_B to be a proportion of tau_A
# *************************************************************************************

fit_IJPP3_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                          rarray, # distance intervals for each point count
                          tarray, # time intervals for each point count
                          
                          X1 = NULL,
                          X2 = NULL,
                          X3 = NULL, 
                          X4 = NULL, 
                          
                          tau_A_inits = NULL,
                          tau_B_inits = NULL,
                          phi_inits = NULL,
                          pA_inits = NULL,
                          
                          method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_params <- colnames(X3)[1] <- "log_phi"
  }
  
  if (!is.null(X4)){
    X4 <- X4[Ykeep, ]
    pA_params <- colnames(X4)
  } else {
    X4 <- matrix(1,nrow = nsurvey,ncol=1)
    pA_params <- colnames(X4)[1] <- "logit_pA"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_inits) != ncol(X3)) phi_inits <- NULL
  if (is.null(phi_inits)) {
    phi_inits <- rep(0, ncol(X3))
    names(phi_inits) <- phi_params
  }
  
  if (length(pA_inits) != ncol(X4)) pA_inits <- NULL
  if (is.null(pA_inits)) {
    pA_inits <- rep(0, ncol(X4))
    names(pA_inits) <- pA_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_inits,pA_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    multB <- binomial("logit")$linkinv(drop(X2 %*% params[(length(tau_A_params) + 1):(length(tau_A_params) + length(tau_B_params))]))
    phi <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params) + length(tau_B_params) + 1):((length(tau_A_params) + length(tau_B_params) + length(phi_params)))]))
    pA <-   binomial("logit")$linkinv(drop(X4 %*% params[(length(tau_A_params) + length(tau_B_params) + length(phi_params) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_A[k] * multB[k]
      phi_k <- phi[k]
      pA_k = pA[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-pA)*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))) ),
                               list(tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    phi = phi_k, 
                                    pA = pA_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
  
}

offsets_IJPP3_fn <- function (Yarray = Yarray,
                              rarray, # distance intervals for each point count
                              tarray, # time intervals for each point count
                              
                              # Estimates for each survey
                              tau_A = NULL,
                              tau_B  = NULL,
                              phi  = NULL,
                              pA  = NULL,
                              ...) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_k <- phi[k]
    pA_k <- pA[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi*tmax)) + (1-pA)*(1-exp(-phi*tmax*exp(-(dmax/tau_B)^2))) ),
                             list(tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  phi = phi_k, 
                                  pA = pA_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}

# *************************************************************************************
# IJ + PP model 4
# *************************************************************************************

fit_IJPP4_fn <- function (Yarray, # Array with dimensions (nsurvey x nrint x ntint)
                          rarray, # distance intervals for each point count
                          tarray, # time intervals for each point count
                          
                          X1 = NULL,
                          X2 = NULL,
                          X3 = NULL, 
                          X4 = NULL, 
                          X5 = NULL, 
                          
                          tau_A_inits = NULL,
                          tau_B_inits = NULL,
                          phi_A_inits = NULL,
                          phi_B_inits = NULL,
                          pA_inits = NULL,
                          
                          method = "Nelder-Mead", ...) {
  
  logdmultinom <- function (x, size, prob) {lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))}
  
  # ----------------------------
  # Only conduct analysis on point counts with non-zero total counts
  # ----------------------------
  
  Ysum <- apply(Yarray,1,sum,na.rm = TRUE)
  Ykeep <- which(Ysum > 0)
  if (length(Ykeep) != length(Ysum)){
    Yarray <- Yarray[Ykeep, , ] %>% array(.,dim=c(length(Ykeep),dim(Yarray)[2],dim(Yarray)[3]))
    rarray<- rarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[2]))
    tarray<- tarray[Ykeep, ] %>% array(.,dim=c(length(Ykeep),ncol = dim(Yarray)[3]))
    Ysum <- Ysum[Ykeep]
  }
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Format parameters and check they are named via column names in design matrix
  if (!is.null(X1)){
    X1 <- X1[Ykeep, ]
    tau_A_params <- colnames(X1)
  } else {
    X1 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_A_params <- colnames(X1)[1] <- "log_tau_A"
  }
  
  if (!is.null(X2)){
    X2 <- X2[Ykeep, ]
    tau_B_params <- colnames(X2)
  } else {
    X2 <- matrix(1,nrow = nsurvey,ncol=1)
    tau_B_params <- colnames(X2)[1] <- "log_tau_B"
  }
  
  if (!is.null(X3)){
    X3 <- X3[Ykeep, ]
    phi_A_params <- colnames(X3)
  } else {
    X3 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_A_params <- colnames(X3)[1] <- "log_phi_A"
  }
  
  if (!is.null(X4)){
    X4 <- X4[Ykeep, ]
    phi_B_params <- colnames(X4)
  } else {
    X4 <- matrix(1,nrow = nsurvey,ncol=1)
    phi_B_params <- colnames(X4)[1] <- "log_phi_B"
  }
  
  if (!is.null(X5)){
    X5 <- X5[Ykeep, ]
    pA_params <- colnames(X5)
  } else {
    X5 <- matrix(1,nrow = nsurvey,ncol=1)
    pA_params <- colnames(X5)[1] <- "logit_pA"
  }
  
  # Initial values
  if (length(tau_A_inits) != ncol(X1)) tau_A_inits <- NULL
  if (is.null(tau_A_inits)) {
    tau_A_inits <- rep(0, ncol(X1))
    names(tau_A_inits) <- tau_A_params
  }
  
  if (length(tau_B_inits) != ncol(X2)) tau_B_inits <- NULL
  if (is.null(tau_B_inits)) {
    tau_B_inits <- rep(0, ncol(X2))
    names(tau_B_inits) <- tau_B_params
  }
  
  if (length(phi_A_inits) != ncol(X3)) phi_A_inits <- NULL
  if (is.null(phi_A_inits)) {
    phi_A_inits <- rep(0, ncol(X3))
    names(phi_A_inits) <- phi_params
  }
  
  if (length(phi_B_inits) != ncol(X4)) phi_B_inits <- NULL
  if (is.null(phi_B_inits)) {
    phi_B_inits <- rep(0, ncol(X4))
    names(phi_B_inits) <- phi_params
  }
  
  if (length(pA_inits) != ncol(X5)) pA_inits <- NULL
  if (is.null(pA_inits)) {
    pA_inits <- rep(0, ncol(X5))
    names(pA_inits) <- pA_params
  }
  
  inits <- c(tau_A_inits,tau_B_inits,phi_A_inits,phi_B_inits,pA_inits)
  
  # Function to calculate multinomial cell probabilities for each point count
  nll.fun <- function(params) {
    
    tau_A <- poisson("log")$linkinv(drop(X1 %*% params[1:length(tau_A_params)]))
    tau_B <- poisson("log")$linkinv(drop(X2 %*% params[(length(tau_A_params) + 1):(length(tau_A_params) + length(tau_B_params))]))
    phi_A <- poisson("log")$linkinv(drop(X3 %*% params[(length(tau_A_params) + length(tau_B_params) + 1):((length(tau_A_params) + length(tau_B_params) + length(phi_A_params)))]))
    phi_B <- poisson("log")$linkinv(drop(X4 %*% params[(length(tau_A_params) + length(tau_B_params) + length(phi_A_params) + 1):((length(tau_A_params) + length(tau_B_params) + length(phi_A_params) + length(phi_B_params)))]))
    pA <-   binomial("logit")$linkinv(drop(X5 %*% params[((length(tau_A_params) + length(tau_B_params) + length(phi_A_params) + length(phi_B_params)) + 1):length(params)]))
    
    nll <- rep(0,nsurvey)
    
    for (k in 1:nsurvey){
      
      tau_A_k <- tau_A[k]
      tau_B_k <- tau_B[k]
      phi_A_k <- phi_A[k]
      phi_B_k <- phi_B[k]
      pA_k = pA[k]
      
      # Calculate CDF and p
      f_d = function(dmax){
        
        integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi_A*tmax)) + (1-pA)*(1-exp(-phi_B*tmax*exp(-(dmax/tau_B)^2))) ),
                               list(tau_A = tau_A_k,
                                    tau_B = tau_B_k,
                                    phi_A = phi_A_k, 
                                    phi_B = phi_B_k, 
                                    pA = pA_k,
                                    tmax = tmax))
        eval(integrand)
      }
      
      # Calculate CDF
      Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
      
      CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
      
      for (j in 1:ntint[k]){
        
        tmax = tarray[k,j] # How many minutes have elapsed so far?
        
        for (i in 1:nrint[k]){
          
          CDF_binned[i,j] = integrate(f_d,lower=0,
                                      upper = rarray[k,i], 
                                      subdivisions = 1000)$value
        }
      }
      
      # Difference across distance bins
      tmp1 = CDF_binned
      if (nrow(tmp1)>1){
        for (i in 2:nrint[k]){
          tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
        }
      }
      
      # Difference across time bins
      p_matrix = tmp1
      if (ncol(p_matrix)>1){
        for (j in 2:ntint[k]){
          p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
        }
      }
      
      # This p_matrix gives us the expected total number of birds detected during the point count
      # if Density = 1, given particular values of phi and tau
      p_matrix
      
      # Normalize the p_matrix to yield the multinomial cell probabilities
      p_matrix = p_matrix/sum(p_matrix)
      
      # Calculate the multinomial log likelihood for this point count
      nll[k] <- logdmultinom(Y, Ysum[k], p_matrix)
      
    } # close loop on k
    
    nll <- -sum(nll)
    
    if (nll %in% c(NA, NaN, Inf, -Inf)) nlimit[2] else nll
    
  }
  
  nlimit <- c(.Machine$double.xmin, .Machine$double.xmax)^(1/3)
  
  res <- optim(inits, nll.fun, method = method, hessian = TRUE)
  
  rval <- list(convergence = res$convergence,
               coefficients = res$par,  
               loglik = -res$value)
  
  rval$results <- res
  rval
  
}

offsets_IJPP4_fn <- function (Yarray = Yarray,
                              rarray, # distance intervals for each point count
                              tarray, # time intervals for each point count
                              
                              # Estimates for each survey
                              tau_A = NULL,
                              tau_B  = NULL,
                              phi_A  = NULL,
                              phi_B  = NULL,
                              pA  = NULL,
                              ...) {
  
  nsurvey <- dim(Yarray)[1] # Number of surveys
  
  if (nsurvey == 1){
    nrint = length(na.omit(rarray))
    ntint = length(na.omit(tarray))
  } else{
    nrint <- apply(rarray,1,function(x)length(na.omit(x))) # Number of distance bins for each point count
    ntint <- apply(tarray,1,function(x)length(na.omit(x))) # Number of time bins for each point count
  }
  
  # Multinomial cell probabilities for each point count
  parray <- Yarray*NA
  offset_vec <- rep(NA,nsurvey)
  
  for (k in 1:nsurvey){
    
    tau_A_k <- tau_A[k]
    tau_B_k <- tau_B[k]
    phi_A_k <- phi_A[k]
    phi_B_k <- phi_B[k]
    pA_k <- pA[k]
    
    # Calculate CDF and p
    f_d = function(dmax){
      
      integrand = substitute(2*pi*dmax * ( pA*exp(-(dmax/tau_A)^2)*(1-exp(-phi_A*tmax)) + (1-pA)*(1-exp(-phi_B*tmax*exp(-(dmax/tau_B)^2))) ),
                             list(tau_A = tau_A_k,
                                  tau_B = tau_B_k,
                                  phi_A = phi_A_k, 
                                  phi_B = phi_B_k, 
                                  pA = pA_k,
                                  tmax = tmax))
      eval(integrand)
    }
    
    # Calculate CDF
    Y <- Yarray[k,1:nrint[k],1:ntint[k]] # Data for this survey
    
    CDF_binned <- matrix(NA,nrow=nrint[k],ncol=ntint[k])
    
    for (j in 1:ntint[k]){
      
      tmax = tarray[k,j] # How many minutes have elapsed so far?
      
      for (i in 1:nrint[k]){
        
        CDF_binned[i,j] = integrate(f_d,lower=0,
                                    upper = rarray[k,i], 
                                    subdivisions = 1000)$value
      }
    }
    
    # Difference across distance bins
    tmp1 = CDF_binned
    if (nrow(tmp1)>1){
      for (i in 2:nrint[k]){
        tmp1[i,] <- CDF_binned[i,] - CDF_binned[i-1,]
      }
    }
    
    # Difference across time bins
    p_matrix = tmp1
    if (ncol(p_matrix)>1){
      for (j in 2:ntint[k]){
        p_matrix[,j] <- tmp1[,j] - tmp1[,j-1]
      }
    }
    
    # This p_matrix gives us the expected total number of birds detected during the point count
    # if Density = 1, given particular values of phi and tau
    p_matrix
    offset_vec[k] <- sum(p_matrix)
    
    # Normalize the p_matrix to yield the multinomial cell probabilities
    p_matrix = p_matrix/sum(p_matrix)
    
    # Calculate the multinomial log likelihood for this point count
    parray[k,,] <- p_matrix
    
  } # close loop on k
  
  return(list(offset_vec = offset_vec,
              parray = parray))
}