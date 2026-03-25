###############################################################################
#
#  l1-Regularized Generalized Least Squares
#  Kaveh Salehzadeh-Nobari & Alex Gibberd
#
#  Monte Carlo simulation code accompanying:
#    "l1-Regularized Generalized Least Squares"
#
#  This script implements the simulation study described in Section 4
#  ("Experimental Results") of the paper. It compares estimation error
#  (in l1, l2, and l-infinity norms) across three estimators:
#
#    1. LASSO        — standard l1-penalised OLS  (Eq. 2)
#    2. GLS-LASSO    — l1-penalised GLS with known rho (Section 3.1, Eq. 3)
#    3. FGLS-LASSO   — feasible variant with estimated rho (Section 3.2)
#
#  The data-generating process follows the linear model in Eq. (1):
#
#         y = X * beta + L * u
#
#  where e = L*u has an AR(1) covariance structure Gamma = sigma_u^2 * L * L'
#  (see Section 2.2, Eq. 6). The whitening step pre-multiplies by R = L^{-1}
#  so that the transformed errors are i.i.d. (Section 3.1).
#
#  Simulation parameters (Section 4.1):
#    - Sample sizes:  n = 50, 100, ..., 500
#    - Dimensions:    p in {128, 256, 512}
#    - Sparsity:      s = p/10  (10% non-zero coefficients)
#    - AR(1) params:  rho in {0, 0.5, 0.9, 0.99}
#    - MC reps:       B = 1000
#    - lambda tuned via 2-fold temporal CV (first/second half split)
#
###############################################################################

library(glmnet)
library(ggplot2)
library(tidyverse)
library(caret)
library(pracma)
library(MASS)

## Replication of Fig 2 of Wainwright's Sharp Thresholds for HD and Noisy Sparsity Recovery 

# 'n_min' is the minimum and initial number of observations for the simulation.
# 'n_max' is the maximum and final number of observations for the simulation.
# 'p_rng' specifies the set of dimension sizes to consider. This argument can be scalar, or a vector of dimensions.
# 'sss' is the sparsity percentage.
# 'step' is the increment size for sequencing from 'n_min' to 'n_max'. e.g., 'n_min = 5', 'n_max = 20', 'step = 10', -> n={5,15}.	
# 'iter' is the number of iterations for the simulations. Set to '1000' unless specified otherwise. 
# 'errStd' is the "noise level", corresponding to the standard deviation of the error terms. Set to '1' unless specifised otherwise.
# 'Intcpt' is the intercept for the DGP. Set to 'FALSE' unless specified otherwise.


###############################################################################
# SparseDGP: Sparse Data-Generating Process
#
# Generates data from the model y = X*beta + e  (Eq. 1 of the paper).
#
# - The coefficient vector beta has s = round(sss * p) non-zero entries,
#   each drawn as +/- betaMin with equal probability (Section 4.1).
# - The design matrix X has i.i.d. rows X_t ~ N(0, Sigma_x), where
#   Sigma_x = (1 - theta)*I + theta*11' is an equi-correlation structure.
#   In the main simulations theta = 0, giving i.i.d. N(0,I) rows.
# - The error vector e is drawn from N(0, Omega) where:
#     * If rho != 0: Omega = (1/(1-rho^2)) * Toeplitz(1, rho, rho^2, ...)
#       corresponding to a stationary AR(1) process (Section 2.2, Eq. 6).
#     * If rho = 0:  Omega = sigma^2 * I  (white noise baseline).
#
# Returns: list with Y_Gen, X_Gen, sparseCoefs, eps, Omega
###############################################################################

SparseDGP <- function(n, p, s, theta = 0, mu_x = rep(0, times = p), betaMin = 0.5, fixDesign = TRUE, noC = TRUE, equiProb = TRUE, rho = 0, epStd = 1, equicor = FALSE, EV = FALSE, BIV = FALSE, GARCH = FALSE, KarlRohe = TRUE, scaleShape = 0){
  
  
  if(KarlRohe == TRUE){
    
    n <- 200
    p <- 1000
    
  }
  
  # Generate a random vector of coefficients.
  # When equiProb = TRUE, non-zero entries are +/- betaMin with equal
  # probability (Rademacher signs), matching the simulation setup in
  # Section 4.1.
  
  if(equiProb == TRUE && noC == FALSE && KarlRohe == FALSE){
    
    w <- rbinom(n = p + 1, size = 1, prob = 0.5)
    coefs <- w * betaMin - (1 - w) * (betaMin)
    
  }else if(equiProb == FALSE && noC == FALSE && KarlRohe == FALSE){
    
    coefs <- rnorm(p + 1, mean = 0, sd = 1)
    
  }else if(equiProb == TRUE && noC == TRUE && KarlRohe == FALSE){
    
    w <- rbinom(n = p, size = 1, prob = 0.5)
    coefs <- w * betaMin - (1 - w) * (betaMin)
    
  }else if(equiProb == FALSE && noC == TRUE && KarlRohe == FALSE){
    
    coefs <- rnorm(p, mean = 0, sd = 1)
    
  }else if(KarlRohe == TRUE || (KarlRohe == TRUE && rho != 0)){
    
    coefs <- c(rep(20, times = 10), rep(0, times = p - 10)) 			
    
  }
  
  
  
  if(KarlRohe == FALSE){
    
    
    # Impose sparsity: randomly zero out (1-s)% of coefficients so that
    # ||beta||_0 = round(s * p), as described in Section 4.1.
    
    ss <- round((1 - s) * p)
    
    toReplace <- sample(p, size = ss)
    sCoefs <- replace(coefs, list = toReplace, values = 0)
    
  }else if(KarlRohe == TRUE){
    
    sCoefs <- coefs
    
  }
  
  
  if(KarlRohe == TRUE || fixDesign == TRUE){
    
    set.seed(42)
    
  }
  
  
  # Generate the random design matrix X with rows X_t ~ N(mu_x, Sigma_x).
  # Sigma_x = (1 - theta)*I + theta*11' is the equi-correlation design
  # (theta = 0 gives i.i.d. sub-Gaussian rows as in Assumption 1 of
  # the paper).
  
  sigma_x <- (1 - theta) * diag(p) + theta * matrix(data = 1, ncol = p, nrow = p)
  x <- mvrnorm(n, mu = mu_x, Sigma = sigma_x)
  
  # If 'noC=FALSE', add a column of ones as the first column of matrix 'X'.
  
  if(noC == FALSE){
    
    X <- cbind(rep(1, times = n), x)
    
  }else{
    
    X <- x
    
  }
  
  
  if(KarlRohe == TRUE || fixDesign == TRUE){
    
    rm(.Random.seed, envir = globalenv())
    
  } 
  
  # Generate error vector e ~ N(0, Omega).
  #
  # The main case of interest is the AR(1) error structure (Section 2.2):
  #   e_t = rho * e_{t-1} + u_t,   u_t ~ N(0, 1)
  # which gives the Toeplitz covariance (Eq. 6):
  #   Omega = (1/(1 - rho^2)) * Toeplitz(1, rho, rho^2, ...)
  #
  # The Cholesky factor L of Omega satisfies e = L*u (Eq. 1), and
  # R = L^{-1} is the whitening matrix used in the GLS step (Section 3.1).
  #
  # Additional error structures (equicorrelated, explosive variance,
  # break-in-variance, GARCH) are included for robustness checks but
  # are not the focus of the main paper results.
  
  if(equicor == TRUE && KarlRohe == FALSE){
    
    Omega <- (1 - rho) * diag(n) + rho * matrix(data = 1, ncol = n, nrow = n)
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
    
  }else if(EV == TRUE && KarlRohe == FALSE){
    
    
    Sqnc <- exp(0.5 * (1:n))
    sqnc <- Sqnc^2
    Omega <- diag(sqnc)
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
    
  }else if(BIV == TRUE && KarlRohe == FALSE){
    
    Sigma <- diag(n)
    t <- round(n/2)
    Sigma[t, t] <- 1000
    Omega <- Sigma
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
  }else if(GARCH == TRUE && KarlRohe == FALSE){
    
    e <- rep(0, times = n)
    e[1] <- rnorm(1)
    Sigma <- rep(0, times = n)
    Sigma[1] <- 1
    
    for(t in 2:n){
      
      Sigma[t] <- 0.00037 + 0.0888 * e[t-1]^2 + 0.9024 * Sigma[t-1]
      e[t] <- rnorm(1, mean = 0, sd = sqrt(Sigma[t]))
      
    }
    
    Omega <- diag(Sigma)
    
  }else if(rho == 0 && epStd == 1 && KarlRohe == FALSE){
    
    # White noise case (rho = 0): baseline where LASSO and GLS should
    # perform comparably (Section 4.1).
    
    sqnc <- rho^seq(0, n, by = 1)
    Sigma <- toeplitz(sqnc[1: n])
    Omega <- (1 / (1 - rho^2)) * Sigma
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
  }else if(rho == 0 && epStd != 1 && KarlRohe == FALSE){
    
    Omega <- diag(epStd^2, n, n)
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
  }else if(rho != 0 || (KarlRohe == TRUE && rho != 0)){
    
    # AR(1) errors with parameter rho (Section 2.2, Eq. 6).
    # Omega = (1/(1-rho^2)) * Toeplitz(1, rho, rho^2, ..., rho^{n-1})
    
    rho <- round(rho, digits = 2)
    sqnc <- rho^seq(0, n, by = 1)
    Sigma <- toeplitz(sqnc[1: n])
    Omega <- (1 / (1 - rho^2)) * Sigma
    e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
    
  }else if(KarlRohe == TRUE && rho == 0){
    
    if(scaleShape == 0){
      
      Omega <- diag(1, n, n)
      e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
      
    }else{
      
      sqnc <- rgamma(n, shape = scaleShape, rate = scaleShape)
      sqncSQ <- sqnc^2 
      Omega <- diag(sqncSQ, n, n)
      e <- mvrnorm(n = 1, mu = rep(0, times = n), Sigma = Omega)
      
    }
  }
  
  # Generate response: y = X*beta + e  (Eq. 1)
  
  y <- X%*%sCoefs+e
  
  # In presence of an intercept, - i.e. 'noC=FALSE', the first column of X constituing of ones is removed for estimation. 
  
  if(noC == FALSE){
    
    X_NI <- X[, -1]
    
  }else{
    
    X_NI <- X
  }
  
  # The outcomes are stored in dataSparse. For instance, to obtain the generated Y vector, write 'dataSparse$Y_Gen'.
  
  data_list <- list(y, X_NI, sCoefs, e, Omega)
  names(data_list) <- c("Y_Gen", "X_Gen", "sparseCoefs", "eps", "Omega")
  return(data_list)
}


###############################################################################
# CV.gls: 2-fold temporal cross-validation for lambda selection
#
# As described in Section 4.1, the regularisation parameter lambda_n is
# tuned using 2-fold CV where the folds are the first and second halves
# of the data, preserving temporal ordering. This is important given the
# AR(1) error structure — random k-fold CV would break the dependence
# structure and produce misleading tuning results.
#
# The function searches over a log-spaced grid of 100 lambda values
# and returns the one minimising average out-of-fold RMSE.
###############################################################################

CV.gls <- function(y, x){
  
  
  lambdas <- logspace(x1 = -5, x2 = 2, n = 100)
  l <- length(lambdas)
  n <- length(y)
  w <- round(n/2)
  
  # Split into first and second halves (temporal folds)
  
  y_1 <- head(y, n = w)
  y_2 <- tail(y, n = n - w)
  
  x_1 <- head(x, n = w)
  x_2 <- tail(x, n = n - w)
  
  
  rmse <- matrix(nrow = l)
  
  for(i in 1:l){
    
    model1 <- glmnet(x = x_1, y = y_1, lambda = lambdas[i])
    model2 <- glmnet(x = x_2, y = y_2, lambda = lambdas[i])
    
    prediction1 <- predict(object = model1, newx = x_2)
    prediction2 <- predict(object = model2, newx = x_1)
    
    RMSE1 <- sqrt(mean((y_2 - prediction1)^2))
    RMSE2 <- sqrt(mean((y_1 - prediction2)^2))
    
    rmse[i] <- mean(RMSE1, RMSE2)
  }
  
  index <- which.min(rmse)
  df<- list(lambdas,rmse,lambdas[index],index)
  names(df) <- c("lambdas","rmse","min_lambda","index")
  return(df)
}


###############################################################################
# Rotation: Whitening transformation via Cholesky factorisation
#
# Implements the GLS rotation described in Section 3.1.
# Given the AR(1) parameter (true or estimated), this function:
#   1. Constructs the Toeplitz covariance Omega (Eq. 6)
#   2. Computes its Cholesky factor L (lower-triangular, such that
#      Omega = L * L')
#   3. Pre-multiplies y and X by R = L^{-1} to obtain the whitened
#      data (tilde{y}, tilde{X}) used in the GLS-LASSO (Eq. 3):
#
#        tilde{beta} = argmin [ (1/2n)||R(y - X*beta)||_2^2
#                                + tilde{lambda}_n ||beta||_1 ]
#
# When parma = true rho, this gives the oracle GLS-LASSO.
# When parma = rhoHat (estimated from first-stage residuals), this
# gives the feasible GLS-LASSO (FGLS, Section 3.2).
###############################################################################

Rotation <- function(Y, X, parma){
  
  sqnc <- parma^seq(0, length(Y), by=1)
  Sigma <- toeplitz(sqnc[1:length(Y)])
  Omega <- (1/(1-parma^2))*Sigma
  
  #Cholesky factorization of the covariance matrix
  U <- chol(Omega)
  L <- t(U)
  
  #Transform the variables: tilde{y} = R*y, tilde{X} = R*X  (Section 3.1)
  yStar <- solve(L)%*%Y
  xStar <- solve(L)%*%X
  
  
  #Create a list of the variables
  data_list <- list(yStar, xStar)
  names(data_list) <- c("yStar", "xStar")
  rotatedVars <<- data_list
  
} 


###############################################################################
# LassoEstError: Main Monte Carlo simulation loop
#
# Implements the full simulation study of Section 4. For each combination
# of (n, p, rho), the function runs B = iter Monte Carlo replications of:
#
#   Step 1 — LASSO (Eq. 2):
#     Fit LASSO on raw (y, X) with lambda chosen by temporal CV.
#
#   Step 2 — GLS-LASSO (oracle, Section 3.1):
#     Whiten (y, X) using the TRUE rho, then fit LASSO on (tilde{y},
#     tilde{X}).
#
#   Step 3 — FGLS-LASSO (feasible, Section 3.2):
#     a) Compute residuals hat{e} = y - X*hat{beta} from Step 1.
#     b) Estimate rhoHat via OLS on the AR(1) representation of the
#        residuals (Eq. 10, using ar.ols).
#     c) Whiten (y, X) using rhoHat, then fit LASSO on the rotated data.
#
#   Estimation errors ||hat{Delta}||_q for q = 1, 2, infinity are recorded
#   (Section 4.2), along with 2.5% and 97.5% quantiles across MC reps.
#
# Arguments:
#   n_min, n_max, step — define the sample size grid (50 to 500 by 50)
#   p_rng              — vector of dimensions, e.g. c(128, 256, 512)
#   sss                — sparsity fraction (0.1 = 10% non-zero)
#   Rho                — true AR(1) parameter for the error process
#   iter               — number of Monte Carlo replications (default 1000)
###############################################################################

LassoEstError <- function(n_min, n_max, p_rng, sss, step, Rho, iter = 1000, errStd = 1, intcpt = FALSE, ...){
  
  m <- length(p_rng)
  
  B <- iter
  
  # 'matRow' is the maximum number of samples considered.
  
  matRow <- length(seq(n_min, n_max, by = step))
  
  
  
  
  # 'n_sizes' lists the samples from 'n_min' to 'n_max' sequenced by increments of size 'step'.
  
  n_sizes <- seq(n_min, n_max, by = step)
  
  
  # Initialise storage matrices for estimation errors.
  # Each matrix is matRow x (m*3), where the 3 columns per dimension
  # store: [2.5% quantile, mean, 97.5% quantile] of ||hat{Delta}||_q.
  
  counter_2 <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_1 <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_infty <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  
  counter_2GLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_1GLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_inftyGLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  
  counter_2FGLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_1FGLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  counter_inftyFGLS <- matrix(data = NA, nrow = matRow, ncol = m * 3)
  
  BetaSim <- array(data = NA, c(max(p_rng), m * matRow, B))
  BetaHatGLSSim <- array(data = NA, c(max(p_rng) + 1, m * matRow, B))
  BetaHatFGLSSim <- array(data = NA, c(max(p_rng) + 1, m * matRow, B))
  BetaHatSim <- array(data = NA, c(max(p_rng) + 1, m * matRow, B))
  
  LambdaFGLSSim <- array(data = NA, c(matRow, m, B))
  LambdaGLSSim <- array(data = NA, c(matRow, m, B))
  LambdaSim <- array(data = NA, c(matRow, m, B))
  
  rhoSim <- array(data = NA, c(matRow, m, B))
  
  errorsFGLS <- array(data = NA, c(n_max, m * matRow, B))
  errorsGLS <- array(data = NA, c(n_max, m * matRow, B))
  errors <- array(data = NA, c(n_max, m * matRow, B))
  
  # XSim <- array(data = NA, c(n_max, max(p_rng) + 1, m * matRow, B))
  SNR <- array(data = NA, c(matRow, m, B))
  SNRGLS <- array(data = NA, c(matRow, m, B))
  SNRFGLS <- array(data = NA, c(matRow, m, B))
  
  
  FTT <- 1
  
  errCounter <- 0
  
  # Outer loop: iterate over dimension sizes p in p_rng
  
  for(t in 1:m){
    
    
    # Inner loop: iterate over sample sizes n in n_sizes
    
    
    for(j in 1:matRow){
      
      
      # Frequent messages that show the progress of the code, while it is running.
      
      
      if(j %% 5 == 0){
        
        cat("Obs", j, "column", t,"\n")
        
      }
      
      
      # For each (n, p) pair, run B Monte Carlo replications
      
      
      l2Error <- matrix(data = 0, nrow = B, ncol = 1)
      l1Error <- matrix(data = 0, nrow = B, ncol = 1)
      lInftyError <- matrix(data = 0, nrow = B, ncol = 1)
      
      l2ErrorGLS <- matrix(data = 0, nrow = B, ncol = 1)
      l1ErrorGLS <- matrix(data = 0, nrow = B, ncol = 1)
      lInftyErrorGLS <- matrix(data = 0, nrow = B, ncol = 1)
      
      l2ErrorFGLS <- matrix(data = 0, nrow = B, ncol = 1)
      l1ErrorFGLS <- matrix(data = 0, nrow = B, ncol = 1)
      lInftyErrorFGLS <- matrix(data = 0, nrow = B, ncol = 1)
      
      
      for(l in 1:B){
        
        d <- p_rng[t]
        n_sim <- n_sizes[j]
        
        if(intcpt == FALSE){
          
          # Generate sparse DGP: y = X*beta + e  (Eq. 1)
          
          DGP_list <- SparseDGP(n = n_sim, p = d, s = sss, rho = Rho, epStd = errStd, ...)
          
        }else{
          
          # Generate sparse DGP with intercept
          
          DGP_list <- SparseDGP(n = n_sim, p = d, s = sss, rho = Rho, epStd = errStd, noC = FALSE, ...)
          
        }
        
        
        # Number of non-zero elements of the random coefficient vector.					
        
        k_p <- round(sss * d)
        
        
        
        # Extract the variables
        
        Y <- DGP_list$Y_Gen
        X <- DGP_list$X_Gen
        Beta <- as.matrix(DGP_list$sparseCoefs)
        
        Bete <- Beta
        if(length(Bete) < max(p_rng)){Bete<-matrix(data=c(Bete, rep(NA,max(p_rng)-length(Bete))),nrow=max(p_rng),ncol=1)}
        BetaSim[, j + errCounter, l] <- Bete
        
        
        ########################################################
        # STEP 1: Standard LASSO on raw data (Eq. 2)
        # hat{beta} = argmin [ (1/2n)||y - X*beta||_2^2
        #                       + lambda_n ||beta||_1 ]
        ########################################################
        
        cvfit <- CV.gls(x = X, y = Y) 
        model <- glmnet(x = X, y = Y, lambda = cvfit$min_lambda)
        Betahat <- coef(model)
        
        
        LambdaSim[j, t, l] <- cvfit$min_lambda 
        
        Betehat <- matrix(Betahat)
        if(length(Betehat) < max(p_rng) + 1){Betehat<-matrix(data=c(Betehat,rep(NA,(max(p_rng)+1)-length(Betehat))),nrow=max(p_rng)+1,ncol=1)}
        BetaHatSim[, j + errCounter, l] <- Betehat
        
        cvfit1 <- CV.gls(x = X, y = Y) 
        model1 <- glmnet(x = X, y = Y, lambda = cvfit1$min_lambda)
        Betahat1 <- coef(model1)
        
        LambdaSim[j, t, l] <- cvfit1$min_lambda 
        X_intercept <- matrix(data=c(rep(1, times = n_sim), X), nrow = n_sim, ncol = d + 1) 
        X_c <- X_intercept
        
        
        errers  <- as.matrix(Y - X_c %*% Betahat1)
        
        XSim <- norm(as.matrix(X_c %*% Betahat1), type = "2")
        errSim <- norm(errers, type = "2")
        SNR[j, t, l] <- log10(XSim / errSim) 
        
        if(length(errers)<n_max){errers<-matrix(data=c(errers,rep(NA,n_max-length(errers))),nrow=n_max,ncol=1)}
        errors[, j + errCounter, l] <- errers
        
        
        ########################################################
        # STEP 3a (FGLS): Estimate AR(1) parameter from
        # first-stage LASSO residuals (Section 3.2, Eq. 10).
        #
        # hat{e}_t = y_t - hat{beta}' * x_t
        # hat{phi} estimated via OLS on the AR(1) representation.
        ########################################################
        
        errHat <- as.numeric(Y - X_c %*% Betahat1)
        
        RhoHat <- ar.ols(errHat, aic = FALSE, order.max = 1, na.action = na.fail, demean = TRUE)
        
        rhoHat <- as.numeric(RhoHat[2])
        
        rhoSim[j, t, l] <- rhoHat 
        
        ########################################################
        # STEP 3b (FGLS): Whiten data using estimated rhoHat
        # and fit second-stage LASSO (Section 3.2).
        #
        # tilde{y} = hat{R} * y,  tilde{X} = hat{R} * X
        # where hat{R} = hat{L}^{-1} is constructed from rhoHat.
        ########################################################
        
        objRotHat <- Rotation(Y = Y, X = X, parma = rhoHat)
        
        xRotHat <- objRotHat$xStar
        yRotHat <- objRotHat$yStar
        
        cvfit2 <- CV.gls(x = xRotHat, y = yRotHat) 
        model2 <- glmnet(x = xRotHat, y = yRotHat, lambda = cvfit2$min_lambda)
        Betahat2 <- coef(model2)
        
        ########################################################
        # STEP 2: Oracle GLS-LASSO using TRUE rho (Section 3.1).
        #
        # Whiten using the known AR(1) parameter Rho, then fit
        # LASSO on the rotated data. This serves as the oracle
        # benchmark for the FGLS estimator.
        ########################################################
        
        objRot <- Rotation(Y = Y, X = X, parma = Rho)
        
        xRot <- objRot$xStar
        yRot <- objRot$yStar
        
        cvfitGLS <- CV.gls(x = xRot, y = yRot) 
        modelGLS <- glmnet(x = xRot, y = yRot, lambda = cvfitGLS$min_lambda)
        Betahat1GLS <- coef(modelGLS)
        Betehat1GLS <- matrix(Betahat1GLS)
        
        if(length(Betehat1GLS) < max(p_rng) + 1){Betehat1GLS <- matrix(data = c(Betehat1GLS, rep(NA, (max(p_rng)+1)-length(Betehat1GLS))), nrow = max(p_rng) + 1, ncol = 1)}
        BetaHatGLSSim[, j + errCounter, l] <- Betehat1GLS
        
        LambdaGLSSim[j,t,l] <- cvfitGLS$min_lambda
        LambdaFGLSSim[j,t,l] <- cvfit2$min_lambda 
        
        Betehat2 <- matrix(Betahat2)
        if(length(Betehat2) < max(p_rng) + 1){Betehat2 <- matrix(data = c(Betehat2, rep(NA, (max(p_rng)+1)-length(Betehat2))), nrow = max(p_rng) + 1, ncol = 1)}
        BetaHatFGLSSim[, j + errCounter, l] <- Betehat2
        
        xRotHat_c <- cbind(rep(1, times = nrow(xRotHat)), xRotHat)
        xRot_c <- cbind(rep(1, times = nrow(xRot)), xRot)
        
        errersFGLS  <- as.matrix(yRotHat - xRotHat_c %*% Betahat2)
        errersGLS <- as.matrix(yRot - xRot_c %*% Betahat1GLS)
        
        XSimFGLS <- norm(as.matrix(xRotHat_c %*% Betahat2), type = "2")
        errSimFGLS <- norm(errersFGLS, type = "2")
        SNRFGLS[j, t, l] <- log10(XSimFGLS / errSimFGLS) 
        
        XSimGLS <- norm(as.matrix(xRot_c %*% Betahat1GLS), type = "2")
        errSimGLS <- norm(errersGLS, type = "2")
        SNRGLS[j, t, l] <- log10(XSimGLS / errSimGLS)
        
        if(length(errersFGLS)<n_max){errersFGLS<-matrix(data=c(errersFGLS,rep(NA,n_max-length(errersFGLS))),nrow=n_max,ncol=1)}
        errorsFGLS[, j + errCounter, l] <- errersFGLS
        
        if(length(errersGLS)<n_max){errersGLS<-matrix(data=c(errersGLS,rep(NA,n_max-length(errersGLS))),nrow=n_max,ncol=1)}
        errorsGLS[, j + errCounter, l] <- errersGLS
        
        ########################################################
        # Compute estimation errors ||hat{Delta}||_q (Section 4.2)
        # where hat{Delta} = hat{beta} - beta_0
        # for q = 2 (l2), q = 1 (l1), q = infinity (l-inf)
        ########################################################
        
        
        l2Error[l] <- norm(Beta - Betahat1[2: length(Betahat1)], type = "2")
        l1Error[l] <- norm(Beta - Betahat1[2: length(Betahat1)], type = "1")
        lInftyError[l] <- norm(Beta - Betahat1[2: length(Betahat1)], type = "I")
        
        l2ErrorGLS[l] <- norm(Beta - Betahat1GLS[2: length(Betahat1GLS)], type = "2")
        l1ErrorGLS[l] <- norm(Beta - Betahat1GLS[2: length(Betahat1GLS)], type = "1")
        lInftyErrorGLS[l] <- norm(Beta - Betahat1GLS[2: length(Betahat1GLS)], type = "I")
        
        l2ErrorFGLS[l] <- norm(Beta - Betahat2[2: length(Betahat2)], type = "2")
        l1ErrorFGLS[l] <- norm(Beta - Betahat2[2: length(Betahat2)], type = "1")
        lInftyErrorFGLS[l] <- norm(Beta - Betahat2[2: length(Betahat2)], type = "I")
        
      }
      
      # Store summary statistics across B replications:
      # [2.5% quantile, mean, 97.5% quantile] for each error norm
      
      counter_2[j, FTT] <- quantile(l2Error, probs = 0.025)	
      counter_2[j, FTT+1] <- mean(l2Error)
      counter_2[j, FTT+2] <- quantile(l2Error, probs = 0.975)
      
      counter_1[j, FTT] <- quantile(l1Error, probs = 0.025)	
      counter_1[j, FTT+1] <- mean(l1Error)
      counter_1[j, FTT+2] <- quantile(l1Error, probs = 0.975)
      
      counter_infty[j, FTT] <- quantile(lInftyError, probs = 0.025)	
      counter_infty[j, FTT+1] <- mean(lInftyError)
      counter_infty[j, FTT+2] <- quantile(lInftyError, probs = 0.975)
      
      counter_2GLS[j, FTT] <- quantile(l2ErrorGLS, probs = 0.025)	
      counter_2GLS[j, FTT+1] <- mean(l2ErrorGLS)
      counter_2GLS[j, FTT+2] <- quantile(l2ErrorGLS, probs = 0.975)
      
      counter_1GLS[j, FTT] <- quantile(l1ErrorGLS, probs = 0.025)	
      counter_1GLS[j, FTT+1] <- mean(l1ErrorGLS)
      counter_1GLS[j, FTT+2] <- quantile(l1ErrorGLS, probs = 0.975)
      
      counter_inftyGLS[j, FTT] <- quantile(lInftyErrorGLS, probs = 0.025)	
      counter_inftyGLS[j, FTT+1] <- mean(lInftyErrorGLS)
      counter_inftyGLS[j, FTT+2] <- quantile(lInftyErrorGLS, probs = 0.975)
      
      
      counter_2FGLS[j, FTT] <- quantile(l2ErrorFGLS, probs = 0.025)	
      counter_2FGLS[j, FTT+1] <- mean(l2ErrorFGLS)
      counter_2FGLS[j, FTT+2] <- quantile(l2ErrorFGLS, probs = 0.975)
      
      counter_1FGLS[j, FTT] <- quantile(l1ErrorFGLS, probs = 0.025)	
      counter_1FGLS[j, FTT+1] <- mean(l1ErrorFGLS)
      counter_1FGLS[j, FTT+2] <- quantile(l1ErrorFGLS, probs = 0.975)
      
      counter_inftyFGLS[j, FTT] <- quantile(lInftyErrorFGLS, probs = 0.025)	
      counter_inftyFGLS[j, FTT+1] <- mean(lInftyErrorFGLS)
      counter_inftyFGLS[j, FTT+2] <- quantile(lInftyErrorFGLS, probs = 0.975)
      
      
      
    }
    
    errCounter <- errCounter + matRow
    
    FTT <- FTT + 3
  }
  
  NormErrors <- list(counter_2FGLS, counter_1FGLS, counter_inftyFGLS,counter_2GLS, counter_1GLS, counter_inftyGLS, counter_2, counter_1, counter_infty, errors, errorsGLS ,errorsFGLS, BetaSim, BetaHatSim, BetaHatGLSSim, BetaHatFGLSSim, XSim, LambdaSim, LambdaGLSSim,LambdaFGLSSim, rhoSim, SNR, SNRGLS,SNRFGLS)
  names(NormErrors) <- c("L2FGLS", "L1FGLS", "LinftyFGLS","L2GLS", "L1GLS", "LinftyGLS","L2", "L1", "Linfty", "errors", "errorsGLS", "errorsFGLS", "BetaSim", "BetaHatSim", "BetaHatGLSSim", "BetaHatFGLSSim","XSim", "LambdaSim", "LambdaGLSSim","LambdaFGLSSim", "rhoSim", "SNR", "SNRGLS","SNRFGLS")
  return(NormErrors)
  
  
}


###############################################################################
# Run the simulation for rho = 0 (white noise baseline)
# See Section 4.1 for full parameter specification.
# Repeat with Rho = 0.5, 0.9, 0.99 for the remaining scenarios.
###############################################################################

l2ErrorTbl_0 <- LassoEstError(n_min=50, n_max=500, p_rng=c(128, 256, 512), Rho = 0, sss = 0.1, step = 50, iter = 1000, errStd = 1, intcpt = FALSE, KarlRohe = FALSE)

save(l2ErrorTbl_0, file = "l2ErrorTbl_0.RData")