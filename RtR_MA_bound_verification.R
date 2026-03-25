###############################################################################
#
#  l1-Regularized Generalized Least Squares
#  Kaveh Salehzadeh-Nobari & Alex Gibberd
#
#  R'R vs MA covariance bound verification accompanying:
#    "l1-Regularized Generalized Least Squares"
#
#  This script verifies the spectral norm bound from Lemma 3 (Eq. 11):
#
#    ||R'R||_2  <=  ||Sigma_{MA(pi)}||_2  <=  (2q+1) ||pi||_2^2
#
#  where R = L^{-1} is the whitening matrix (Section 3.1) and
#  Sigma_{MA(pi)} is the autocovariance of the MA(q) process induced
#  by pi = (1, -phi_1, ..., -phi_q).
#
#  The script:
#    1. Generates Figure 3 (labelled "fig:RTR approx" in the paper):
#       heatmaps comparing R'R and Sigma_{MA(pi)} (and their inverses)
#       for an AR(2) model near the stationarity boundary.
#    2. Sweeps over the AR(2) parameter space (phi_1, phi_2) to check
#       where the bound holds, computing Frobenius norms, spectral norms,
#       and the ratio of eigenvalues.
#    3. Produces contour/image plots of diagnostic functions f(phi_1, phi_2)
#       and g(phi_1, phi_2) used in the Appendix proofs.
#
###############################################################################

library(ltsa)
library(fields)
library(latex2exp)
library(gridExtra)
library(colorspace)

output_dir <- "~/Data/"


###############################################################################
# Helper: rotate matrix for correct image() orientation
###############################################################################

rotIm <- function(X){
  n <- dim(X)[1]
  X1 <- apply(X, 2, rev)
  return(t(X1))
}


###############################################################################
# 1. Sweep over AR(2) parameter space (phi_1, phi_2)
#
# For each (phi_1, phi_2) pair, compute:
#   - The MA(q) autocovariance Sigma_{MA(pi)} induced by the AR polynomial
#   - The Cholesky-based precision R'R where Gamma = LL' (Section 3.1)
#   - Spectral norm ratio ||Sigma_{MA}||_2 / ||R'R||_2  (Lemma 3)
#   - Upper bound ratio (2q+1)||pi||_2^2 / ||R'R||_2
#   - Frobenius norm ratio ||Sigma_{MA}||_F^2 / ||R'R||_F^2
#
# Points where the AR model is non-stationary are flagged via tryCatch.
###############################################################################

n <- 10
nbig <- 100
q <- 2
n1 <- 100
n2 <- 50
phis1 <- seq(-2, 2, length.out = n1)
phis2 <- seq(-1, 1, length.out = n2)

stationary <- array(1, dim = c(n1, n2))
state <- array(0, dim = c(n1, n2))
L2 <- Fnorm2 <- stateF <- stateL2 <- stateUB <- state
gammasum <- phi1 <- phi2 <- array(NA, dim = c(n1, n2))

for (i in 1:n1){
  for (j in 1:n2){
    
    # Compute MA autocovariance for this (phi_1, phi_2) pair
    phis <- c(phis1[i], phis2[j])
    ma <- tacvfARMA(phi = numeric(0), theta = phis, maxLag = n-1, sigma2 = 1)
    ma2 <- tacvfARMA(phi = numeric(0), theta = phis, maxLag = nbig-1, sigma2 = 1)
    
    # Check row-sum dominance condition (Appendix)
    r2 <- rev(ma[1:(q+1)])
    r1 <- c(-rev(phis[1:(q-1)]), r2[q+1] - phis[q]^2)
    if(sum(abs(r2)) >= sum(abs(r1))){
      state[i,j] <- 1
    }else{
      state[i,j] <- 0
    }
    
    # Attempt Cholesky factorisation; non-stationary models will fail
    res <- tryCatch({
      artest <- tacvfARMA(phi = phis, theta = numeric(0), maxLag = n-1, sigma2 = 1)
      artest2 <- tacvfARMA(phi = phis, theta = numeric(0), maxLag = nbig-1, sigma2 = 1)
      Gamma2 <- toeplitz(artest2)
      L2 <- t(chol(Gamma2))
      
      Gamma <- toeplitz(artest)
      L <- t(chol(Gamma))
      f <- 0
      
    }, error = function(err) { 
      f <- 1
      return(f)
    })
    
    if(res == 1){
      # Non-stationary: Cholesky failed
      stationary[i,j] <- 0
    }else{
      
      # Construct R'R and Sigma_{MA(pi)} (Lemma 3)
      GammaR <- toeplitz(ma)
      GammaR2 <- toeplitz(ma2)
      R <- solve(L)
      R2 <- solve(L2)
      
      RtR <- t(R) %*% R
      R2tR2 <- t(R2) %*% R2
      
      # Spectral norms and upper bound (Eq. 11)
      evd <- eigen(RtR)
      evd2 <- eigen(GammaR)
      ub <- (2*q+1) * (1 + sum(phis^2))     # (2q+1)||pi||_2^2
      
      L2[i,j] <- max(evd$values)              # ||R'R||_2
      Fnorm2[i,j] <- sum(RtR^2)               # ||R'R||_F^2
      stateF[i,j] <- sum(GammaR^2) / Fnorm2[i,j]   # Frobenius ratio
      stateL2[i,j] <- max(evd2$values) / L2[i,j]    # Spectral ratio
      stateUB[i,j] <- ub / L2[i,j]            # Upper bound ratio
      
      # Row-sum diagnostics for Appendix conditions
      gammasum[i,j] <- sum(abs(ma[-1]))
      phi1[i,j] <- phis[1]
      phi2[i,j] <- phis[2]
    }
  }
}


###############################################################################
# 2. Figure 3 ("fig:RTR approx"): R'R vs Sigma_{MA(pi)} heatmaps
#
# Compares R'R and Sigma_{MA(pi)} (and their inverses) for an AR(2)
# model with phi = (1.96, -0.97), near the stationarity boundary.
#
# Top row:    R'R, (R'R)^{-1} at n=10, (R'R)^{-1} at n=100
# Bottom row: Sigma_{MA(pi)}, Sigma_{MA}^{-1} at n=10, Sigma_{MA}^{-1} at n=100
#
# At n=10 there is visible discrepancy between R'R and Sigma_{MA},
# while at n=100 they converge (as noted in the paper caption).
###############################################################################

image2 <- function(X, x, y, ...){
  n <- dim(X)[1]
  X1 <- apply(X, 2, rev)
  image(x, y, X1, ...)
}

cols <- two.colors(n = 256, start = "blue", end = "red", middle = "white",
                   alpha = 1.0)

pdf(file = file.path(output_dir, "RtR.pdf"),
    width = 6,
    height = 3.8)

par(mfrow = c(2,3), family = "serif")

par(mai = c(0.2, 0.2, 0.4, 0.6))
t1 <- TeX("$R^TR$")
t2 <- TeX("$(R^TR)^{-1}$ , n=10")
t3 <- TeX("$(R^TR)^{-1}$ , n=100")

image.plot(rotIm(RtR), main = t1, col = cols, zlim = c(-6,6), axes = FALSE)
image.plot(rotIm(solve(RtR)), main = t2, col = cols, axes = FALSE)
image.plot(rotIm(solve(R2tR2)), main = t3, col = cols, axes = FALSE)

par(mai = c(0.2, 0.2, 0.4, 0.6))
t4 <- TeX("$ \\Sigma_{MA(\\pi)}$")
t5 <- TeX("$ \\Sigma_{MA(\\pi)}^{-1}$ , n=10")
t6 <- TeX("$ \\Sigma_{MA(\\pi)}^{-1}$ , n=100")

image.plot(rotIm(GammaR), main = t4, col = cols, zlim = c(-6,6), axes = FALSE)
image.plot(rotIm(solve(GammaR)), main = t5, col = cols, axes = FALSE)
image.plot(rotIm(solve(GammaR2)), main = t6, col = cols, axes = FALSE)

dev.off()


###############################################################################
# 3. Contour plots: diagnostic functions over (phi_1, phi_2) space
#
# These plots verify the conditions used in the Appendix proofs for
# the bound in Lemma 3 (Eq. 11). They show:
#   - f(phi_1, phi_2): row-sum dominance condition
#   - g(phi_1, phi_2): relationship between gamma sums and phi norms
#   - Spectral norm ratio ||Sigma_{MA}||_2 / ||R'R||_2
###############################################################################

val <- val2 <- val3 <- array(NA, dim = c(n1, n2))
for (i in 1:n1){
  for (j in 1:n2){
    val[i,j]  <- phis2[j]^2 + abs(phis1[i] * (1 - phis2[j])) - abs(phis1[i])
    val2[i,j] <- gammasum[i,j] - phi2[i,j]^2 - abs(phi1[i,j])
    val3[i,j] <- abs(phi1[i,j] * (phi2[i,j] - 1)) + abs(phi2[i,j]) - phi2[i,j]^2 - abs(phi1[i,j])
  }
}

cols <- two.colors(n = 256, start = "blue", end = "red", alpha = 1.0)

par(mfrow = c(1,3), family = "serif")
ylab <- TeX("$ \\phi_2$")
xlab <- TeX("$ \\phi_1$")
t1 <- TeX("$ f(\\phi_1,\\phi_2)$")
t2 <- TeX("$ g(\\phi_1,\\phi_2)$")
t3 <- TeX("$  \\|\\Sigma_{MA}\\|_2 / \\| R^TR \\|_2$")

par(mai = c(0.2, 0.2, 0.4, 0.6))
contour(phis1, phis2, val, xlab = xlab, ylab = ylab, main = t1, nlevels = 20)
contour(phis1, phis2, val2, xlab = xlab, ylab = ylab, main = t2)
image.plot(phis1, phis2, val3, main = t2, xlab = xlab, ylab = ylab, col = cols)