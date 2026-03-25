###############################################################################
#
#  l1-Regularized Generalized Least Squares
#  Kaveh Salehzadeh-Nobari & Alex Gibberd
#
#  Covariance and Cholesky factor heatmaps accompanying:
#    "l1-Regularized Generalized Least Squares"
#
#  This script generates Figure 1 of the paper, which visualises
#  the normalised autocovariance matrix Gamma/v^2 (top row) and the
#  lower-triangular Cholesky factor Psi_q/v (bottom row) for three
#  AR error structures (Section 2.2):
#
#    - AR(1)  with phi = 0.9          (persistent, slow decay)
#    - AR(2)  with phi = (1.96, -0.97) (oscillatory)
#    - AR(10) with phi_10 = 0.9       (seasonal/long-lag dependence)
#
#  The Cholesky factor L satisfies Gamma = L*L' (Eq. 1), and its
#  inverse R = L^{-1} is the whitening matrix used in the GLS step
#  (Section 3.1). The heatmaps illustrate how different AR structures
#  produce qualitatively different sparsity patterns in L, motivating
#  the discussion in Section 2.2 (see also Figure 2).
#
###############################################################################

library(fields)
library(ltsa)
library(latex2exp)
library(gridExtra)
library(colorspace)


###############################################################################
# Helper functions for image orientation
###############################################################################

# rotIm: Rotate a matrix for correct orientation in image()
rotIm <- function(X){
  n <- dim(X)[1]
  X1 <- apply(X, 2, rev)
  return(t(X1))
}

# image2: Wrapper for image() with rotated matrix
image2 <- function(X, ...){
  n <- dim(X)[1]
  X1 <- apply(X, 2, rev)
  image(1:n, 1:n, t(X1), ...)
}


###############################################################################
# Generate Figure 1: Heatmaps of Gamma/v^2 and Psi_q/v
#
# Top row:    Normalised autocovariance Gamma/v^2 for each AR model
# Bottom row: Normalised Cholesky factor L/v (= Psi_q/v in the paper)
#
# The dimension is n = 40. The colour scale is symmetric about zero
# (blue-white-red) to highlight positive and negative correlations.
###############################################################################

output_dir <- "~/Data/"

pdf(file = file.path(output_dir, "arexample2.pdf"),
    width = 6,
    height = 4.4)

par(mfrow = c(2,3), family = "serif")

# Define the three AR models (Section 2.2)
phis <- labs <- list()
phis[[1]] <- 0.9                    # AR(1): single persistent root
phis[[2]] <- c(1.96, -0.97)        # AR(2): oscillatory behaviour
phis[[3]] <- c(rep(0, 9), 0.9)     # AR(10): long-lag dependence

labs[[1]] <- "AR(1)   v="
labs[[2]] <- "AR(2)   v="
labs[[3]] <- "AR(10)   v="

cols <- two.colors(n = 256, start = "blue", end = "red", middle = "white",
                   alpha = 1.0)

# Row i=1: Gamma/v^2 (autocovariance)
# Row i=2: L/v (Cholesky factor)
for (i in 1:2){
  for(j in 1:3){
    
    n <- 40
    phi <- phis[[j]]
    lab <- labs[[j]]
    
    # Compute theoretical autocovariance function and Toeplitz matrix
    # Gamma = Toeplitz(gamma(0), gamma(1), ..., gamma(n-1))
    # (Section 2.2, Eq. 6 for the AR(1) special case)
    gamma3 <- tacvfARMA(phi = phi, theta = numeric(0), maxLag = n-1, sigma2 = 1)
    Gamma <- toeplitz(gamma3)
    v <- sqrt(gamma3[1])
    lab <- paste0(lab, round(v, 2))
    
    # Cholesky decomposition: Gamma = L * L'
    # R = L^{-1} is the whitening matrix (Section 3.1)
    L <- t(chol(Gamma))
    R <- solve(L)
    
    # Adjust margins for panel layout
    if(j == 1){
      par(mai = c(0.4, 0.45, 0.3, 0.0))
    }else if(j == 2){
      par(mai = c(0.4, 0.2, 0.3, 0.2))
    }else{
      par(mai = c(0.4, 0.0, 0.3, 0.6))
    }
    
    if(i == 1){
      # Top row: normalised autocovariance Gamma/v^2
      if(j == 1){
        ylab <- TeX("$ \\Gamma/v^2$")
      }else{
        ylab <- ""
      }
      if(j == 3){
        # Rightmost panel includes colour bar
        image.plot(rotIm(Gamma/(v^2)), ylab = "", xlab = "", main = lab,
                   font.main = 1, axes = FALSE, zlim = c(-1,1), col = cols)
      }else{
        image(rotIm(Gamma/(v^2)), ylab = "", xlab = "", main = lab,
              font.main = 1, axes = FALSE, zlim = c(-1,1), col = cols)
      }
      title(ylab = ylab, line = 1, cex.lab = 1)
      
    }else{
      # Bottom row: normalised Cholesky factor L/v (= Psi_q/v)
      if(j == 1){
        ylab <- TeX("$ \\Psi_{q}/v^2$")
      }else{
        ylab <- ""
      }
      if(j == 3){
        image.plot(rotIm(L/v), main = "", ylab = "", xlab = "",
                   font.main = 1, axes = FALSE, zlim = c(-1,1), col = cols)
      }else{
        image(rotIm(L/v), main = "", ylab = "", xlab = "",
              font.main = 1, axes = FALSE, zlim = c(-1,1), col = cols)
      }
      title(ylab = ylab, line = 1, cex.lab = 1)
    }
  }
}

dev.off()