###############################################################################
#
#  l1-Regularized Generalized Least Squares
#  Kaveh Salehzadeh-Nobari & Alex Gibberd
#
#  Sign recovery analysis accompanying:
#    "l1-Regularized Generalized Least Squares"
#
#  This script computes and plots the empirical sign recovery probability
#  P[sign(hat{beta}) = sign(beta_0)] across Monte Carlo replications,
#  as described in Section 4.2 of the paper (see also Figure 4).
#
#  Sign recovery requires exact agreement for all p components:
#    sign(hat{beta}_i) = sign(beta_{0,i})  for all i = 1, ..., p
#
#  The script loads simulation output from HPC_estimation_error_simulation.R
#  and compares recovery rates for the three estimators:
#    1. LASSO       — standard l1-penalised OLS (Eq. 2)
#    2. GLS-LASSO   — oracle GLS with known rho (Section 3.1)
#    3. FGLS-LASSO  — feasible GLS with estimated rho (Section 3.2)
#
#  Recovery probabilities are plotted against sample size n for each
#  dimension p in {128, 256, 512}.
#
###############################################################################

library(latex2exp)
library(ggplot2)

###############################################################################
# Setup: Set the working directory to wherever the .RData files are stored.
# Modify this path as needed for your local machine.
###############################################################################

data_dir <- "~/Data/"

###############################################################################
# loadRData: Helper to load an .RData file and return the object inside it,
# regardless of the variable name used when saving.
###############################################################################

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load simulation results for a given rho.
# Change the filename to match the desired scenario:
#   l2ErrorTbl_0.RData    -> rho = 0    (white noise, Section 4.1)
#   l2ErrorTbl_0_5.RData  -> rho = 0.5
#   l2ErrorTbl_0.9.RData  -> rho = 0.9
#   l2ErrorTbl_0.99.RData -> rho = 0.99 (near unit root)

data <- loadRData(file.path(data_dir, "l2ErrorTbl_0_5.RData"))

###############################################################################
# Compute sign recovery for each (n, p) combination across B replications.
#
# For each MC replication i and each (n, p) index j, we check whether
#   sign(round(hat{beta})) == sign(round(beta_0))
# holds for ALL p components simultaneously (Section 4.2).
#
# This is done separately for each estimator:
#   results_Hat  -> LASSO
#   results_GLS  -> GLS-LASSO (oracle)
#   results_FGLS -> FGLS-LASSO (feasible)
###############################################################################

results_Hat <- array(data = NA, c(dim(data$BetaSim)[2],dim(data$BetaSim)[3]))
results_GLS <- array(data = NA, c(dim(data$BetaSim)[2],dim(data$BetaSim)[3]))
results_FGLS <- array(data = NA, c(dim(data$BetaSim)[2],dim(data$BetaSim)[3]))


for(i in 1:dim(data$BetaSim)[3]){
  
  for(j in 1:dim(data$BetaSim)[2]){
    
    z_Hat <- all(sign(round(data$BetaSim[,j,i], digits = 0)) == sign(round(data$BetaHatSim[2:length(data$BetaHatSim[,1,1]),j,i], digits = 0)), na.rm = TRUE)
    z_GLS <- all(sign(round(data$BetaSim[,j,i], digits = 0)) == sign(round(data$BetaHatGLSSim[2:length(data$BetaHatSim[,1,1]),j,i], digits = 0)), na.rm = TRUE)
    z_FGLS <- all(sign(round(data$BetaSim[,j,i], digits = 0)) == sign(round(data$BetaHatFGLSSim[2:length(data$BetaHatSim[,1,1]),j,i], digits = 0)), na.rm = TRUE)
    
    results_Hat[j,i] <- z_Hat
    results_GLS[j,i] <- z_GLS
    results_FGLS[j,i] <- z_FGLS 
    
  }
  
  RecoveryProb <- list(results_Hat, results_GLS, results_FGLS)
  names(RecoveryProb) <- c("Hat","GLS","FGLS")
}

###############################################################################
# Compute mean recovery probability hat{P}[sign(hat{beta}) = sign(beta_0)]
# for each sample size n, separately for each dimension p = 128, 256, 512.
#
# The simulation stores results for all three dimensions sequentially:
#   rows 1-10:  p = 128  (n = 50, 100, ..., 500)
#   rows 11-20: p = 256
#   rows 21-30: p = 512
###############################################################################

mean_128 <- rep(NA, times = 10)
mean_256 <- rep(NA, times = 10)
mean_512 <- rep(NA, times = 10)

mean_128_GLS <- rep(NA, times = 10)
mean_256_GLS <- rep(NA, times = 10)
mean_512_GLS <- rep(NA, times = 10)

mean_128_FGLS <- rep(NA, times = 10)
mean_256_FGLS <- rep(NA, times = 10)
mean_512_FGLS <- rep(NA, times = 10)


for(j in 1:dim(RecoveryProb$Hat)[1]){
  
  if(j<= 10){
    
    # p = 128 block
    
    if(j == 1){counter = 1}
    
    mean_128[counter] <- mean(RecoveryProb$Hat[j, ]) 
    mean_128_GLS[counter] <- mean(RecoveryProb$GLS[j, ]) 
    mean_128_FGLS[counter] <- mean(RecoveryProb$FGLS[j, ]) 
    
    
    counter = counter + 1
    
  }else if(j>10 && j<=20){
    
    # p = 256 block
    
    if(j == 11){counter = 1}
    
    mean_256[counter] <- mean(RecoveryProb$Hat[j, ]) 
    mean_256_GLS[counter] <- mean(RecoveryProb$GLS[j, ]) 
    mean_256_FGLS[counter] <- mean(RecoveryProb$FGLS[j, ]) 
    
    counter = counter + 1
    
  }else if(j>20){
    
    # p = 512 block
    
    if(j == 21){counter = 1}
    
    mean_512[counter] <- mean(RecoveryProb$Hat[j, ]) 
    mean_512_GLS[counter] <- mean(RecoveryProb$GLS[j, ]) 
    mean_512_FGLS[counter] <- mean(RecoveryProb$FGLS[j, ]) 
    
    counter = counter + 1
    
  }
  
}

###############################################################################
# Plot sign recovery probability vs sample size (Figure 4 of the paper).
#
# Each panel shows one estimator (FGLS, GLS, LASSO), with separate curves
# for p = 128 (solid black), p = 256 (dashed green), p = 512 (twodash red).
#
# As discussed in Section 4.2, the GLS and FGLS estimators achieve higher
# recovery probabilities than the unadjusted LASSO when rho is large,
# reflecting the theoretical advantages of the whitening procedure.
###############################################################################

x_axis = seq(50,500,by = 50)

par(mfrow=c(1,3))

plot(x = x_axis, y = mean_128_FGLS, type="l", ylim=c(0,1),
     ylab = "P(Recovery)",xlab = "n", font.main=1, main="FGLS")
lines(x = x_axis, y = mean_256_FGLS, lty="dashed",col="green")
lines(x = x_axis, y = mean_512_FGLS, lty="twodash",col="red")

plot(x = x_axis, y = mean_128_GLS, type="l", ylim=c(0,1),
     ylab = "P(Recovery)" ,xlab = "n", font.main=1, main="GLS")
lines(x = x_axis, y = mean_256_GLS, lty="dashed",col="green")
lines(x = x_axis, y = mean_512_GLS, lty="twodash",col="red")


plot(x = x_axis, y = mean_128, type="l", ylim=c(0,1),
     ylab = "P(Recovery)", xlab = "n", font.main=1, main="Lasso")
lines(x = x_axis, y = mean_256, lty="dashed",col="green")
lines(x = x_axis, y = mean_512, lty="twodash",col="red")