###############################################################################
#
#  l1-Regularized Generalized Least Squares
#  Kaveh Salehzadeh-Nobari & Alex Gibberd
#
#  Paper figure generation accompanying:
#    "l1-Regularized Generalized Least Squares"
#
#  This script produces the estimation error figures from Section 4.2:
#
#    - l2_512.pdf     (Figure 3 top):    ||hat{Delta}||_2     vs n
#    - linfty_512.pdf (Figure 3 bottom): ||hat{Delta}||_infty vs n
#
#  Layout: 3 columns (GLS, FGLS, LASSO) x 3 rows (rho = 0, 0.9, 0.99),
#  each panel showing mean error (solid) with 95% empirical confidence
#  bands (dashed), for p = 512.
#
###############################################################################

library(latex2exp)

data_dir <- "~/Data/"

###############################################################################
# 1. Load simulation results
###############################################################################

rho_files <- c("l2ErrorTbl_0.RData",
               "l2ErrorTbl_0.9.RData",
               "l2ErrorTbl_0.99.RData")

rho_labels <- c("0", "0.9", "0.99")

data_list <- list()
for (i in seq_along(rho_files)){
  env <- new.env()
  nm <- load(file.path(data_dir, rho_files[i]), env)
  data_list[[i]] <- env[[nm]]
}


###############################################################################
# 2. Plotting helper
#
# For p = 512 the summary columns are 7 (2.5%), 8 (mean), 9 (97.5%).
# Each panel plots one estimator at one rho: solid mean, dashed bands.
###############################################################################

x_axis <- seq(50, 500, by = 50)

# Column indices for p = 512 (third dimension block)
col_lower <- 7
col_mean  <- 8
col_upper <- 9

plot_panel <- function(mat, ylab_text, title_text, ylim_range){
  plot(x_axis, mat[, col_mean], type = "l", lwd = 1.5,
       ylim = ylim_range, xlab = "n", ylab = ylab_text,
       main = title_text, font.main = 1)
  lines(x_axis, mat[, col_upper], lty = 2, lwd = 1)
  lines(x_axis, mat[, col_lower], lty = 2, lwd = 1)
}


###############################################################################
# 3. Figure 3 top: ||hat{Delta}||_2 for p = 512 (l2_512.pdf)
###############################################################################

pdf(file = file.path(data_dir, "l2_512.pdf"), width = 8, height = 8)

par(mfrow = c(3, 3), mar = c(4, 4.5, 2, 1), family = "serif")

for (i in seq_along(data_list)){
  
  df <- data_list[[i]]
  rho <- rho_labels[i]
  ylab <- bquote(group("||", hat(Delta), "||")[2] ~ " (" * rho == .(rho) * ")")
  
  # Common y-axis range across all three estimators for this rho
  all_vals <- c(df$L2GLS[, col_upper], df$L2FGLS[, col_upper], df$L2[, col_upper])
  ylim_range <- c(0, max(all_vals, na.rm = TRUE) * 1.05)
  
  # GLS (Section 3.1) | FGLS (Section 3.2) | LASSO (Eq. 2)
  plot_panel(df$L2GLS,  ylab, "GLS",   ylim_range)
  plot_panel(df$L2FGLS, ylab, "FGLS",  ylim_range)
  plot_panel(df$L2,     ylab, "Lasso", ylim_range)
}

dev.off()


###############################################################################
# 4. Figure 3 bottom: ||hat{Delta}||_infty for p = 512 (linfty_512.pdf)
###############################################################################

pdf(file = file.path(data_dir, "linfty_512.pdf"), width = 8, height = 8)

par(mfrow = c(3, 3), mar = c(4, 4.5, 2, 1), family = "serif")

for (i in seq_along(data_list)){
  
  df <- data_list[[i]]
  rho <- rho_labels[i]
  ylab <- bquote(group("||", hat(Delta), "||")[infinity] ~ " (" * rho == .(rho) * ")")
  
  all_vals <- c(df$LinftyGLS[, col_upper], df$LinftyFGLS[, col_upper], df$Linfty[, col_upper])
  ylim_range <- c(0, max(all_vals, na.rm = TRUE) * 1.05)
  
  plot_panel(df$LinftyGLS,  ylab, "GLS",   ylim_range)
  plot_panel(df$LinftyFGLS, ylab, "FGLS",  ylim_range)
  plot_panel(df$Linfty,     ylab, "Lasso", ylim_range)
}

dev.off()
