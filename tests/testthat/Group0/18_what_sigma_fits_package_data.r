#!/usr/bin/env Rscript
# Key test: What sigma value best fits PACKAGE-generated data?

cat("\n\n-------------- What Sigma Fits Package Data? ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)
library(ggdmcPrior)

Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

# TRUE parameters used for simulation
true_means <- c(0.5, 0.2)
true_sigma <- 0.2
true_guess <- c(.1, .2, .3, .4, .5)
true_slip <- c(.2, .4, .6, .8, .9)

true_p_vector <- c(
  guess1 = true_guess[1], guess2 = true_guess[2], guess3 = true_guess[3],
  guess4 = true_guess[4], guess5 = true_guess[5],
  mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
  slip1 = true_slip[1], slip2 = true_slip[2], slip3 = true_slip[3],
  slip4 = true_slip[4], slip5 = true_slip[5]
)

cat("Data generated with TRUE parameters:\n")
cat("  sigma =", true_sigma, "\n")
cat("  means =", true_means, "\n\n")

# Build model
model <- BuildModel(
  p_map = list(
    guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
    mean1 = "1", mean2 = "1", sigma = "1",
    slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
  ),
  factors = NULL,
  constants = NULL,
  match_map = NULL,
  accumulators = Q,
  type = "cdm",
  verbose = FALSE
)

sub_model <- setCDM(model,
  q_matrix = model@cdm_info$q_matrix,
  profile_probability = model@cdm_info$profile_probability,
  rule = "DINA",
  use_mvn = TRUE
)

# Simulate data with PACKAGE
cat("Simulating data with PACKAGE...\n")
set.seed(123)
N <- 10000

dat <- simulate(sub_model,
  nsim = N,
  parameter_vector = true_p_vector,
  nschool = 1,
  debug = FALSE
)

cat("Data simulated. N =", N, "\n\n")

# Build DMI
sub_dmis <- BuildDMI(dat$responses, model,
  q_matrix = model@cdm_info$q_matrix,
  profile_probability = model@cdm_info$profile_probability,
  rule = "DINA",
  use_mvn = TRUE
)

# Load helper
source("Group0/00_mle_helper.r")

# Test fine grid around true value
cat(rep("=", 70), "\n", sep = "")
cat("Testing likelihood at different sigma values\n")
cat("(All other parameters FIXED at TRUE values)\n")
cat(rep("=", 70), "\n\n", sep = "")

sigma_grid <- seq(0.05, 0.50, by = 0.01)

results <- data.frame(
  sigma = sigma_grid,
  loglik = NA
)

for (i in seq_along(sigma_grid)) {
  sigma <- sigma_grid[i]

  test_params <- true_p_vector
  test_params["sigma"] <- sigma

  loglik <- sll_from_p(test_params, sub_dmis[[1]])
  results$loglik[i] <- loglik

  if (sigma %in% c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)) {
    cat(sprintf("sigma = %.2f: loglik = %.2f\n", sigma, loglik))
  }
}

# Find maximum
max_idx <- which.max(results$loglik)
ml_sigma <- results$sigma[max_idx]

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("RESULT\n")
cat(rep("=", 70), "\n\n", sep = "")

cat(sprintf("TRUE sigma used to generate data: %.2f\n", true_sigma))
cat(sprintf("ML sigma that best fits the data: %.2f\n", ml_sigma))
cat(sprintf("Difference: %.2f\n\n", ml_sigma - true_sigma))

if (abs(ml_sigma - true_sigma) < 0.02) {
  cat("SUCCESS: Likelihood correctly recovers true sigma!\n")
  cat("The simulation and likelihood are BOTH correct.\n")
} else if (ml_sigma > true_sigma + 0.05) {
  cat("*** BUG: Likelihood prefers HIGHER sigma! ***\n")
  cat("Data appears to have been generated with HIGHER correlation.\n")
} else if (ml_sigma < true_sigma - 0.05) {
  cat("*** BUG: Likelihood prefers LOWER sigma! ***\n")
  cat("Data appears to have been generated with LOWER correlation.\n")
} else {
  cat("Marginal: Small difference, within sampling variability.\n")
}

# Create plot
cat("\nCreating likelihood profile plot...\n")
pdf("/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat/Group0/package_data_sigma_profile.pdf",
  width = 10, height = 6
)

plot(results$sigma, results$loglik,
  type = "l", lwd = 2, col = "blue",
  xlab = "Sigma (correlation parameter)",
  ylab = "Log-Likelihood",
  main = sprintf("Likelihood Profile for Package-Generated Data (N=%d)", N),
  las = 1
)

grid()
abline(v = true_sigma, col = "red", lwd = 2, lty = 2)
abline(v = ml_sigma, col = "darkgreen", lwd = 2, lty = 3)

legend("bottomleft",
  legend = c(
    sprintf("True sigma = %.2f", true_sigma),
    sprintf("ML sigma = %.2f", ml_sigma)
  ),
  col = c("red", "darkgreen"),
  lty = c(2, 3),
  lwd = 2
)

dev.off()

cat("Plot saved.\n\n")
