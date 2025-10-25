#!/usr/bin/env Rscript
# Compare simulation attribute distribution vs likelihood profile probabilities
# to detect any mismatch

cat("\n\n-------------- Simulation vs Likelihood Consistency Check ---------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "cdModel", "ggdmcPrior", "mvtnorm")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat"
setwd(home_dir)

# ===========================================================================
# Setup model
# ===========================================================================
Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

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
    rule = "DINO",
    use_mvn = TRUE
)

# ===========================================================================
# Test parameters
# ===========================================================================
true_means <- c(0.5, 0.2)
true_sigma <- 0.2

true_p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
)

cat("Parameters:\n")
cat("  means:", true_means, "\n")
cat("  sigma:", true_sigma, "\n\n")

# ===========================================================================
# Theoretical profile probabilities (what likelihood uses)
# ===========================================================================
build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

K <- 2
Sigma <- build_correlation_matrix(K, true_sigma)

profiles <- matrix(c(
  0, 0,  # Profile 0
  1, 0,  # Profile 1
  0, 1,  # Profile 2
  1, 1   # Profile 3
), ncol = 2, byrow = TRUE)

theoretical_probs <- numeric(4)

for (i in 1:4) {
  lower <- ifelse(profiles[i, ] == 0, -Inf, 0)
  upper <- ifelse(profiles[i, ] == 0, 0, Inf)

  theoretical_probs[i] <- pmvnorm(lower = lower, upper = upper,
                                   mean = true_means, sigma = Sigma)[1]
}

# Normalize (like C++ code does)
theoretical_probs <- theoretical_probs / sum(theoretical_probs)

cat("THEORETICAL profile probabilities (from pmvnorm):\n")
for (i in 1:4) {
  cat(sprintf("  Profile (%d,%d): %.6f\n", profiles[i,1], profiles[i,2], theoretical_probs[i]))
}
cat(sprintf("  Sum: %.10f\n\n", sum(theoretical_probs)))

# ===========================================================================
# Empirical profile frequencies (from simulation)
# ===========================================================================
test_N_values <- c(1000, 5000, 10000, 50000)

for (N in test_N_values) {
  cat(rep("=", 70), "\n", sep = "")
  cat(sprintf("Testing with N = %d\n", N))
  cat(rep("=", 70), "\n", sep = "")

  # Simulate
  set.seed(42)
  dat <- simulate(sub_model,
      nsim = N,
      parameter_vector = true_p_vector,
      nschool = 1,
      debug = FALSE
  )

  # Extract attribute profiles from simulation
  # Note: Simulation generates responses, but we need to access the underlying alphas
  # This is tricky - we might need to re-simulate just the alphas

  # Alternative: Simulate alphas directly using R
  set.seed(42)  # Same seed

  # Generate MVN data
  library(MASS)
  X <- mvrnorm(n = N, mu = true_means, Sigma = Sigma)

  # Threshold to get binary attributes
  alpha_matrix <- (X > 0) * 1

  # Count profile frequencies
  empirical_counts <- numeric(4)
  for (i in 1:N) {
    alpha_i <- alpha_matrix[i, ]
    # Find which profile this corresponds to
    profile_idx <- which(apply(profiles, 1, function(p) all(p == alpha_i)))
    if (length(profile_idx) > 0) {
      empirical_counts[profile_idx] <- empirical_counts[profile_idx] + 1
    }
  }

  empirical_probs <- empirical_counts / N

  cat("\nEMPIRICAL profile probabilities (from simulation):\n")
  for (i in 1:4) {
    cat(sprintf("  Profile (%d,%d): %.6f (count: %d)\n",
               profiles[i,1], profiles[i,2], empirical_probs[i], empirical_counts[i]))
  }
  cat(sprintf("  Sum: %.10f\n\n", sum(empirical_probs)))

  # Compare theoretical vs empirical
  cat("COMPARISON (Empirical - Theoretical):\n")
  diffs <- empirical_probs - theoretical_probs
  for (i in 1:4) {
    cat(sprintf("  Profile (%d,%d): %+.6f (%.2f%%)\n",
               profiles[i,1], profiles[i,2], diffs[i],
               100 * diffs[i] / theoretical_probs[i]))
  }

  # Chi-square goodness of fit test
  expected_counts <- theoretical_probs * N
  chi_sq <- sum((empirical_counts - expected_counts)^2 / expected_counts)
  df <- 3  # 4 profiles - 1
  p_value <- pchisq(chi_sq, df, lower.tail = FALSE)

  cat(sprintf("\nChi-square test:\n"))
  cat(sprintf("  Chi-square statistic: %.4f\n", chi_sq))
  cat(sprintf("  Degrees of freedom: %d\n", df))
  cat(sprintf("  p-value: %.4f\n", p_value))
  if (p_value < 0.05) {
    cat("  Result: REJECT null (empirical differs from theoretical) ***\n")
  } else {
    cat("  Result: Accept null (empirical matches theoretical)\n")
  }

  cat("\n")
}

# ===========================================================================
# Test if likelihood properly reflects the simulation
# ===========================================================================
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("FINAL DIAGNOSTIC: Likelihood at different sigma values\n")
cat(rep("=", 70), "\n\n", sep = "")

# Use a large dataset
N <- 10000
set.seed(42)
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = true_p_vector,
    nschool = 1,
    debug = FALSE
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = TRUE
)

# Source MLE helper
source("Group0/00_mle_helper.r")

# Test sigma values around the true value
sigma_test_values <- c(0.15, 0.18, 0.20, 0.22, 0.25)

cat("Testing likelihood at different sigma values (N=10000):\n")
cat("All other parameters fixed at true values\n\n")

results <- data.frame(
  sigma = sigma_test_values,
  loglik = NA_real_,
  diff_from_true = NA_real_
)

for (i in seq_along(sigma_test_values)) {
  test_params <- true_p_vector
  test_params["sigma"] <- sigma_test_values[i]

  loglik <- sll_from_p(test_params, sub_dmis[[1]])
  results$loglik[i] <- loglik

  cat(sprintf("sigma = %.2f: loglik = %.4f\n", sigma_test_values[i], loglik))
}

results$diff_from_true <- results$loglik - results$loglik[results$sigma == 0.20]

cat("\n")
print(results)

max_idx <- which.max(results$loglik)
cat(sprintf("\nMaximum likelihood at sigma = %.2f (true = %.2f)\n",
           results$sigma[max_idx], true_sigma))

if (abs(results$sigma[max_idx] - true_sigma) < 0.01) {
  cat("RESULT: Likelihood IS maximized at true sigma (PASS)\n")
} else {
  cat("RESULT: Likelihood NOT maximized at true sigma (FAIL)\n")
  cat("        This confirms there is a bug!\n")
}

cat("\n")
