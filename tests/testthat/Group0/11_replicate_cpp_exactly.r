#!/usr/bin/env Rscript
# Replicate EXACTLY what the C++ code does, including the normalization bug

cat("\n\n-------------- Exact C++ Replication Test ---------------\n")
rm(list = ls())

library(mvtnorm)
library(MASS)

# C++ style normalization (WITH THE BUG)
cpp_normalize <- function(probs) {
  L <- length(probs)
  # Line 908: tmp_prob = tmp_prob / arma::accu(tmp_prob + 1e-10);
  denominator <- sum(probs + 1e-10)  # = sum(probs) + L * 1e-10
  normalized <- probs / denominator

  # Line 909-910: clamp to [1e-10, 1 - 1e-10]
  clamped <- pmax(pmin(normalized, 1.0 - 1e-10), 1e-10)

  return(clamped)
}

# Correct normalization
correct_normalize <- function(probs) {
  normalized <- probs / sum(probs)
  return(normalized)
}

# Manual DINO implementation
dino_item_probability <- function(alpha, Q_row, guess, slip) {
  eta <- as.numeric(sum(alpha * Q_row) >= 1)
  p <- eta * (1 - slip) + (1 - eta) * guess
  return(p)
}

compute_dino_probabilities <- function(profiles, Q, guess, slip) {
  L <- nrow(profiles)
  J <- nrow(Q)
  P <- matrix(0, nrow = L, ncol = J)
  for (l in 1:L) {
    for (j in 1:J) {
      P[l, j] <- dino_item_probability(profiles[l,], Q[j,], guess[j], slip[j])
    }
  }
  return(P)
}

compute_log_likelihood <- function(Y_data, P_LJ, profile_probs) {
  N <- nrow(Y_data)
  J <- ncol(Y_data)
  L <- nrow(P_LJ)

  total_log_lik <- 0

  for (i in 1:N) {
    log_lik_profiles <- numeric(L)

    for (l in 1:L) {
      log_lik_l <- 0
      for (j in 1:J) {
        if (Y_data[i,j] == 1) {
          log_lik_l <- log_lik_l + log(P_LJ[l, j])
        } else {
          log_lik_l <- log_lik_l + log(1 - P_LJ[l, j])
        }
      }
      log_lik_profiles[l] <- log_lik_l
    }

    # Log-sum-exp with profile probabilities
    max_log_lik <- max(log_lik_profiles)
    exp_shifted <- exp(log_lik_profiles - max_log_lik)
    weighted_sum <- sum(exp_shifted * profile_probs)
    log_likelihood_i <- max_log_lik + log(weighted_sum)

    total_log_lik <- total_log_lik + log_likelihood_i
  }

  return(total_log_lik)
}

# Setup
Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)

profiles <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1
), ncol = 2, byrow = TRUE)

# True parameters
true_means <- c(0.5, 0.2)
true_sigma <- 0.2
true_guess <- c(.1, .2, .3, .4, .5)
true_slip <- c(.2, .4, .6, .8, .9)

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

compute_profile_probs_raw <- function(means, sigma) {
  Sigma <- build_correlation_matrix(2, sigma)
  probs <- numeric(4)
  for (i in 1:4) {
    lower <- ifelse(profiles[i,] == 0, -Inf, 0)
    upper <- ifelse(profiles[i,] == 0, 0, Inf)
    probs[i] <- pmvnorm(lower=lower, upper=upper, mean=means, sigma=Sigma)[1]
  }
  return(probs)
}

# Simulate data
cat("Simulating data...\n")
set.seed(123)
N <- 10000

Sigma_true <- build_correlation_matrix(2, true_sigma)
X <- mvrnorm(n = N, mu = true_means, Sigma = Sigma_true)
alpha_data <- (X > 0) * 1

Y_data <- matrix(0, nrow = N, ncol = 5)
for (i in 1:N) {
  for (j in 1:5) {
    p_ij <- dino_item_probability(alpha_data[i,], Q[j,], true_guess[j], true_slip[j])
    Y_data[i, j] <- rbinom(1, 1, p_ij)
  }
}

cat("Data simulated.\n\n")

# Test both normalization methods
cat(rep("=", 70), "\n", sep = "")
cat("Comparing CORRECT vs C++ (BUGGY) normalization\n")
cat(rep("=", 70), "\n\n", sep = "")

sigma_grid <- seq(0.05, 0.50, by = 0.05)
P_LJ <- compute_dino_probabilities(profiles, Q, true_guess, true_slip)

results_correct <- data.frame(
  sigma = sigma_grid,
  log_lik = NA
)

results_cpp <- data.frame(
  sigma = sigma_grid,
  log_lik = NA
)

for (i in seq_along(sigma_grid)) {
  sigma <- sigma_grid[i]

  # Get raw probabilities from pmvnorm
  raw_probs <- compute_profile_probs_raw(true_means, sigma)

  # CORRECT normalization
  profile_probs_correct <- correct_normalize(raw_probs)
  log_lik_correct <- compute_log_likelihood(Y_data, P_LJ, profile_probs_correct)
  results_correct$log_lik[i] <- log_lik_correct

  # C++ (BUGGY) normalization
  profile_probs_cpp <- cpp_normalize(raw_probs)
  log_lik_cpp <- compute_log_likelihood(Y_data, P_LJ, profile_probs_cpp)
  results_cpp$log_lik[i] <- log_lik_cpp
}

# Find maxima
max_idx_correct <- which.max(results_correct$log_lik)
max_idx_cpp <- which.max(results_cpp$log_lik)

cat("CORRECT normalization:\n")
cat(sprintf("  Maximum at sigma = %.2f (true = %.2f)\n",
           results_correct$sigma[max_idx_correct], true_sigma))
cat(sprintf("  Log-lik at true = %.2f\n", results_correct$log_lik[results_correct$sigma == true_sigma]))
cat(sprintf("  Log-lik at max = %.2f\n\n", max(results_correct$log_lik)))

cat("C++ (BUGGY) normalization:\n")
cat(sprintf("  Maximum at sigma = %.2f (true = %.2f)\n",
           results_cpp$sigma[max_idx_cpp], true_sigma))
cat(sprintf("  Log-lik at true = %.2f\n", results_cpp$log_lik[results_cpp$sigma == true_sigma]))
cat(sprintf("  Log-lik at max = %.2f\n\n", max(results_cpp$log_lik)))

# Compare side by side
comparison <- data.frame(
  sigma = sigma_grid,
  loglik_correct = results_correct$log_lik,
  loglik_cpp = results_cpp$log_lik,
  difference = results_cpp$log_lik - results_correct$log_lik
)

cat("Comparison table:\n")
print(comparison)

cat("\n")
cat("Maximum absolute difference:", max(abs(comparison$difference)), "\n")

if (max_idx_correct != max_idx_cpp) {
  cat("\n*** DIFFERENT MAXIMA! The normalization bug DOES affect the result! ***\n")
} else {
  cat("\nSame maximum - the tiny normalization error doesn't affect the result.\n")
  cat("The bug must be elsewhere.\n")
}

cat("\n")
