#!/usr/bin/env Rscript
# Trace through the likelihood calculation step by step to find the bug

cat("\n\n-------------- Likelihood Calculation Trace ---------------\n")
rm(list = ls())

library(mvtnorm)

# Manual implementation of DINO CDM likelihood
dino_item_probability <- function(alpha, Q_row, guess, slip) {
  # eta = 1 if sum(alpha * Q_row) >= 1 (at least one required skill)
  eta <- as.numeric(sum(alpha * Q_row) >= 1)
  # P(Y=1|alpha) = eta * (1-slip) + (1-eta) * guess
  p <- eta * (1 - slip) + (1 - eta) * guess
  return(p)
}

compute_dino_probabilities <- function(profiles, Q, guess, slip) {
  # Returns L x J matrix of P(Y=1|alpha_l) for each profile l and item j
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

compute_likelihood_one_student <- function(Y, P_LJ, profile_probs) {
  # Y: 1 x J response vector
  # P_LJ: L x J probability matrix
  # profile_probs: length L vector

  L <- nrow(P_LJ)
  J <- ncol(P_LJ)

  # For each profile, compute P(Y | alpha_l)
  log_lik_profiles <- numeric(L)

  for (l in 1:L) {
    log_lik_l <- 0
    for (j in 1:J) {
      if (Y[j] == 1) {
        log_lik_l <- log_lik_l + log(P_LJ[l, j])
      } else {
        log_lik_l <- log_lik_l + log(1 - P_LJ[l, j])
      }
    }
    log_lik_profiles[l] <- log_lik_l
  }

  # Log-sum-exp with weights
  max_log_lik <- max(log_lik_profiles)
  exp_shifted <- exp(log_lik_profiles - max_log_lik)
  weighted_sum <- sum(exp_shifted * profile_probs)
  log_likelihood <- max_log_lik + log(weighted_sum)

  return(list(
    log_lik_profiles = log_lik_profiles,
    log_likelihood = log_likelihood,
    likelihood = exp(log_likelihood),
    profile_contributions = exp_shifted * profile_probs / weighted_sum
  ))
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
  0, 0,  # Profile 0
  1, 0,  # Profile 1
  0, 1,  # Profile 2
  1, 1   # Profile 3
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

compute_profile_probs <- function(means, sigma) {
  Sigma <- build_correlation_matrix(2, sigma)
  probs <- numeric(4)
  for (i in 1:4) {
    lower <- ifelse(profiles[i,] == 0, -Inf, 0)
    upper <- ifelse(profiles[i,] == 0, 0, Inf)
    probs[i] <- pmvnorm(lower=lower, upper=upper, mean=means, sigma=Sigma)[1]
  }
  return(probs)
}

# Test with a simple response pattern
Y_test <- c(1, 1, 0, 1, 0)  # Example response

cat("\nTest response pattern Y =", Y_test, "\n\n")

# Compute DINO probabilities
P_LJ <- compute_dino_probabilities(profiles, Q, true_guess, true_slip)

cat("P(Y_j=1 | alpha_l) matrix (L x J):\n")
print(round(P_LJ, 4))
cat("\n")

cat("For each profile, check eta values:\n")
for (l in 1:4) {
  cat(sprintf("  Profile (%d,%d): ", profiles[l,1], profiles[l,2]))
  for (j in 1:5) {
    eta <- sum(profiles[l,] * Q[j,]) >= 1
    cat(sprintf("eta_%d=%d ", j, eta))
  }
  cat("\n")
}
cat("\n")

# Test different sigma values
sigma_values <- c(0.10, 0.15, 0.20, 0.25, 0.30)

cat(rep("=", 70), "\n", sep = "")
cat("Testing likelihood at different sigma values\n")
cat(rep("=", 70), "\n\n", sep = "")

results <- data.frame(
  sigma = sigma_values,
  sum_profile_probs = NA,
  log_lik = NA,
  likelihood = NA
)

for (i in seq_along(sigma_values)) {
  sigma <- sigma_values[i]

  # Compute profile probabilities
  profile_probs <- compute_profile_probs(true_means, sigma)

  cat(sprintf("sigma = %.2f\n", sigma))
  cat("Profile probabilities:\n")
  for (l in 1:4) {
    cat(sprintf("  P(%d,%d) = %.6f\n", profiles[l,1], profiles[l,2], profile_probs[l]))
  }
  cat(sprintf("  Sum = %.10f\n", sum(profile_probs)))

  # Compute likelihood for this response
  lik_result <- compute_likelihood_one_student(Y_test, P_LJ, profile_probs)

  cat(sprintf("  Log-likelihood = %.6f\n", lik_result$log_likelihood))
  cat(sprintf("  Likelihood = %.6f\n\n", lik_result$likelihood))

  results$sum_profile_probs[i] <- sum(profile_probs)
  results$log_lik[i] <- lik_result$log_likelihood
  results$likelihood[i] <- lik_result$likelihood

  # Show profile contributions
  cat("  Profile posterior contributions:\n")
  for (l in 1:4) {
    cat(sprintf("    Profile (%d,%d): %.6f (log_lik=%.2f, prior=%.4f)\n",
               profiles[l,1], profiles[l,2],
               lik_result$profile_contributions[l],
               lik_result$log_lik_profiles[l],
               profile_probs[l]))
  }
  cat("\n")
}

cat(rep("=", 70), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
print(results)

max_idx <- which.max(results$log_lik)
cat(sprintf("\nMaximum log-likelihood at sigma = %.2f (true = %.2f)\n",
           results$sigma[max_idx], true_sigma))

# Check if there's a monotonic trend
cat("\nLog-likelihood differences from true sigma:\n")
results$diff_from_true <- results$log_lik - results$log_lik[results$sigma == true_sigma]
print(results[, c("sigma", "log_lik", "diff_from_true")])

cat("\n")

# CRITICAL: Check derivative
cat(rep("=", 70), "\n", sep = "")
cat("DERIVATIVE ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("If likelihood is correctly specified, the derivative w.r.t. sigma\n")
cat("should be approximately zero at the true value (for large N).\n\n")

# Numerical derivative at true sigma
idx_true <- which(results$sigma == true_sigma)
if (idx_true > 1 && idx_true < nrow(results)) {
  deriv <- (results$log_lik[idx_true + 1] - results$log_lik[idx_true - 1]) /
           (results$sigma[idx_true + 1] - results$sigma[idx_true - 1])
  cat(sprintf("Numerical derivative at sigma=%.2f: %.6f\n", true_sigma, deriv))

  if (abs(deriv) < 0.01) {
    cat("Derivative is close to zero - GOOD\n")
  } else if (deriv > 0) {
    cat("Derivative is POSITIVE - likelihood still increasing!\n")
    cat("This suggests sigma should be HIGHER than true value.\n")
  } else {
    cat("Derivative is NEGATIVE - likelihood decreasing\n")
  }
}

cat("\n")
