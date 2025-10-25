#!/usr/bin/env Rscript
# Test likelihood with actual simulated data to see if bias persists

cat("\n\n-------------- Real Data Likelihood Test ---------------\n")
rm(list = ls())

library(mvtnorm)
library(MASS)

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

# Simulate data with TRUE parameters
cat("Simulating data with true parameters...\n")
cat("  N = 10000\n")
cat("  true_means =", true_means, "\n")
cat("  true_sigma =", true_sigma, "\n\n")

set.seed(123)
N <- 10000

# Generate alphas from MVN
Sigma_true <- build_correlation_matrix(2, true_sigma)
X <- mvrnorm(n = N, mu = true_means, Sigma = Sigma_true)
alpha_data <- (X > 0) * 1

# Check empirical profile frequencies
cat("Empirical alpha profile frequencies:\n")
for (i in 1:4) {
  count <- sum(alpha_data[,1] == profiles[i,1] & alpha_data[,2] == profiles[i,2])
  freq <- count / N
  cat(sprintf("  Profile (%d,%d): %d (%.4f)\n", profiles[i,1], profiles[i,2], count, freq))
}

# Theoretical frequencies
theory_probs <- compute_profile_probs(true_means, true_sigma)
cat("\nTheoretical profile probabilities:\n")
for (i in 1:4) {
  cat(sprintf("  Profile (%d,%d): %.4f\n", profiles[i,1], profiles[i,2], theory_probs[i]))
}
cat("\n")

# Generate responses from alpha using DINO rule
Y_data <- matrix(0, nrow = N, ncol = 5)

for (i in 1:N) {
  for (j in 1:5) {
    p_ij <- dino_item_probability(alpha_data[i,], Q[j,], true_guess[j], true_slip[j])
    Y_data[i, j] <- rbinom(1, 1, p_ij)
  }
}

cat("Data simulated successfully.\n")
cat("Sample responses (first 10):\n")
print(head(Y_data, 10))
cat("\n")

# Now compute likelihood at different sigma values
cat(rep("=", 70), "\n", sep = "")
cat("Computing likelihood at different sigma values\n")
cat("(all other parameters fixed at TRUE values)\n")
cat(rep("=", 70), "\n\n", sep = "")

sigma_grid <- seq(0.05, 0.50, by = 0.05)
P_LJ <- compute_dino_probabilities(profiles, Q, true_guess, true_slip)

results <- data.frame(
  sigma = sigma_grid,
  log_lik = NA,
  diff_from_true = NA
)

for (i in seq_along(sigma_grid)) {
  sigma <- sigma_grid[i]
  profile_probs <- compute_profile_probs(true_means, sigma)

  log_lik <- compute_log_likelihood(Y_data, P_LJ, profile_probs)

  results$log_lik[i] <- log_lik

  if (i %% 2 == 0) {
    cat(sprintf("sigma = %.2f: log-lik = %.2f\n", sigma, log_lik))
  }
}

results$diff_from_true <- results$log_lik - results$log_lik[results$sigma == true_sigma]

cat("\n")
print(results)

# Find maximum
max_idx <- which.max(results$log_lik)
cat(sprintf("\nMaximum at sigma = %.2f (true = %.2f)\n",
           results$sigma[max_idx], true_sigma))
cat(sprintf("Difference = %.2f\n", results$sigma[max_idx] - true_sigma))

# Plot profile
cat("\nCreating plot...\n")
pdf("/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat/Group0/likelihood_profile_manual.pdf",
    width = 10, height = 6)
par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))

plot(results$sigma, results$log_lik,
     type = "b", pch = 19, col = "blue",
     xlab = "Sigma (correlation parameter)",
     ylab = "Log-Likelihood",
     main = sprintf("Likelihood Profile (N=%d, true sigma=%.2f)", N, true_sigma),
     las = 1, lwd = 2)

grid()
abline(v = true_sigma, col = "red", lwd = 2, lty = 2)
abline(v = results$sigma[max_idx], col = "darkgreen", lwd = 2, lty = 3)

legend("bottomright",
       legend = c(
         sprintf("True sigma = %.2f", true_sigma),
         sprintf("ML sigma = %.2f", results$sigma[max_idx])
       ),
       col = c("red", "darkgreen"),
       lty = c(2, 3),
       lwd = 2)

dev.off()

cat("Plot saved.\n\n")

# Check derivative
idx_true <- which(results$sigma == true_sigma)
if (idx_true > 1 && idx_true < nrow(results)) {
  deriv <- (results$log_lik[idx_true + 1] - results$log_lik[idx_true - 1]) /
           (results$sigma[idx_true + 1] - results$sigma[idx_true - 1])

  cat(sprintf("Numerical derivative at true sigma: %.4f\n", deriv))

  if (deriv > 0.1) {
    cat("*** STRONG POSITIVE DERIVATIVE ***\n")
    cat("The likelihood is INCREASING at the true value!\n")
    cat("This confirms the bug.\n")
  } else if (abs(deriv) < 0.1) {
    cat("Derivative close to zero - likelihood appears correct.\n")
  }
}

cat("\n")
