#!/usr/bin/env Rscript
# Fixed comparison with proper data reshaping

cat("\n\n-------------- FIXED Package vs Manual Comparison ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)
library(ggdmcPrior)
library(mvtnorm)

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

compute_log_likelihood_manual <- function(Y_data, P_LJ, profile_probs) {
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

    max_log_lik <- max(log_lik_profiles)
    exp_shifted <- exp(log_lik_profiles - max_log_lik)
    weighted_sum <- sum(exp_shifted * profile_probs)
    log_likelihood_i <- max_log_lik + log(weighted_sum)

    total_log_lik <- total_log_lik + log_likelihood_i
  }

  return(total_log_lik)
}

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

compute_profile_probs <- function(means, sigma, profiles) {
  Sigma <- build_correlation_matrix(2, sigma)
  L <- nrow(profiles)
  probs <- numeric(L)
  for (i in 1:L) {
    lower <- ifelse(profiles[i,] == 0, -Inf, 0)
    upper <- ifelse(profiles[i,] == 0, 0, Inf)
    probs[i] <- pmvnorm(lower=lower, upper=upper, mean=means, sigma=Sigma)[1]
  }
  return(probs / sum(probs))
}

# Setup
Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

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

true_p_vector <- c(
  guess1 = true_guess[1], guess2 = true_guess[2], guess3 = true_guess[3],
  guess4 = true_guess[4], guess5 = true_guess[5],
  mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
  slip1 = true_slip[1], slip2 = true_slip[2], slip3 = true_slip[3],
  slip4 = true_slip[4], slip5 = true_slip[5]
)

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
    rule = "DINO",
    use_mvn = TRUE
)

# Simulate
cat("Simulating data...\n")
set.seed(123)
N <- 1000

dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = true_p_vector,
    nschool = 1,
    debug = FALSE
)

cat("Data simulated. N =", N, "\n\n")

# Reshape data from long to wide format
# dat$responses has columns: student, item, C, s
# We need N x J matrix where Y[i,j] = C value for student i, item j

cat("Reshaping data from long to wide format...\n")
J <- 5  # number of items
Y_matrix <- matrix(NA, nrow = N, ncol = J)

for (i in 1:nrow(dat$responses)) {
  student_id <- dat$responses$student[i]
  item_id <- dat$responses$item[i]
  response <- dat$responses$C[i]

  Y_matrix[student_id, item_id] <- response
}

cat("Y_matrix dimensions:", dim(Y_matrix), "\n")
cat("First 10 rows:\n")
print(head(Y_matrix, 10))
cat("\n")

# Build DMI
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = TRUE
)

# Load MLE helper
source("Group0/00_mle_helper.r")

# Compare likelihoods
cat(rep("=", 70), "\n", sep = "")
cat("Comparing PACKAGE vs MANUAL likelihood\n")
cat(rep("=", 70), "\n\n", sep = "")

sigma_test <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)

results <- data.frame(
  sigma = sigma_test,
  loglik_package = NA,
  loglik_manual = NA,
  difference = NA
)

P_LJ <- compute_dino_probabilities(profiles, Q, true_guess, true_slip)

for (i in seq_along(sigma_test)) {
  sigma <- sigma_test[i]

  # PACKAGE likelihood
  test_params <- true_p_vector
  test_params["sigma"] <- sigma
  loglik_pkg <- sll_from_p(test_params, sub_dmis[[1]])

  # MANUAL likelihood
  profile_probs <- compute_profile_probs(true_means, sigma, profiles)
  loglik_manual <- compute_log_likelihood_manual(Y_matrix, P_LJ, profile_probs)

  results$loglik_package[i] <- loglik_pkg
  results$loglik_manual[i] <- loglik_manual
  results$difference[i] <- loglik_pkg - loglik_manual

  cat(sprintf("sigma = %.2f: Pkg=%.2f, Manual=%.2f, Diff=%.4f\n",
             sigma, loglik_pkg, loglik_manual, results$difference[i]))
}

cat("\n")
print(results)

max_pkg <- which.max(results$loglik_package)
max_manual <- which.max(results$loglik_manual)

cat(sprintf("\nPackage maximum at sigma = %.2f (true = %.2f)\n", results$sigma[max_pkg], true_sigma))
cat(sprintf("Manual maximum at sigma = %.2f (true = %.2f)\n\n", results$sigma[max_manual], true_sigma))

if (max_pkg == which(results$sigma == true_sigma) && max_manual == which(results$sigma == true_sigma)) {
  cat("*** BOTH CORRECT! ***\n")
} else if (max_pkg == which(results$sigma == true_sigma)) {
  cat("Package CORRECT, Manual WRONG\n")
} else if (max_manual == which(results$sigma == true_sigma)) {
  cat("Manual CORRECT, Package WRONG\n")
} else {
  cat("*** BOTH WRONG! ***\n")
}

# Check difference pattern
cat("\nDifference pattern:\n")
cat("  Mean:", mean(results$difference), "\n")
cat("  SD:", sd(results$difference), "\n")
cat("  Max abs:", max(abs(results$difference)), "\n")

if (max(abs(results$difference)) < 0.01) {
  cat("\nDifferences are tiny - implementations match!\n")
} else if (sd(results$difference) < 0.01) {
  cat("\nConstant offset - simple bug.\n")
} else {
  cat("\nVarying differences - parameter-dependent bug!\n")
}

cat("\n")
