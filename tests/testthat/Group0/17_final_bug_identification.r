#!/usr/bin/env Rscript
# Final bug identification: comprehensive test

cat("\n\n-------------- FINAL BUG IDENTIFICATION ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)
library(mvtnorm)
library(MASS)

# Setup model
Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

true_means <- c(0.5, 0.2)
true_sigma <- 0.2

true_p_vector <- c(
  guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
  mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
  slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
)

cat("Test parameters:\n")
cat("  mean1 =", true_means[1], "\n")
cat("  mean2 =", true_means[2], "\n")
cat("  sigma =", true_sigma, "\n\n")

# Build package model
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

# Test different sample sizes
N_values <- c(1000, 5000, 10000)

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

# Compute theoretical profile probabilities
Sigma <- build_correlation_matrix(2, true_sigma)
profiles <- matrix(c(0,0, 1,0, 0,1, 1,1), ncol=2, byrow=TRUE)

theory_probs <- numeric(4)
for (i in 1:4) {
  lower <- ifelse(profiles[i,] == 0, -Inf, 0)
  upper <- ifelse(profiles[i,] == 0, 0, Inf)
  theory_probs[i] <- pmvnorm(lower=lower, upper=upper, mean=true_means, sigma=Sigma)[1]
}

cat("Theoretical profile probabilities:\n")
for (i in 1:4) {
  cat(sprintf("  P(%d,%d) = %.6f\n", profiles[i,1], profiles[i,2], theory_probs[i]))
}
cat("\n")

summary_table <- data.frame(
  N = integer(),
  emp_cor = numeric(),
  pkg_pvalue = numeric(),
  manual_pvalue = numeric(),
  pkg_reject = logical(),
  manual_reject = logical()
)

for (N in N_values) {
  cat(rep("=", 70), "\n", sep = "")
  cat(sprintf("Testing with N = %d\n", N))
  cat(rep("=", 70), "\n\n", sep = "")

  # PACKAGE simulation
  cat("1. PACKAGE simulation:\n")
  set.seed(123)
  dat_pkg <- simulate(sub_model,
      nsim = N,
      parameter_vector = true_p_vector,
      nschool = 1,
      debug = FALSE
  )

  # Extract alphas
  K <- 2
  alpha_pkg <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:nrow(dat_pkg$alpha)) {
    student_id <- dat_pkg$alpha$student[i]
    skill_id <- dat_pkg$alpha$skill[i]
    alpha_val <- dat_pkg$alpha$C[i]
    alpha_pkg[student_id, skill_id] <- alpha_val
  }

  # Check empirical profile frequencies
  emp_freqs_pkg <- numeric(4)
  for (i in 1:4) {
    count <- sum(alpha_pkg[,1] == profiles[i,1] & alpha_pkg[,2] == profiles[i,2])
    emp_freqs_pkg[i] <- count / N
  }

  # Chi-square test
  expected_counts <- theory_probs * N
  observed_counts_pkg <- emp_freqs_pkg * N
  chi_sq_pkg <- sum((observed_counts_pkg - expected_counts)^2 / expected_counts)
  p_value_pkg <- pchisq(chi_sq_pkg, df=3, lower.tail=FALSE)

  # Correlation
  cor_pkg <- cor(alpha_pkg)[1,2]

  cat(sprintf("  Empirical correlation: %.4f (true: %.2f)\n", cor_pkg, true_sigma))
  cat(sprintf("  Chi-square: %.4f, p-value: %.4f\n", chi_sq_pkg, p_value_pkg))
  if (p_value_pkg < 0.05) {
    cat("  *** REJECT: Profile frequencies don't match theory! ***\n")
  } else {
    cat("  Accept: Profile frequencies match theory\n")
  }
  cat("\n")

  # MANUAL R simulation (correct implementation)
  cat("2. MANUAL R simulation:\n")
  set.seed(456)  # Different seed since RNG streams differ
  Z <- matrix(rnorm(N * 2), nrow = N, ncol = 2)
  L <- t(chol(Sigma))
  X <- Z %*% t(L)
  X <- sweep(X, 2, true_means, "+")
  alpha_manual <- (X > 0) * 1

  # Check empirical profile frequencies
  emp_freqs_manual <- numeric(4)
  for (i in 1:4) {
    count <- sum(alpha_manual[,1] == profiles[i,1] & alpha_manual[,2] == profiles[i,2])
    emp_freqs_manual[i] <- count / N
  }

  # Chi-square test
  observed_counts_manual <- emp_freqs_manual * N
  chi_sq_manual <- sum((observed_counts_manual - expected_counts)^2 / expected_counts)
  p_value_manual <- pchisq(chi_sq_manual, df=3, lower.tail=FALSE)

  # Correlation
  cor_manual <- cor(alpha_manual)[1,2]

  cat(sprintf("  Empirical correlation: %.4f (true: %.2f)\n", cor_manual, true_sigma))
  cat(sprintf("  Chi-square: %.4f, p-value: %.4f\n", chi_sq_manual, p_value_manual))
  if (p_value_manual < 0.05) {
    cat("  REJECT: Profile frequencies don't match theory\n")
  } else {
    cat("  Accept: Profile frequencies match theory\n")
  }
  cat("\n")

  # Add to summary table
  summary_table <- rbind(summary_table, data.frame(
    N = N,
    emp_cor_pkg = cor_pkg,
    emp_cor_manual = cor_manual,
    pkg_pvalue = p_value_pkg,
    manual_pvalue = p_value_manual,
    pkg_reject = p_value_pkg < 0.05,
    manual_reject = p_value_manual < 0.05
  ))
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("SUMMARY TABLE\n")
cat(rep("=", 70), "\n\n", sep = "")

print(summary_table)

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("CONCLUSION\n")
cat(rep("=", 70), "\n\n", sep = "")

if (any(summary_table$pkg_reject) && !any(summary_table$manual_reject)) {
  cat("*** BUG CONFIRMED IN PACKAGE! ***\n\n")
  cat("The C++ package simulation generates alphas that DON'T match the MVN theory,\n")
  cat("while the correct R implementation DOES match theory.\n\n")
  cat("The bug is in the C++ simulation code in cdm.h\n")
  cat("Likely location: lines 850-866 (MVN generation)\n\n")

  cat("Evidence:\n")
  cat(sprintf("- Package empirical correlations: %.4f to %.4f\n",
             min(summary_table$emp_cor_pkg), max(summary_table$emp_cor_pkg)))
  cat(sprintf("- Manual empirical correlations: %.4f to %.4f\n",
             min(summary_table$emp_cor_manual), max(summary_table$emp_cor_manual)))
  cat(sprintf("- True continuous correlation: %.2f\n", true_sigma))
} else if (all(summary_table$pkg_reject) && all(summary_table$manual_reject)) {
  cat("BOTH implementations reject - the theory prediction might be wrong\n")
  cat("(This is actually expected for binary threshold of continuous MVN)\n")
} else {
  cat("Results are inconclusive. Need more investigation.\n")
}

cat("\n")
