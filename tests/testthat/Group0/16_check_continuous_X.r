#!/usr/bin/env Rscript
# Check if the continuous X values (before thresholding) have correct correlation

cat("\n\n-------------- Continuous X Correlation Check ---------------\n")
rm(list = ls())

library(MASS)

# Replicate the C++ simulation exactly
true_means <- c(0.5, 0.2)
true_sigma <- 0.2

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

Sigma <- build_correlation_matrix(2, true_sigma)

cat("Parameters:\n")
cat("  means:", true_means, "\n")
cat("  sigma:", true_sigma, "\n")
cat("  Correlation matrix:\n")
print(Sigma)
cat("\n")

# Generate data the same way as C++ code
set.seed(123)
N <- 10000

# Step 1: Z ~ N(0, I)
Z <- matrix(rnorm(N * 2), nrow = N, ncol = 2)

# Step 2: Cholesky decomposition
L <- t(chol(Sigma))  # R's chol() returns upper, we want lower

cat("Cholesky decomposition L:\n")
print(L)
cat("\nVerify L %*% t(L) = Sigma:\n")
print(L %*% t(L))
cat("\n")

# Step 3: X = Z * L^T + mean
X <- Z %*% t(L)  # N x K

# Step 4: Add mean
X <- sweep(X, 2, true_means, "+")

# Check X statistics BEFORE thresholding
cat("Continuous X statistics (BEFORE thresholding):\n")
cat("  Mean of X[,1]:", mean(X[,1]), "(should be", true_means[1], ")\n")
cat("  Mean of X[,2]:", mean(X[,2]), "(should be", true_means[2], ")\n")
cat("  SD of X[,1]:", sd(X[,1]), "(should be ~1)\n")
cat("  SD of X[,2]:", sd(X[,2]), "(should be ~1)\n\n")

cor_X <- cor(X)
cat("Correlation matrix of X:\n")
print(cor_X)
cat("\nOff-diagonal correlation:", cor_X[1,2], "\n")
cat("True correlation:", true_sigma, "\n")
cat("Difference:", cor_X[1,2] - true_sigma, "\n\n")

# Step 5: Threshold to get alpha
alpha <- (X > 0) * 1

# Check alpha statistics AFTER thresholding
cat("Binary alpha statistics (AFTER thresholding):\n")
cat("  P(alpha[,1] = 1):", mean(alpha[,1]), "\n")
cat("  P(alpha[,2] = 1):", mean(alpha[,2]), "\n\n")

cor_alpha <- cor(alpha)
cat("Pearson correlation of binary alpha:\n")
print(cor_alpha)
cat("\nOff-diagonal correlation:", cor_alpha[1,2], "\n")
cat("True continuous correlation:", true_sigma, "\n")
cat("Difference:", cor_alpha[1,2] - true_sigma, "\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("CONCLUSION\n")
cat(rep("=", 70), "\n\n", sep = "")

if (abs(cor_X[1,2] - true_sigma) < 0.01) {
  cat("Continuous X has CORRECT correlation (~", round(cor_X[1,2], 3), ")\n")
  cat("The thresholding reduces the Pearson correlation of binary alphas.\n")
  cat("This is EXPECTED behavior (tetrachoric correlation).\n\n")

  cat("The simulation is CORRECT.\n")
  cat("The likelihood calculation must properly account for this.\n")
} else {
  cat("*** BUG FOUND! ***\n")
  cat("Continuous X has WRONG correlation!\n")
  cat("Expected:", true_sigma, "\n")
  cat("Got:", cor_X[1,2], "\n")
  cat("The simulation MVN generation is broken.\n")
}

cat("\n")
