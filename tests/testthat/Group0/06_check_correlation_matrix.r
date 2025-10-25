#!/usr/bin/env Rscript
# Test to verify correlation matrix construction and Cholesky decomposition

cat("\n\n-------------- Correlation Matrix Verification ---------------\n")

# Test the correlation matrix that should be built in C++
build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

# Test case 1: K=2, sigma=0.2 (from diagnostic)
K <- 2
sigma <- 0.2
Sigma <- build_correlation_matrix(K, sigma)

cat("\nTest 1: K=2, sigma=0.2\n")
cat("Correlation matrix:\n")
print(Sigma)

cat("\nEigenvalues:\n")
eig <- eigen(Sigma)
print(eig$values)
cat("Is positive definite?", all(eig$values > 0), "\n")

cat("\nCholesky decomposition:\n")
L <- t(chol(Sigma))  # R's chol() returns upper tri, we want lower
cat("L (lower triangular):\n")
print(L)

cat("\nVerify L %*% t(L) == Sigma:\n")
print(L %*% t(L))
cat("Match?", all.equal(L %*% t(L), Sigma), "\n")

# Test case 2: Simulate and check profile probabilities
cat("\n\n-------------- Profile Probability Check ---------------\n")

library(mvtnorm)

mean_vec <- c(0.5, 0.2)
sigma_val <- 0.2
Sigma <- build_correlation_matrix(2, sigma_val)

cat("\nMean:", mean_vec, "\n")
cat("Sigma value:", sigma_val, "\n")
cat("Correlation matrix:\n")
print(Sigma)

# Calculate profile probabilities for all 4 profiles (2^2)
# Profile (0,0): both skills not mastered
# Profile (1,0): skill 1 mastered, skill 2 not
# Profile (0,1): skill 1 not, skill 2 mastered
# Profile (1,1): both skills mastered

profiles <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1
), ncol = 2, byrow = TRUE)

profile_probs <- numeric(4)

for (i in 1:4) {
  lower <- ifelse(profiles[i, ] == 0, -Inf, 0)
  upper <- ifelse(profiles[i, ] == 0, 0, Inf)

  profile_probs[i] <- pmvnorm(lower = lower, upper = upper,
                               mean = mean_vec, sigma = Sigma)[1]
}

cat("\nProfile probabilities:\n")
for (i in 1:4) {
  cat(sprintf("  Profile (%d,%d): %.6f\n", profiles[i,1], profiles[i,2], profile_probs[i]))
}

cat("\nSum of probabilities:", sum(profile_probs), "\n")
cat("Error from 1.0:", sum(profile_probs) - 1.0, "\n")

# Marginal probabilities
P_A1 <- pnorm(0, mean = mean_vec[1], sd = 1, lower.tail = FALSE)
P_A2 <- pnorm(0, mean = mean_vec[2], sd = 1, lower.tail = FALSE)

cat("\nMarginal probabilities (if independent):\n")
cat("  P(A1=1) =", P_A1, "\n")
cat("  P(A2=1) =", P_A2, "\n")

cat("\nIf independent, P(1,1) would be:", P_A1 * P_A2, "\n")
cat("With correlation sigma=", sigma_val, ", P(1,1) =", profile_probs[4], "\n")
cat("Difference:", profile_probs[4] - P_A1 * P_A2, "\n")

# Test case 3: Check what happens with different sigma values
cat("\n\n-------------- Sigma Sensitivity Check ---------------\n")

sigma_grid <- c(0.01, 0.1, 0.2, 0.3, 0.5, 0.8)
results <- data.frame(
  sigma = sigma_grid,
  P_00 = NA,
  P_10 = NA,
  P_01 = NA,
  P_11 = NA,
  sum = NA
)

for (j in seq_along(sigma_grid)) {
  Sigma_test <- build_correlation_matrix(2, sigma_grid[j])

  for (i in 1:4) {
    lower <- ifelse(profiles[i, ] == 0, -Inf, 0)
    upper <- ifelse(profiles[i, ] == 0, 0, Inf)

    prob <- pmvnorm(lower = lower, upper = upper,
                     mean = mean_vec, sigma = Sigma_test)[1]
    results[j, i + 1] <- prob
  }
  results$sum[j] <- sum(results[j, 2:5])
}

print(results)

cat("\nObservations:\n")
cat("- All sums should be very close to 1.0\n")
cat("- As sigma increases (more positive correlation), P(both mastered) should increase\n")
cat("- As sigma increases, P(both not mastered) should also increase\n")
cat("- P(one mastered, one not) should decrease\n")

# Verify these patterns
cat("\nP(1,1) trend as sigma increases:\n")
print(results[, c("sigma", "P_11")])
cat("Is P(1,1) increasing? ", all(diff(results$P_11) > -1e-10), "\n")

cat("\n")
