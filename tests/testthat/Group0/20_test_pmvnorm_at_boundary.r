#!/usr/bin/env Rscript
# Test what pmvnorm returns at boundary values

cat("\n\n-------------- pmvnorm Boundary Test ---------------\n")

library(mvtnorm)

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

means <- c(0.5, 0.2)
profiles <- matrix(c(0,0, 1,0, 0,1, 1,1), ncol=2, byrow=TRUE)

# Test different sigma values
sigma_values <- c(0.90, 0.95, 0.99, 0.999, 0.9999, 1.0, 1.0001)

cat("Testing pmvnorm at different sigma values:\n\n")

for (sigma in sigma_values) {
  cat(sprintf("sigma = %.5f:\n", sigma))

  Sigma <- build_correlation_matrix(2, sigma)

  # Check if matrix is positive definite
  eig_vals <- eigen(Sigma, only.values = TRUE)$values
  cat(sprintf("  Eigenvalues: %.6f, %.6f\n", eig_vals[1], eig_vals[2]))
  cat(sprintf("  Min eigenvalue: %.10f\n", min(eig_vals)))

  if (min(eig_vals) <= 0) {
    cat("  *** NOT POSITIVE DEFINITE! ***\n\n")
    next
  }

  # Try to compute pmvnorm for each profile
  probs <- numeric(4)
  valid <- TRUE

  for (i in 1:4) {
    lower <- ifelse(profiles[i,] == 0, -Inf, 0)
    upper <- ifelse(profiles[i,] == 0, 0, Inf)

    tryCatch({
      result <- pmvnorm(lower=lower, upper=upper, mean=means, sigma=Sigma)
      probs[i] <- result[1]

      if (!is.finite(probs[i]) || probs[i] < 0 || probs[i] > 1) {
        cat(sprintf("  Profile (%d,%d): %.10f [INVALID!]\n",
                   profiles[i,1], profiles[i,2], probs[i]))
        valid <- FALSE
      }
    }, error = function(e) {
      cat(sprintf("  Profile (%d,%d): ERROR - %s\n",
                 profiles[i,1], profiles[i,2], e$message))
      probs[i] <<- NA
      valid <<- FALSE
    })
  }

  if (valid && all(is.finite(probs))) {
    cat(sprintf("  Profile probabilities: %.6f, %.6f, %.6f, %.6f\n",
               probs[1], probs[2], probs[3], probs[4]))
    cat(sprintf("  Sum: %.10f\n", sum(probs)))

    if (abs(sum(probs) - 1.0) > 1e-6) {
      cat("  *** SUM != 1.0 ***\n")
    }
  }

  cat("\n")
}

cat(rep("=", 70), "\n", sep = "")
cat("FINDINGS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("If sigma = 1.0 causes:\n")
cat("1. Non-positive definite matrix -> chol() will fail\n")
cat("2. pmvnorm() returns NaN/Inf -> profile probs become invalid\n")
cat("3. Profile probs don't sum to 1 -> normalization fails\n\n")

cat("Any of these will cause catastrophic likelihood calculations.\n\n")

# Test what happens when we try to use these probabilities
cat("Testing what happens with sigma = 1.0 probabilities:\n\n")

sigma_test <- 1.0
Sigma_test <- build_correlation_matrix(2, sigma_test)

tryCatch({
  probs_test <- numeric(4)
  for (i in 1:4) {
    lower <- ifelse(profiles[i,] == 0, -Inf, 0)
    upper <- ifelse(profiles[i,] == 0, 0, Inf)
    probs_test[i] <- pmvnorm(lower=lower, upper=upper, mean=means, sigma=Sigma_test)[1]
  }

  cat("Profile probs:", probs_test, "\n")
  cat("Are finite?", all(is.finite(probs_test)), "\n")
  cat("Sum:", sum(probs_test), "\n")

  # Try to normalize (what C++ code does)
  L <- length(probs_test)
  denominator <- sum(probs_test + 1e-10)
  normalized <- probs_test / denominator
  clamped <- pmax(pmin(normalized, 1.0 - 1e-10), 1e-10)

  cat("\nAfter normalization and clamping:\n")
  cat("  Normalized:", normalized, "\n")
  cat("  Clamped:", clamped, "\n")
  cat("  Sum of clamped:", sum(clamped), "\n")

  # What likelihood would this give?
  # If profile probs are all 0.25 (uniform), log(0.25) = -1.386
  # For 100 observations: 100 * (-1.386) = -138.6
  # If they're much smaller: log(1e-10) = -23.03, 100 * -23.03 = -2303

  cat("\nIf all profile probs = 1e-10 (due to clamping):\n")
  cat("  Per-observation log-lik: log(1e-10) =", log(1e-10), "\n")
  cat("  For N=100: 100 *", log(1e-10), "=", 100 * log(1e-10), "\n")

  if (any(clamped == 1e-10)) {
    cat("\n*** Some probs hit lower clamp bound! This causes terrible likelihood! ***\n")
  }

}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})

cat("\n")
