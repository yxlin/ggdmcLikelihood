#!/usr/bin/env Rscript
# Test if the normalization + clamping causes probabilities to not sum to 1

cat("\n\n-------------- Normalization Bug Test ---------------\n")

library(mvtnorm)

# Simulate what the C++ code does at lines 900-910 in cdm.h
normalize_and_clamp <- function(probs) {
  # Line 908: tmp_prob = tmp_prob / arma::accu(tmp_prob + 1e-10);
  # This should be: tmp_prob / sum(tmp_prob + 1e-10)
  # But wait - in Armadillo, arma::accu(vec) sums all elements
  # So arma::accu(tmp_prob + 1e-10) = sum(tmp_prob) + length(tmp_prob)*1e-10

  L <- length(probs)
  normalized <- probs / (sum(probs) + L * 1e-10)

  # Line 909-910: clamp to [1e-10, 1-1e-10]
  clamped <- pmax(pmin(normalized, 1.0 - 1e-10), 1e-10)

  return(list(
    original = probs,
    normalized = normalized,
    clamped = clamped,
    sum_original = sum(probs),
    sum_normalized = sum(normalized),
    sum_clamped = sum(clamped)
  ))
}

# Test case 1: Normal probabilities
cat("\nTest 1: Normal probabilities (all positive)\n")
cat(rep("-", 50), "\n", sep = "")
probs1 <- c(0.157766, 0.262975, 0.150772, 0.428488)
result1 <- normalize_and_clamp(probs1)
print(result1)

cat("\nDoes clamped sum to 1? ", abs(result1$sum_clamped - 1.0) < 1e-10, "\n")
cat("Error:", result1$sum_clamped - 1.0, "\n")

# Test case 2: Include some very small probabilities
cat("\n\nTest 2: With very small probabilities\n")
cat(rep("-", 50), "\n", sep = "")
probs2 <- c(0.001, 0.299, 0.200, 0.500)
result2 <- normalize_and_clamp(probs2)
print(result2)

cat("\nDoes clamped sum to 1? ", abs(result2$sum_clamped - 1.0) < 1e-10, "\n")
cat("Error:", result2$sum_clamped - 1.0, "\n")

# Test case 3: Including zeros
cat("\n\nTest 3: With exact zeros\n")
cat(rep("-", 50), "\n", sep = "")
probs3 <- c(0.0, 0.3, 0.2, 0.5)
result3 <- normalize_and_clamp(probs3)
print(result3)

cat("\nDoes clamped sum to 1? ", abs(result3$sum_clamped - 1.0) < 1e-10, "\n")
cat("Error:", result3$sum_clamped - 1.0, "\n")

# Test case 4: What pmvnorm actually returns
cat("\n\nTest 4: Actual pmvnorm output\n")
cat(rep("-", 50), "\n", sep = "")

build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

mean_vec <- c(0.5, 0.2)
Sigma <- build_correlation_matrix(2, 0.2)

profiles <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1
), ncol = 2, byrow = TRUE)

probs4 <- numeric(4)
for (i in 1:4) {
  lower <- ifelse(profiles[i, ] == 0, -Inf, 0)
  upper <- ifelse(profiles[i, ] == 0, 0, Inf)
  probs4[i] <- pmvnorm(lower = lower, upper = upper,
                       mean = mean_vec, sigma = Sigma)[1]
}

result4 <- normalize_and_clamp(probs4)
print(result4)

cat("\nDoes clamped sum to 1? ", abs(result4$sum_clamped - 1.0) < 1e-10, "\n")
cat("Error:", result4$sum_clamped - 1.0, "\n")

# CRITICAL TEST: What if normalization is WRONG in C++?
cat("\n\n", rep("=", 70), "\n", sep = "")
cat("CRITICAL BUG CHECK\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nThe C++ code has:\n")
cat("  tmp_prob = tmp_prob / arma::accu(tmp_prob + 1e-10);\n\n")

cat("This divides by: sum(tmp_prob + 1e-10)\n")
cat("Which equals: sum(tmp_prob) + L*1e-10\n\n")

cat("For our test case with L=4:\n")
cat("  sum(tmp_prob) = ", sum(probs4), "\n")
cat("  L*1e-10 = ", 4*1e-10, "\n")
cat("  Denominator = ", sum(probs4) + 4*1e-10, "\n\n")

cat("After division:\n")
cat("  sum(normalized) = sum(tmp_prob) / (sum(tmp_prob) + L*1e-10)\n")
cat("                  = ", sum(probs4), " / ", sum(probs4) + 4*1e-10, "\n")
cat("                  = ", sum(probs4) / (sum(probs4) + 4*1e-10), "\n\n")

cat("This is LESS than 1.0 by:", 1.0 - sum(probs4) / (sum(probs4) + 4*1e-10), "\n\n")

cat("After clamping (assuming no values hit bounds):\n")
cat("  sum(clamped) = ", result4$sum_normalized, "\n\n")

if (abs(result4$sum_normalized - 1.0) > 1e-15) {
  cat("*** BUG FOUND! ***\n")
  cat("The normalization REDUCES the sum below 1.0!\n")
  cat("This means the likelihood uses profile probabilities that don't sum to 1.\n")
  cat("This could cause the likelihood to be biased!\n")
} else {
  cat("No bug - probabilities still sum to 1 within machine precision.\n")
}

cat("\n")
