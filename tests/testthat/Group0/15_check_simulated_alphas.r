#!/usr/bin/env Rscript
# Check if simulated alphas match the specified MVN parameters

cat("\n\n-------------- Simulated Alpha Check ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)
library(mvtnorm)

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

# Simulate with large N
cat("Simulating data with N = 10000...\n")
set.seed(123)
N <- 10000

dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = true_p_vector,
    nschool = 1,
    debug = FALSE
)

cat("Data simulated.\n\n")

# Extract alphas
cat("Checking dat$alpha:\n")
cat("  class:", class(dat$alpha), "\n")
cat("  dim:", dim(dat$alpha), "\n")
cat("  colnames:", colnames(dat$alpha), "\n\n")

print(head(dat$alpha, 20))

# Reshape alpha from long to wide format
# dat$alpha has columns: student, skill, C, s
# C contains the alpha values
K <- 2  # number of skills
alpha_matrix <- matrix(NA, nrow = N, ncol = K)

for (i in 1:nrow(dat$alpha)) {
  student_id <- dat$alpha$student[i]
  skill_id <- dat$alpha$skill[i]
  alpha_val <- dat$alpha$C[i]

  alpha_matrix[student_id, skill_id] <- alpha_val
}

cat("Alpha matrix dimensions:", dim(alpha_matrix), "\n")
cat("First 10 rows:\n")
print(head(alpha_matrix, 10))
cat("\n")

cat("\n\nEmpirical alpha statistics:\n")
cat("P(A1=1) =", mean(alpha_matrix[,1]), "\n")
cat("P(A2=1) =", mean(alpha_matrix[,2]), "\n\n")

# Theoretical marginal probabilities
P_A1_theory <- pnorm(0, mean = true_means[1], sd = 1, lower.tail = FALSE)
P_A2_theory <- pnorm(0, mean = true_means[2], sd = 1, lower.tail = FALSE)

cat("Theoretical marginal probabilities:\n")
cat("P(A1=1) =", P_A1_theory, "\n")
cat("P(A2=1) =", P_A2_theory, "\n\n")

cat("Differences:\n")
cat("P(A1=1):", mean(alpha_matrix[,1]) - P_A1_theory, "\n")
cat("P(A2=1):", mean(alpha_matrix[,2]) - P_A2_theory, "\n\n")

# Check empirical correlation
emp_cor <- cor(alpha_matrix)
cat("Empirical correlation matrix:\n")
print(emp_cor)

cat("\nOff-diagonal correlation:", emp_cor[1,2], "\n")
cat("True correlation (sigma):", true_sigma, "\n")
cat("Difference:", emp_cor[1,2] - true_sigma, "\n\n")

# Check joint probabilities
profiles <- matrix(c(
  0, 0,
  1, 0,
  0, 1,
  1, 1
), ncol = 2, byrow = TRUE)

cat("Profile frequencies:\n")
empirical_freqs <- numeric(4)
for (i in 1:4) {
  count <- sum(alpha_matrix[,1] == profiles[i,1] & alpha_matrix[,2] == profiles[i,2])
  empirical_freqs[i] <- count / N
  cat(sprintf("  Profile (%d,%d): %.6f (n=%d)\n",
             profiles[i,1], profiles[i,2], empirical_freqs[i], count))
}
cat(sprintf("  Sum: %.10f\n\n", sum(empirical_freqs)))

# Theoretical probabilities
build_correlation_matrix <- function(K, sigma) {
  Sigma <- matrix(sigma, nrow = K, ncol = K)
  diag(Sigma) <- 1
  return(Sigma)
}

Sigma <- build_correlation_matrix(2, true_sigma)
theory_probs <- numeric(4)
for (i in 1:4) {
  lower <- ifelse(profiles[i,] == 0, -Inf, 0)
  upper <- ifelse(profiles[i,] == 0, 0, Inf)
  theory_probs[i] <- pmvnorm(lower=lower, upper=upper, mean=true_means, sigma=Sigma)[1]
}

cat("Theoretical profile probabilities:\n")
for (i in 1:4) {
  cat(sprintf("  Profile (%d,%d): %.6f\n", profiles[i,1], profiles[i,2], theory_probs[i]))
}
cat(sprintf("  Sum: %.10f\n\n", sum(theory_probs)))

cat("Differences (Empirical - Theoretical):\n")
for (i in 1:4) {
  diff <- empirical_freqs[i] - theory_probs[i]
  cat(sprintf("  Profile (%d,%d): %+.6f (%.2f%%)\n",
             profiles[i,1], profiles[i,2], diff, 100 * diff / theory_probs[i]))
}

# Chi-square test
expected_counts <- theory_probs * N
observed_counts <- empirical_freqs * N
chi_sq <- sum((observed_counts - expected_counts)^2 / expected_counts)
df <- 3
p_value <- pchisq(chi_sq, df, lower.tail = FALSE)

cat(sprintf("\nChi-square test:\n"))
cat(sprintf("  Statistic: %.4f\n", chi_sq))
cat(sprintf("  df: %d\n", df))
cat(sprintf("  p-value: %.4f\n", p_value))

if (p_value < 0.05) {
  cat("  *** REJECT null: Simulated alphas DO NOT match MVN! ***\n")
  cat("  This would explain the likelihood bug!\n")
} else {
  cat("  Accept null: Simulated alphas match MVN distribution.\n")
  cat("  The bug is not in the simulation.\n")
}

cat("\n")
