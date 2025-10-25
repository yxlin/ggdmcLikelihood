#!/usr/bin/env Rscript
# Diagnose why MCMC chains are still flat

cat("\n\n-------------- MCMC Flat Chain Diagnosis ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)
library(ggdmcPrior)

Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

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

cat("True parameters:\n")
print(true_p_vector)
cat("\n")

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

# Simulate data
cat("Simulating data...\n")
set.seed(123)
N <- 100  # Small for faster testing

dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = true_p_vector,
    nschool = 1,
    debug = FALSE
)

cat("Data simulated. N =", N, "\n\n")

# Build DMI
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = TRUE
)

# Load helper
source("Group0/00_mle_helper.r")

# Test 1: Check if likelihood works at various parameter values
cat(rep("=", 70), "\n", sep = "")
cat("TEST 1: Likelihood at different sigma values\n")
cat(rep("=", 70), "\n\n", sep = "")

sigma_test <- c(0.01, 0.10, 0.20, 0.30, 0.50, 0.80, 0.95, 0.99)

for (sigma in sigma_test) {
  test_params <- true_p_vector
  test_params["sigma"] <- sigma

  tryCatch({
    loglik <- sll_from_p(test_params, sub_dmis[[1]])
    cat(sprintf("sigma = %.2f: loglik = %.4f\n", sigma, loglik))
  }, error = function(e) {
    cat(sprintf("sigma = %.2f: ERROR - %s\n", sigma, e$message))
  })
}

# Test 2: Check if extreme guess/slip values work
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TEST 2: Likelihood with extreme guess/slip values\n")
cat(rep("=", 70), "\n\n", sep = "")

extreme_tests <- list(
  list(name = "All guess = 0.01", params = c(
    guess1 = .01, guess2 = .01, guess3 = .01, guess4 = .01, guess5 = .01,
    mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
  )),
  list(name = "All guess = 0.99", params = c(
    guess1 = .99, guess2 = .99, guess3 = .99, guess4 = .99, guess5 = .99,
    mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
  )),
  list(name = "All slip = 0.01", params = c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
    slip1 = .01, slip2 = .01, slip3 = .01, slip4 = .01, slip5 = .01
  )),
  list(name = "All slip = 0.99", params = c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
    slip1 = .99, slip2 = .99, slip3 = .99, slip4 = .99, slip5 = .99
  ))
)

for (test in extreme_tests) {
  tryCatch({
    loglik <- sll_from_p(test$params, sub_dmis[[1]])
    cat(sprintf("%-20s: loglik = %.4f\n", test$name, loglik))
  }, error = function(e) {
    cat(sprintf("%-20s: ERROR - %s\n", test$name, e$message))
  })
}

# Test 3: Test parameter validation bounds
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TEST 3: Parameter validation at boundaries\n")
cat(rep("=", 70), "\n\n", sep = "")

boundary_tests <- list(
  list(name = "sigma = 0.999", sigma = 0.999),
  list(name = "sigma = 1.000", sigma = 1.000),
  list(name = "sigma = -0.499", sigma = -0.499),  # For K=2, min = -1/(2-1) = -1
  list(name = "sigma = -0.999", sigma = -0.999),
  list(name = "guess = slip", params = c(
    guess1 = .5, guess2 = .5, guess3 = .5, guess4 = .5, guess5 = .5,
    mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
    slip1 = .5, slip2 = .5, slip3 = .5, slip4 = .5, slip5 = .5
  ))
)

for (test in boundary_tests) {
  if (!is.null(test$sigma)) {
    test_params <- true_p_vector
    test_params["sigma"] <- test$sigma
    test_name <- test$name
  } else {
    test_params <- test$params
    test_name <- test$name
  }

  tryCatch({
    loglik <- sll_from_p(test_params, sub_dmis[[1]])
    cat(sprintf("%-20s: loglik = %.4f [ACCEPTED]\n", test_name, loglik))
  }, error = function(e) {
    cat(sprintf("%-20s: ERROR - %s\n", test_name, e$message))
  })
}

# Test 4: Check if small perturbations change likelihood
cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("TEST 4: Likelihood sensitivity to small changes\n")
cat(rep("=", 70), "\n\n", sep = "")

baseline_loglik <- sll_from_p(true_p_vector, sub_dmis[[1]])
cat(sprintf("Baseline (true params): %.6f\n\n", baseline_loglik))

perturbations <- c(0.001, 0.01, 0.05, 0.10)

cat("Perturbations in sigma:\n")
for (delta in perturbations) {
  test_params <- true_p_vector
  test_params["sigma"] <- true_sigma + delta

  loglik <- sll_from_p(test_params, sub_dmis[[1]])
  diff <- loglik - baseline_loglik
  cat(sprintf("  sigma + %.3f: loglik = %.6f (diff = %+.6f)\n",
             delta, loglik, diff))
}

cat("\nPerturbations in mean1:\n")
for (delta in perturbations) {
  test_params <- true_p_vector
  test_params["mean1"] <- true_means[1] + delta

  loglik <- sll_from_p(test_params, sub_dmis[[1]])
  diff <- loglik - baseline_loglik
  cat(sprintf("  mean1 + %.3f: loglik = %.6f (diff = %+.6f)\n",
             delta, loglik, diff))
}

cat("\nPerturbations in guess1:\n")
for (delta in perturbations) {
  test_params <- true_p_vector
  test_params["guess1"] <- min(0.99, true_guess[1] + delta)

  loglik <- sll_from_p(test_params, sub_dmis[[1]])
  diff <- loglik - baseline_loglik
  cat(sprintf("  guess1 + %.3f: loglik = %.6f (diff = %+.6f)\n",
             delta, loglik, diff))
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("If TEST 1 shows errors at extreme sigma (0.99, etc.), the validation\n")
cat("bounds might be too strict and rejecting valid MCMC proposals.\n\n")

cat("If TEST 4 shows very small likelihood differences (<< 1), the likelihood\n")
cat("surface might be too flat, making MCMC unable to explore effectively.\n\n")

cat("If all tests pass but MCMC is still flat, the issue is in the MCMC\n")
cat("sampler itself (DE-MCMC algorithm, proposal generation, etc.).\n\n")
