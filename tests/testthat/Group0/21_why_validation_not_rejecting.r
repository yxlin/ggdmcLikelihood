#!/usr/bin/env Rscript
# Test why validation is not rejecting sigma=1.0

cat("\n\n-------------- Validation Debug Test ---------------\n")

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

# Simulate minimal data
set.seed(123)
dat <- simulate(sub_model,
    nsim = 10,
    parameter_vector = c(
      guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
      mean1 = 0.5, mean2 = 0.2, sigma = 0.2,
      slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
    ),
    nschool = 1,
    debug = FALSE
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = TRUE
)

source("Group0/00_mle_helper.r")

# Test with debug=TRUE to see if validation warnings appear
cat("Testing sigma values near boundary with debug=TRUE:\n\n")

sigma_tests <- c(0.95, 0.99, 0.999, 1.0, 1.001)

for (sigma in sigma_tests) {
  cat(rep("-", 70), "\n", sep = "")
  cat(sprintf("Testing sigma = %.4f:\n\n", sigma))

  test_params <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = 0.5, mean2 = 0.2, sigma = sigma,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
  )

  # Call with debug to see warnings
  lik_vec <- ggdmcLikelihood::compute_subject_likelihood(
    sub_dmis[[1]],
    test_params,
    debug = TRUE
  )

  loglik <- sum(log(pmax(lik_vec[[1]], .Machine$double.xmin)))

  cat(sprintf("\nResult: loglik = %.4f\n", loglik))
  cat(sprintf("Min likelihood value: %.10e\n", min(lik_vec[[1]])))
  cat(sprintf("Max likelihood value: %.10e\n", max(lik_vec[[1]])))

  if (any(lik_vec[[1]] < 1e-8)) {
    cat("*** WARNING: Some likelihood values < 1e-8! ***\n")
  }

  cat("\n")
}

cat("\n")
cat(rep("=", 70), "\n", sep = "")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("If no 'Parameter validation failed' warnings appear for sigma=1.0,\n")
cat("then validation is NOT being called or is passing incorrectly.\n\n")

cat("If warnings appear but likelihood is still calculated, then the\n")
cat("validation failure is being handled by setting lik=1e-10, which\n")
cat("still causes terrible log-likelihood.\n\n")
