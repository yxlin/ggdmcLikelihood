#!/usr/bin/env Rscript
cat("\n\n-------------- Testing CDM MLE with MVN (Multiple N's) ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "cdModel", "ggdmcPrior")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))
home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat"
mle_fun <- file.path(home_dir, "Group0/00_mle_helper.r")
setwd(home_dir)
cat("\nWorking directory: ", getwd(), "\n")
source(mle_fun)

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
), ncol = 2, byrow = TRUE)

colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3", "Item4", "Item5")

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
    verbose = TRUE
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINO",
    use_mvn = TRUE
)
means <- c(0.5, 0.2)
sigma <- 0.2

p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = means[1], mean2 = means[2], sigma = sigma,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
)


# Test multiple sample sizes
N_values <- c(50, 500, 5000, 50000, 500000)

# Set display options
options(digits = 4)
options(scipen = 999)

# Store results for comparison
all_results <- list()

cat("\n=======================================================================\n")
cat("Testing ggdmcLikelihood::compute_subject_likelihood with different N's\n")
cat("=======================================================================\n\n")

for (i in seq_along(N_values)) {
    N <- N_values[i]

    cat("\n-----------------------------------------------------------------------\n")
    cat("Test ", i, "/", length(N_values), ": N = ", N, "\n")
    cat("-----------------------------------------------------------------------\n")

    # Simulate data
    set.seed(123) # Same seed for reproducibility
    dat <- simulate(sub_model,
        nsim = N,
        parameter_vector = p_vector,
        nschool = 1,
        debug = FALSE
    )

    # Build DMI
    sub_dmis <- BuildDMI(dat$responses, model,
        q_matrix = model@cdm_info$q_matrix,
        profile_probability = model@cdm_info$profile_probability,
        rule = "DINO",
        use_mvn = TRUE
    )

    # Jittered starting values (reproducible)
    set.seed(456)
    start <- p_vector * runif(length(p_vector), 0.8, 1.2)

    # Fit model
    cat("\nFitting model with optim (L-BFGS-B)...\n")
    fit <- optim(
        par = start,
        fn = nll,
        method = "L-BFGS-B",
        lower = rep(1e-6, length(p_vector)),
        upper = rep(1 - 1e-6, length(p_vector))
    )

    # Compute results
    res <- rbind(
        True = p_vector,
        Estimate = fit$par,
        Diff = fit$par - p_vector,
        AbsDiff = abs(fit$par - p_vector)
    )

    # Store for later comparison
    all_results[[i]] <- list(
        N = N,
        estimates = fit$par,
        loglik = -fit$value,
        convergence = fit$convergence
    )

    # Display results
    cat("\nResults:\n")
    print(res)
    cat("\nMax Log-Likelihood: ", -fit$value, "\n")
    cat("Convergence: ", fit$convergence, " (0 = success)\n")
    cat("Mean Absolute Error: ", mean(abs(fit$par - p_vector)), "\n")

    # Quick sanity check with compute_subject_likelihood
    test_lik <- ggdmcLikelihood::compute_subject_likelihood(
        sub_dmis[[1]],
        p_vector,
        debug = FALSE
    )
    cat(
        "Test likelihood (true params): sum(log(lik)) = ",
        sum(log(pmax(test_lik[[1]], .Machine$double.xmin))), "\n"
    )
}

cat("\n=======================================================================\n")
cat("Summary Across All Sample Sizes\n")
cat("=======================================================================\n\n")

summary_mat <- matrix(NA, nrow = length(N_values), ncol = 4)
colnames(summary_mat) <- c("N", "LogLik", "MAE", "Convergence")
for (i in seq_along(all_results)) {
    summary_mat[i, 1] <- all_results[[i]]$N
    summary_mat[i, 2] <- all_results[[i]]$loglik
    summary_mat[i, 3] <- mean(abs(all_results[[i]]$estimates - p_vector))
    summary_mat[i, 4] <- all_results[[i]]$convergence
}
print(summary_mat)

cat("\n=======================================================================\n")
cat("Test complete!\n")
cat("=======================================================================\n\n")

# =======================================================================
# Summary Across All Sample Sizes
# =======================================================================

#           N     LogLik      MAE Convergence
# [1,]     50     -152.4 0.242307           0
# [2,]    500    -1571.7 0.197135           0
# [3,]   5000   -15737.1 0.014724           0
# [4,]  50000  -157897.8 0.011158           0
# [5,] 500000 -1578972.8 0.007057           0
