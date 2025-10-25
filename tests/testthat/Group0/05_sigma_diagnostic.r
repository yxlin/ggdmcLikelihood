#!/usr/bin/env Rscript
# Diagnostic script to test if likelihood is maximized at true sigma value
# Tests the hypothesis that fixing all parameters except sigma, the likelihood
# might be larger at non-true sigma values

cat("\n\n-------------- Sigma Likelihood Profile Diagnostic ---------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "cdModel", "ggdmcPrior")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat"
mle_fun <- file.path(home_dir, "Group0/00_mle_helper.r")
setwd(home_dir)
cat("Working directory: ", getwd(), "\n\n")
source(mle_fun)

# ============================================================================
# Model setup
# ============================================================================
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
    rule = "DINA",
    use_mvn = TRUE
)

# ============================================================================
# True parameter values
# ============================================================================
true_means <- c(0.5, 0.2)
true_sigma <- 0.2 # TRUE VALUE we want to recover

true_p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
)

cat("True parameter values:\n")
print(true_p_vector)
cat("\n")

# ============================================================================
# Function to evaluate likelihood profile for sigma
# ============================================================================
evaluate_sigma_profile <- function(dmi, true_params, sigma_grid) {
    results <- data.frame(
        sigma = sigma_grid,
        loglik = NA_real_,
        lik = NA_real_
    )

    for (i in seq_along(sigma_grid)) {
        # Create parameter vector with all true values except sigma
        test_params <- true_params
        test_params["sigma"] <- sigma_grid[i]

        # Compute likelihood
        lik_vec <- ggdmcLikelihood::compute_subject_likelihood(
            dmi, test_params,
            debug = FALSE
        )[[1]]

        # Sum log-likelihood (with safeguard against zeros)
        loglik <- sum(log(pmax(lik_vec, .Machine$double.xmin)))

        results$loglik[i] <- loglik
        results$lik[i] <- exp(loglik)

        if (i %% 10 == 0) {
            cat(sprintf(
                "  Progress: %d/%d (sigma=%.3f, loglik=%.2f)\n",
                i, length(sigma_grid), sigma_grid[i], loglik
            ))
        }
    }

    return(results)
}

# ============================================================================
# Function to summarize and test results
# ============================================================================
summarize_profile <- function(results, true_sigma, N) {
    # Find maximum likelihood
    max_idx <- which.max(results$loglik)
    ml_sigma <- results$sigma[max_idx]
    ml_loglik <- results$loglik[max_idx]

    # Find likelihood at true sigma
    true_idx <- which.min(abs(results$sigma - true_sigma))
    true_loglik <- results$loglik[true_idx]

    # Calculate difference
    loglik_diff <- ml_loglik - true_loglik

    cat("\n", rep("=", 70), "\n", sep = "")
    cat(sprintf("SUMMARY for N = %d\n", N))
    cat(rep("=", 70), "\n", sep = "")
    cat(sprintf("True sigma:              %.4f\n", true_sigma))
    cat(sprintf("ML estimate (sigma):     %.4f\n", ml_sigma))
    cat(sprintf("Difference:              %.4f\n", ml_sigma - true_sigma))
    cat(sprintf("\n"))
    cat(sprintf("LogLik at true sigma:    %.4f\n", true_loglik))
    cat(sprintf("LogLik at ML sigma:      %.4f\n", ml_loglik))
    cat(sprintf("Difference:              %.4f\n", loglik_diff))
    cat(sprintf("\n"))

    if (abs(loglik_diff) < 0.01) {
        cat("RESULT: Likelihood is maximized at or very near true sigma (PASS)\n")
        status <- "PASS"
    } else if (loglik_diff > 0.01) {
        cat("WARNING: Likelihood is HIGHER at non-true sigma value (FAIL)\n")
        cat(sprintf("         This suggests a potential issue with the likelihood function!\n"))
        status <- "FAIL"
    } else {
        cat("NOTE: Likelihood is slightly lower at ML estimate (check needed)\n")
        status <- "CHECK"
    }
    cat(rep("=", 70), "\n\n", sep = "")

    return(list(
        ml_sigma = ml_sigma,
        true_sigma = true_sigma,
        ml_loglik = ml_loglik,
        true_loglik = true_loglik,
        diff = loglik_diff,
        status = status
    ))
}

# ============================================================================
# Test with multiple sample sizes
# ============================================================================
N_values <- c(100, 500, 2000, 5000)
# Create sigma grid: focused around true value with fine resolution
sigma_grid <- seq(0.01, 0.95, by = 0.02)

cat("\n", rep("=", 70), "\n", sep = "")
cat("Testing likelihood profile across different sample sizes\n")
cat("Sigma grid: ", length(sigma_grid), " values from ",
    min(sigma_grid), " to ", max(sigma_grid), "\n",
    sep = ""
)
cat(rep("=", 70), "\n\n", sep = "")

all_results <- list()
all_summaries <- list()

for (i in seq_along(N_values)) {
    N <- N_values[i]

    cat("\n", rep("-", 70), "\n", sep = "")
    cat(sprintf("TEST %d/%d: N = %d\n", i, length(N_values), N))
    cat(rep("-", 70), "\n", sep = "")

    # Generate data with true parameters
    set.seed(123 + i) # Different seed for each N
    dat <- simulate(sub_model,
        nsim = N,
        parameter_vector = true_p_vector,
        nschool = 1,
        debug = FALSE
    )

    cat(sprintf("Generated %d observations\n", N))

    # Build DMI
    sub_dmis <- BuildDMI(dat$responses, model,
        q_matrix = model@cdm_info$q_matrix,
        profile_probability = model@cdm_info$profile_probability,
        rule = "DINA",
        use_mvn = TRUE
    )

    cat("Evaluating likelihood profile...\n")

    # Evaluate likelihood profile
    results <- evaluate_sigma_profile(sub_dmis[[1]], true_p_vector, sigma_grid)
    results$N <- N
    all_results[[i]] <- results

    # Summarize and test
    summary <- summarize_profile(results, true_sigma, N)
    all_summaries[[i]] <- summary
}

# ============================================================================
# Create combined plot
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("Creating diagnostic plots...\n")
cat(rep("=", 70), "\n\n", sep = "")

# Combine all results
combined_results <- do.call(rbind, all_results)

# Create plot filename
plot_file <- file.path(home_dir, "Group0/sigma_likelihood_profile.pdf")
pdf(plot_file, width = 12, height = 8)

# Set up 2x2 plot layout
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

for (i in seq_along(N_values)) {
    N <- N_values[i]
    data_i <- all_results[[i]]
    summary_i <- all_summaries[[i]]

    # Plot likelihood profile
    plot(data_i$sigma, data_i$loglik,
        type = "l", lwd = 2, col = "blue",
        xlab = "Sigma", ylab = "Log-Likelihood",
        main = sprintf("N = %d (Status: %s)", N, summary_i$status),
        las = 1
    )

    # Add grid
    grid(col = "gray80")

    # Mark true sigma
    abline(v = true_sigma, col = "red", lwd = 2, lty = 2)

    # Mark ML estimate
    abline(v = summary_i$ml_sigma, col = "darkgreen", lwd = 2, lty = 3)

    # Add horizontal line at true loglik
    abline(h = summary_i$true_loglik, col = "red", lwd = 1, lty = 2)

    # Add legend
    legend("bottomright",
        legend = c(
            sprintf("True sigma = %.3f", true_sigma),
            sprintf("ML sigma = %.3f", summary_i$ml_sigma),
            sprintf("Diff = %.4f", summary_i$diff)
        ),
        col = c("red", "darkgreen", "black"),
        lty = c(2, 3, 0),
        lwd = c(2, 2, 0),
        bty = "n",
        cex = 0.8
    )
}

dev.off()
cat(sprintf("Plot saved to: %s\n\n", plot_file))

# ============================================================================
# Final summary table
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("FINAL SUMMARY TABLE\n")
cat(rep("=", 70), "\n", sep = "")

summary_table <- data.frame(
    N = N_values,
    True_Sigma = sapply(all_summaries, function(x) x$true_sigma),
    ML_Sigma = sapply(all_summaries, function(x) x$ml_sigma),
    Sigma_Error = sapply(all_summaries, function(x) x$ml_sigma - x$true_sigma),
    LogLik_True = sapply(all_summaries, function(x) x$true_loglik),
    LogLik_ML = sapply(all_summaries, function(x) x$ml_loglik),
    LogLik_Diff = sapply(all_summaries, function(x) x$diff),
    Status = sapply(all_summaries, function(x) x$status)
)

print(summary_table, row.names = FALSE)
cat("\n")

# Count failures
n_fail <- sum(summary_table$Status == "FAIL")
n_pass <- sum(summary_table$Status == "PASS")

cat(rep("=", 70), "\n", sep = "")
cat(sprintf("OVERALL RESULT: %d PASS, %d FAIL\n", n_pass, n_fail))
cat(rep("=", 70), "\n", sep = "")

if (n_fail > 0) {
    cat("\nWARNING: Some tests failed!\n")
    cat("The likelihood function may not be correctly identifying the true sigma.\n")
    cat("This suggests a potential bug in the MVN CDM implementation.\n")
} else {
    cat("\nAll tests passed! Likelihood appears to be correctly specified.\n")
}

cat("\n")

# ============================================================================
# Additional diagnostic: Check curvature at true sigma
# ============================================================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("CURVATURE DIAGNOSTIC (at true sigma)\n")
cat(rep("=", 70), "\n", sep = "")

for (i in seq_along(N_values)) {
    N <- N_values[i]
    data_i <- all_results[[i]]

    # Find neighborhood around true sigma
    true_idx <- which.min(abs(data_i$sigma - true_sigma))
    neighborhood <- (true_idx - 2):(true_idx + 2)
    neighborhood <- neighborhood[neighborhood > 0 & neighborhood <= nrow(data_i)]

    cat(sprintf("\nN = %d:\n", N))
    cat("  Sigma values near true:\n")
    for (idx in neighborhood) {
        marker <- if (abs(data_i$sigma[idx] - true_sigma) < 0.01) " <-- TRUE" else ""
        cat(sprintf(
            "    sigma = %.3f, loglik = %.4f%s\n",
            data_i$sigma[idx], data_i$loglik[idx], marker
        ))
    }

    # Check if true sigma is a local maximum
    if (true_idx > 1 && true_idx < nrow(data_i)) {
        is_local_max <- (data_i$loglik[true_idx] >= data_i$loglik[true_idx - 1]) &&
            (data_i$loglik[true_idx] >= data_i$loglik[true_idx + 1])
        cat(sprintf(
            "  Is true sigma a local maximum? %s\n",
            ifelse(is_local_max, "YES", "NO")
        ))
    }
}

cat("\n", rep("=", 70), "\n\n", sep = "")
cat("Diagnostic complete!\n\n")
