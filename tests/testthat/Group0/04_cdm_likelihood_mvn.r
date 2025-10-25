#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n-------------- Testing CDM likelihood ---------------")
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
    verbose = TRUE
)
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINO",
    use_mvn = TRUE
)
means <- c(0.5, 0.2)
# sigma <- 0.2
sigma <- .1

p_vector <- c(
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .2,
    mean1 = means[1], mean2 = means[2], sigma = sigma,
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)


N <- 100
dat <- simulate(sub_model,
    nsim = N, parameter_vector = p_vector,
    nschool = 1,
    debug = FALSE, seed = 123
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINO",
    use_mvn = TRUE
)
sub_dmis[[1]]@use_mvn

sigma <- 1
p_vector <- c(
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .5,
    mean1 = means[1], mean2 = means[2], sigma = sigma,
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)
res <- ggdmcLikelihood::compute_subject_likelihood(
    sub_dmis[[1]],
    p_vector,
    debug = FALSE
)

(result <- sum(log(res[[1]])))

# expected <- -286.5494
# testthat::expect_equal(result, expected, tolerance = 1e-4)
# -286.5494 - (-287.0394)
# -286.5494 - (-287.3693)
# -286.5494 - (-287.6752)
# -286.5494 - (-286.5193)
# -286.5494 - (-290.079)
