# q(save = "no")
cat("\n\n-------------- Testing CDM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "cdModel", "ggdmcPrior")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmcLikelihood/tests/testthat"
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
    rule = "DINA",
    use_mvn = TRUE
)


cat("sub_model = ", sub_model@use_mvn, "\n")

sigma <- 0
means <- c(0, 0)
Sigma <- matrix(c(1, sigma, sigma, 1), ncol = 2)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)


sim_p_vector <- c(
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .2,
    mean1 = means[1], mean2 = means[2], sigma = sigma,
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)

N <- 50
dat <- simulate(sub_model,
    nsim = N, parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE, seed = 123
)



model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = c(pi1 = 0.25, pi2 = 0.25, pi3 = 0.25),
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)


true_p_vector <- c(
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .2,
    p1 = skill_probs$probability[1], p2 = skill_probs$probability[2],
    p3 = skill_probs$probability[3],
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)

res <- ggdmcLikelihood::compute_subject_likelihood(
    sub_dmis[[1]],
    true_p_vector,
    debug = TRUE
)

expected <- -150.9091
print(sum(log(res[[1]])))
testthat::expect_equal(sum(log(res[[1]])), expected, tolerance = 1e-6)
