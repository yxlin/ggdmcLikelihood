# q(save = "no")
cat("\n\n-------------- Testing CDM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcLikelihood", "cdModel", "ggdmcPrior", "ggplot2")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

home_dir <- "/media/yslin/Tui/01_Project/ggdmcLikelihood/tests/testthat/"
# save_path <- paste0(home_dir, "cdm_data0.rda")

# Negative log-likelihood for optim
nll <- function(par_vec, dmi = sub_dmis[[1]]) {
    names(par_vec) <- names(p_vector) # keep parameter names
    -sll_from_p(par_vec, dmi)
}

# ---- Core: sum log-likelihood for a parameter vector ----
sll_from_p <- function(p_vec, dmi = sub_dmis[[1]]) {
    lik <- compute_subject_likelihood(dmi, p_vec, FALSE)[[1]]
    # Guard against zeros to avoid -Inf
    sum(log(pmax(lik, .Machine$double.xmin)))
}

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = NULL,
    type = "cdm",
    verbose = TRUE
)

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")

n_item <- nrow(Q)
n_skill <- ncol(Q)
n_profile <- 2^(n_skill)
pi_uniform <- rep(1 / n_profile, n_profile)

# ---- True item parameters (guessing g_j, slipping s_j)
pop_mean <- c(
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .2,
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)
pop_scale <- c(
    guess1 = .01, guess2 = .01, guess3 = .01, guess4 = .01, guess5 = .01,
    slip1 = .01, slip2 = .01, slip3 = .01, slip4 = .01, slip5 = .01
)

pop_dist <- BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

sub_model <- setCDM(model, q_matrix = Q, prior_pi = pi_uniform, rule = "DINO")
pop_model <- setCDM(model,
    population_distribution = pop_dist, q_matrix = Q,
    prior_pi = pi_uniform, rule = "DINO"
)


p_vector <- c(
    guess1 = .2, guess2 = .5, guess3 = .2, guess4 = .2, guess5 = .2,
    slip1 = .1, slip2 = .8, slip3 = .1, slip4 = .1, slip5 = .1
)


p0 <- rep(0, model@npar)
names(p0) <- model@pnames

p_prior <- BuildPrior(
    p0 = p0,
    p1 = rep(1, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
sub_priors <- set_priors(p_prior = p_prior)


# N <- 1000000
N <- 10
dat <- simulate(sub_model,
    nsim = N, parameter_vector = p_vector,
    nschool = 1,
    debug = F, seed = 123
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = Q, prior_pi = pi_uniform,
    rule = "DINO"
)


start <- p_vector * runif(length(p_vector), 0.8, 1.2) # jittered start

fit <- optim(
    par = start,
    fn = nll,
    method = "L-BFGS-B",
    lower = rep(1e-6, length(p_vector)),
    upper = rep(1 - 1e-6, length(p_vector))
)

# fit$par # estimated parameters

res <- rbind(fit$par, p_vector, fit$par - p_vector)
row.names(res) <- c("Estimate", "True", "Diff")
options(digits = 3)
print(res)
print(-fit$value)

res <- compute_subject_likelihood(
    sub_dmis[[1]],
    p_vector
)
sum(log(res[[1]]))

# 100
#          guess1 guess2  guess3  guess4  guess5 slip1  slip2  slip3  slip4
# Estimate  0.257 0.5972 0.20761  0.1307  0.0635 0.203 0.8135 0.1695  1e-06
# True      0.200 0.5000 0.20000  0.2000  0.2000 0.100 0.8000 0.1000  1e-01
# Diff      0.057 0.0972 0.00761 -0.0693 -0.1365 0.103 0.0135 0.0695 -1e-01
#             slip5
# Estimate  0.00678
# True      0.10000
# Diff     -0.09322


# 3000
#          guess1  guess2  guess3 guess4  guess5   slip1   slip2   slip3   slip4
# Estimate 0.2143  0.4729  0.1925 0.2267 0.20706  0.0807  0.7893  0.0596 0.10355
# True     0.2000  0.5000  0.2000 0.2000 0.20000  0.1000  0.8000  0.1000 0.10000
# Diff     0.0143 -0.0271 -0.0075 0.0267 0.00706 -0.0193 -0.0107 -0.0404 0.00355
#             slip5
# Estimate 0.100749
# True     0.100000
# Diff     0.000749

# 30000
#           guess1  guess2  guess3  guess4   guess5    slip1   slip2    slip3
# Estimate 0.20557  0.4897 0.20271 0.20381  0.19882  0.09432 0.80177  0.09831
# True     0.20000  0.5000 0.20000 0.20000  0.20000  0.10000 0.80000  0.10000
# Diff     0.00557 -0.0103 0.00271 0.00381 -0.00118 -0.00568 0.00177 -0.00169
#             slip4   slip5
# Estimate  0.09895 0.10614
# True      0.10000 0.10000
# Diff     -0.00105 0.00614
# [1] -93580

# 80000
#   guess1   guess2  guess3   guess4   guess5    slip1   slip2   slip3
# Estimate  0.19578  0.49669  0.1987  0.19832  0.19694  0.09572 0.80291 0.10313
# True      0.20000  0.50000  0.2000  0.20000  0.20000  0.10000 0.80000 0.10000
# Diff     -0.00422 -0.00331 -0.0013 -0.00168 -0.00306 -0.00428 0.00291 0.00313
#             slip4     slip5
# Estimate 0.100192  0.099628
# True     0.100000  0.100000
# Diff     0.000192 -0.000372
# [1] -248946

# 100000
#       guess1    guess2 guess3   guess4    guess5    slip1   slip2
# Estimate  0.19887  0.499869  2e-01  0.19979  0.199714 1.00e-01 0.80102
# True      0.20000  0.500000  2e-01  0.20000  0.200000 1.00e-01 0.80000
# Diff     -0.00113 -0.000131  4e-04 -0.00021 -0.000286 8.77e-05 0.00102
#              slip3    slip4     slip5
# Estimate  0.099445 0.100929  0.099886
# True      0.100000 0.100000  0.100000
# Diff     -0.000555 0.000929 -0.000114
