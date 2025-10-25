q(save = "no")
cat("\n\n-------------- Testing RRUM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcLikelihood", "cdModel", "ggdmcPrior")

sapply(pkg, require, character.only = TRUE)
cat("\nWorking directory: ", getwd(), "\n")
home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat/Group0"

Q <- matrix(c(
    1, 0, 0,
    0, 1, 1,
    1, 1, 0
), ncol = 3, byrow = TRUE)
colnames(Q) <- c("A1", "A2", "A3")
rownames(Q) <- c("Item1", "Item2", "Item3")

model <- BuildModel(
    p_map = list(
        beta1__1 = "1", beta1__A1 = "1",
        beta2__1 = "1", beta2__A2 = "1", beta2__A3 = "1",
        beta2__A2xA3 = "1",
        beta3__1 = "1",
        beta3__A1 = "1",
        beta3__A2 = "1",
        beta3__A1xA2 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1", pi4 = "1",
        pi5 = "1", pi6 = "1", pi7 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL, # Must enter NULL
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)


# ---- True item parameters (guessing g_j, slipping s_j)

pop_mean <- c(
    beta1__1 = .1, beta1__A1 = .2,
    beta2__1 = .3, beta2__A2 = .4, beta2__A3 = .5,
    beta2__A2xA3 = .6,
    beta3__1 = .1,
    beta3__A1 = .2,
    beta3__A2 = .3,
    beta3__A1xA2 = .4,
    pi1 = .05, pi2 = .1, pi3 = .06, pi4 = .2,
    pi5 = .07, pi6 = .3, pi7 = .08
)
# sort(names(p_vector))
pop_scale <- c(
    beta1__1 = .01, beta1__A1 = .02,
    beta2__1 = .03, beta2__A2 = .04, beta2__A3 = .05,
    beta2__A2xA3 = .06,
    beta3__1 = .01,
    beta3__A1 = .02,
    beta3__A2 = .03,
    beta3__A1xA2 = .04,
    pi1 = .05, pi2 = .01, pi3 = .06, pi4 = .02,
    pi5 = .07, pi6 = .03, pi7 = .08
)


pop_dist <- BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "LCDM"
)


pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "LCDM"
)


p_vector <- c(
    beta1__1 = .1, beta1__A1 = .2,
    beta2__1 = .3, beta2__A2 = .4, beta2__A3 = .5,
    beta2__A2xA3 = .6,
    beta3__1 = .1,
    beta3__A1 = .2,
    beta3__A2 = .3,
    beta3__A1xA2 = .4,
    pi1 = .05, pi2 = .1, pi3 = .06, pi4 = .2,
    pi5 = .07, pi6 = .3, pi7 = .08
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
# plot_prior(p_prior)

# safest: use stringi to force C collation just for this call
# ascii_sorted <- stringi::stri_sort(names(p_vector), locale = "C")
# ncell <- length(sub_model@model@cell_names)
# name_sorted_p_vector <- p_vector[ascii_sorted]

# nschool <- 1
# param_matrix <- t(sapply(seq_len(nschool), function(i) {
#     matrix(name_sorted_p_vector, nrow = 1, ncol = sub_model@model@npar)
# }))

N_total_student <- 1000
dat <- simulate(sub_model,
    nsim = N_total_student, parameter_vector = p_vector,
    nschool = 1,
    seed = 123,
    debug = F
)


sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "LCDM"
)



# Negative log-likelihood for optim
# ---- Core: sum log-likelihood for a parameter vector ----
sll_from_p <- function(p_vec, dmi = sub_dmis[[1]]) {
    lik <- compute_subject_likelihood(dmi, p_vec, FALSE)[[1]]
    # Guard against zeros to avoid -Inf
    sum(log(pmax(lik, .Machine$double.xmin)))
}

nll <- function(par_vec, dmi = sub_dmis[[1]]) {
    names(par_vec) <- names(name_sorted_p_vector) # keep parameter names
    -sll_from_p(par_vec, dmi)
}


npar <- model@npar
# npar <- length(name_sorted_p_vector)
start <- p_vector * runif(npar, 0.8, 1.2) # jittered start
round(start, 2)

lik <- compute_subject_likelihood(sub_dmis[[1]], start, TRUE)

fit <- optim(
    par = start,
    fn = nll,
    method = "L-BFGS-B",
    lower = rep(1e-6, npar),
    upper = rep(1 - 1e-6, npar)
)


# fit$par # estimated parameters
# -fit$value # max log-likelihood

options(digits = 3)
rbind(fit$par, name_sorted_p_vector)
#                        pi1   pi2   pi3 pi_item1 pi_item2 pi_item3 pi_item4
#                      0.236 0.315 0.449    0.227     0.34    0.302     0.38
# name_sorted_p_vector 0.100 0.200 0.300    0.100     0.20    0.300     0.40
#                      pi_item5   r11   r22   r31   r32   r41   r52
#                         0.427 0.263 0.556 0.246 0.448 0.379 0.575
# name_sorted_p_vector    0.500 0.100 0.400 0.300 0.600 0.400 0.660
# cdm_est <- c(0.6166, 0.1627, 0.1813, 0.0394)
# sum(cdm_est)

res1 <- compute_subject_likelihood(sub_dmis[[1]], p_vector, TRUE)
sum(log(res1[[1]]))

head(dat$responses)

long <- tibble::as_tibble(dat$responses)
wide <- long2wide(long)


mod1 <- CDM::gdina(data = data.frame(wide[, -1]), q.matrix = Q, rule = "RRUM")
summary(mod1)
x <- mod1
# x <- CDM::gdina(..., rule = "RRUM")
tab <- x$itempar

by_item <- split(tab, tab$itemno)

pi_rrum <- sapply(by_item, function(df) {
    delta0 <- df$est[df$partype == 0]
    deltas <- df$est[df$partype > 0]
    exp(delta0 + sum(deltas))
})

r_rrum <- lapply(by_item, function(df) {
    use <- df$partype > 0
    setNames(exp(-df$est[use]), df$partype.attr[use])
})
