q(save = "no")
cat("\n\n-------------- Testing RRUM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcLikelihood", "cdModel", "ggdmcPrior")

sapply(pkg, require, character.only = TRUE)
cat("\nWorking directory: ", getwd(), "\n")
long2wide <- function(df) {
    require(dplyr)
    require(tidyr)
    require(stringr)
    wide <- df %>%
        # keep only what we need
        select(student, item, C) %>%
        # make sure ids and items are clean, ordered, and item names look like E1..E28
        mutate(
            id   = as.integer(as.character(student)),
            item = paste0("E", as.integer(as.character(item))),
            C    = as.integer(C)
        ) %>%
        select(id, item, C) %>%
        # if there could be duplicates per (id,item), pick one; change to mean/max if needed
        distinct(id, item, .keep_all = TRUE) %>%
        # ensure a full rectangular layout (so missing item administrations show up as NA)
        complete(id, item) %>%
        # pivot to wide
        pivot_wider(names_from = item, values_from = C, values_fill = NA_integer_) %>%
        arrange(id)

    wide
}

home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat/Group0"

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
        pi1 = "1", pi2 = "1", pi3 = "1",
        pi_item1 = "1", pi_item2 = "1", pi_item3 = "1", pi_item4 = "1",
        pi_item5 = "1",
        r11 = "1", r22 = "1",
        r31 = "1", r32 = "1",
        r41 = "1",
        r52 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)


# ---- True item parameters (guessing g_j, slipping s_j)
pop_mean <- c(
    pi1 = .1, pi2 = .2, pi3 = .3,
    pi_item1 = .1, pi_item2 = .2, pi_item3 = .3, pi_item4 = .4,
    pi_item5 = .5,
    r11 = .1,
    r22 = .4,
    r31 = .3, r32 = .6,
    r41 = .4,
    r52 = .66
)
pop_scale <- c(
    pi1 = .01, pi2 = .02, pi3 = .03,
    pi_item1 = .01, pi_item2 = .02, pi_item3 = .03,
    pi_item4 = .01, pi_item5 = .02,
    r11 = .01,
    r22 = .04,
    r31 = .03, r32 = .06,
    r41 = .04,
    r52 = .06
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
    rule = "RRUM"
)


pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "RRUM"
)

p_vector <- c(
    pi1 = .1, pi2 = .2, pi3 = .3,
    pi_item1 = .1, pi_item2 = .2, pi_item3 = .3, pi_item4 = .4,
    pi_item5 = .5,
    r11 = .1,
    r22 = .4,
    r31 = .3, r32 = .6,
    r41 = .4,
    r52 = .66
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
ascii_sorted <- stringi::stri_sort(names(p_vector), locale = "C")
ncell <- length(sub_model@model@cell_names)
name_sorted_p_vector <- p_vector[ascii_sorted]

nschool <- 1
param_matrix <- t(sapply(seq_len(nschool), function(i) {
    matrix(name_sorted_p_vector, nrow = 1, ncol = sub_model@model@npar)
}))



dat <- simulate(sub_model,
    nsim = 10000, parameter_vector = name_sorted_p_vector,
    nschool = 1,
    debug = F, seed = 123
)


sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "RRUM"
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



npar <- length(name_sorted_p_vector)
start <- name_sorted_p_vector * runif(npar, 0.8, 1.2) # jittered start
round(start, 2)

# lik <- compute_subject_likelihood(sub_dmis[[1]], start, TRUE)

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
