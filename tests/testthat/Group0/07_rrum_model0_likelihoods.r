# q(save = "no")
cat("\n\n-------------- Testing RRUM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcLikelihood", "cdModel", "ggdmcPrior")

sapply(pkg, require, character.only = TRUE)
cat("\nWorking directory: ", getwd(), "\n")
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmcLikelihood/tests/testthat/Group0"

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
), ncol = 2, byrow = TRUE)
K <- nrow(Q)
colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- paste0("Item", seq_len(K))
p_map <- cdModel::generate_rrum_pmap(Q)

#  "pi_00"    "pi_01"    "pi_10"    "pi_item1" "pi_item2" "pi_item3"
#  "pi_item4" "pi_item5" "r11"      "r22"      "r31"      "r32"
#  "r41"      "r52"

# print_rrum_pmap(p_map, Q)
model <- BuildModel(
    p_map = p_map,
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

# Set RRUM rule without multivariate normal
use_mvn <- FALSE
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "RRUM",
    use_mvn = use_mvn
)

p_vector <- c(
    pi_00 = .1, pi_01 = .2, pi_10 = .3,
    pi_item1 = .1, pi_item2 = .2, pi_item3 = .3, pi_item4 = .4,
    pi_item5 = .5,
    r11 = .1,
    r22 = .4,
    r31 = .3, r32 = .6,
    r41 = .4,
    r52 = .66
)

# use_mvn <- FALSE
dat <- simulate(sub_model,
    nsim = 10000, parameter_vector = p_vector,
    nschool = 1,
    debug = FALSE, seed = 123
)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "RRUM"
    # use_mvn = use_mvn
)

# ---- Core: sum log-likelihood for a parameter vector ----
sll_from_p <- function(p_vec, dmi = sub_dmis[[1]]) {
    lik <- compute_subject_likelihood(dmi, p_vec, FALSE)[[1]]
    # Guard against zeros to avoid -Inf
    sum(log(pmax(lik, .Machine$double.xmin)))
}

nll <- function(par_vec, dmi = sub_dmis[[1]]) {
    names(par_vec) <- names(p_vector) # keep parameter names
    -sll_from_p(par_vec, dmi)
}


npar <- length(p_vector)
start <- p_vector * runif(npar, 0.8, 1.2) # jittered start
round(start, 2)



# print(sub_dmis[[1]]@use_mvn)
lik <- compute_subject_likelihood(sub_dmis[[1]], start, TRUE)
sum(log(lik[[1]]))

fit <- optim(
    par = start,
    fn = nll,
    method = "L-BFGS-B",
    lower = rep(1e-6, npar),
    upper = rep(1 - 1e-6, npar)
)


fit$par # estimated parameters
-fit$value # max log-likelihood

options(digits = 3)
rbind(fit$par, p_vector)


# mod1 <- CDM::gdina(data = data.frame(wide[, -1]), q.matrix = Q, rule = "RRUM")
# summary(mod1)
# x <- mod1
# # x <- CDM::gdina(..., rule = "RRUM")
# tab <- x$itempar

# by_item <- split(tab, tab$itemno)

# pi_rrum <- sapply(by_item, function(df) {
#     delta0 <- df$est[df$partype == 0]
#     deltas <- df$est[df$partype > 0]
#     exp(delta0 + sum(deltas))
# })

# r_rrum <- lapply(by_item, function(df) {
#     use <- df$partype > 0
#     setNames(exp(-df$est[use]), df$partype.attr[use])
# })



# # ---- True item parameters (guessing g_j, slipping s_j)
# pop_mean <- c(
#     pi1 = .1, pi2 = .2, pi3 = .3,
#     pi_item1 = .1, pi_item2 = .2, pi_item3 = .3, pi_item4 = .4,
#     pi_item5 = .5,
#     r11 = .1,
#     r22 = .4,
#     r31 = .3, r32 = .6,
#     r41 = .4,
#     r52 = .66
# )
# pop_scale <- c(
#     pi1 = .01, pi2 = .02, pi3 = .03,
#     pi_item1 = .01, pi_item2 = .02, pi_item3 = .03,
#     pi_item4 = .01, pi_item5 = .02,
#     r11 = .01,
#     r22 = .04,
#     r31 = .03, r32 = .06,
#     r41 = .04,
#     r52 = .06
# )


# pop_dist <- BuildPrior(
#     p0 = pop_mean,
#     p1 = pop_scale,
#     lower = rep(0, model@npar),
#     upper = rep(NA, model@npar),
#     dists = rep("tnorm", model@npar),
#     log_p = rep(F, model@npar)
# )



# pop_model <- setCDM(model,
#     population_distribution = pop_dist,
#     q_matrix = model@cdm_info$q_matrix,
#     profile_probability = model@cdm_info$profile_probability,
#     rule = "RRUM"
# )

# p_vector <- c(
#     pi1 = .1, pi2 = .2, pi3 = .3,
#     pi_item1 = .1, pi_item2 = .2, pi_item3 = .3, pi_item4 = .4,
#     pi_item5 = .5,
#     r11 = .1,
#     r22 = .4,
#     r31 = .3, r32 = .6,
#     r41 = .4,
#     r52 = .66
# )


# p0 <- rep(0, model@npar)
# names(p0) <- model@pnames

# p_prior <- BuildPrior(
#     p0 = p0,
#     p1 = rep(1, model@npar),
#     lower = rep(0, model@npar),
#     upper = rep(NA, model@npar),
#     dist = rep("unif", model@npar),
#     log_p = rep(TRUE, model@npar)
# )

# sub_priors <- set_priors(p_prior = p_prior)
# # plot_prior(p_prior)

# # safest: use stringi to force C collation just for this call
# ascii_sorted <- stringi::stri_sort(names(p_vector), locale = "C")
# ncell <- length(sub_model@model@cell_names)
# name_sorted_p_vector <- p_vector[ascii_sorted]

# nschool <- 1
# param_matrix <- t(sapply(seq_len(nschool), function(i) {
#     matrix(name_sorted_p_vector, nrow = 1, ncol = sub_model@model@npar)
# }))




# # Negative log-likelihood for optim
