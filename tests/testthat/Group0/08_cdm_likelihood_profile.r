q(save = "no")
cat("\n\n-------------- Testing CDM likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcLikelihood", "cdModel", "ggdmcPrior", "ggplot2")

sapply(pkg, require, character.only = TRUE)
cat("\nWorking directory: ", getwd(), "\n")

home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood/tests/testthat/Group0"
# save_path <- paste0(wkdir, "cdm_data0.rda")

loglik_dina <- function(Y, Q, guess, slip, pi_class) {
    N <- nrow(Y)
    J <- ncol(Y)
    K <- ncol(Q)

    # guess = g_true
    # slip = s_true
    # pi_class = pi_uniform

    stopifnot(length(guess) == J, length(slip) == J)

    # Build all profiles (L × K) and eta (L × J)
    A <- as.matrix(expand.grid(rep(list(0:1), K))) # L × K

    need <- rowSums(Q)
    eta <- (A %*% t(Q)) == matrix(need, nrow = nrow(A), ncol = nrow(Q), byrow = TRUE) # L × J
    eta <- 1L * eta
    L <- nrow(A)
    stopifnot(length(pi_class) == L)

    # Class-conditional p_j(ℓ)
    p <- sweep(eta, 2, 1 - slip, `*`) + sweep(1 - eta, 2, guess, `*`) # L × J
    logp <- log(p) # L x J
    log1mp <- log(1 - p)

    # For each person i, each class ℓ: log P(y_i | ℓ)
    # Compute Y %*% logp^T + (1-Y) %*% log(1-p)^T
    # Shapes: (N×J) %*% (J×L) = N×L
    li_mat <- tcrossprod(Y, logp) + tcrossprod(1 - Y, log1mp) # N × L

    # log-sum-exp over classes with weights pi_class
    # ll_i = log(sum_ℓ pi_ℓ * exp(li_mat[i,ℓ]))
    maxrow <- apply(li_mat, 1, max)
    ll_i <- maxrow + log(rowSums(exp(li_mat - maxrow) * rep(pi_class, each = N)))

    return(ll_i) # length N
}

# ---- Core: sum log-likelihood for a parameter vector ----
sll_from_p <- function(p_vec, dmi = sub_dmis[[1]]) {
    lik <- compute_subject_likelihood(dmi, p_vec, FALSE)[[1]]
    # Guard against zeros to avoid -Inf
    sum(log(pmax(lik, .Machine$double.xmin)))
}

# Utility to make a sequence around a base value, clipped to (0,1)
around <- function(base, width = 0.15, n = 60) {
    lo <- max(1e-6, base - width)
    hi <- min(1 - 1e-6, base + width)
    seq(lo, hi, length.out = n)
}

# ---- Make a 2D surface for a (guess_i, slip_i) pair ----
make_surface <- function(par_guess, par_slip,
                         w_guess = 0.20, w_slip = 0.12, n = 60) {
    xs <- around(p_vector[[par_guess]], width = w_guess, n = n) # guess axis
    ys <- around(p_vector[[par_slip]], width = w_slip, n = n) # slip axis
    grid <- expand.grid(x = xs, y = ys)

    grid$sll <- vapply(seq_len(nrow(grid)), function(i) {
        p_try <- p_vector
        p_try[[par_guess]] <- grid$x[i]
        p_try[[par_slip]] <- grid$y[i]
        sll_from_p(p_try)
    }, numeric(1))

    grid$par_guess <- par_guess
    grid$par_slip <- par_slip
    grid$true_x <- p_vector[[par_guess]]
    grid$true_y <- p_vector[[par_slip]]
    grid$item <- sub("^guess", "", par_guess) # "1","2",...
    grid
}

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
        pi1 = "1", pi2 = "1", pi3 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
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
    guess1 = .2, guess2 = .2, guess3 = .2, guess4 = .2, guess5 = .2,
    pi1 = 0.1, pi2 = 0.2, pi3 = 0.5,
    slip1 = .1, slip2 = .1, slip3 = .1, slip4 = .1, slip5 = .1
)
pop_scale <- c(
    guess1 = .01, guess2 = .01, guess3 = .01, guess4 = .01, guess5 = .01,
    pi1 = 0.01, pi2 = 0.01, pi3 = 0.01,
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

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability
)


pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability
)


p_vector <- c(
    guess1 = .2, guess2 = .5, guess3 = .2, guess4 = .2, guess5 = .2,
    pi1 = 0.1, pi2 = 0.2, pi3 = 0.5,
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


# plot_prior(p_prior)
sub_priors <- set_priors(p_prior = p_prior)
ncell <- length(sub_model@model@cell_names)

slotNames(sub_model)

name_sorted_p_vector <- p_vector[sort(names(p_vector))]
name_sorted_p_vector
nschool <- 1
param_matrix <- t(sapply(seq_len(nschool), function(i) {
    matrix(name_sorted_p_vector, nrow = 1, ncol = sub_model@model@npar)
}))


dat <- simulate(sub_model,
    nsim = 2000, parameter_vector = p_vector,
    nschool = 1,
    debug = F, seed = 123
)
# head(dat$responses)
# names(dat)
# head(dat$alpha)


sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINA"
)

# tibble::as_tibble(dat$responses)
# head(dat$responses$C)

# Negative log-likelihood for optim
nll <- function(par_vec, dmi = sub_dmis[[1]]) {
    names(par_vec) <- names(p_vector) # keep parameter names
    -sll_from_p(par_vec, dmi)
}

slotNames(sub_dmis[[1]])

start <- p_vector * runif(length(p_vector), 0.8, 1.2) # jittered start
start

fit <- optim(
    par = start,
    fn = nll,
    method = "L-BFGS-B",
    lower = rep(1e-6, length(p_vector)),
    upper = rep(1 - 1e-6, length(p_vector))
)

fit$par # estimated parameters
-fit$value # max log-likelihood

options(digits = 3)
rbind(fit$par, p_vector)

res1 <- compute_subject_likelihood(sub_dmis[[1]], p_vector, TRUE)
sum(log(res1[[1]]))


# Build all five pairs: (guess1, slip1), ..., (guess5, slip5)
pairs_gs <- lapply(1:5, function(i) c(sprintf("guess%d", i), sprintf("slip%d", i)))
# Compute surfaces (adjust n for speed/precision)
surface_list <- lapply(pairs_gs, function(ps) make_surface(ps[1], ps[2], n = 50))
surface_df <- do.call(rbind, surface_list)





# 1) Relative SLL per item + a nice facet label
library(dplyr)

surface_df2 <- surface_df %>%
    group_by(item) %>%
    mutate(rel_sll = sll - max(sll)) %>%
    ungroup() %>%
    mutate(panel = paste0("Item ", item, " (", par_guess, " vs ", par_slip, ")"))

# 2) Best (max) point per item (break ties deterministically)
best_points <- surface_df2 %>%
    group_by(panel) %>%
    slice_max(rel_sll, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(panel, x, y)

# 3) One row per item with the true pair
true_points <- surface_df2 %>%
    distinct(panel, true_x, true_y)


# 4) Plot — map fill only on the raster layer; avoid inheriting fill elsewhere
p0 <- ggplot(surface_df2, aes(x, y)) +
    geom_raster(aes(fill = rel_sll), interpolate = TRUE) +
    geom_contour(aes(z = rel_sll), color = "white", alpha = 0.8, bins = 12) +
    scale_fill_gradientn(
        colors = c("red", "orange", "white", "skyblue", "blue"),
        guide = "colorbar"
    ) +
    geom_point(
        data = best_points, inherit.aes = FALSE,
        aes(x = x, y = y),
        color = "black", shape = 3, size = 3
    ) +
    geom_point(
        data = true_points, inherit.aes = FALSE,
        aes(x = true_x, y = true_y),
        color = "black", shape = 4, size = 3, stroke = 1.2
    ) +
    facet_wrap(~panel, scales = "free", ncol = 2) +
    labs(
        title = "Relative SLL surfaces per item: guess vs slip",
        x = "guess_i", y = "slip_i", fill = "Δ SLL (vs max)"
    ) +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_size = 16)

fn <- file.path(home_dir, "profile1.png")
png(filename = fn)
print(p0)
dev.off()


d <- data.table::data.table(surface_df)
d[x <= 0.22 & x > 0.19 & y < 0.15 & y > 0.09 & item == 1]





# ---- Profile one partner by maximizing the SLL over a grid ----
profile_pair <- function(par_guess, par_slip,
                         w_guess = 0.20, w_slip = 0.12,
                         n_g = 60, n_s = 60) {
    # Grid for guess and slip
    xs <- around(p_vector[[par_guess]], width = w_guess, n = n_g)
    ys <- around(p_vector[[par_slip]], width = w_slip, n = n_s)

    # For each guess x, profile (maximize) over slip
    sll_profile_over_slip <- vapply(xs, function(x) {
        slls_y <- vapply(ys, function(y) {
            p_try <- p_vector
            p_try[[par_guess]] <- x
            p_try[[par_slip]] <- y
            sll_from_p(p_try)
        }, numeric(1))
        max(slls_y)
    }, numeric(1))

    # For each slip y, profile (maximize) over guess
    sll_profile_over_guess <- vapply(ys, function(y) {
        slls_x <- vapply(xs, function(x) {
            p_try <- p_vector
            p_try[[par_guess]] <- x
            p_try[[par_slip]] <- y
            sll_from_p(p_try)
        }, numeric(1))
        max(slls_x)
    }, numeric(1))

    rbind(
        data.frame(
            item = sub("^guess", "", par_guess),
            varying = par_guess, value = xs,
            sll = sll_profile_over_slip,
            type = "profile_over_slip",
            base = p_vector[[par_guess]]
        ),
        data.frame(
            item = sub("^guess", "", par_guess),
            varying = par_slip, value = ys,
            sll = sll_profile_over_guess,
            type = "profile_over_guess",
            base = p_vector[[par_slip]]
        )
    )
}

profiles_list <- lapply(pairs_gs, function(ps) profile_pair(ps[1], ps[2]))
profiles_df <- do.call(rbind, profiles_list)

# Plot profiled curves, one facet per item, color by which partner was profiled out
p1 <- ggplot(profiles_df, aes(value, sll, color = type)) +
    geom_line() +
    geom_vline(aes(xintercept = base), linetype = 2) +
    facet_wrap(~ paste0("Item ", item, " — varying: ", varying),
        scales = "free_x", ncol = 2
    ) +
    labs(
        title = "Profiled SLL: varying one of (guess_i, slip_i), profiling out the other",
        x = "Parameter value", y = "Sum log-likelihood", color = "Profile type"
    ) +
    theme_bw(base_size = 16)

# Parameters to profile (edit as you wish)
params_to_profile <- names(p_vector) # e.g., c("guess1","slip1","guess2",...)

library(ggplot2)

profile_df_list <- lapply(params_to_profile, function(par) {
    base <- p_vector[par]
    xs <- around(base, width = if (grepl("^guess", par)) 0.20 else 0.12, n = 60)

    slls <- vapply(xs, function(x) {
        p_try <- p_vector
        p_try[par] <- x
        sll_from_p(p_try)
    }, numeric(1))

    data.frame(parameter = par, value = xs, sll = slls, base = as.numeric(base))
})

profile_df <- do.call(rbind, profile_df_list)



# Plot: one line per parameter, faceted
d <- data.table::data.table(profile_df)


p2 <- ggplot(profile_df, aes(value, sll)) +
    geom_line() +
    geom_vline(aes(xintercept = base), linetype = 2) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(
        title = "Profile log-likelihoods (one parameter at a time)",
        x = "Parameter value", y = "Sum log-likelihood"
    ) +
    theme_bw(base_size = 18)


fn <- file.path(home_dir, "profile2.png")
png(filename = fn)
print(p2)
dev.off()

par1 <- "guess1"
par2 <- "slip1"

# Grids around the true values (adjust widths if needed)
xs <- around(p_vector[[par1]], width = 0.20, n = 60)
ys <- around(p_vector[[par2]], width = 0.12, n = 60)

grid <- expand.grid(x = xs, y = ys)

# Evaluate SLL on the grid
grid$sll <- vapply(seq_len(nrow(grid)), function(i) {
    p_try <- p_vector
    p_try[[par1]] <- grid$x[i]
    p_try[[par2]] <- grid$y[i]
    sll_from_p(p_try)
}, numeric(1))

# Heatmap + contours
p3 <- ggplot(grid, aes(x, y, fill = sll)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = sll), color = "white", alpha = 0.6, bins = 12) +
    geom_point(aes(x = p_vector[[par1]], y = p_vector[[par2]]),
        color = "black", size = 2, shape = 4, stroke = 1.2
    ) +
    labs(
        title = sprintf("Likelihood surface: %s vs %s", par1, par2),
        x = par1, y = par2, fill = "SLL"
    ) +
    coord_cartesian(expand = FALSE)


pairs <- list(c("guess1", "slip1"), c("guess2", "slip2")) # add more pairs

make_surface <- function(par1, par2, w1 = 0.20, w2 = 0.12, n = 60) {
    xs <- around(p_vector[[par1]], width = w1, n = n)
    ys <- around(p_vector[[par2]], width = w2, n = n)
    grid <- expand.grid(x = xs, y = ys)
    grid$sll <- vapply(seq_len(nrow(grid)), function(i) {
        p_try <- p_vector
        p_try[[par1]] <- grid$x[i]
        p_try[[par2]] <- grid$y[i]
        sll_from_p(p_try)
    }, numeric(1))
    grid$par1 <- par1
    grid$par2 <- par2
    grid
}

surface_list <- lapply(pairs, function(p) make_surface(p[1], p[2]))
surface_df <- do.call(rbind, surface_list)

# Add true values to the surface_df
surface_df$true_x <- p_vector[surface_df$par1]
surface_df$true_y <- p_vector[surface_df$par2]

p4 <- ggplot(surface_df, aes(x, y, fill = sll)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = sll), color = "white", alpha = 0.6, bins = 10) +
    geom_point(aes(x = true_x, y = true_y),
        color = "black", shape = 4, size = 3, stroke = 1.2
    ) + # shape=4 = X mark
    facet_wrap(~ par1 + par2, scales = "free", ncol = 2) +
    labs(
        title = "Likelihood surfaces across parameter pairs",
        x = "Parameter 1", y = "Parameter 2", fill = "SLL"
    ) +
    coord_cartesian(expand = FALSE)


p5 <- ggplot(surface_df, aes(x, y, fill = sll)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = sll), color = "white", alpha = 0.6, bins = 10) +
    facet_wrap(~ par1 + par2, scales = "free", ncol = 2) +
    labs(
        title = "Likelihood surfaces across parameter pairs",
        x = "Parameter 1", y = "Parameter 2", fill = "SLL"
    ) +
    coord_cartesian(expand = FALSE)
