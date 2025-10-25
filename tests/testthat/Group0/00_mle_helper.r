# Negative log-likelihood for optim
nll <- function(par_vec, dmi = sub_dmis[[1]]) {
    names(par_vec) <- names(p_vector) # keep parameter names
    -sll_from_p(par_vec, dmi)
}

# ---- Core: sum log-likelihood for a parameter vector ----
sll_from_p <- function(p_vec, dmi = sub_dmis[[1]]) {
    lik <- ggdmcLikelihood::compute_subject_likelihood(dmi, p_vec, FALSE)[[1]]
    # Guard against zeros to avoid -Inf
    sum(log(pmax(lik, .Machine$double.xmin)))
}
