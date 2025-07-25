This package helps researchers analyse how people make decisions over time by estimating models of choice and reaction time. It supports a variety of response time models, including the well-known Decision Diffusion Model, and works with data from both individuals and groups. Designed for speed and scalability, it makes it easier to fit complex models—especially when using advanced methods like genetic sampling with many Markov chains. The package is especially useful for experimental psychologists and behavioral scientists who work with large datasets and need fast, accurate model estimation.

# Getting Started
This package is mainly intended to work alongside 'ggdmc', and it integrates most smoothly when used as part of the ggdmc workflow. You can use it on its own, 
but it’s optimised for that ecosystem.

```
model <- ggdmcModel::BuildModel(
    p_map = list(
        A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1",
        t0 = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)

dmis <- ggdmcModel::BuildDMI(hdat, model)

nsubject <- length(unique(hdat$s))
parameters <- list()
for (i in seq_len(nsubject)) {
    new_p_vector <- p_vector[model@pnames]
    parameters[[i]] <- new_p_vector
}

result <- compute_subject_likelihood(sub_dmis[[1]], parameters[[1]], F)
result

sll <- sum(sapply(result, function(x) {
    sum(log(x))
}))
print(sll)

result <- compute_likelihood(pop_dmis, parameters, F)
n_subject <- length(pop_dmis)
for (i in seq_len(n_subject)) {
    res <- result[[i]]
    sll <- sum(sapply(res, function(x) {
        sum(log(x))
    }))
    cat("Subject ", i, "results in summed log likelihood = ", sll, "\n")
}


```

# Prerequisites
R (>= 3.5.0), Rcpp (>= 1.0.7), RcppArmadillo (>= 0.10.7.5.0), and ggdmcHeaders.

# Installation

From CRAN:
```
install.packages("ggdmcLikelihood")
```