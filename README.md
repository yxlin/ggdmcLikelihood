<<<<<<< HEAD
# ggdmcLikelihood
ggdmcLikelihood provides likelihood computation functions for the ggdmc package. It supports both single-subject and multi-subject likelihood computations.

# Getting Started
The package is mainly to support ggdmc, so you can use it together with other ggdmc supporting packages.
=======
This package provides a user interface for efficiently computing likelihoods in design-based choice response time models. It is built to integrate with the 'ggdmc' package, enabling rapid initialisation of starting samples for genetic algorithms, which typically use a large number of Markov chains. The likelihood computation supports both single-subject and multi-subject models. Implemented methods include fast evaluation of likelihoods across trial-level data, support for various response time distributions (e.g., diffusion models). Ideal for cognitive scientists and decision modellers, this tool enhances the usability and computational performance of complex model fitting tasks in experimental psychology and behavioral data research.

# Getting Started
The package is mainly to support the 'ggdmc' package. Although you could use it in other situations, it works most seamlessly with the 'ggdmc' ecosystem.
>>>>>>> 1804e22 (reconnect to the hub, after moving the folder location)

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