# ggdmcLikelihood

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ggdmcLikelihood)](https://cran.r-project.org/package=ggdmcLikelihood)
[![Downloads](https://cranlogs.r-pkg.org/badges/ggdmcLikelihood)](https://cran.r-project.org/package=ggdmcLikelihood)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ggdmcLikelihood/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ggdmcLikelihood/actions/workflows/R-CMD-check.yaml)

**ggdmcLikelihood 0.2.9.1 (development)** provides efficient likelihood computation for choice response time models, supporting both accuracy and response time analysis. These models, collectively known as *choice response time models*, include established frameworks such as the Diffusion Decision Model (DDM) and the Linear Ballistic Accumulation Model (LBA). The development version also extends functionality to the Cognitive Diagnostic Model (CDM), enabling applications in skill diagnosis and educational measurement. Designed for both individual- and group-level data, the package is optimised for speed and scalability. It is particularly suited for experimental psychologists and behavioural scientists analysing large datasets who require fast and accurate model estimation.
 
---

## 📦 Prerequisites
- **R** (≥ 3.5.0)  
- **Rcpp** (≥ 1.0.7)  
- **RcppArmadillo** (≥ 0.10.7.5.0)  
- **ggdmcHeaders**

---

## 📥 Installation

### Development version from GitHub

```r
# install.packages("remotes")
remotes::install_github("yxlin/ggdmcHeaders", ref = "dev")
```


## 🚀 Getting Started
This package is primarily designed to work with ggdmc and integrates seamlessly into the ggdmc workflow.

### Example: LBA Model Likelihood Computation
```r
library(ggdmcModel)
library(ggdmcLikelihood)

# Build a minimal LBA model
model <- BuildModel(
    p_map = list(
        A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1", t0 = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)

# Build DMI object from data and model
dmis <- BuildDMI(hdat, model)

# Prepare parameter sets for each subject
n_subject <- length(unique(hdat$s))
parameters <- vector("list", n_subject)
for (i in seq_len(n_subject)) {
    new_p_vector <- p_vector[model@pnames]
    parameters[[i]] <- new_p_vector
}

# Compute likelihood for a single subject
result <- compute_subject_likelihood(sub_dmis[[1]], parameters[[1]], FALSE)
single_sll <- sum(sapply(result, function(x) sum(log(x))))
print(single_sll)

# Compute likelihood for all subjects
result <- compute_likelihood(pop_dmis, parameters, FALSE)
for (i in seq_len(length(pop_dmis))) {
    sll <- sum(sapply(result[[i]], function(x) sum(log(x))))
    cat("Subject", i, "summed log-likelihood =", sll, "\n")
}
```

### Example: CDM Likelihood Computation

```r
model <- ggdmcModel::BuildModel(
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

p_vector <- c(
    guess1 = .2, guess2 = .5, guess3 = .2, guess4 = .2, guess5 = .2,
    slip1 = .1, slip2 = .8, slip3 = .1, slip4 = .1, slip5 = .1
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


sub_model <- cdModel::setCDM(model, q_matrix = Q, prior_pi = pi_uniform, rule = "DINO")

N <- 10000
dat <- cdModel::simulate(sub_model,
    nsim = N, parameter_vector = p_vector,
    nschool = 1,
    debug = FALSE, seed = 123
)


sub_dmis <- ggdmcModel::BuildDMI(dat$responses, model,
    q_matrix = Q, prior_pi = pi_uniform,
    rule = "DINO"
)


res <- compute_subject_likelihood(
    sub_dmis[[1]],
    p_vector
)
sum(log(res[[1]]))


```
## 📄 License
GPL (≥ 3)

