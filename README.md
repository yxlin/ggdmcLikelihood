# ggdmcLikelihood

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ggdmcLikelihood)](https://cran.r-project.org/package=ggdmcLikelihood)
[![Downloads](https://cranlogs.r-pkg.org/badges/ggdmcLikelihood)](https://cran.r-project.org/package=ggdmcLikelihood)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ggdmcLikelihood/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ggdmcLikelihood/actions/workflows/R-CMD-check.yaml)


**ggdmcLikelihood** provides efficient likelihood computation for choice response time models, supporting both accuracy and response time analysis.  
This class of models, known as *choice response time modelling*, includes well-known frameworks such as the Diffusion Decision Model (DDM) and the Linear Ballistic Accumulator (LBA).  

The package works with data from individuals or groups, and is optimised for speed and scalability.  
It is particularly suited for experimental psychologists and behavioural scientists who analyse large datasets and require fast, accurate model estimation.

---

## 📦 Prerequisites
- **R** (≥ 3.5.0)  
- **Rcpp** (≥ 1.0.7)  
- **RcppArmadillo** (≥ 0.10.7.5.0)  
- **ggdmcHeaders**

---

## 📥 Installation

From CRAN:
```r
install.packages("ggdmcLikelihood")
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

## 📄 License
GPL (≥ 3)

