# q(save = "no")
cat("\n\n-------------- Testing Build new likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmcModel", "ggdmcLikelihood")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data0.rda"
load(fn)

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


sll <- sum(sapply(result, function(x) {
    sum(log(x))
}))

testthat::expect_true(all.equal(sll, -137.7996, tolerance = 1e-6))


result <- compute_likelihood(pop_dmis, parameters, F)
n_subject <- length(pop_dmis)
for (i in seq_len(n_subject)) {
    res <- result[[i]]
    sll <- sum(sapply(res, function(x) {
        sum(log(x))
    }))
    cat("Subject ", i, "results in summed log likelihood = ", sll, "\n")
}

# testthat::expect_true(all.equal(sll, -144.1973, tolerance = 1e-6))
