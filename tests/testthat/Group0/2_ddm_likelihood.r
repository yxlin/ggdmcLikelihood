q(save = "no")
cat("\n\n-------------- Testing Build new likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmcModel", "ggdmcLikelihood")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ddm_data0.rda"
load(fn)

# DDM
model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1", t0 = "1",
        st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "fastdm"
)

sub_model <- ddModel::setDDM(model)

p_vector <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.2, z = .38)
dat <- ddModel::simulate(sub_model, nsim = 32, parameter_vector = p_vector, n_subject = 1)
sub_dmis <- ggdmcModel::BuildDMI(dat, model)

dat
options(digits = 5)
p_vector <- c(a = 1.5, sz = 0.25, t0 = 0.15, v = 2.2, z = .38)
p_vector <- c(a = 1.5, sz = 0.25, t0 = 0.15, v = 3.2, z = .38)
p_vector <- sub_samples@theta[, 3, 1]
p_vector
result <- compute_subject_likelihood(sub_dmis[[1]], p_vector, FALSE)
sum(sapply(result, function(x) {
    sum(log(x))
}))

nmc <- 3
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 846671, verbose = T)
sub_samples@theta[, , 1]
sub_samples@log_likelihoods[, 1]
