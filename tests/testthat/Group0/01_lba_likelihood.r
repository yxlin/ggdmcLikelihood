# q(save = "no")
cat("\n\n-------------- Testing Build new lba likelihood ---------------")
rm(list = ls())
pkg <- c("ggdmcModel", "ggdmcLikelihood", "lbaModel")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
home_dir <- "/media/yslin/Tui/01_Projects/ggdmcLikelihood"


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


sub_model <- setLBA(model)
p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)
dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)
sub_dmis <- ggdmcModel::BuildDMI(dat, model)

result <- compute_subject_likelihood(sub_dmis[[1]], p_vector)

sll <- sum(sapply(result, function(x) {
    sum(log(x))
}))
print(sll)
# result
