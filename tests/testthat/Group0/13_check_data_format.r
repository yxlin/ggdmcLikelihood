#!/usr/bin/env Rscript
# Check what format the package data is in

cat("\n\n-------------- Data Format Check ---------------\n")
rm(list = ls())

library(ggdmc)
library(ggdmcModel)
library(cdModel)

Q <- matrix(c(
  1, 0,
  0, 1,
  1, 1,
  1, 0,
  0, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")
rownames(Q) <- paste0("Item", 1:5)

true_means <- c(0.5, 0.2)
true_sigma <- 0.2
true_guess <- c(.1, .2, .3, .4, .5)
true_slip <- c(.2, .4, .6, .8, .9)

true_p_vector <- c(
  guess1 = true_guess[1], guess2 = true_guess[2], guess3 = true_guess[3],
  guess4 = true_guess[4], guess5 = true_guess[5],
  mean1 = true_means[1], mean2 = true_means[2], sigma = true_sigma,
  slip1 = true_slip[1], slip2 = true_slip[2], slip3 = true_slip[3],
  slip4 = true_slip[4], slip5 = true_slip[5]
)

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        mean1 = "1", mean2 = "1", sigma = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = FALSE
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = TRUE
)

# Simulate
set.seed(123)
N <- 10

dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = true_p_vector,
    nschool = 1,
    debug = FALSE
)

cat("Data structure:\n")
cat("  class(dat):", class(dat), "\n")
cat("  names(dat):", names(dat), "\n\n")

cat("dat$responses:\n")
cat("  class:", class(dat$responses), "\n")
cat("  dim:", dim(dat$responses), "\n")
cat("  colnames:", colnames(dat$responses), "\n\n")

cat("First 10 rows:\n")
print(head(dat$responses, 10))

cat("\n\ndat$data_matrix:\n")
if ("data_matrix" %in% names(dat)) {
  cat("  class:", class(dat$data_matrix), "\n")
  cat("  dim:", dim(dat$data_matrix), "\n")
  print(head(dat$data_matrix, 10))
} else {
  cat("  (does not exist)\n")
}

cat("\n\nOther components:\n")
for (name in names(dat)) {
  if (name != "responses" && name != "data_matrix") {
    cat(sprintf("  %s: class=%s, length/dim=%s\n",
               name, class(dat[[name]])[1],
               paste(dim(dat[[name]]), collapse="x")))
  }
}

# Check if responses are 0/1 or some other format
cat("\n\nUnique values in responses:\n")
print(table(as.vector(dat$responses)))

# Convert to matrix
Y_matrix <- as.matrix(dat$responses)
cat("\n\nConverted to matrix:\n")
cat("  class:", class(Y_matrix), "\n")
cat("  dim:", dim(Y_matrix), "\n")
print(head(Y_matrix, 10))

cat("\n")
