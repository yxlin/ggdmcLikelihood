Sys.setenv("R_TESTS" = "")
## Workaround for the error,
## "cannot open file 'startup.Rs': No such file or directory" in Windows 10

library(testthat)
library(ggdmcLikelihood)
library(ggdmcModel)

cat("\nRunning testthat in the directory: ")
cat(getwd(), "\n")


cat("\n================= Group 0 tests =======================\n\n")
test_file(path = "Group0/1_compute_likelihood.r")
