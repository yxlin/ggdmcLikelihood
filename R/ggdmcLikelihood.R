#' Likelihood Tools for the \code{ggdmc} Package
#'
#' Efficient computation of likelihoods in design-based choice response time
#' (RT) and cognitive diagnostic models. The package enables rapid evaluation of
#' likelihood functions for both single- and multi-subject models across
#' trial-level data for the choice RT model. For the cognitive diagnostic model,
#' the pacakge, for now, treats the school as the upper and the student as the
#' lower level. The package also offers fast initialisation of starting
#' parameters for genetic sampling with many Markov chains, facilitating
#' estimation in complex models typically found in experimental psychology and
#' behavioural science. These optimisations help reduce computational overhead
#' in large-scale model fitting tasks.
#'
#' @name ggdmcLikelihood
#' @keywords internal
#' @author  Yi-Shin Lin <yishinlin001@gmail.com>
#' @importFrom Rcpp evalCpp
#' @useDynLib ggdmcLikelihood
"_PACKAGE"
NULL
