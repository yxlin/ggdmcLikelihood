#' Likelihood tools for the ggdmc package
#'
#' This package provides a user interface for efficiently computing 
#' likelihoods in design-based choice response time models. It is built to
#' integrate with the 'ggdmc' package, enabling rapid initialisation of 
#' starting samples for genetic algorithms, which typically use a large 
#' number of Markov chains. The likelihood computation supports both 
#' single-subject and multi-subject models. Implemented methods include 
#' fast evaluation of likelihoods across trial-level data, support for 
#' various response time distributions (e.g., diffusion models). Ideal 
#' for cognitive scientists and decision modellers, this tool enhances 
#' the usability and computational performance of complex model fitting
#' tasks in experimental psychology and behavioral data research.
#' 
#' @name ggdmcLikelihood
#' @keywords internal
#' @author  Yi-Shin Lin <yishinlin001@gmail.com>
#' @importFrom Rcpp evalCpp
#' @useDynLib ggdmcLikelihood
"_PACKAGE"
NULL
