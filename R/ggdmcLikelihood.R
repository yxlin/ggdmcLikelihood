#' Likelihood Tools for the \code{ggdmc} Package
#'
#' Efficient computation of likelihoods in design-based choice response time
#' models, including the Decision Diffusion Model, is supported. The package 
#' enables rapid evaluation of likelihood functions for both 
#' single- and multi-subject models across trial-level data. It also offers 
#' fast initilisation of starting parameters for genetic sampling with many 
#' Markov chains, facilitating estimation in complex models typically found 
#' in experimental psychology and behavioural science. These optimisations 
#' help reduce computational overhead in large-scale model fitting tasks.
#'
#' @name ggdmcLikelihood
#' @keywords internal
#' @author  Yi-Shin Lin <yishinlin001@gmail.com>
#' @importFrom Rcpp evalCpp
#' @useDynLib ggdmcLikelihood
"_PACKAGE"
NULL
