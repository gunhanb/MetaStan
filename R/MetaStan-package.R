#' The 'MetaStan' package.
#'
#' @description Fitting Bayesian meta-analysis models via Rstan.
#'
#' @details To fit meta-analysis models using frequentist methods, there are many R packages
#' available including `metafor`. On the other hand, Bayesian estimation methods such as
#' Markov chain Monte Carlo (MCMC) are very attractive for meta-analysis, especially
#' because they can be used to fit more complicated models. These include binomial-normal
#' hierarchical models and beta-binomial models which are based on the exact distributional
#' assumptions unlike (commonly used) normal-normal hierarchical model. Another advantage of
#' Bayesian methods to be able to use informative prior distributions for example to
#' regularize heterogeneity estimates in case of low number of studies. Thus, we developed
#' `MetaStan` which uses Stan (a modern MCMC engine) to fit several pairwise meta-analysis
#' models including binomial-normal hierarchical model and beta-binomial model. This package is
#' also the accompanying package of Guenhan et al (2020). Another important functionality of the
#' package is the model-based meta-analysis models.
# '
#' @docType package
#' @author Burak Kuersad Guenhan <burak.gunhan@med.uni-goettingen.de>
#' @name MetaStan-package
#' @aliases MetaStan
#' @useDynLib MetaStan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom forestplot forestplot fpDrawPointCI fpDrawSummaryCI fpDrawBarCI fpColors
#' @importFrom rstan sampling As.mcmc.list
#' @importFrom loo loo_compare
#' @importFrom metafor escalc
#' @importFrom stats qnorm
#' @importFrom grid gpar unit
#' @importFrom coda mcmc.list
#' @importFrom HDInterval hdi
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.3.
#' http://mc-stan.org
#'
#' Günhan, B and Röver, C and Friede, T (2020). Random-effects meta-analysis of few studies
#' involving rare events. Research Synthesis Methods. doi = 10.1002/jrsm.1370.
NULL
