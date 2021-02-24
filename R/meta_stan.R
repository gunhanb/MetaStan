#' Fitting a meta-analysis model using Stan
#'
#' `meta_stan` fits a meta-analysis model using Stan.
#'
#' @export
#' @param data Data frame created by `convert_data_arm`
#' @param mu_prior A numerical vector specifying the parameter of the normal prior
#' density for baseline risks, first value is parameter for mean, second is for variance.
#' Default is c(0, 10).
#' @param theta_prior A numerical vector specifying the parameter of the normal prior
#' density for treatment effect estimate, first value is parameter for mean, second
#' is for variance. Default is NULL.
#' @param delta A numerical value specifying the upper bound of the a priori interval for
#' treatment effect on odds ratio scale (\emph{Guenhan et al (2020)}). This is used to calculate
#' a normal weakly informative prior.
#' for theta. Thus when this argument is specified, `theta` should be left empty. Default is NULL.
#' @param tau_prior A numerical value specifying the standard dev. of the prior density
#' for heterogeneity stdev. Default is 0.5.
#' @param tau_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param beta_prior A numerical vector specifying the parameter of the normal prior
#' density for beta coefficients in a meta-regression model, first value is parameter for mean, second
#' is for variance. Default is c(0, 100).
#' @param likelihood A string specifying the likelihood function defining the statistical
#' model. Options include  "normal", "binomial", and "Poisson".
#' @param re A string specifying whether random-effects are included to the model. When `FALSE`, the
#' model corresponds to a fixed-effects model. The default is `TRUE`.
#' @param ncp A string specifying whether to use a non-centered parametrization.
#' The default is `TRUE`.
#' @param mreg A string specifying whether to fit a meta-regression model.
#' The default is `FALSE`.
#' @param cov A numeric matrix specifying trial-level covariates in each row.
#' This is need when `mreg = TRUE`.
#' @param adapt_delta A numerical value specifying the target average proposal acceptance
#' probability for adaptation. See Stan manual for details. Default is 0.95. In general
#' you should not need to change adapt_delta unless you see a warning message about
#' divergent transitions, in which case you can increase adapt_delta from the
#' default to a value closer to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc).
#' @param iter A positive integer specifying the number of iterations for each chain
#' (including warmup). The default is 2000.
#' @param warmup A positive integer specifying the number of warmup (aka burnin)
#' iterations per chain. The default is 1000.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#' @param ... Further arguments passed to or from other methods.
#' @return an object of class `MetaStan`.
#' @references Guenhan BK, Roever C, Friede T. Random-effects meta-analysis of few studies involving
#' rare events \emph{Resarch Synthesis Methods} 2020; doi:10.1002/jrsm.1370.
#' @references Jackson D, Law M, Stijnen T, Viechtbauer W, White IR. A comparison of 7
#' random-effects models for meta-analyses that estimate the summary odds ratio.
#' \emph{Stat Med} 2018;37:1059--1085.
#' @references Kuss O.Statistical methods for meta-analyses including information
#' from studies without any events-add nothing to nothing and succeed nevertheless,
#' \emph{Stat Med}, 2015; 4; 1097--1116, doi: 10.1002/sim.6383.
#'
#' @examples
#' \dontrun{
#' data('dat.Crins2014', package = "MetaStan")
#' dat_MetaStan <- convert_data_arm(dat.Crins2014$exp.total, dat.Crins2014$cont.total,
#'                              dat.Crins2014$exp.AR.events, dat.Crins2014$cont.AR.events)
#' bnhm.Crins  <- meta_stan(data = dat_MetaStan, likelihood = "binomial",
#'                          mu_prior = c(0, 10), theta_prior = c(0, 100),
#'                          tau_prior =  0.5)
#' print(bnhm.Crins)
#'
#' forestplot(bnhm.Crins)
#'
#' ## TB dataset
#' data('dat.Berkey1995', package = "MetaStan")
#' ## Fitting a Binomial-Normal Hierarchical model using WIP priors
#' dat_MetaStan <- convert_data_arm(dat.Berkey1995$nt, dat.Berkey1995$nc,
#'                             dat.Berkey1995$rt, dat.Berkey1995$rc)
#'
#' meta.reg.stan  <- meta_stan(data = dat_MetaStan,
#'                            likelihood = "binomial",
#'                            mu_prior = c(0, 10),
#'                            theta_prior = c(0, 100),
#'                            tau_prior = 0.5,
#'                            tau_prior_dist = "half-normal",
#'                            mreg = TRUE,
#'                            cov = matrix(dat.Berkey1995$Latitude, nrow = 1))
#'
#' print(meta.reg.stan)
#' }
#'
meta_stan = function(data = NULL,
                     likelihood = NULL,
                     mu_prior = c(0, 10),
                     theta_prior = NULL,
                     tau_prior = 0.5,
                     tau_prior_dist = "half-normal",
                     beta_prior = c(0, 100),
                     delta = NULL,
                     re = TRUE,
                     ncp = TRUE,
                     mreg = FALSE,
                     cov = NULL,
                     chains = 4,
                     iter = 2000,
                     warmup = 1000,
                     adapt_delta = 0.95,
                     ...) {
  ################ check likelihood used
  if(is.null(likelihood) == TRUE){
    stop("Function argument \"likelihood\" must be specified !!!")
  }

  if (likelihood %in% c("binomial", "normal", "poisson") == FALSE) {
    stop("Function argument \"likelihood\" must be equal to \"binomial\" or \"normal\" or \"poisson\"!!!")
  }


  ################ check prior for treatment effect parameter
  if(is.null(delta) == TRUE & is.null(theta_prior) == TRUE){
    stop("Function argument \"delta\" or \"theta_prior\" must be specified !!!")
  }

  if(is.null(delta) == FALSE & is.null(theta_prior) == FALSE){
    stop("Either \"delta\" OR \"theta_prior\" can be specified, but not both !!!")
  }

  if(mreg == TRUE & is.null(cov) == TRUE){
    stop("Function argument \"cov\" argument must be specified, when \"mreg\" is TRUE !!!")
  }

  ## Calculate standard deviation of theta_prior from delta
  if(is.null(delta) == FALSE){
    theta_prior[1] = 0
    theta_prior[2] = (log(delta) - log(1/delta)) / (2 * 1.96)
  }

  ################ check prior for heterogeneity parameter
  if(is.null(tau_prior_dist) == TRUE){
    stop("Function argument \"half-normal\" or \"uniform\" or \"half-cauchy\" must be specified !!!")
  }

  if(tau_prior_dist == "half-normal") { tau_prior_dist_num = 1 }
  if(tau_prior_dist == "uniform")     { tau_prior_dist_num = 2 }
  if(tau_prior_dist == "half-cauchy") { tau_prior_dist_num = 3 }

  if(likelihood == "normal")   { link = 1 }
  if(likelihood == "binomial") { link = 2 }
  if(likelihood == "poisson")  { link = 3 }


  ################ check data
  data_wide = data$data_wide
  data = data$data_long

  if (!is.data.frame(data)) { stop("Data MUST be a data frame!!!") }

  ################ prepare data
  Nobs = nrow(data)

  y <- array(rep(0,Nobs))
  y_se <- array(rep(0,Nobs))
  r   <- array(rep(0,Nobs))
  n <- array(rep(1,Nobs))
  count <- array(rep(0, Nobs))
  exposure <- array(rep(0, Nobs))

  if(likelihood == "binomial") {
    r = data$r
    n = data$n
  }
  if(likelihood == "normal") {
    y = data$y
    y_se = data$y_se
  }
  if(likelihood == "poisson") {
    count = data$count
    exposure = data$exposure
  }

  if(mreg == FALSE){
    cov = matrix(rep(0, max(as.numeric(as.vector(data$mu)))), nrow = 1)
  }
  ncov = nrow(cov)

  if(re == TRUE)   {re = 1}
  if(re == FALSE)  {re = 0}
  if(ncp == TRUE)  {ncp = 1}
  if(ncp == FALSE) {ncp = 0}
  if(mreg == TRUE)   {mreg = 1}
  if(mreg == FALSE)  {mreg = 0}

  ## Create a list to be used with Stan
  stanDat <- list(Nobs = Nobs,
                  t = data$theta,
                  link = link,
                  st = as.numeric(as.vector(data$mu)),
                  y = y,
                  y_se = y_se,
                  r = r,
                  n = n,
                  count = count,
                  exposure = exposure,
                  mu_prior = mu_prior,
                  theta_prior = theta_prior,
                  tau_prior = tau_prior,
                  tau_prior_dist = tau_prior_dist_num,
                  beta_prior = beta_prior,
                  re = re,
                  ncp = ncp,
                  mreg = mreg,
                  cov = cov,
                  ncov = ncov)



  ## Fitting the model
  fit = rstan::sampling(stanmodels$SMA,
                        data = stanDat,
                        chains = chains,
                        iter = iter,
                        warmup = warmup,
                        control = list(adapt_delta = adapt_delta),
                        ...)


  ## MODEL FINISHED
  fit_sum <- rstan::summary(fit)$summary
  Rhat.max <- max(fit_sum[,"Rhat"], na.rm=TRUE)
  N_EFF.min <- min(fit_sum[,"n_eff"], na.rm=TRUE)

  if(Rhat.max > 1.1)
    warning("Maximal Rhat > 1.1. Consider increasing meta_stan MCMC parameters.")

  ## finally include a check if the Stan NuTS sample had any
  ## divergence in the sampling phase, these are not supposed to
  ## happen and can often be avoided by increasing adapt_delta
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  n_divergent <- sum(sapply(sampler_params, function(x) sum(x[,'divergent__'])) )
  if(n_divergent > 0) {
    warning(paste("In total", n_divergent, "divergent transitions occured during the sampling
                  phase.\nPlease consider increasing adapt_delta closer to 1."))
  }


  out = list(fit = fit,
             fit_sum = fit_sum,
             data_wide = data_wide,
             data_long = data,
             stanDat = stanDat,
             Rhat.max = Rhat.max,
             N_EFF.min = N_EFF.min,
             tau_prior_dist = tau_prior_dist)

  class(out) <- "meta_stan"

  return(out)

}

