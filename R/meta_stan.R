#' Fitting a meta-analysis model using Stan
#'
#' `meta_stan` fits a meta-analysis model using Stan.
#'
#' @export
#' @param ntrt Number of subjects in treatment arm
#' @param nctrl Number of subjects in control arm
#' @param rtrt Number of events in treatment arm
#' @param rctrl Number of events in contrl arm
#' @param mu_prior A numerical vector specifying the parameter of the normal prior
#' density for baseline risks, first value is parameter for mean, second is for variance.
#' Default is c(0, 10).
#' @param theta_prior A numerical vector specifying the parameter of the normal prior
#' density for treatment effect estimate, first value is parameter for mean, second
#' is for variance. Default is NULL.
#' @param delta_u A numerical value specifying the upper bound of the a priori inerval for
#' treatment effect on odds ratio scale. This is used to cacluate a normal weakly informative prior
#' for theta. Thus when this argument is pecified, `theta` should be left empty. Default is NULL.
#' @param tau_prior A numerical value specifying the standard dev. of the prior density
#' for heterogenety stdev. Default is 0.5.
#' @param tau_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is 'half-normal' for half-normal prior.
#' @param adapt_delta A numerical value specfying the target average proposal acceptance
#' probability for adaptation. See Stan manual for details. Default is 0.95. In general
#' you should not need to change adapt_delta unless you see a warning message about
#' divergent transitions, in which case you can increase adapt_delta from the
#' default to a value closer to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc).
#' @param iter A positive integer specifying the number of iterations for each chain
#' (including warmup). The default is 2000.
#' @param warmup A positive integer specifying the number of warmup (aka burnin)
#' iterations per chain.
#' @param chains A positive integer specifying the number of Markov chains.
#' The default is 4.
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' \dontrun{
#' data('dat.Crins2014', package = "MetaStan")
#' ## Subset of dataset where PTLD outcomes available
#' dat.Crins2014.PTLD = subset(dat.Crins2014, is.finite(exp.PTLD.events))
#' ## Fitting a Binomial-Normal Hierarchial model
#' bnhm.vague.PTLD.stan  <- meta_stan(ntrt = dat.Crins2014.PTLD$exp.total,
#'                                    nctrl = dat.Crins2014.PTLD$cont.total,
#'                                    rtrt = dat.Crins2014.PTLD$exp.PTLD.events,
#'                                    rctrl = dat.Crins2014.PTLD$cont.PTLD.event,¨
#'                                    mu_prior = c(0, 10), theta_prior = c(0, 100),
#'                                    tau_prior =  0.5)
#' }
#'
meta_stan = function(ntrt = NULL,
                     nctrl = NULL,
                     rtrt = NULL,
                     rctrl = NULL,
                     mu_prior = c(0, 10),
                     theta_prior = NULL,
                     tau_prior = 0.5,
                     tau_prior_dist = "half-normal",
                     delta_u = NULL,
                     chains = 4,
                     iter = 2000,
                     warmup = 1000,
                     adapt_delta = 0.95) {

  ################ check prior for treatment effect parameter
  if(is.null(delta_u) == TRUE & is.null(theta_prior) == TRUE){
    stop("Function argument \"delta_u\" or \"theta_prior\" must be specified !!!")
  }

  if(is.null(delta_u) == FALSE & is.null(theta_prior) == FALSE){
    stop("Either \"delta_u\" OR \"theta_prior\" can be specified, but not both !!!")
  }

  ## Calculate standard devation of theta_prior from delta_u
  if(is.null(delta_u) == FALSE){
    theta_prior[1] = 0
    theta_prior[2] = (log(delta_u) - log(1/delta_u)) / (2 * 1.96)
  }


  if(tau_prior_dist == "half-normal") { tau_prior_dist_num = 1 }
  ## Create a list to be used with Stan
  stanDat <- list(N = length(ntrt),
                  ntrt = ntrt,
                  nctrl = nctrl,
                  rtrt = rtrt,
                  rctrl = rctrl,
                  mu_prior = mu_prior,
                  theta_prior = theta_prior,
                  tau_prior = tau_prior,
                  tau_prior_dist = tau_prior_dist_num)
  ## Ftiing the model
  fit = rstan::sampling(stanmodels$BNHM,
                        data = stanDat,
                        chains = chains,
                        iter = iter,
                        warmup = warmup,
                        control = list(adapt_delta = adapt_delta))
  ## MODEL FINISHED
  fit_sum <- rstan::summary(fit)$summary


  Rhat.max <- max(fit_sum[,"Rhat"], na.rm = TRUE)

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
             fit_sum = fit_sum)

  return(out)
}

