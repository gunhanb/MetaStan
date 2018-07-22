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
#' @param theta_prior A numerical vector specifying the parameter of the normal prior
#' density for treatment effect estimate, first value is parameter for mean, second
#' is for variance.
#' @param tau_prior A numerical vector specifying the parameter of the prior density
#' for heterogenety stdev.
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @examples
#' \dontrun{
#' data('dat.Crins2014', package = "MetaStan")
#' ## Subset of dataset where PTLD outcomes available
#' dat.Crins2014.PTLD = subset(dat.Crins2014, is.na(exp.PTLD.events) == FALSE)
#' ## Fitting a Binomial-Normal Hierarchial model
#' bnhm.vague.PTLD.stan  <- meta_stan(ntrt = dat.Crins2014.PTLD$exp.total,
#'                                    nctrl = dat.Crins2014.PTLD$cont.total,
#'                                    rtrt = dat.Crins2014.PTLD$exp.PTLD.events,
#'                                    rctrl = dat.Crins2014.PTLD$cont.PTLD.event,¨
#'                                    mu_prior = c(0, 10), theta_prior = c(0, 100),
#'                                    tau_prior = c(0, 0.5))
#' }
#'
meta_stan = function(ntrt = NULL, nctrl = NULL, rtrt = NULL, rctrl = NULL,
                     mu_prior = c(0, 10), theta_prior = c(0, 100),
                     tau_prior = c(0, 0.5), tau_prior_dist = "half-normal") {

  ## Create a list to be used with Stan
  stanDat <- list(N = length(ntrt),
                  ntrt = ntrt,
                  nctrl = nctrl,
                  rtrt = rtrt,
                  rctrl = rctrl,
                  mupar = mupar,
                  dpar = dpar,
                  taupar = taupar)
  ## Ftiing the model
  out = rstan::sampling(stanmodels$BNHM, data = stanDat)
  return(out)
}

