#' Fitting a model-based meta-analysis model using Stan
#'
#' `MBMA_stan` fits a model-based meta-analysis model using Stan.
#'
#' @export
#' @param data An object of `create_MBMA_dat`.
#' @param Pred_doses A numerical vector specifying the doses which prediction will be made.
#' @param mu_prior A numerical vector specifying the parameter of the normal prior
#' density for baseline risks, first value is parameter for mean, second is for variance.
#' Default is c(0, 10).
#' @param alpha_prior A numerical vector specifying the parameter of the normal prior
#' density for the alpha parameter, first value is parameter for mean, second is for variance.
#' Default is c(0, 10). Needed for linear and linear log-dose models.
#' @param Emax_prior A numerical vector specifying the parameter of the normal prior
#' density for Emax parameter, first value is parameter for mean, second
#' is for standard deviation. Default is c(0, 10). Needed for emax and sigmoidal models.
#' @param ED50_prior A numerical vector specifying the parameter of the normal prior
#' density for ED50 parameter, first value is parameter for mean, second
#' is for standard deviation. Default is c(0, 10). Needed for emax and sigmoidal models.
#' @param ED50_prior_dist A string specifying the prior density for the ED50 parameter,
#' `functional` is for a functional uniform prior, `half-normal` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param tau_prior A numerical value specifying the standard dev. of the prior density
#' for heterogenety stdev. Default is 0.5.
#' @param tau_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param gamma_prior A numerical vector specifying the parameter of the normal prior
#' density for gamma parameter, first value is parameter for mean, second
#' is for standard deviation. Default is c(1, 2). Needed for sigmoidal model.
#' @param dose_response A string specifying the function defining the dose-response model.
#' Options include "linear", "log-linear", "emax", and "sigmoidal".
#' @param likelihood A string specifying the likelihood of distributions defining the statistical
#' model. Options include "normal", "binomial", and "Poisson".
#' @param re A string specifying whether random-effects are included to the model. When `FALSE`, the
#' model corresponds to a fixed-effects model. The default is `TRUE`.
#' @param ncp A string specifying whether to use a non-centered parametrization.
#' The default is `TRUE`.
#' @param adapt_delta A numerical value specfying the target average proposal acceptance
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
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @references Boucher M, et al. The many flavors of model-based meta-analysis:
#' Part I-Introduction and landmark data. \emph{CPT: Pharmacometrics and Systems Pharmacology}.
#' 2016;5:54-64.
#' @references Guenhan BK, Roever C, Friede T. MetaStan: An R package for meta-analysis
#' and model-based meta-analysis using Stan. In preparation.
#' @references Mawdsley D, et al. Model-based network meta-analysis: A
#' framework for evidence synthesis of clinical trial data.
#' \emph{CPT: Pharmacometrics and Systems Pharmacology}. 2016;5:393-401.
#' @references  Zhang J, et al. (2014). Network meta-analysis of randomized
#' clinical trials: Reporting the proper summaries. \emph{Clinical Trials}.
#' 11(2), 246–262.
#' @references Dias S, et al. Absolute or relative effects?
#' Arm-based synthesis of trial data. \emph{Research Synthesis Methods}.
#' 2016;7:23–28.
#'
#' @examples
#' \dontrun{
#'## Load the dataset
#' data('dat.Eletriptan', package = "MetaStan")
## Fitting a MBMA model
#' datMBMA = create_MetaStan_dat(dat = dat.Eletriptan,
#'                               armVars = c(dose = "d",
#'                                           responders = "r",
#'                                           sampleSize = "n"),
#'                               nArmsVar = "nd")
#'
#' MBMA.Emax  <- MBMA_stan(data = datMBMA,
#'                         likelihood = "binomial",
#'                         dose_response = "emax",
#'                         Pred_doses = seq(0, 80, length.out = 11),
#'                         mu_prior = c(0, 100),
#'                         Emax_prior = c(0, 100),
#'                         tau_prior_dist = "half-normal",
#'                         tau_prior = 0.5)
#' plot(MBMA.Emax) + ggplot2::xlab("Doses (mg)") + ggplot2::ylab("response probabilities")
#'
#' }
#'
MBMA_stan = function(data = NULL,
                     likelihood = NULL,
                     dose_response = "emax",
                     mu_prior = c(0, 10),
                     Emax_prior = c(0, 100),
                     alpha_prior = c(0, 100),
                     tau_prior = 0.5,
                     tau_prior_dist = "half-normal",
                     ED50_prior = c(-2.5, 1.8),
                     ED50_prior_dist = "functional",
                     gamma_prior = c(1, 2),
                     Pred_doses,
                     re = TRUE,
                     ncp = TRUE,
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

  ################ check model used
  if (dose_response %in% c("linear", "log-linear", "emax", "sigmoidal") == FALSE) {
    stop("Function argument \"dose-response\" must be equal to \"linear\" or \"log-linear\" or \"emax\" or
         \"sigmoidal\" !!!")
  }

  ################ check data
  data_wide = data$data_wide
  data = data$data_long

  if (!is.data.frame(data)) { stop("Data MUST be a data frame!!!") }


  if (missing(Pred_doses)) {
  Pred_doses = seq(from = 0, to = max(as.numeric(as.character(data$dose)),
                                      na.rm = TRUE), length.out =  30)

  }


  ################ check prior for heterogeneity parameter
  if(is.null(tau_prior_dist) == TRUE){
    stop("Function argument \"half-normal\" or \"uniform\" or \"half-cauchy\" must be specified !!!")
  }

  ################ prior for heterogeneity
  if(tau_prior_dist == "half-normal") { tau_prior_dist_num = 1 }
  if(tau_prior_dist == "uniform")     { tau_prior_dist_num = 2 }
  if(tau_prior_dist == "half-cauchy") { tau_prior_dist_num = 3 }

  if(ED50_prior_dist == "functional") { ED50_prior_dist_num = 1 }
  if(ED50_prior_dist == "half-normal") { ED50_prior_dist_num = 2 }


  if(likelihood == "normal")   { link = 1 }
  if(likelihood == "binomial") { link = 2 }
  if(likelihood == "poisson")  { link = 3 }

  if(dose_response == "linear")     { dose_response = 1 }
  if(dose_response == "log-linear") { dose_response = 2 }
  if(dose_response == "emax")       { dose_response = 3 }
  if(dose_response == "sigmoidal")  { dose_response = 4 }


  ################ Create a list to be used with Stan
  Nobs = nrow(data)

  y <- array(rep(0,Nobs))
  y_se <- array(rep(0,Nobs))
  r   <- array(rep(0,Nobs))
  n <- array(rep(1,Nobs))
  count <- array(rep(0, Nobs))
  exposure <- array(rep(0, Nobs))

  if(likelihood == "binomial") {
    r = data$responders
    n = data$sampleSize
  }
  if(likelihood == "normal") {
    y = data$mean
    y_se = data$std.err
  }
  if(likelihood == "poisson") {
    count = data$count
    exposure = data$exposure
  }

  if(re == TRUE)   {re = 1}
  if(re == FALSE)  {re = 0}
  if(ncp == TRUE)  {ncp = 1}
  if(ncp == FALSE) {ncp = 0}
  ## Create a list to be used with Stan
  data = cbind(ID = 1:nrow(data), data)
  b_ndx = data[data$dose == 0, ]$ID
  t_ndx = data[data$dose != 0, ]$ID

  if(dose_response == 1 || dose_response == 2) {alpha_Emax_prior = alpha_prior}
  if(dose_response == 3 || dose_response == 4)    {alpha_Emax_prior = Emax_prior}


  stan_dat_long <- list(Nobs = Nobs,
                        link = link,
                        st = data$study,
                        Nst = max(data$study),
                        dose = as.numeric(as.character(data$dose)),
                        ndose = data_wide$nd,
                        Npred = length(Pred_doses),
                        Pred_doses = Pred_doses,
                        y = y,
                        y_se = y_se,
                        r = r,
                        n = n,
                        count = count,
                        exposure = exposure,
                        maxdose = max(as.numeric(data$dose)),
                        mu_prior = mu_prior,
                        alpha_prior = alpha_Emax_prior,
                        ED50_prior = ED50_prior,
                        gamma_prior = gamma_prior,
                        ED50_prior_dist = ED50_prior_dist_num,
                        tau_prior = tau_prior,
                        tau_prior_dist = tau_prior_dist_num,
                        dose_response = dose_response,
                        re = re,
                        ncp = ncp,
                        b_ndx = b_ndx,
                        t_ndx = t_ndx)



  ## Fitting the model
  fit = rstan::sampling(stanmodels$MBMA,
                        data = stan_dat_long,
                        chains = chains,
                        iter = iter,
                        warmup = warmup,
                        control = list(adapt_delta = adapt_delta),
                        ...)


  ## MODEL FINISHED
  fit_sum   <- rstan::summary(fit)$summary
  Rhat.max  <- max(fit_sum[,"Rhat"], na.rm=TRUE)
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
             data = stan_dat_long,
             data_wide = data_wide,
             data_long = data,
             Rhat.max = Rhat.max,
             N_EFF.min = N_EFF.min,
             tau_prior_dist = tau_prior_dist,
             ED50_prior_dist = ED50_prior_dist)

  class(out) <- "MBMA_stan"

  return(out)

}

