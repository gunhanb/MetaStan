#' Fitting a model-based meta-analysis model using Stan
#'
#' `MBMA_stan` fits a model-based meta-analysis model using Stan.
#'
#' @export
#' @param data An object of `create_MBMA_dat`.
#' @param Pred_doses A numerical vector specfying the doses which prediction will be made.
#' Default is NULL.
#' @param mu_prior A numerical vector specifying the parameter of the normal prior
#' density for baseline risks, first value is parameter for mean, second is for variance.
#' Default is c(0, 10).
#' @param Emax_prior A numerical vector specifying the parameter of the normal prior
#' density for Emax parameter, first value is parameter for mean, second
#' is for standard devation. Default is c(0, 10).
#' @param ED50_prior A numerical vector specifying the parameter of the normal prior
#' density for ED50 parameter, first value is parameter for mean, second
#' is for standard devation. Default is c(0, 10).
#' @param tau_prior A numerical value specifying the standard dev. of the prior density
#' for heterogenety stdev. Default is 0.5.
#' @param tau_prior_dist A string specifying the prior density for the heterogeneity standard deviation,
#' option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
#' half-cauchy prior.
#' @param model A string specifying the model used. Available options are `Baseline_Emax`
#' (Baseline random effects model (\emph{Boucher and Bennets, 2016)}), `CB_Emax`
#' (contrast-based random effects model (\emph{Mawdsley et al ,2016)}), `AB_Emax`
#' (arm-based random effects model, adapted from (\emph{Zhang et al ,2017)}), `CBPlusBaseline_Emax`
#' (contrast-based plus baseline random effects model, adapted from (\emph{Dias et al ,2013)}).
#' Default is `AB_Emax`.
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
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @references Boucher M, et al. The many flavors of model-based meta-analysis:
#' Part I-Introduction and landmark data. \emph{CPT: Pharmacometrics & Systems Pharmacology}.
#' 2016;5:54–64.
#' @references Mawdsley D, et al. Model-based network meta-analysis: A
#' framework for evidence synthesis of clinical trial data.
#' \emph{CPT: Pharmacometrics & Systems Pharmacology}. 2016;5:393-401.
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
#'data('dat.Eletriptan', package = "MetaStan")
#'## Fitting a Binomial-Normal Hierarchial model using WIP priors
#'datMBMA = create_MBMA_dat(dat = dat.Eletriptan,
#'                          armVars = c(dose = "d", responders = "r",
#'                                      sampleSize = "n"),
#'                          nArmsVar = "nd")
#'
#'MBMA.AB.Emax  <- MBMA_stan(data = datMBMA,
#'                           model = "AB_Emax",
#'                           Pred_doses = seq(0, 80, length.out = 11),
#'                           Emax_prior = c(0, 10),
#'                           tau_prior_dist = "half-normal",
#'                           tau_prior = 0.5)
#' }
#'
MBMA_stan = function(data = NULL,
                     Pred_doses = NULL,
                     model = "AB_Emax",
                     mu_prior = c(0, 10),
                     Emax_prior = c(0,10),
                     ED50_prior = c(0,10),
                     tau_prior = 0.5,
                     tau_prior_dist = "half-normal",
                     chains = 4,
                     iter = 2000,
                     warmup = 1000,
                     adapt_delta = 0.95) {
  ################ check model used
  if (model %in% c("Baseline_Emax", "CB_Emax", "AB_Emax", "CBPlusBaseline_Emax") == FALSE) {
    stop("Function argument \"model\" must be equal to \"Baseline_Emax\" or \"CB_Emax\" or \"AB_Emax\" or
         \"CBPlusBaseline_Emax\" !!!")
  }

  ################ check data
  if (!is.data.frame(data)) { stop("Data MUST be a data frame!!!") }

  ################ prior for heterogeneity
  if(tau_prior_dist == "half-normal") { tau_prior_dist_num = 1 }
  if(tau_prior_dist == "uniform")     { tau_prior_dist_num = 2 }
  if(tau_prior_dist == "half-cauchy") { tau_prior_dist_num = 3 }


  ################ Create a list to be used with Stan
  ## Ftiing the model
  if(model == "Baseline_Emax") {
    stanDat <- list(Nobs = nrow(data),
                    Nst = max(data$study),
                    dose = as.numeric(as.character(data$dose)),
                    r = data$responders,
                    n = data$sampleSize,
                    ndose = subset(data, dose == 0)$na,
                    NPred = length(Pred_doses),
                    Pred_doses = Pred_doses)

    fit = rstan::sampling(stanmodels$Baseline_Emax,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }


  if(model == "CB_Emax") {
    data = cbind("ID" = 1:nrow(data), data)
    b_ndx = data[data$dose == 0,]$ID
    ## Nonbaseline index
    t_ndx = data[data$dose != 0,]$ID
    st = data$study


    stanDat <- list(Nobs = nrow(data),
                    Nst = max(data$study),
                    dose = as.numeric(as.character(data$dose)),
                    r = data$responders,
                    n = data$sampleSize,
                    ndose = subset(data, dose == 0)$na,
                    NPred = length(Pred_doses),
                    Pred_doses = Pred_doses,
                    st = st,
                    b_ndx = b_ndx,
                    t_ndx = t_ndx)

    fit = rstan::sampling(stanmodels$CB_Emax,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }


  if(model == "AB_Emax") {
    stanDat <- list(Nobs = nrow(data),
                    Nst = max(data$study),
                    dose = as.numeric(as.character(data$dose)),
                    r = data$responders,
                    n = data$sampleSize,
                    ndose = subset(data, dose == 0)$na,
                    NPred = length(Pred_doses),
                    Pred_doses = Pred_doses)

    fit = rstan::sampling(stanmodels$AB_Emax,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }


  if(model == "CBPlusBaseline_Emax") {
    data = cbind("ID" = 1:nrow(data), data)
    b_ndx = data[data$dose == 0,]$ID
    ## Nonbaseline index
    t_ndx = data[data$dose != 0,]$ID

    stanDat <- list(Nobs = nrow(data),
                    Nst = max(data$study),
                    dose = as.numeric(as.character(data$dose)),
                    r = data$responders,
                    n = data$sampleSize,
                    ndose = subset(data, dose == 0)$na,
                    NPred = length(Pred_doses),
                    Pred_doses = Pred_doses,
                    b_ndx = b_ndx,
                    t_ndx = t_ndx)
    fit = rstan::sampling(stanmodels$CBPlusBaseline_Emax,
                          data = stanDat,
                          chains = chains,
                          iter = iter,
                          warmup = warmup,
                          control = list(adapt_delta = adapt_delta))
  }



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
             fit_sum = fit_sum,
             model = model,
             data = stanDat)
  class(out) <- "MBMA_stan"

  return(out)

}

