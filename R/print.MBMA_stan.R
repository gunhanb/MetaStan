#' Print MBMA object
#'
#' Takes an \code{MBMA_stan} object which is obtained by function \code{MBMA_stan} and print
#' the model and data information such as model type used in the model.
#'
#' @param x A \code{MBMA_stan} object.
#' @param digits An integer indicating the number of decimal places.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}
#' @export
print.MBMA_stan <- function(x, digits = 2, ...) {
  if (!is.element("MBMA_stan", class(x)))
    stop("Argument 'x' must be an object of class \"MBMA_stan\".")
  cat("Model-based meta-analysis using MetaStan \n\n")
  cat("Maximum Rhat:", signif(x$Rhat.max, digits=digits),"\n")
  cat("Minimum Effective Sample Size:", signif(x$N_EFF.min, digits=digits),"\n\n")

  cat("mu prior: Normal")
  cat("(")
  cat(x$data$mu_prior[1])
  cat(",")
  cat(x$data$mu_prior[2])
  cat(")")
  cat("\n")
  cat("alpha prior: Normal")
  cat("(")
  cat(x$data$alpha_prior[1])
  cat(",")
  cat(round(x$data$alpha_prior[2], 2))
  cat(")")
  cat("\n")

  if (x$data$dose_response == 3 || x$data$dose_response == 4){
    cat("ED50 prior:")
    cat(x$ED50_prior_dist)
    cat("(")
    cat(x$data$ED50_prior[1])
    cat(",")
    cat(round(x$data$ED50_prior[2], 2))
    cat(")")
    cat("\n\n")
  }

  if (x$data$dose_response == 4){
    cat("gamma prior: Normal")
    cat("(")
    cat(x$data$gamma_prior[1])
    cat(",")
    cat(round(x$data$gamma_prior[2], 2))
    cat(")")
    cat("\n")
  }




  if (x$data$re){
    cat("tau prior:")
    cat(x$tau_prior_dist)
    cat("(")
    cat(x$data$tau_prior)
    cat(")")
    cat("\n\n")
  }



  if (x$data$dose_response == 1) {
    cat("Dose-response function = linear \n\n")
    cat("alpha estimates\n")
    print(round(x$fit_sum['alpha', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")

  }
  if (x$data$dose_response == 2) {
    cat("Dose-response function = log-linear \n\n")
    cat("alpha estimates\n")
    print(round(x$fit_sum['alpha', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")

  }
  if (x$data$dose_response == 3) {
    cat("Dose-response function = emax \n\n")
    cat("Emax estimates\n")
    print(round(x$fit_sum['alpha', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")
    cat("ED50 estimates\n")
    print(round(x$fit_sum['ED50[1]', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")


  }
  if (x$data$dose_response == 4) {
    cat("Dose-response function = sigmoidal emax \n\n")
    cat("Emax estimates\n")
    print(round(x$fit_sum['alpha', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")
    cat("ED50 estimates\n")
    print(round(x$fit_sum['ED50[1]', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")
    cat("gamma estimates\n")
    print(round(x$fit_sum['gamma[1]', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")

  }


  if (x$data$re == TRUE){
 #   cat("Heterogeneity stdev (tau)\n")
 #   print(round(x$fit_sum['tau[1]', -c(2, 3, 5, 7, 9, 10)], digits))
    mcmc = coda::mcmc.list(rstan::As.mcmc.list(x$fit, pars = c("tau[1]")))
    tau_int = HDInterval::hdi(mcmc, credMass = 0.95)
    print(round(c(tau_int[1], x$fit_sum['tau[1]', "50%"], tau_int[2]), digits))

  }

  return(invisible())
}
