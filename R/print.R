#' Print meta_stan object
#'
#' Takes an \code{meta_stan} object which is obtained by function \code{meta_stan} and print
#' the model and data information such as model type used in the model.
#'
#' The resulting data.frame can be used as data argument in \code{meta_stan}.
#'
#' @param x A \code{meta_stan} object.
#' @param digits An integer indicating the number of decimal places.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}
#' @export
print.meta_stan <- function(x, digits = 2, ...) {
  if (!is.element("meta_stan", class(x)))
    stop("Argument 'x' must be an object of class \"meta_stan\".")
  cat("Meta-analysis using MetaStan \n\n")
  cat("Maximum Rhat:", signif(x$Rhat.max, digits=3),"\n")
  cat("Minimum Effective Sample Size:", signif(x$N_EFF.min, digits=digits),"\n\n")

  cat("mu prior: Normal")
  cat("(")
  cat(x$data$mu_prior[1])
  cat(",")
  cat(x$data$mu_prior[2])
  cat(")")
  cat("\n")
  cat("theta prior: Normal")
  cat("(")
  cat(x$data$theta_prior[1])
  cat(",")
  cat(round(x$data$theta_prior[2], 2))
  cat(")")
  cat("\n")

  if (x$data$mreg){
    cat("beta prior: Normal")
    cat("(")
    cat(x$data$beta_prior[1])
    cat(",")
    cat(x$data$beta_prior[2])
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

  cat("Treatment effect (theta) estimates\n")
  print(round(x$fit_sum['theta', -c(2, 3, 5, 7, 9, 10)], digits))
  cat("\n")

  if (x$data$mreg){
    cat("Beta coeffients\n")
    print(round(x$fit_sum['beta[1,1]', -c(2, 3, 5, 7, 9, 10)], digits))
    cat("\n")
  }

  if (x$data$re){
    cat("Heterogeneity stdev (tau)\n")
    print(round(x$fit_sum['tau[1]', -c(2, 3, 5, 7, 9, 10)], digits))
  }


  return(invisible())
}
