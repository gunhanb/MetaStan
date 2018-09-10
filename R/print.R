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
  results = x$fit_sum
  cat("Meta-analysis using MetaStan\n")
  cat("Mean treatment effect\n")
  print(round(results['theta', -c(2, 3, 5, 7, 9, 10)], digits))

  if (x$model %in% c("BNHM1", "BNHM2") == TRUE){
    cat("Heterogeneity stdev\n")
    print(round(results['tau', -c(2, 3, 5, 7, 9, 10)], digits))
  }

  return(invisible())
}
