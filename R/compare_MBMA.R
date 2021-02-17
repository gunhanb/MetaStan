#' Compare MBMA fits using LOO-IC
#'
#' Takes a vector of \code{MBMA_stan} fits and give
#' the model comparison results based on LOO-IC criteria. This
#' is useful to compare different dose-response models. The function
#' depends on \code{loo_compare} function from \code{loo} package.
#'
#'
#' @param x A vector of \code{MBMA_stan} object.
#' @param digits An integer indicating the number of decimal places.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}
#' @references Vehtari, A, A Gelman, and J Gabry (Sept. 2017).
#' "Practical Bayesian model evaluation using leave one-out
#' cross-validation and WAIC." In: Statistics and Computing 27.5, pp. 1413–1432.
#' @export
compare_MBMA <- function(model_list, digits = 2, ...) {

    log_lik_list <- lapply(model_list, extract_log_lik)
    # optional but recommended
    r_eff_list <- lapply(model_list, function(x) {
      ll_array <- extract_log_lik(x, merge_chains = FALSE)
      relative_eff(exp(ll_array))
    })
    # can also pass a list of psis_loo objects to avoid recomputing loo
    loo_list <- lapply(1:length(log_lik_list), function(j) {
      loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
    })
    res = loo_compare(loo_list, ...)[,7:8]

    return(round(res, digits))

}
