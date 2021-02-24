#' Plot a dose-response plot
#'
#' Takes a \code{MBMA_stan} object which is obtained by function \code{MBMA_stan} and plot
#' a dose-response plot, showing observed event probabilities and the estimated dose-response
#'function with pointwise 95% credible intervals.
#'
#' @param x A \code{MBMA_stan} object.
#' @param ... Further arguments passed to ggplot.
#' @return The return value is invisible \code{NULL}.
#' @examples
#' \dontrun{
#'## Load the dataset
#'data('dat.Eletriptan', package = "MetaStan")
#'## Fitting a Binomial-Normal Hierarchial model using WIP priors
#'datMBMA = create_MBMA_dat(dat = dat.Eletriptan,
#'                          armVars = c(dose = "d", r = "r",
#'                                      n = "n"),
#'                          nArmsVar = "nd")
#'
#' MBMA.Emax  <- MBMA_stan(data = datMBMA,
#'                           likelihood = "binomial",
#'                           dose_response = "emax",
#'                           Emax_prior = c(0, 10),
#'                           ED50_prior = "functional",
#'                           tau_prior_dist = "half-normal",
#'                           tau_prior = 0.5)
#' print(MBMA.Emax)
#' plot(MBMA.Emax)
#' # use ggplot2 package for further adjustments
#' Library(ggplot2)
#' theme_set(theme_classic())
#' plot(MBMA.Emax) + xlab("Doses (mg)") + ylab("response probabilities")
#'
#' }
#' @source This function uses \code{ggplot} function from \code{ggplot2}
#' R package.
#' @author Christian Roever and Burak Kuersad Guenhan
#' @seealso \code{ggplot2::ggplot}
#' @export
plot.MBMA_stan = function(x = MBMA.stan,...) {



  ## Plot the observed probabilities
  r = as.vector(rbind(x$data_wide$r1, x$data_wide$r2,
                      x$data_wide$r3, x$data_wide$r4))
  dose = as.numeric(as.vector(rbind(x$data_wide$d1, x$data_wide$d2,
                                    x$data_wide$d3, x$data_wide$d4)))
  n = as.vector(rbind(x$data_wide$n1, x$data_wide$n2,
                      x$data_wide$n3, x$data_wide$n4))
  obs.probs = r/n

  mydata = data.frame(dose = dose,
                      y.obs = obs.probs,
                      y.samplesizes = n,
                      Study = gl(nrow(x$data_wide), 4,
                                 labels = letters[1:nrow(x$data_wide)]))

  Pred_probs_args = rep(NA, times = x$data$Npred)

  for(i in 1:x$data$Npred) {
    Pred_probs_args[i] = sprintf("Pred_probs[%i%s", i, "]")
  }

  pred.probs = x$fit_sum[Pred_probs_args, c(4, 5, 6, 7, 8)]

  x_data = data.frame(x$fit_sum)

  Preddata = data.frame(Pred_doses = x$data$Pred_doses,
                        y.pred =  pred.probs[,3],
                        y.lo =  pred.probs[,1],
                        y.lo25 =  pred.probs[,2],
                        y.hi =  pred.probs[,5],
                        y.hi75 =  pred.probs[,4])

  dose.plot = ggplot2::ggplot(mydata, ggplot2::aes(x = x_data, y = y)) +
    ggplot2::geom_point(ggplot2::aes(x = dose, y = y.obs, size = y.samplesizes)) +
    ggplot2::geom_line(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.pred)) +
    ggplot2::geom_line(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.lo), linetype = "dashed") +
    ggplot2::geom_line(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.hi), linetype = "dashed") +
    ggplot2::geom_line(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.lo25), linetype = "dashed") +
    ggplot2::geom_line(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.hi75), linetype = "dashed") +
    ggplot2::geom_ribbon(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.pred, ymin=y.lo25,ymax=y.hi75),
                fill="blue", alpha=0.75) +
    ggplot2::geom_ribbon(data = Preddata, ggplot2::aes(x = Pred_doses, y = y.pred, ymin=y.lo,ymax=y.hi),
                fill="blue", alpha=0.25) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::theme(legend.position = "none")

 dose.plot
}
