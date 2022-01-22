#' Plot a forest plot
#'
#' Takes a \code{meta_stan} object which is obtained by function \code{meta_stan} and plot
#' a forestplot, showing individual estimates along with their 95 percent credible intervals,
#' resulting effect estimate and prediction interval.
#'
#' @param x A \code{meta_stan} object.
#' @param labels Optional vector with labels for the studies (publication author/year).
#' @param digits A numerical value specifying the number of significant digits to be shown.
#' Default is 2.
#' @param boxsize A numerical value specifying the box size. Default is 0.3.
#' @param heterogeneity A logical value to include heterogeneity estimates (DEFAULT = TRUE)
#' @param col A function specifying the colors. See \code{forestplot::fpColors} for details.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}.
#' @examples
#' \dontrun{
#' data('dat.Crins2014', package = "MetaStan")
#' dat_long <- create_MetaStan_dat(dat = dat.Crins2014,
#'                                     armVars = c(responders = "r", sampleSize = "n"))
#' bnhm.Crins  <- meta_stan(data = dat_long, likelihood = "binomial",
#'                          mu_prior = c(0, 10), theta_prior = c(0, 100),
#'                          tau_prior =  0.5)
#' forest_plot(bnhm.Crins, xlab = "log-OR", labels = dat.Crins2014$publication)
#'
#' }
#' @source This function is based \code{foresplot} function from \code{foresplot}
#' R package.
#' @author Christian Roever and Burak Kuersad Guenhan
#' @seealso \code{foresplot::foresplot}
#' @export
forest_plot = function(x = NULL,
                       labels = NULL,
                       digits = 2,
                       boxsize = 0.3,
                       heterogeneity = TRUE,
                       col,
                       ...) {

  if (!requireNamespace("forestplot", quietly=TRUE))
    stop("required 'forestplot' package not available!")
  if (utils::packageVersion("forestplot") < "1.6")
    warning("you may need to update 'forestplot' to a more recent version (>=1.6).")

  # some sanity checks for the provided arguments:
  stopifnot(is.element("meta_stan", class(x)),
            length(digits)==1, digits==round(digits), digits>=0)

  ## number of studies
  data_wide = x$data_wide

  Nstud = nrow(data_wide)
  if(is.null(labels)) {labels = paste("Study ID:", 1:Nstud)}
  data_wide$labels = labels

  raw.data <- metafor::escalc(measure = "OR",
                              ai = r2, n1i = n2,
                              ci = r1, n2i = n1,
                              slab = labels, data = data_wide)


  # summary estimate
  summary.est <- cbind.data.frame("method" = c("Summary", "Prediction"),
                                  "y"      = x$fit_sum[c("theta", "theta_pred[1]"), "50%"],
                                  "lower"  = x$fit_sum[c("theta", "theta_pred[1]"), "2.5%"],
                                  "upper"  = x$fit_sum[c("theta", "theta_pred[1]"), "97.5%"],
                                  stringsAsFactors = FALSE)

  # summary estimate
  if(heterogeneity == TRUE) {
    mcmc = coda::mcmc.list(rstan::As.mcmc.list(x$fit, pars = c("tau[1]")))
    tau_int = HDInterval::hdi(mcmc, credMass = 0.95)
    temp = matrix(round(c(x$fit_sum['tau[1]', 1], tau_int[1],
                          x$fit_sum['tau[1]', "50%"], tau_int[2]), digits), ncol = 4)

    summary.est <- cbind.data.frame("method" = c("Summary", "Prediction", "Heterogeneity"),
                                    "y"      = c(x$fit_sum[c("theta", "theta_pred[1]"), "50%"], temp[, 3]),
                                    "lower"  = c(x$fit_sum[c("theta", "theta_pred[1]"), "2.5%"], temp[, 2]),
                                    "upper"  = c(x$fit_sum[c("theta", "theta_pred[1]"), "97.5%"], temp[, 4]),
                                    stringsAsFactors = FALSE)
  }



  ################################################################################
  # specify particular summary plotting functions:
  # estimates' colours:
  if (missing(col)) {
    col <- forestplot::fpColors(box=c("black", "grey45"),
                                lines=c("black","grey45"),
                                summary="blue")
    fpDrawPred <- function(lower_limit, estimate, upper_limit, size, col, y.offset = 0.5,...)
    {
      forestplot::fpDrawBarCI(lower=lower_limit, esti=estimate, upper=upper_limit,
                              col = "purple",
                              size=size,
                              y.offset=y.offset,...)
    }

  } else {
    fpDrawPred <- function(lower_limit, estimate, upper_limit, size, col, y.offset = 0.5,...)
    {
      forestplot::fpDrawBarCI(lower=lower_limit, esti=estimate, upper=upper_limit,
                              size=size,
                              y.offset=y.offset,...)
    }


  }


  fpDrawSummary <- function(lower_limit, estimate, upper_limit, size, y.offset = 0.5,...)
  {
    forestplot::fpDrawSummaryCI(lower=lower_limit, esti=estimate, upper=upper_limit,
                                size=size,
                                y.offset=y.offset,...)
  }



  q975 <- stats::qnorm(0.975) # (normal quantile for CIs)


  # need to specify the "tabular" part of the forest plot
  # in terms of a (character) matrix.
  # Here we need 3 columns, 7 rows (including title row):
  format <- "%.2f"     # (common format for "sprintf()" calls)

  mlabtext <- matrix("", nrow = Nstud + 3, ncol = 3)
  mlabtext[1,] <- c("Study", "Estimate", "95% CI")
  mlabtext[,1] <- c("Study", raw.data$labels, "Summary", "Prediction")

  mlabtext[2:(Nstud + 1),2] <- sprintf(format, raw.data$yi)
  mlabtext[2:(Nstud + 1),3] <- sprintf(paste0("[",format,", ",format,"]"),
                                       raw.data$yi-q975*sqrt(raw.data$vi),
                                       raw.data$yi+q975*sqrt(raw.data$vi))

  mlabtext[(Nstud + 2):(Nstud + 3),2] <- sprintf(format, summary.est$y[1:2])
  mlabtext[(Nstud + 2):(Nstud + 3),3] <- sprintf(paste0("[",format,", ",format,"]"),
                                                 summary.est$lower[1:2], summary.est$upper[1:2])


  horizl <- list(grid::gpar(col="black"), grid::gpar(col="darkgrey"), grid::gpar(col="black"))
  names(horizl) <- as.character(c(2,(Nstud + 2),(Nstud + 4)))


  fn.ci_sum <- list(NULL)
  fn.ci_sum[[Nstud + 2]] <- fpDrawSummary
  fn.ci_sum[[Nstud + 3]] <- fpDrawPred
  for(i in 1:(Nstud + 1)) {
    fn.ci_sum[[i]] = forestplot::fpDrawSummaryCI
  }


  res = forestplot::forestplot(labeltext = mlabtext,
                               hrzl = horizl,
                               boxsize = boxsize,
                               new_page = TRUE,
                               align = c("l", "c", "c"),
                               is.summary=c(TRUE, rep(FALSE,Nstud), TRUE, TRUE),
                               mean  = c(NA, raw.data$yi, summary.est[1:2,"y"]),
                               lower = c(NA, raw.data$yi - q975 * sqrt(raw.data$vi),
                                         summary.est[1:2,"lower"]),
                               upper = c(NA, raw.data$yi + q975 * sqrt(raw.data$vi),
                                         summary.est[1:2,"upper"]),
                               graphwidth = grid::unit(80, 'mm'),
                               txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex=1),
                                                            xlab = grid::gpar(cex=0.9)),
                               fn.ci_sum = fn.ci_sum,
                               col = col,
                               ...)

  plot(res)

  # add heterogeneity phrase at bottom left:
  if (heterogeneity) {
    tauFigures <- summary.est[3,c("y", "lower", "upper")]
    tauFigures <- tauFigures[tauFigures > 0]
    formatstring <- paste0("%.", digits, "f")
    tauphrase <- sprintf(paste0("Heterogeneity (tau): ",formatstring,
                                " [",formatstring,", ",formatstring,"]"),
                         summary.est[3,"y"],
                         summary.est[3,"lower"],
                         summary.est[3,"upper"])

    tvp <- grid::viewport(x=grid::unit(0.0, "npc"),
                          y=grid::unit(0.0, "npc"),
                          width=grid::stringWidth(tauphrase),
                          height=grid::unit(2, "lines"),
                          just=c("left","bottom"),
                          name="heterogeneityEstimate")
    grid::pushViewport(tvp)
    grid::grid.text(tauphrase,
                    x=grid::unit(0.0, "npc"),
                    y=grid::unit(0.5, "npc"),
                    just=c("left", "centre"),
                    gp=grid::gpar(fontface="oblique"))
  }



}
