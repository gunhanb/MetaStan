#' Plot a forest plot
#'
#' Takes a \code{meta_stan} object which is obtained by function \code{meta_stan} and plot
#' a forestplot, showing individual estimates along with their 95 percent credible intervals,
#' resulting effect estimate and prediction interval.
#'
#' @param x A \code{meta_stan} object.
#' @param data_wide Dataset in a wide format, one-row-per-study.
#' @param digits A numerical value specifying the number of significant digits to be shown.
#' Default is 2.
#' @param boxsize A numerical value specifying the box size. Default is 0.3.
#' @param col A function specifying the colors. See \code{forestplot::fpColors} for details.
#' @param ... Further arguments passed to or from other methods.
#' @return The return value is invisible \code{NULL}.
#' #' @examples
#' \dontrun{
#' data('dat.Crins2014', package = "MetaStan")
#' dat_long <- convert_data_arm(dat.Crins2014$exp.total, dat.Crins2014$cont.total,
#'                              dat.Crins2014$exp.AR.events, dat.Crins2014$cont.AR.events)
#' bnhm.Crins  <- meta_stan(data = dat_long, family = "binomial",
#'                          mu_prior = c(0, 10), theta_prior = c(0, 100),
#'                          tau_prior =  0.5)
#' forestplot(bnhm.Crins)
#'
#' }
#' @source This function is based \code{foresplot} function from \code{foresplot}
#' R package.
#' @author Christian Roever and Burak Kuersad Guenhan
#' @seealso \code{foresplot::foresplot}
#' @export
forestplot.meta_stan = function(x = NULL,
                               data_wide = NULL,
                               digits = 2,
                               boxsize = 0.3,
                               col,
                               ...) {

  if (!requireNamespace("forestplot", quietly=TRUE))
    stop("required 'forestplot' package not available!")
  if (utils::packageVersion("forestplot") < "1.6")
    warning("you may need to update 'forestplot' to a more recent version (>=1.6).")
  # auxiliary function:
  decplaces <- function(x, signifdigits=digits)

    # some sanity checks for the provided arguments:
    stopifnot(is.element("meta_stan", class(x)),
              length(digits)==1, digits==round(digits), digits>=0)

  ## number of studies
  Nstud = nrow(data_wide)
  raw.data <- escalc(measure="OR",
                     ai=exp.AR.events,  n1i=exp.total,
                     ci=cont.AR.events, n2i=cont.total,
                     slab=publication, data= data_wide)

  # summary estimate
  summary.est <- cbind.data.frame("method" = c("Summary", "Prediction"),
                                  "y"      = summary(x$fit)$summary[c("d", "d_pred"), "50%"],
                                  "lower"  = summary(x$fit)$summary[c("d", "d_pred"), "2.5%"],
                                  "upper"  = summary(x$fit)$summary[c("d", "d_pred"), "97.5%"],
                                  stringsAsFactors=FALSE)


  ################################################################################
  # specify particular summary plotting functions:
  # estimates' colours:
  if (missing(col)) {
    col <- forestplot::fpColors(box=c("black", "grey45"),
                                lines=c("black","grey45"),
                                summary="blue")
    fpDrawPred <- function(lower_limit, estimate, upper_limit, size, col, y.offset = 0.5,...)
    {
      forestplot::fpDrawSummaryCI(lower=lower_limit, esti=estimate, upper=upper_limit, col = "purple",
                      size=size,
                      y.offset=y.offset,...)
    }

  } else {
    fpDrawPred <- function(lower_limit, estimate, upper_limit, size, y.offset = 0.5,...)
    {
      forestplot::fpDrawSummaryCI(lower=lower_limit, esti=estimate, upper=upper_limit, size=size,
                      y.offset=y.offset,...)
    }


  }


  fpDrawSummary <- function(lower_limit, estimate, upper_limit, size, y.offset = 0.5,...)
  {
    forestplot::fpDrawSummaryCI(lower=lower_limit, esti=estimate, upper=upper_limit, size=size,
                    y.offset=y.offset,...)
  }



  q975 <- qnorm(0.975) # (normal quantile for CIs)


  # need to specify the "tabular" part of the forest plot
  # in terms of a (character) matrix.
  # Here we need 3 columns, 7 rows (including title row):
  format <- "%.2f"     # (common format for "sprintf()" calls)

  mlabtext <- matrix("", nrow=Nstud + 3, ncol = 3)
  mlabtext[1,] <- c("Study", "Estimate", "95% CI")
  mlabtext[,1] <- c("Study", raw.data$publication, "Summary", "Prediction")

  mlabtext[2:(Nstud + 1),2] <- sprintf(format, raw.data$yi)
  mlabtext[2:(Nstud + 1),3] <- sprintf(paste0("[",format,", ",format,"]"),
                                       raw.data$yi-q975*sqrt(raw.data$vi), raw.data$yi+q975*sqrt(raw.data$vi))

  mlabtext[(Nstud + 2):(Nstud + 3),2] <- sprintf(format, summary.est$y)
  mlabtext[(Nstud + 2):(Nstud + 3),3] <- sprintf(paste0("[",format,", ",format,"]"),
                                                 summary.est$lower, summary.est$upper)


  horizl <- list(grid::gpar(col="black"),grid::gpar(col="darkgrey"),grid::gpar(col="black"))
  names(horizl) <- as.character(c(2,(Nstud + 2),(Nstud + 4)))


  fn.ci_sum <- list(NULL)
  fn.ci_sum[[Nstud + 2]] <- fpDrawSummary
  fn.ci_sum[[Nstud + 3]] <- fpDrawPred
  for(i in 1:(Nstud + 1)) {
    fn.ci_sum[[i]] = fpDrawSummaryCI
  }


  forestplot::forestplot(labeltext = mlabtext,
             hrzl = horizl,
             boxsize = boxsize,
             new_page = TRUE, #  <-- may need to set this argument for PDF or EPS plots
             align = c("l", "c", "c"),
             is.summary=c(TRUE, rep(FALSE,Nstud), TRUE, TRUE), # first line (title) is a "summary"
             # (printed in bold),
             mean  = c(NA, raw.data$yi, summary.est[,"y"]),
             lower = c(NA, raw.data$yi-q975*sqrt(raw.data$vi), summary.est[,"lower"]),
             upper = c(NA, raw.data$yi+q975*sqrt(raw.data$vi), summary.est[,"upper"]),
             xlab="log-OR",                                   # (x-axis label)
             graphwidth = unit(80, 'mm'),
             txt_gp = fpTxtGp(ticks = gpar(cex=1), xlab = gpar(cex=0.9)), # (increase axis label size a bit)
             fn.ci_sum=fn.ci_sum,
             col = col,
             ...)                        # specify x-axis ticks (if necssary)

}
