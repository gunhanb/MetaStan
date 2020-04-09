#' Convert contrast-based dataset to arm-based dataset
#'
#' \code{convert_data_arm} creates a dataframe to fit a Binomial-Normal
#' Hierarchical model using \code{glmer} function.
#'
#' @param nt Number of subjects in treatment arm
#' @param nc Number of subjects in control arm
#' @param pt Number of events in treatment arm
#' @param pc Number of events in treatment arm
#' @return A dataframe object
#' @examples
#' data('dat.Crins2014', package = "MetaStan")
#' ## Subset of dataset where PTLD outcomes available
#' dat.Crins2014.PTLD = subset(dat.Crins2014, is.finite(exp.PTLD.events))
#' ## Create arm-based dataset
#' dat.Crins2014.PTLD.arm <- convert_data_arm(dat.Crins2014.PTLD$exp.total,
#' dat.Crins2014.PTLD$cont.total,dat.Crins2014.PTLD$exp.PTLD.events,
#' dat.Crins2014.PTLD$cont.PTLD.events)
#'
#'
#' @export
convert_data_arm <- function(nt, nc, pt, pc) {
    data <- NULL
    data$pt <- pt
    data$nt <- nt
    data$pc <- pc
    data$nc <- nc
    data <- data.frame(data)
    N <- nrow(data)
    r <- as.vector(rbind(data$pc, data$pt))  # number of events
    sampleSize <- as.vector(rbind(data$nc, data$nt))  # number of all patients
    theta <- rep(0:1, times = N)
    theta12 <- rep((-0.5):0.5, times = N)             # Needed for Smith model
    het <- as.vector(rbind(rep(NA, times = N), 1:N))  # ID for random effects
    # Dataset for arm-level meta-analysis
    data.arm <- data.frame(cbind(r, sampleSize, theta, theta12, het))
    data.arm$mu <- as.factor(as.numeric(gl(n = N, k = 2)))
    return(data.arm)
}
