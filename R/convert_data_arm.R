#' Convert contrast-based dataset to arm-based dataset
#'
#' \code{convert_data_arm} creates a dataframe to fit a meta-analysis
#' model using \code{meta_Stan} function.
#'
#' @param nt Number of subjects in treatment arm
#' @param nc Number of subjects in control arm
#' @param pt Number of events in treatment arm
#' @param pc Number of events in treatment arm
#' @param pub The corresponding publication
#' @param data Optional data frame containing the variables given to the arguments above.
#' @return A dataframe object
#' @examples
#' data('dat.Crins2014', package = "MetaStan")
#' ## Subset of dataset where PTLD outcomes available
#' dat.Crins2014.PTLD = subset(dat.Crins2014, is.finite(exp.PTLD.events))
#' ## Create arm-based dataset
#' data('dat.Crins2014', package = "MetaStan")
#' dat_converted <- convert_data_arm(exp.total, cont.total,
#'                                   exp.AR.events, cont.AR.events,
#'                                   publication, data = dat.Crins2014)
#'
#'
#' @export
convert_data_arm <- function(nt, nc, pt, pc, pub, data = NULL) {

    ################ check data argument
    if (is.null(data))
        data <- sys.frame(sys.parent())
    data_wide <- match.call()
    data_wide$data <- NULL
    data_wide[[1]] <- as.name("data.frame")
    data_wide <- eval(data_wide,data)

    N     <- nrow(data_wide)
    r     <- as.vector(rbind(data_wide$pc, data_wide$pt))  # number of events
    n     <- as.vector(rbind(data_wide$nc, data_wide$nt))  # number of all patients
    theta <- rep(0:1, times = N)
    # Dataset for arm-level meta-analysis
    data_long    <- data.frame(cbind(r, n, theta))
    data_long$mu <- as.factor(as.numeric(gl(n = N, k = 2)))

    final = list(data_long = data_long,
                 data_wide = data_wide)

    return(final)
}
