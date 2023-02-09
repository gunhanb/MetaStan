#' Prepare meta-analysis dataset for meta_stan function
#'
#' \code{metastan_data} creates datasets suitable for meta_stan function
#'
#' The resulting data.frame can be used as data argument in \code{meta_stan}.
#'
#' @param id           # a vector of study IDs (labels)
#' @param treatment    # (optional) vector indicating treatment groups
#' @param x            # (optional) covariable/regressor vector or matrix
#' @param type         # type of outcome
#' @param y            # estimates (for normal outcomes)
#' @param se           # associated standard errors (normal outcomes)
#' @param v            # variances / squared standard errors (normal outcomes)
#' @param count        # event count (for binomial or Poisson outcomes)
#' @param total        # sample size (for binomial outcomes)
#' @param exposure     # exposure (for Poisson outcomes)
#' @param labels       # (optional) vector of row labels
#' @param data         # a data frame from which (some of) above variables may be taken
#' @param tidy         # if TRUE (the default), output rows will be sorted (keeping entries corresponding to the same study ("id") next to each other, and sorting by "treatment" and "x")
#' @return The returned object is a "list" of class "metastan_data"
#'  containing the following elements:
#'    * "id"        :  study identifier (a vector of type factor)
#'    * "treatment" :  (optional) a treatment identifier (a vector of type factor, or NULL)
#'    * "x"         :  (optional) covariable(s), (a numeric matrix, or NULL)
#'    * "type"      :  a flag identifying the outcome type ("normal", "binomial", or "poisson")
#'    * "outcome"   :  the actual outcome data (a two-column numeric matrix)
#'    * "call"      :  the original function call
#' @details   NB: arguments "id", "treatment", "x", "y", "se", "v", "count",
#' "total", "exposure", "labels"
#'       may be taken from local environment OR from "data" (data frame) argument.
#'
#' @examples
#' # 3 studies, binomial endpoint, one covariable (named "dose"):
#' msd <- metastan_data(id=c("Smith","Smith", "Taylor","Taylor",
#'                           "Jones", "Jones","Taylor"),
#'                      treatment=c("placebo","verum","placebo","verum",
#'                                  "verum","placebo","verum"),
#'                      type="binomial",
#'                      count=c(3, 2, 7, 5, 4, 5, 6),
#'                      total=c(10, 10, 25, 27, 15, 16, 23),
#'                      x=cbind("dose"=c(0, 50, 0, 50, 30, 0, 20)))
#' msd
#' print(msd, n=7)     # (show all 7 lines of data)
#' print.default(msd)  # default view of returned object
#' @export
metastan_data <- function(id,
                          treatment,
                          x,
                          type = c("normal", "binomial", "poisson"),
                          y, se, v,     # for normal outcomes
                          count, total, # for binomial outcomes (event count, sample size)
                          exposure,     # for Poisson outcomes (count, exposure)
                          labels,       # (optional) vector of row labels
                          data,         # a data frame from which (some of) above variables may be taken
                          tidy=TRUE,
                          checkForConflicts=TRUE)
{
  ##############################
  # gather data from arguments;
  # check whether a "data" argument was supplied:
  if (missing(data))
    data <- NULL
  if (is.null(data)) {
    dataArg <- FALSE
    data <- sys.frame(sys.parent())
  } else {
    dataArg <- TRUE
    if (!is.data.frame(data))
      data <- data.frame(data)
  }
  mc <- match.call()
  # gather arguments from environment _OR_ from "data":
  mc.id <- mc[[match("id", names(mc))]]
  id <- eval(mc.id, data, enclos=sys.frame(sys.parent()))
  mc.treatment <- mc[[match("treatment", names(mc))]]
  treatment <- eval(mc.treatment, data, enclos=sys.frame(sys.parent()))
  mc.x <- mc[[match("x", names(mc))]]
  x <- eval(mc.x, data, enclos=sys.frame(sys.parent()))
  mc.y <- mc[[match("y", names(mc))]]
  y <- eval(mc.y, data, enclos=sys.frame(sys.parent()))
  mc.se <- mc[[match("se", names(mc))]]
  se <- eval(mc.se, data, enclos=sys.frame(sys.parent()))
  mc.v <- mc[[match("v", names(mc))]]
  v <- eval(mc.v, data, enclos=sys.frame(sys.parent()))
  mc.count <- mc[[match("count", names(mc))]]
  count <- eval(mc.count, data, enclos=sys.frame(sys.parent()))
  mc.total <- mc[[match("total", names(mc))]]
  total <- eval(mc.total, data, enclos=sys.frame(sys.parent()))
  mc.exposure <- mc[[match("exposure", names(mc))]]
  exposure <- eval(mc.exposure, data, enclos=sys.frame(sys.parent()))
  mc.labels <- mc[[match("labels", names(mc))]]
  labels <- eval(mc.labels, data, enclos=sys.frame(sys.parent()))
  ###################################
  # check for (potential) conflicts:
  if (dataArg & checkForConflicts) {
    dataName <- deparse(substitute(data))
    id2 <- try(eval(mc.id), silent=TRUE)
    if ((!is.element("try-error", class(id2))) && (!identical(id, id2)))
      warning(paste0("potential conflict -- interpreted \"id=", deparse(substitute(mc.id)),
                     "\" argument within data frame \"", dataName, "\""))
    treatment2 <- try(eval(mc.treatment), silent=TRUE)
    if ((!is.element("try-error", class(treatment2))) && (!identical(treatment, treatment2)))
      warning(paste0("potential conflict -- interpreted \"id=", deparse(substitute(mc.treatment)),
                     "\" argument within data frame \"", dataName, "\""))
    x2 <- try(eval(mc.x), silent=TRUE)
    if ((!is.element("try-error", class(x2))) && (!identical(x, x2)))
      warning(paste0("potential conflict -- interpreted \"x=", deparse(substitute(mc.x)),
                     "\" argument within data frame \"", dataName, "\""))
    y2 <- try(eval(mc.y), silent=TRUE)
    if ((!is.element("try-error", class(y2))) && (!identical(y, y2)))
      warning(paste0("potential conflict -- interpreted \"y=", deparse(substitute(mc.y)),
                     "\" argument within data frame \"", dataName, "\""))
    se2 <- try(eval(mc.se), silent=TRUE)
    if ((!is.element("try-error", class(se2))) && (!identical(se, se2)))
      warning(paste0("potential conflict -- interpreted \"se=", deparse(substitute(mc.se)),
                     "\" argument within data frame \"", dataName, "\""))
    v2 <- try(eval(mc.v), silent=TRUE)
    if ((!is.element("try-error", class(v2))) && (!identical(v, v2)))
      warning(paste0("potential conflict -- interpreted \"v=", deparse(substitute(mc.v)),
                     "\" argument within data frame \"", dataName, "\""))
    count2 <- try(eval(mc.count), silent=TRUE)
    if ((!is.element("try-error", class(count2))) && (!identical(count, count2)))
      warning(paste0("potential conflict -- interpreted \"count=", deparse(substitute(mc.count)),
                     "\" argument within data frame \"", dataName, "\""))
    total2 <- try(eval(mc.total), silent=TRUE)
    if ((!is.element("try-error", class(total2))) && (!identical(total, total2)))
      warning(paste0("potential conflict -- interpreted \"total=", deparse(substitute(mc.total)),
                     "\" argument within data frame \"", dataName, "\""))
    exposure2 <- try(eval(mc.exposure), silent=TRUE)
    if ((!is.element("try-error", class(exposure2))) && (!identical(exposure, exposure2)))
      warning(paste0("potential conflict -- interpreted \"exposure=", deparse(substitute(mc.exposure)),
                     "\" argument within data frame \"", dataName, "\""))
    labels2 <- try(eval(mc.labels), silent=TRUE)
    if ((!is.element("try-error", class(labels2))) && (!identical(labels, labels2)))
      warning(paste0("potential conflict -- interpreted \"labels=", deparse(substitute(mc.labels)),
                     "\" argument within data frame \"", dataName, "\""))

    rm(list=c("id2", "y2", "se2", "v2", "count2",
              "total2", "exposure2", "x2", "labels2"))
  }
  rm(list=c("mc.id", "mc.treatment", "mc.x", "mc.y", "mc.se", "mc.v",
            "mc.count", "mc.total", "mc.exposure", "mc.labels"))
  ##########
  type <- match.arg(type)
  # some sanity checks:
  stopifnot(!is.null(id), is.vector(id), length(id)>0)
  if (!is.null(treatment))
    stopifnot(is.vector(treatment), (length(treatment)==length(id)))
  stopifnot(is.vector(tidy), is.logical(tidy), length(tidy)==1, !is.na(tidy))
  # process outcome data:
  if (type=="binomial") { # (proper binomial counts & totals required)
    stopifnot(!is.null(count) & !is.null(total),
              is.vector(count), is.numeric(count), all(is.finite(count)),
              all(count>=0), length(count)==length(id),
              is.vector(total), is.numeric(total), all(is.finite(total)),
              all(total>0), all(total>=count),
              length(total)==length(id))
    outcome <- matrix(NA_integer_, ncol=2, nrow=length(id),
                      dimnames=list(NULL, c("count","total")))
    outcome[,"count"] <- count
    outcome[,"total"] <- total
  }
  if (type=="poisson") { # (proper Poisson counts & exposures required)
    stopifnot(!is.null(count) & !is.null(exposure),
              is.vector(count), is.numeric(count), all(is.finite(count)),
              all(count>=0), length(count)==length(id),
              is.vector(exposure), is.numeric(exposure), all(is.finite(exposure)),
              all(exposure>0), length(exposure)==length(id))
    outcome <- matrix(NA_integer_, ncol=2, nrow=length(id),
                      dimnames=list(NULL, c("count","exposure")))
    outcome[,"count"]    <- count
    outcome[,"exposure"] <- exposure
  }
  if (type=="normal") {  # (estimates and standard errors / variances required)
    stopifnot(!is.null(y) & (!is.null(v) | !is.null(se)),
              is.vector(y), is.numeric(y), all(is.finite(y)),
              length(y)==length(id))
    if (!is.null(v)) {  # check variances, convert to standard errors:
      stopifnot(is.vector(v), is.numeric(v), all(is.finite(v)),
                all(v>0), length(v)==length(id))
      se <- sqrt(v)
    }                    # check standard errors:
    stopifnot(is.vector(se), is.numeric(se), all(is.finite(se)),
              all(se>0), length(se)==length(id))
    outcome <- matrix(NA_integer_, ncol=2, nrow=length(id),
                      dimnames=list(NULL, c("y","se")))
    outcome[,"y"]  <- y
    outcome[,"se"] <- se
  }
  # check "id" vector, convert to factor (if necessary):
  if (!is.factor(id)) {
    id <- as.factor(id)
  } else {  # check whether all levels are actually used:
    if (length(unique(id)) != length(levels(id)))
      warning("\"id\" argument contains un-used factor levels.")
  }
  # check "treatment" vector, convert to factor (if necessary):
  if (!is.null(treatment)) {
    if (!is.factor(treatment)) {
      treatment <- as.factor(treatment)
    } else {  # check whether all levels are actually used:
      if (length(unique(treatment)) != length(levels(treatment)))
        warning("\"treatment\" argument contains un-used factor levels.")
    }
  }
  # set up regressor matrix
  # (NB: always a matrix, but commonly just a single column):
  if (!is.null(x)) {
    stopifnot(is.vector(x) | is.matrix(x))
    if (is.vector(x)) {
      stopifnot(is.numeric(x), all(is.finite(x)), length(x)==length(id))
      x <- matrix(x, nrow=length(id), ncol=1,
                  dimnames=list(NULL, "x"))
    } else if (is.matrix(x)) {
      stopifnot(is.numeric(x), all(is.finite(x)), nrow(x)==length(id))
      matcolnames <- colnames(x)
      colnames(x) <- make.names(matcolnames, unique=TRUE)
    }
  } else { # (the default)
    x <- NULL
  }
  if (tidy) { # sort output rows sensibly:
    # create a "sorted ID" vector (in order to keep studies "in order of appearance"):
    sortedID <- factor(id, levels=as.character(unique(id)))
    if (!is.null(treatment)) {
      if (!is.null(x)) {
        rowOrder <- order(sortedID, treatment, x[,1])
      } else {
        rowOrder <- order(sortedID, treatment)
      }
    } else {
      if (!is.null(x)) {
        rowOrder <- order(sortedID, x[,1])
      } else {
        rowOrder <- order(sortedID)
      }
    }
    # permute/sort rows of output elements:
    id <- id[rowOrder]
    if (!is.null(treatment)) treatment <- treatment[rowOrder]
    if (!is.null(x)) x <- x[rowOrder,,drop=FALSE]
    outcome <- outcome[rowOrder,,drop=FALSE]
    if (!is.null(labels)) labels <- labels[rowOrder]
  }
  # generate a sensible row label vector (essentially for "internal use" only):
  if (!is.null(labels)) {
    stopifnot(is.vector(labels), is.character(labels),
              !any(duplicated(labels)))
  } else {
    labels <- make.names(as.character(id))
    countvec <- rep(NA_integer_, length(id))
    for (i in 1:nlevels(id)) {
      idx <- which(id == levels(id)[i])
      if (length(idx)>0) {
        countvec[idx] <- 1:length(idx)
        # check for duplicate regressor entries:
        if (length(idx)>1) {
          if (is.null(treatment) & is.null(x)) {
            warning("Duplicate entries in \"id\" argument (id=\"",levels(id)[i],"\")!")
          } else if (any(duplicated(cbind(as.numeric(treatment), x)[idx,]))) {
            # (NB: "duplicated()" works for vectors and matrices/data.frames (row-wise)!)
            warning("Duplicate covariable entries (id=\"",levels(id)[i],"\")!")
          }
        }
      }
    }
    if (max(table(id)) > 1) {
      labels <- paste(labels, sprintf("%.02d",countvec), sep=".")
    }
  }
  # apply row labels throughout:
  rownames(outcome) <- labels
  if (!is.null(treatment)) names(treatment) <- labels
  if (!is.null(x)) rownames(x) <- labels
  # set up object to hold results:
  result <- list("id"        = id,        # study ID
                 "treatment" = treatment, # (optional) treatment ID
                 "x"         = x,         # (optional) covariable(s)
                 "type"      = type,      # type of endpoint
                 "outcome"   = outcome,   # outcome data
                 "call"      = match.call(expand.dots=FALSE))
  class(result) <- "metastan_data"
  return(result)
}



################################################################################
# "print()" method for a "metastan_data" object:
#

print.metastan_data <- function(x, n=6, ...)
{
  cat("\"metastan_data\" object.\n")
  cat(paste0("outcome type         : \"", x$type, "\"\n"))
  cat(paste0("number of studies    : ", length(unique(x$id)), "\n"))
  cat(paste0("number of treatments : ",
             ifelse(is.null(x$treatment), 1, length(unique(x$treatment))),
             "\n"))
  tab <- table(x$id)
  if (length(unique(tab))==1)
    cat(paste0("outcomes per study   : ", tab[1], "\n"))
  else
    cat(paste0("outcomes per study   : ", min(tab), "-", max(tab), "\n"))
  if (!is.null(x$x)) { # (only printed if present)
    cat(paste0("covariable(s)        : ", as.character(paste0("\"",colnames(x$x),
                                                              "\"",collapse=", ")), "\n"))
  }
  cat("top of data set:\n")
  dat <- cbind.data.frame("id"=as.character(x$id))
  if (!is.null(x$treatment)) {
    dat <- cbind.data.frame(dat, "treatment"=x$treatment)
  }
  if (!is.null(x$x)) {
    dat <- cbind.data.frame(dat, x$x)
  }
  dat <- cbind.data.frame(dat, x$outcome)
  rownames(dat) <- NULL
  print(head(dat, n=n))
  invisible(x)
}
