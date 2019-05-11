#' Prepare model-based meta-analysis dataset for Stan.
#'
#' \code{create_MBMA_dat} converts datasets in the one-study-per-row
#' format to one-arm-per-row format,
#'
#' The resulting data.frame can be used as data argument in \code{MBMA_stan}.
#'
#' @param dat Data in one-study-per-row format.
#' @param armVars Vector of per-arm variables
#'  The name of each component will be the column name in the resulting dataset.
#' @param nArmsVar Variable holding the number of arms for each study.
#' @return A data frame with the generated coloumns.
#' @author Burak Kuersad Guenhan, \email{burak.gunhan@med.uni-goettingen.de}
#' and Gert van Valkenhoef
#' @seealso \code{gemtc::mtc.data.studyrow} and \code{nmaINLA::create_INLA_dat}
#' @examples
#' data('dat.Eletriptan')
#' ## Create the dataset suitable for MBMA_stan
#' EletriptanDat <- create_MBMA_dat(dat = dat.Eletriptan,
#' armVars = c("dose" = "d", "r" = "r", "n" = "n"), nArmsVar = 'nd')
#' ## Check that the data are correct
#' print(EletriptanDat)
#' @export
create_MBMA_dat <- function(dat = dat,
                            armVars = c(dose = "t", responders = "r", sampleSize = "n"),
                            nArmsVar = "nd") {
  ######################################################## THIS CODE IS COPY PASTE FROM gemtc::mtc.data.studyrow
  studyNames = 1:nrow(dat)
  patterns = c("%s", "%s%d")
  studyVars = c()
  treatmentNames = NA
  colsOrNA <- function(row, cols) {
    rval <- rep(NA, length(cols))
    sel <- cols %in% colnames(row)
    rval[sel] <- row[cols[sel]]
    rval
  }

  nArmsCol <- sprintf(patterns[1], nArmsVar)
  studyCols <- sprintf(patterns[1], studyVars)

  out <- do.call(rbind, lapply(1:nrow(dat), function(i) {
    row <- dat[i, ]
    na <- row[nArmsCol]
    studyEntries <- row[studyCols]
    names(studyEntries) <- names(studyVars)
    do.call(rbind, lapply(1:unlist(na), function(j) {
      armCols <- sprintf(patterns[2], armVars, j)
      armEntries <- colsOrNA(row, armCols)
      names(armEntries) <- names(armVars)
      c(list(study = i), studyEntries, armEntries)
    }))
  }))

  colNames <- colnames(out)
  out <- lapply(colnames(out), function(col) {
    unlist(out[, col])
  })
  names(out) <- colNames

  out[["study"]] <- studyNames[out[["study"]]]
  if (all(!is.na(treatmentNames))) {
    out[["treatment"]] <- treatmentNames[out[["treatment"]]]
  }
  datMBMA <- do.call(data.frame, out)
  #####################      mtc.data.studyrow is finished
  # Adding indicator variable needed for INLA
  datMBMA$na <- rep(dat[[paste(nArmsVar)]], times = dat[[paste(nArmsVar)]])
  N <- nrow(datMBMA)

  return(datMBMA)
}

