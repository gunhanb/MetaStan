#' Prepare model-based meta-analysis dataset for Stan.
#'
#' \code{create_MetaStan_dat} converts datasets in the one-study-per-row
#' format to one-arm-per-row format,
#'
#' The resulting data.frame can be used as data argument in \code{MBMA_stan}.
#'
#' @param dat Data in one-study-per-row format.
#' @param armVars Vector of per-arm variables
#'  The name of each component will be the column name in the resulting dataset.
#' @param nArmsVar Variable holding the number of arms for each study.
#' @return A data frame with the generated columns.
#' @author Burak Kuersad Guenhan, \email{burak.gunhan@med.uni-goettingen.de}
#' and Gert van Valkenhoef
#' @seealso \code{gemtc::mtc.data.studyrow} and \code{nmaINLA::create_INLA_dat}
#' @examples
#' \dontrun{
#' data('dat.Eletriptan')
#' ## Create the dataset suitable for MBMA_stan
#' EletriptanDat <- create_MetaStan_dat(dat = dat.Eletriptan,
#'                                      armVars = c(dose = "d",
#'                                                  responders = "r",
#'                                                  sampleSize = "n"),
#'                                      nArmsVar = 'nd')
#' ## Check that the data are correct
#' print(EletriptanDat)
#' }
#' @export
create_MetaStan_dat <- function(dat = NULL,
                                armVars = c(dose = "d",
                                            responders = "r",
                                            sampleSize = "n"),
                                nArmsVar = "nd") {

  if(is.null(dat$nd)) {dat$nd = 2}

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
  datStan <- do.call(data.frame, out)
  #####################      mtc.data.studyrow is finished
  # Adding indicator variable
  datStan$na <- rep(dat[[paste(nArmsVar)]], times = dat[[paste(nArmsVar)]])

  final = list(data_long = datStan,
               data_wide = dat)


  return(final)
}

