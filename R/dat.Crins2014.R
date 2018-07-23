#' Pediatric liver transplant example data
#'
#' Numbers of cases and events (PTLDs or deaths)
#' in experimental and control groups of six studies.
#'
#' @format A data frame with following coloumns
#' \describe{
#'   \item{publication}{publication identifier (first author and publication year)}
#'   \item{year}{publication year}
#'   \item{randomized}{randomization status (y/n)}
#'   \item{control.type}{type of control group ("concurrent" or "historical")}
#'   \item{comparison}{type of comparison ("IL-2RA only", "delayed CNI", or "no/low steroids")}
#'   \item{followup}{t	 follow-up time in months}
#'   \item{exp.PTLD.events}{number of PTLD events in experimental group}
#'   \item{exp.death.events}{number of deaths in experimental group}
#'   \item{exp.total}{number of patients in experimental group}
#'   \item{cont.PTLD.events}{number of PTLD events in control group}
#'   \item{cont.death.events}{number of deaths in control group}
#'   \item{cont.total}{number of patients in control group}
#'  }
#' @source N.D. Crins, C. Roever, A.D. Goralczyk, T. Friede. Interleukin-2
#' receptor antagonists for pediatric liver transplant recipients:
#' A systematic review and meta-analysis of controlled studies. Pediatric
#' Transplantation, 18(8):839-850, 2014.
"dat.Crins2014"

