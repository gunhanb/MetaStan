#' Trials investigating effectiveness of the BCG vaccine against TB
#'
#' A dataset containing the results from 13 trials examining the efficacy
#' of Bacillus Calmette-Guerin (BCG) vaccine against tuberculosis (TB).
#'
#' @format A data frame with following coloumns
#' \describe{
#'   \item{Trial}{Trial number}
#'   \item{rt}{number of TB events in treatment arm}
#'   \item{nt}{number of subjects in treatment arm}
#'   \item{rc}{number of TB events in control arm}
#'   \item{nc}{number of subjects in control arm}
#'   \item{Latitude}{absolute latitude of the study location}
#' }
#' @source Berkey, C.S., Hoaglin, D.C., Mosteller, F. and Colditz,
#' G.A., 1995. A random-effects regression model for meta-analysis. 
#' Statistics in medicine, 14(4), pp.395-411
"dat.Berkey1995"