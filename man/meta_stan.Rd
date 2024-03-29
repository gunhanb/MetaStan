% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_stan.R
\name{meta_stan}
\alias{meta_stan}
\title{Fitting a meta-analysis model using Stan}
\usage{
meta_stan(
  data = NULL,
  likelihood = NULL,
  mu_prior = c(0, 10),
  theta_prior = NULL,
  tau_prior = 0.5,
  tau_prior_dist = "half-normal",
  beta_prior = c(0, 100),
  delta = NULL,
  param = "Smith",
  re = TRUE,
  ncp = TRUE,
  interval.type = "shortest",
  mreg = FALSE,
  cov = NULL,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  adapt_delta = 0.95,
  ...
)
}
\arguments{
\item{data}{Data frame created by `create_MetaStan_dat`}

\item{likelihood}{A string specifying the likelihood function defining the statistical
model. Options include  `normal`, `binomial`, and `Poisson`.}

\item{mu_prior}{A numerical vector specifying the parameter of the normal prior
density for baseline risks, first value is parameter for mean, second is for variance.
Default is c(0, 10).}

\item{theta_prior}{A numerical vector specifying the parameter of the normal prior
density for treatment effect estimate, first value is parameter for mean, second
is for variance. Default is NULL.}

\item{tau_prior}{A numerical value specifying the standard dev. of the prior density
for heterogeneity stdev. Default is 0.5.}

\item{tau_prior_dist}{A string specifying the prior density for the heterogeneity standard deviation,
option is `half-normal` for half-normal prior, `uniform` for uniform prior, `half-cauchy` for
half-cauchy prior.}

\item{beta_prior}{A numerical vector specifying the parameter of the normal prior
density for beta coefficients in a meta-regression model, first value is parameter for mean, second
is for variance. Default is c(0, 100).}

\item{delta}{A numerical value specifying the upper bound of the a priori interval for
treatment effect on odds ratio scale (\emph{Guenhan et al (2020)}). This is used to calculate
a normal weakly informative prior.
for theta. Thus when this argument is specified, `theta` should be left empty. Default is NULL.}

\item{param}{Paramteriztaion used. The default is the `Smith` model suggested by
Smith et al (1995). The alternative is `Higgins` is the common meta-analysis
model (Simmonds and Higgins, 2014).}

\item{re}{A string specifying whether random-effects are included to the model. When `FALSE`, the
model corresponds to a fixed-effects model. The default is `TRUE`.}

\item{ncp}{A string specifying whether to use a non-centered parametrization.
The default is `TRUE`.}

\item{interval.type}{A string specifying the type of interval estimate. Options
 include
shortest credible interval `shortest` (default) and qui-tailed credible interval
`central`.}

\item{mreg}{A string specifying whether to fit a meta-regression model.
The default is `FALSE`.}

\item{cov}{A numeric vector or matrix specifying trial-level covariates (in each row).
This is needed when `mreg = TRUE`.}

\item{chains}{A positive integer specifying the number of Markov chains.
The default is 4.}

\item{iter}{A positive integer specifying the number of iterations for each chain
(including warmup). The default is 2000.}

\item{warmup}{A positive integer specifying the number of warmup (aka burnin)
iterations per chain. The default is 1000.}

\item{adapt_delta}{A numerical value specifying the target average proposal acceptance
probability for adaptation. See Stan manual for details. Default is 0.95. In general
you should not need to change adapt_delta unless you see a warning message about
divergent transitions, in which case you can increase adapt_delta from the
default to a value closer to 1 (e.g. from 0.95 to 0.99, or from 0.99 to 0.999, etc).}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
an object of class `MetaStan`.
}
\description{
`meta_stan` fits a meta-analysis model using Stan.
}
\examples{
\dontrun{

## TB dataset
data('dat.Berkey1995', package = "MetaStan")
## Fitting a Binomial-Normal Hierarchical model using WIP priors
dat_MetaStan <- create_MetaStan_dat(dat = dat.Berkey1995,
                                    armVars = c(responders = "r", sampleSize = "n"))

 ma.stan  <- meta_stan(data = dat_MetaStan,
                           likelihood = "binomial",
                           mu_prior = c(0, 10),
                           theta_prior = c(0, 100),
                           tau_prior = 0.5,
                           tau_prior_dist = "half-normal")
print(ma.stan)
forest_plot(ma.stan)


meta.reg.stan  <- meta_stan(data = dat_MetaStan,
                           likelihood = "binomial",
                           mu_prior = c(0, 10),
                           theta_prior = c(0, 100),
                           tau_prior = 0.5,
                           tau_prior_dist = "half-normal",
                           mreg = TRUE,
                           cov = dat.Berkey1995$Latitude)

print(meta.reg.stan)
}

}
\references{
Guenhan BK, Roever C, Friede T. MetaStan: An R package for meta-analysis
and model-based meta-analysis using Stan. In preparation.

Guenhan BK, Roever C, Friede T. Random-effects meta-analysis of
few studies involving
rare events \emph{Resarch Synthesis Methods} 2020; doi:10.1002/jrsm.1370.

Jackson D, Law M, Stijnen T, Viechtbauer W, White IR. A comparison of 7
random-effects models for meta-analyses that estimate the summary odds ratio.
\emph{Stat Med} 2018;37:1059--1085.

Kuss O. Statistical methods for meta-analyses including information
from studies without any events-add nothing to nothing and succeed nevertheless,
\emph{Stat Med}, 2015; 4; 1097--1116, doi: 10.1002/sim.6383.
}
