
context("Checking meta-analysis example: Crins dataset")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model
test_that("Results are correct for fitting binomial normal hierarchical model using WIP priors.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Berkey1995', package = "MetaStan")
  ## Fitting a Binomial-Normal Hierarchical model using WIP priors

  dat_long = create_MetaStan_dat(dat = dat.Berkey1995,
                                 armVars = c(responders = "r",
                                             sampleSize = "n"),
                                 nArmsVar = "nd")

  bnhm.wip.bnhm1.stan  <- meta_stan(data = dat_long,
                                    likelihood = "binomial",
                                    mu_prior = c(0, 10),
                                    delta = 250,
                                    tau_prior = 0.5,
                                    tau_prior_dist = "half-normal")
  ### compare with results
  results = bnhm.wip.bnhm1.stan$fit_sum

  expect_equivalent(round(results['theta', '50%'], 2), -0.75, tolerance = 0.1)

})

## Fitting a meta-regression model
test_that("Results are correct for a meta-regression model.", {
  skip_on_cran()

  set.seed(11112)
  ## Load the dataset
  data('dat.Berkey1995', package = "MetaStan")
  ## Fitting a Binomial-Normal Hierarchical model using WIP priors
  data_converted = create_MetaStan_dat(dat = dat.Berkey1995,
                                       armVars = c(responders = "r", sampleSize = "n"))

  meta.reg.stan  <- meta_stan(data = data_converted,
                              likelihood = "binomial",
                              mu_prior = c(0, 10),
                              theta_prior = c(0, 100),
                              tau_prior = 0.5,
                              tau_prior_dist = "half-normal",
                              mreg = TRUE,
                              cov = dat.Berkey1995$Latitude)
  ### compare with results
  results = meta.reg.stan$fit_sum

  expect_equivalent(round(results['beta[1,1]', '50%'], 2), -0.03, tolerance = 0.1)

})


