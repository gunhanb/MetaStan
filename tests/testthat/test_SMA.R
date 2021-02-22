
context("Checking meta-analysis example: Crins dataset")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model
test_that("Results are correct for fitting binomial normal hierarchical model using WIP priors.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Berkey1995', package = "MetaStan")
  ## Fitting a Binomial-Normal Hierarchical model using WIP priors
  dat_long <- convert_data_arm(dat.Berkey1995$nt, dat.Berkey1995$nc,
                               dat.Berkey1995$rt, dat.Berkey1995$rc,
                               dat.Berkey1995$publication)

bnhm.wip.bnhm1.stan  <- meta_stan(data = dat_long,
                                    family = "binomial",
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
  data_converted <- convert_data_arm(dat.Berkey1995$nt, dat.Berkey1995$nc,
                               dat.Berkey1995$rt, dat.Berkey1995$rc,
                               dat.Berkey1995$publication)

  meta.reg.stan  <- meta_stan(data = data_converted,
                              family = "binomial",
                              mu_prior = c(0, 10),
                              theta_prior = c(0, 100),
                              tau_prior = 0.5,
                              tau_prior_dist = "half-normal",
                              mreg = TRUE,
                              cov = matrix(dat.Berkey1995$Latitude, nrow = 1))
  ### compare with results
  results = meta.reg.stan$fit_sum

  expect_equivalent(round(results['beta[1,1]', '50%'], 2), -0.03, tolerance = 0.1)

})


