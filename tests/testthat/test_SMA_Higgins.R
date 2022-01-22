
context("Checking meta-analysis example: BCG dataset")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model
test_that("Results are correct for fitting binomial normal hierarchical model using Higgins parametrization", {
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
                                    param = "Higgins",
                                    tau_prior_dist = "half-normal")
  ### compare with results
  results = bnhm.wip.bnhm1.stan$fit_sum

  expect_equivalent(round(results['theta', '50%'], 2), -0.75, tolerance = 0.1)

})

