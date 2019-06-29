
context("Checking model-based meta-analysis example: Eletriptans")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model with Smith parametrization
test_that("Results are correct for fitting model-based meta-analysis using Arm-based Emax model.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Eletriptan', package = "MetaStan")
  ## Fitting a Binomial-Normal Hierarchial model using WIP priors
  datMBMA = create_MBMA_dat(dat = dat.Eletriptan,
                  armVars = c(dose = "d", responders = "r",
                              sampleSize = "n"),
                  nArmsVar = "nd")

  MBMA.AB.Emax  <- MBMA_stan(data = datMBMA,
                             model = "AB_Emax",
                             Pred_doses = seq(0, 80, length.out = 11),
                             Emax_prior = c(0, 10),
                             tau_prior_dist = "half-normal",
                             tau_prior = 0.5)
  ### compare with results
  results = MBMA.AB.Emax$fit_sum

  expect_equivalent(round(results['theta_1', '50%'], 2), 2.26, tolerance = 0.1)


  MBMA.CB.Emax  <- MBMA_stan(data = datMBMA,
                             model = "CB_Emax",
                             Pred_doses = seq(0, 80, length.out = 11),
                             Emax_prior = c(0, 10),
                             tau_prior_dist = "half-normal",
                             tau_prior = 0.5)


})


