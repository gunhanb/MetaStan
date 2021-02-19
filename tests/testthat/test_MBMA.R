
context("Checking model-based meta-analysis example")
suppressPackageStartupMessages(library(rstan))

## Fitting a model-based meta-analysis model
test_that("Results are correct for fitting model-based meta-analysis Eleptriptan dataset.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Eletriptan', package = "MetaStan")
  ## Fitting a MBMA model
  datMBMA = create_MBMA_dat(dat = dat.Eletriptan,
                            armVars = c(dose = "d", r = "r",
                                        n = "n"),
                            nArmsVar = "nd")

  MBMA.Emax  <- MBMA_stan(data = datMBMA,
                          family = "binomial",
                          dose_response = "emax",
                          Pred_doses = seq(0, 80, length.out = 11),
                          mu_prior = c(0, 100),
                          Emax_prior = c(0, 100),
                          tau_prior_dist = "half-normal",
                          tau_prior = 0.5)
  ### compare with results
  results = MBMA.Emax$fit_sum

  expect_equivalent(round(results['alpha', '50%'], 2), 2.44, tolerance = 0.1)



})

## Fitting a model-based meta-analysis model
test_that("Results are correct for fitting model-based meta-analysis Paresthesia dataset.", {
  skip_on_cran()

  set.seed(23344)
  ## Load the dataset
  data('dat.Paresthesia', package = "MetaStan")
  ## Fitting a MBMA model
  datMBMA = create_MBMA_dat(dat = dat.Paresthesia,
                            armVars = c(dose = "d", r = "r",
                                        n = "n"),
                            nArmsVar = "nd")

  MBMA.Emax  <- MBMA_stan(data = datMBMA,
                          family = "binomial",
                          dose_response = "emax",
                          mu_prior = c(0, 100),
                          Emax_prior = c(0, 100),
                          tau_prior_dist = "half-normal",
                          tau_prior = 0.5)
  ### compare with results
  results = MBMA.Emax$fit_sum

  expect_equivalent(round(results['alpha', '50%'], 2), 2.91, tolerance = 0.1)



})




