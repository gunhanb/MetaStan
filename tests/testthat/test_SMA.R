
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



## Relative risk and risk difference transformations
test_that("Relative risk and risk difference fits agree with metafor.", {
  skip_on_cran()

  set.seed(4321)
  data('dat.Berkey1995', package = "MetaStan")
  dat_long = create_MetaStan_dat(dat = dat.Berkey1995,
                                 armVars = c(responders = "r",
                                             sampleSize = "n"))

  ref_rr = metafor::rma(measure = "RR",
                        ai = dat.Berkey1995$r2, n1i = dat.Berkey1995$n2,
                        ci = dat.Berkey1995$r1, n2i = dat.Berkey1995$n1)
  rr.stan = meta_stan(data = dat_long,
                      likelihood = "binomial",
                      transformation = "relative_risk",
                      mu_prior = c(0, 10),
                      theta_prior = c(0, 10),
                      tau_prior = 0.5,
                      refresh = 0)
  expect_equal(rr.stan$transformation, "relative_risk")
  expect_equal(rr.stan$stanDat$link, 4)
  expect_equivalent(rr.stan$fit_sum['theta', '50%'], as.numeric(ref_rr$b), tolerance = 0.1)

  ref_rd = metafor::rma(measure = "RD",
                        ai = dat.Berkey1995$r2, n1i = dat.Berkey1995$n2,
                        ci = dat.Berkey1995$r1, n2i = dat.Berkey1995$n1)
  rd.stan = meta_stan(data = dat_long,
                      likelihood = "binomial",
                      transformation = "risk_difference",
                      mu_prior = c(0, 10),
                      theta_prior = c(0, 1),
                      tau_prior = 0.1,
                      refresh = 0)
  expect_equal(rd.stan$stanDat$link, 5)
  expect_equivalent(rd.stan$fit_sum['theta', '50%'], as.numeric(ref_rd$b), tolerance = 0.02)
  expect_false(any(is.na(rd.stan$fit_sum[grep("^p\\[", rownames(rd.stan$fit_sum)), "mean"])))
})

test_that("Invalid transformation arguments are rejected.", {
  data('dat.Berkey1995', package = "MetaStan")
  dat_long = create_MetaStan_dat(dat = dat.Berkey1995,
                                 armVars = c(responders = "r",
                                             sampleSize = "n"))

  expect_error(meta_stan(data = dat_long, likelihood = "binomial",
                         transformation = "log", theta_prior = c(0, 1)),
               "transformation")
  expect_error(meta_stan(data = dat_long, likelihood = "normal",
                         transformation = "relative_risk", theta_prior = c(0, 1)),
               "only valid for")
  expect_error(meta_stan(data = dat_long, likelihood = "binomial",
                         transformation = "relative_risk", param = "Higgins",
                         theta_prior = c(0, 1)),
               "Smith")
  expect_error(meta_stan(data = dat_long, likelihood = "binomial",
                         transformation = "risk_difference", delta = 250),
               "ratio scale")
})
