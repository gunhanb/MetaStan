
context("Checking forest plot: Crins dataset")
suppressPackageStartupMessages(library(rstan))

## Fitting a binomial-normal hierachical model
test_that("Plot is correct for Crins dataset.", {
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
                                    theta_prior  = c(0, 10),
                                    tau_prior = 0.5,
                                    tau_prior_dist = "half-normal")
  plotforest = forest_plot(x = bnhm.wip.bnhm1.stan,
                           xlab = "log-OR",
                           labels = dat.Berkey1995$publication)
  ### compare with results
  vdiffr::expect_doppelganger("Forestplot: Crins data with publication names",
                              plotforest)

})


