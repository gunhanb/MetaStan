
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
  ## the MCMC draws differ by toolchain, so check the fitted summaries with a
  ## tolerance and pin them to fixed values for the snapshot. The tau interval
  ## is computed from the draws inside forest_plot, so it is left off the plot
  expect_equivalent(bnhm.wip.bnhm1.stan$fit_sum['theta', '50%'], -0.75, tolerance = 0.1)
  expect_equivalent(bnhm.wip.bnhm1.stan$fit_sum['tau[1]', '50%'], 0.58, tolerance = 0.1)

  pinned = bnhm.wip.bnhm1.stan
  pinned$fit_sum['theta', c('2.5%', '50%', '97.5%')] = c(-1.17, -0.75, -0.36)
  pinned$fit_sum['theta_pred[1]', c('2.5%', '50%', '97.5%')] = c(-2.04, -0.74, 0.55)

  ## forest_plot draws to the device rather than returning a plot, so vdiffr
  ## needs a function that draws it
  plotforest = function() {
    forest_plot(x = pinned,
                xlab = "log-OR",
                labels = dat.Berkey1995$publication,
                heterogeneity = FALSE)
  }
  ### compare with results
  vdiffr::expect_doppelganger("Forestplot: Crins data with publication names",
                              plotforest)

})


