
context("Checking meta-analysis example: Berkey (1995) Tuberclosis dataset")

## Fitting a pairwise random effects meta-analysis model
test_that("Results are correct for fitting binomial normal hierarchical model using WIP priors.", {
  skip_on_cran()
  ## Load the dataset
  data('dat.Berkey1995', package = "MetaStan")
  ## Fitting a Binomial-Normal Hierarchial model using WIP priors
  bnhm.wip.stan  <- meta_stan(ntrt = dat.Berkey1995$nt,
                              nctrl = dat.Berkey1995$nc,
                              rtrt = dat.Berkey1995$rt,
                              rctrl = dat.Berkey1995$rc,
                              mu_prior = c(0, 10),
                              delta_u = 250,
                              tau_prior = c(0, 0.5),
                              tau_prior_dist = "half-normal")
  ### compare with results
  results = summary(bnhm.wip.stan)[[1]]
  #print(results[14,1])
  expect_equivalent(round(results[14,1], 2), -0.76)

})

