/* A Binomial-Normal hierarchical model for
   pairwise meta-analysis
   Model 4 from Jackson et al (2018)
*/


data {
  int<lower=1> N;                           // num studies
  int<lower=0> rctrl[N];                    // num events, control
  int<lower=1> nctrl[N];                    // num patients, control
  int<lower=0> rtrt[N];                     // num events, treatment
  int<lower=1> ntrt[N];                     // num patients, treatment
  vector[2] mu_prior;                         // Prior parameters for mu
  vector[2] theta_prior;                          // Prior parameters for d
  real tau_prior;                             // Prior parameters for tau
  int tau_prior_dist;                        // Indicator for distribution of prior for tau
}

parameters {
  vector[N] mu;                             // baseline risks (log odds)
  real theta;                               // relative treatment effect (log odds ratio)
  vector[N] zeta;                          // individual treatment effects
  real<lower=0> tau;                        // heterogeneity stdev.
}

transformed parameters {
  real pctrl[N];
  real ptrt[N];

  for(i in 1:N) {
    pctrl[i] = inv_logit(mu[i] - theta * 0.5 - zeta[i] * tau * 0.5);
    ptrt[i]  = inv_logit(mu[i] + theta * 0.5 + zeta[i] * tau * 0.5);
  }
}

model {
  // latent variable (random effects)
  zeta ~ normal(0, 1);
  // prior distributions
  mu ~ normal(mu_prior[1], mu_prior[2]);
  theta ~ normal(theta_prior[1], theta_prior[2]);
  if(tau_prior_dist == 1)  tau ~ normal(0, tau_prior)T[0,];
  if(tau_prior_dist == 2)  tau ~ uniform(0, tau_prior);
  if(tau_prior_dist == 3)  tau ~ cauchy(0, tau_prior)T[0,];
  // likelihood
  rctrl~ binomial(nctrl, pctrl);                // cntrl
  rtrt ~ binomial(ntrt, ptrt);                  // trt
}
