/* A fixed effect model for
   pairwise meta-analysis
   using Binomial likelihood
*/


data {
  int<lower=1> N;                           // num studies
  int<lower=0> rctrl[N];                    // num events, control
  int<lower=1> nctrl[N];                    // num patients, control
  int<lower=0> rtrt[N];                     // num events, treatment
  int<lower=1> ntrt[N];                     // num patients, treatment
  vector[2] mu_prior;                         // Prior parameters for mu
  vector[2] theta_prior;                          // Prior parameters for d
}

parameters {
  vector[N] mu;                             // baseline risks (log odds)
  real theta;                               // relative treatment effect (log odds ratio)
}

transformed parameters {
  real pctrl[N];
  real ptrt[N];

  for(i in 1:N) {
    pctrl[i] = inv_logit(mu[i] - theta * 0.5);
    ptrt[i] =  inv_logit(mu[i] + theta * 0.5);
  }
}

model {
  // prior distributions
  mu ~ normal(mu_prior[1], mu_prior[2]);
  theta ~ normal(theta_prior[1], theta_prior[2]);
  // likelihood
  rctrl~ binomial(nctrl, pctrl);                // cntrl
  rtrt ~ binomial(ntrt, ptrt);                  // trt
}
