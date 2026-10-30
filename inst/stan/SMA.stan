/*
   Author: Burak Kuersad Guenhan
   A standard Meta-analysis model
   Dataset need is one-row-per-arm format
*/


data {
  // num of observations
  int<lower=1> Nobs;

  // treatment arm indicator (0=control, 1=experimental)
  vector[Nobs] t;

  // Study number for each observation
  int<lower=1> st[Nobs];

  // link function
  // 1 = continuous (MD or SMD)
  // 2 = binary logit (log OR)
  // 3 = count log (rate ratio)
  // 4 = binary log (log RR)
  // 5 = binary identity (RD)
  int<lower=1,upper=5> link;

  // normal data, link=identity=1
  vector[Nobs] y;
  vector[Nobs] y_se;

  // binomial data, link=logit=2
  int<lower=0>  r[Nobs];
  int<lower=1>  n[Nobs];


  // count data, link=log=3
  int<lower=0> count[Nobs];
  vector[Nobs]    exposure;

  // Priors
  vector[2] mu_prior;                         // Prior parameters for mu
  vector[2] theta_prior;                      // Prior parameters for theta
  real tau_prior;                             // Prior parameter for tau
  int tau_prior_dist;                        // Indicator for distribution of prior for tau
  vector[2] beta_prior;                         // Prior parameters for beta

  // Fixed-effects or Random-effects model: (0: fe, 1: re)
  int<lower=0,upper=1> re;

  // non-centered parametrization: (0: NO, 1: YES)
  int<lower=0,upper=1> ncp;

  // met-regression: (0: NO, 1: YES)
  int<lower=0,upper=1> mreg;

  // number of covariates
  int<lower=0> ncov;

  // trial-level covariate information
  matrix[ncov,max(st)] cov;

}

parameters {
  vector[Nobs] mu;                             // baseline risks (log odds)
  real theta;                               // relative treatment effect (log odds ratio)
  vector[Nobs] u[re];                          // individual treatment effects
  real<lower=0> tau[re];                        // heterogeneity stdev.
  vector[ncov] beta[mreg];                      // beta coeffients in meta-regression
}

transformed parameters {
  vector[Nobs] d;
  vector[max(st)] temp;
  vector[max(st)] beta_cov;

  vector[Nobs] p; // event probabilities (binomial links only)

  if(mreg) {
    for(i in 1:Nobs)
      temp[st[i]] = 0.5 * beta[1,1] * cov[1,st[i]];

    if(ncov == 1) {
      for(i in 1:Nobs)
      beta_cov[st[i]] = temp[st[i]];

    } else {

      for(j in 1:ncov) {
        for(i in 1:Nobs) {
          beta_cov[st[i]] = temp[st[i]] + 0.5 * beta[1,j] * cov[j,st[i]];
          temp[st[i]] = beta_cov[st[i]];
        }
      }
    }
  } else {
    for(i in 1:Nobs) {
     beta_cov[st[i]] = 0;
     temp[st[i]] = 0;
    }
    };



  if(re) {

  if (ncp) {
    for(i in 1:Nobs) {
      if (t[i] == 0) {
        d[i] =  mu[st[i]]  - 0.5 * theta - 0.5 * u[1,st[i]] * tau[1] - beta_cov[st[i]];
       } else {
        d[i] = mu[st[i]] + 0.5 * theta + 0.5 * u[1,st[i]] * tau[1] + beta_cov[st[i]];
       }
    }

   } else {
      for(i in 1:Nobs) {
        if (t[i] == 0) {
          d[i] = mu[st[i]] - 0.5 * u[1,st[i]] - beta_cov[st[i]];
        }   else {
          d[i] = mu[st[i]] + 0.5 * u[1,st[i]] + beta_cov[st[i]];
        }
      }
  }
} else {
        for(i in 1:Nobs) {
          if (t[i] == 0) {
            d[i] = mu[st[i]] - 0.5 * theta - beta_cov[st[i]];
          }   else {
            d[i] = mu[st[i]] + 0.5 * theta + beta_cov[st[i]];
        }
      }
}

  // event probabilities, computed once d is defined
  for (i in 1:Nobs) {
    if (link == 2) {
      // logit link (OR)
      p[i] = inv_logit(d[i]);
    } else if (link == 4) {
      // log link (RR)
      p[i] = exp(d[i]);
    } else if (link == 5) {
      // identity link (RD)
      p[i] = d[i];
    } else {
      // not used for normal or poisson likelihoods
      p[i] = 0;
    }
  }

}

model {

    if(re) {
      if (ncp) {
      // latent variable (random effects)
      u[1] ~ normal(0, 1);
      } else {
      // random effect distribution
      u[1] ~ normal(d, tau[1]);
      }
    }
  // prior distributions
  mu ~ normal(mu_prior[1], mu_prior[2]);
  theta ~ normal(theta_prior[1], theta_prior[2]);

  if(re) {
    // prior for tau
    if(tau_prior_dist == 1)  tau[1] ~ normal(0, tau_prior)T[0,];
    if(tau_prior_dist == 2)  tau[1] ~ uniform(0, tau_prior);
    if(tau_prior_dist == 3)  tau[1] ~ cauchy(0, tau_prior)T[0,];
  }

  if(mreg)  beta[1] ~ normal(beta_prior[1], beta_prior[2]);

  // likelihood
  if (link == 1) {
    // MD or SMD
    y ~ normal(d, y_se);
  }

  if (link == 2) {
    // logit link (OR)
    r ~ binomial_logit(n, d);
  }

  if (link == 3) {
    // log rate ratio
    count ~ poisson_log(exposure + d);
  }

  if (link == 4) {
    // log link (RR)
    for (i in 1:Nobs)
      r[i] ~ binomial(n[i], p[i]);
  }

  if (link == 5) {
    // identity link (RD)
    for (i in 1:Nobs)
      r[i] ~ binomial(n[i], p[i]);
  }

}

generated quantities {
  vector[Nobs] log_lik;                // pointwise log-likelihood contribution
  real theta_pred[re];                 // predicted log-odds ratio for the new study


  for (s in 1:Nobs) {
    if (link == 1)
    log_lik[s] = normal_lpdf(y[s] | d[s], y_se[s]);


    if (link == 2)
    log_lik[s] = binomial_logit_lpmf(r[s] | n[s], d[s]);


    if (link == 3)
    log_lik[s] = poisson_log_lpmf(count[s] | exposure[s] + d[s]);


    if (link == 4)
    log_lik[s] = binomial_lpmf(r[s] | n[s], p[s]);


    if (link == 5)
    log_lik[s] = binomial_lpmf(r[s] | n[s], p[s]);
  }

  // Prediction for new study
   if(re) { theta_pred[1] = normal_rng(theta, tau[1]);}

}
