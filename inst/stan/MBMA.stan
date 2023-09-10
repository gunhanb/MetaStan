/*
   Author: Burak Kuersad Guenhan
   Model-based Meta-analysis
   Dataset needed is one-row-per-arm format
*/

data {
  // num of arm
  int<lower=1> Nobs;

  // num of studies
  int<lower=1> Nst;

  // Study number for each arm
  array[Nobs] int<lower=1> st;

  // Dose amount for each arm
  array[Nobs] real<lower=0> dose;


  // num doses in each trial
  array[Nst] int<lower=1> ndose;

  // link function (1=normal, 2=binary, 3=poisson)
  int<lower=1,upper=3> link;

  // dose-response function (1=Linear, 2=Log-linear, 3=Emax, 4=sigmoidal emax)
  int<lower=1,upper=4> dose_response;

  // normal data, link=identity=1
  vector[Nobs] y;
  vector[Nobs] y_se;

  // binomial data, link=logit=2
  array[Nobs] int<lower=0>  r;
  array[Nobs] int<lower=1>  n;

  // count data, link=log=3
  array[Nobs] int<lower=0> count;
  vector[Nobs]    exposure;

  // Priors
  vector[2] mu_prior;
  vector[2] alpha_prior;
  vector[2] ED50_prior;
  real tau_prior;
  int tau_prior_dist;                        // Indicator for distribution of prior for tau
  int ED50_prior_dist;                        // Indicator for distribution of prior for ED50
  vector[2] gamma_prior;


  // Number of prediction points in the dose-response curve
  int<lower=1> Npred;

  // Prediction points in the dose-response curve
  array[Npred] real Pred_doses;

  //  Indicator for placebo arm
  array[Nst] int<lower=1> b_ndx;

  // Indicator for non-placebo arm
  array[Nobs-Nst] int<lower=1> t_ndx;

  // Fixed-effects or Random-effects model: (0: fe, 1: re)
  int<lower=0,upper=1> re;

  // non-centered parametrization: (0: NO, 1: YES)
  int<lower=0,upper=1> ncp;

}
transformed data{
  // Is it an Emax model: (0: Linear or Log-linear, 1: Emax or Sigmoidal emax)
  int<lower=0,upper=1> emax;
  // Is it a sigmoidal Emax model: (0: Linear or Log-linear or 0: Emax, 1: Sigmoidal emax)
  int<lower=0,upper=1> hill;

  real<lower=0> maxdose;

  if(dose_response == 1) {
    emax = 0;
    hill = 0;
  }
  if(dose_response == 2) {
    emax = 0;
    hill = 0;
  }
  if(dose_response == 3) {
    emax = 1;
    hill = 0;
  }
  if(dose_response == 4) {
    emax = 1;
    hill = 1;
  }

  maxdose = max(dose);

}


parameters {
  vector[Nst] mu;                             // baseline risks
  real alpha;
  array[emax] real<lower=0, upper=1.5> ED50_raw;
  array[hill] real<lower=0.5, upper=10> gamma;              // Hill parameter
  vector[Nobs - Nst] u;
  array[re] real<lower=0> tau;
}

transformed parameters{
  vector[Nobs - Nst] delta_param;
  cov_matrix[max(ndose)] Sigma;
  matrix[max(ndose), max(ndose)] Sigma_chol;
  vector[Nobs] md;
  vector[Nobs] delta;
  vector[Nobs] theta;
  array[emax] real<lower=0> ED50;

  // used for the prior distribution of ED50 (ED50) parameter
  if(emax == 1) {ED50[1] = ED50_raw[1] * maxdose;}

  // Dose-response: Emax model
  for(i in 1:Nobs) {
    if(dose_response == 1) {md[i] = (alpha * dose[i]);}                        // linear
    if(dose_response == 2) {md[i] = (alpha * log(dose[i] + 1));}               // log-linear
    if(dose_response == 3) {md[i] = (alpha * dose[i])/(ED50[1] + dose[i]);} // Emax
    if(dose_response == 4) {md[i] = (alpha * dose[i]^gamma[1])/
    (ED50[1]^gamma[1] + dose[i]^gamma[1]);}                                           // Sigmoidal Emax
  }
  // Covariance matrix
  for(i in 1:max(ndose)){
    Sigma[i,i] = tau[1]^2;
    for(j in i+1:max(ndose)) {
      Sigma[i,j] = (tau[1]^2)/2;
      Sigma[j,i] = Sigma[i,j];
    }
  }

  Sigma_chol = cholesky_decompose(Sigma);

  {
      int pos_delta;
      int pos_md;
      pos_delta = 1;
      pos_md = 2;
      for (j in 1:Nst) {
        if (ndose[j] == 2) { // two arm trials
          delta_param[pos_delta] = md[pos_md] + tau[1] * u[pos_delta];
        } else {           // multi-arm trials
            delta_param[pos_delta:(pos_delta+ndose[j]-2)] = md[pos_md:(pos_md+ndose[j]-2)] +
            Sigma_chol[1:(ndose[j]-1), 1:(ndose[j]-1)] * u[pos_delta:(pos_delta+ndose[j]-2)];
        }
        pos_delta = pos_delta + ndose[j] - 1;
        pos_md    = pos_md + ndose[j];

      }

  }
  for (i in 1:Nst)
    delta[b_ndx[i]] = 0;
  for (i in 1:(Nobs-Nst))
    delta[t_ndx[i]] = delta_param[i];

  if(re) {
    if (ncp) {
      for(i in 1:Nobs)
        theta[i] = mu[st[i]] + delta[i];
    }
  }

}

model {

    {
      int pos_delta;
      pos_delta = 1;
      for (j in 1:Nst) {
        if (ndose[j] == 2) { // two arm trials
          u[pos_delta] ~ normal(0,1);
        } else {           // multi-arm trials
            u[pos_delta:(pos_delta+ndose[j]-2)] ~ normal(0,1);
        }
        pos_delta = pos_delta + ndose[j] - 1;

      }
    }


  // prior distributions
  mu    ~ normal(mu_prior[1], mu_prior[2]);
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);

  if(emax == 1) {
      // approximation to the functional uniform prior
      if(ED50_prior_dist == 1) ED50_raw[1]  ~ lognormal(ED50_prior[1], ED50_prior[2]);
      // a half normal distribution
      if(ED50_prior_dist == 2) ED50_raw[1]  ~ normal(ED50_prior[1], ED50_prior[2])T[0,];
       // a uniform distribution
      if(ED50_prior_dist == 3) ED50_raw[1]  ~ uniform(0, ED50_prior[1]);;

  }

  if(hill == 1) {gamma[1]  ~ normal(gamma_prior[1], gamma_prior[2])T[0,];}

  // prior for tau
  if(tau_prior_dist == 1)  tau[1] ~ normal(0, tau_prior)T[0,];
  if(tau_prior_dist == 2)  tau[1] ~ uniform(0, tau_prior);
  if(tau_prior_dist == 3)  tau[1] ~ cauchy(0, tau_prior)T[0,];

  // likelihood
  if(link == 1) y ~ normal(theta, y_se);
  if(link == 2) r ~ binomial_logit(n, theta);
  if(link == 3) count ~ poisson_log(exposure + theta);
}


generated quantities {
  real mean_mu;
  array[Npred] real md_pred;
  array[Npred] real delta_pred;                                      // Predicted delta
  array[Npred] real<lower=0, upper=1> Pred_probs;                    // Predicted probabolities
  vector[Nobs] log_lik;                                        // pointwise log-likelihood contribution


  mean_mu = mean(mu);

  Pred_probs[1] = inv_logit(mean(mu));
  md_pred[1] = 0;
  delta_pred[1] = 0;
  for(i in 2:num_elements(Pred_doses)) {
    if(dose_response == 1) {md_pred[i] = (alpha * Pred_doses[i]);}                        // linear
    if(dose_response == 2) {md_pred[i] = (alpha * log(Pred_doses[i] + 1));}               // log-linear
    if(dose_response == 3) {md_pred[i] = (alpha * Pred_doses[i])/(ED50[1] + Pred_doses[i]);} // Emax
    if(dose_response == 4) {md_pred[i] = (alpha * Pred_doses[i]^gamma[1])/
    (ED50[1]^gamma[1] + Pred_doses[i]^gamma[1]);}                                           // Sigmoidal Emax

    delta_pred[i] = normal_rng(mean_mu + md_pred[i], tau[1]);
    Pred_probs[i] = inv_logit(delta_pred[i]);
  }

  for (s in 1:Nobs) {
    if(link == 1)  log_lik[s] = normal_lpdf(y | theta[s], y_se[s]);
    if(link == 2)  log_lik[s] = binomial_logit_lpmf(r[s] | n[s], theta[s]);
    if(link == 3)  log_lik[s] = poisson_log_lpmf(count[s] | exposure[s] + theta[s]);
  }


}


