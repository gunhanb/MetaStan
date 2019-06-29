/*
  Model based meta-analysis
  Contrast-based model
  Emax
  Random-effects
*/


data {
  int<lower=1> Nobs;                           // num observations
  int<lower=1> Nst;                           // num studies
  int<lower=1> NPred;                           // num predicted doses
  int<lower=0> r[Nobs];                       // num events
  int<lower=1> n[Nobs];                       // num patients
  real<lower=0> dose[Nobs];
  int<lower=1> st[Nobs];                    // Study for data point
  int<lower=1> ndose[Nst];                       // num doses in each trial
  real Pred_doses[NPred];                         // predicted doses
  int<lower=1> b_ndx[Nst];
  int<lower=1> t_ndx[Nobs-Nst];
}

parameters {
  real theta_1;                               //
  real<lower=0> theta_2;                      //
  vector[Nobs - Nst] delta_param_raw;
  real<lower=0> tau_1;
  real<lower=0> tau_2;
  real mu;                             // baseline risks
  vector[Nst] zeta;                          // individual treatment effects

}

transformed parameters{
  vector[Nobs - Nst] delta_param;
  cov_matrix[max(ndose)] Sigma;
  matrix[max(ndose), max(ndose)] Sigma_chol;
  vector[Nobs] md;
  vector[Nobs] delta;


  // Dose-response: Emax model
  for(i in 1:Nobs)
    md[i] = (theta_1 * dose[i])/(theta_2 + dose[i]);


  // Covariance matrix
  for(i in 1:max(ndose)){
    Sigma[i,i] = tau_1^2;
    for(j in i+1:max(ndose)) {
      Sigma[i,j] = (tau_1^2)/2;
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
          delta_param[pos_delta] = md[pos_md] + tau_1 * delta_param_raw[pos_delta];
        } else {           // multi-arm trials
            delta_param[pos_delta:(pos_delta+ndose[j]-2)] = md[pos_md:(pos_md+ndose[j]-2)] + Sigma_chol[1:(ndose[j]-1), 1:(ndose[j]-1)] * delta_param_raw[pos_delta:(pos_delta+ndose[j]-2)];
        }
        pos_delta = pos_delta + ndose[j] - 1;
        pos_md    = pos_md + ndose[j];

      }
    }

  for (i in 1:Nst)
    delta[b_ndx[i]] = 0;
  for (i in 1:(Nobs-Nst))
    delta[t_ndx[i]] = delta_param[i];


}

model {
  // latent variable (random effects)
  zeta ~ normal(0, 1);




    {
      int pos_delta;
      int pos_md;
      pos_delta = 1;
      pos_md = 2;
      for (j in 1:Nst) {
        if (ndose[j] == 2) { // two arm trials
          delta_param_raw[pos_delta] ~ normal(0,1);
        } else {           // multi-arm trials
            delta_param_raw[pos_delta:(pos_delta+ndose[j]-2)] ~ normal(0,1);
        }
        pos_delta = pos_delta + ndose[j] - 1;
        pos_md    = pos_md + ndose[j];

      }
    }


  // likelihood
  for (i in 1:Nobs)
    r[i] ~ binomial_logit(n[i], mu + zeta[st[i]] * tau_2 + delta[i]);

  theta_1 ~ normal(0, 10);
  theta_2 ~ normal(0, 10);
  tau_1 ~ normal(0, 0.5)T[0,];
  tau_2 ~ normal(0, 0.5)T[0,];

}


generated quantities {
  real delta_pred[NPred]; // Predicted delta
  real<lower=0, upper=1> Pred_probs[NPred]; // Predicted probabolities
  vector[NPred] zeta_pred; // Predicted delta

  for(i in 1:num_elements(Pred_doses)) {
    zeta_pred[i] = normal_rng(0, 1);
    delta_pred[i] = normal_rng(mu + zeta_pred[i] * tau_2 +  (theta_1 * Pred_doses[i])/(theta_2 + Pred_doses[i]), tau_1);
    Pred_probs[i] = inv_logit(delta_pred[i]);
  }

}

