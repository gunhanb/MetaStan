/*
   Author: Burak Kuersad Guenhan
   A Model-based Meta-analysis model
   Model taken from Boucher and Bennets (2016)
   "Baseline random-effects model"
   Dose-response: Emax model (Random effects model)

*/

data {
  int<lower=1> Nstud;                            // num studies
  int<lower=1> Ndose;                            // num doses (including)
  int<lower=1> NPred;                           // num predicted doses
  real<lower=0> dose[Nstud,Ndose];               // matrix of doses
  real<lower=0> dose_na[Nstud,Ndose];            // locating NA values
  int<lower=0> r[Nstud,Ndose];                   // num events
  int<lower=1> n[Nstud,Ndose];                   // num patients
  real Pred_doses[NPred];                         // predicted doses

}

parameters {
  real mu;
  real theta_1;
  real<lower=0> theta_2;
  vector[Nstud] zeta;                          // individual treatment effects
  real<lower=0> tau;                        // heterogeneity stdev.

}

transformed parameters {
  real<lower=0, upper=1> p[Nstud,Ndose];

  for(i in 1:Nstud) {
    for(k in 1:Ndose) {
      if(dose_na[i,k] == 0) {  // if dose is not missing
        p[i,k]  = inv_logit(mu + ((theta_1 * dose[i,k])/(theta_2 + dose[i,k])) + zeta[i] * tau);
      }
      if(dose_na[i,k] == 1) {  // if dose is missing
        p[i,k] = 0.999;
      }
  }
}
}


model {

  // latent variable (random effects)
  zeta ~ normal(0, 1);
  // prior distributions
  mu ~ normal(0, 10); // baseline risks
  theta_1 ~ normal(0, 10);
  theta_2 ~ normal(0, 10);
  tau ~ normal(0, 0.5)T[0,];

  // likelihood
  for(i in 1:Nstud) {
    for(k in 1:Ndose) {
      if(dose_na[i,k] == 0) {  // if dose is not missing
        target += binomial_lpmf(r[i,k] | n[i,k], p[i,k]);
      }
    }
  }
}

// Posterior predictive checking
// and Model comparison
generated quantities {
  vector[NPred] zeta_pred; // Predicted delta
  real<lower=0, upper=1> Pred_probs[NPred]; // Predicted probabolities

  // placebo arm (dose = 0)
  Pred_probs[1] = inv_logit(mu);
  zeta_pred[1] = 0;

  for(i in 2:NPred) {
    zeta_pred[i] = normal_rng(0, 1);
    Pred_probs[i] = inv_logit(mu + (theta_1 * Pred_doses[i])/(theta_2 + Pred_doses[i]) + zeta_pred[i] * tau);
}





}





