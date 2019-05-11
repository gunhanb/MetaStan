/* 
   Author: Burak Kuersad Guenhan
   A Model-based Meta-analysis model
   Model inspired by NMA literature
   "Arm-based random-effects model"
   Dose-response: Emax model (Random effects model)
   Dataset need is one-row-per-arm format
*/

data {
  int<lower=1> Nobs;                            // num observations
  int<lower=1> Nst;                             // num studies
  int<lower=1> NPred;                           // num predicted doses 
  int<lower=0> r[Nobs];                         // num events
  int<lower=1> n[Nobs];                         // num patients
  real<lower=0> dose[Nobs];     
  int<lower=1> ndose[Nst];                      // num doses in each study
  real Pred_doses[NPred];                       // predicted doses
}

parameters {
  real mu;                                    // fixed baseline risk parameter
  real theta_1;                               // Emax parameter
  real<lower=0> theta_2;                      // ED50 parameter
  real<lower=0> tau;                          // heterogeneity
  vector[Nobs] delta_raw;

}

transformed parameters{
  cov_matrix[max(ndose)] Sigma;
  matrix[max(ndose), max(ndose)] Sigma_chol;
  vector[Nobs] md;
  vector[Nobs] delta;
  
  // Dose-response: Emax model
  for (i in 1:Nobs) 
    md[i] = mu + (theta_1 * dose[i])/(theta_2 + dose[i]);


  // Homogeneous variance-covariance matrix
  for(i in 1:max(ndose)){
    Sigma[i,i] = tau^2;
    for(j in i+1:max(ndose)) {
      Sigma[i,j] = (tau^2)/2;
      Sigma[j,i] = Sigma[i,j];
    }
  }
  
  Sigma_chol = cholesky_decompose(Sigma);
  
      {
      int pos_delta;
      int pos_md;
      
      pos_delta = 1;
      pos_md = 1;
      for (j in 1:Nst) {
        
        delta[pos_delta:(pos_delta+ndose[j]-1)] = md[pos_md:(pos_md+ndose[j]-1)] + Sigma_chol[1:(ndose[j]), 1:(ndose[j])] * delta_raw[pos_delta:(pos_delta+ndose[j]-1)];
        
        pos_delta = pos_delta + ndose[j];
        pos_md    = pos_md + ndose[j];


      }
    }

}

model {


    {
      int pos_delta;
      int pos_md;
      
      pos_delta = 1;
      pos_md = 1;
      for (j in 1:Nst) {
        
        delta_raw[pos_delta:(pos_delta+ndose[j]-1)] ~ normal(0,1);
        
        pos_delta = pos_delta + ndose[j];
        pos_md    = pos_md + ndose[j];


      }
    }
    

  // likelihood
  r ~ binomial_logit(n, delta);

  // prior distributions
  mu ~ normal(0, 10);
  theta_1 ~ normal(0, 10);
  theta_2 ~ normal(0, 10);
  tau ~ normal(0, 0.5)T[0,];
}


generated quantities {
  vector[NPred] delta_pred_mean;                                     // Predicted delta mean
  matrix[NPred, NPred] delta_pred_Sigma;
  vector[NPred] delta_pred;                                          // Predicted delta mean
  real<lower=0, upper=1> Pred_probs[NPred];                          // Predicted probabolities


  for(i in 1:NPred) {
    delta_pred_mean[i] = (theta_1 * Pred_doses[i])/(theta_2 + Pred_doses[i]);
  }
    // Covariance matrix
  for(i in 1:NPred){
    delta_pred_Sigma[i,i] = tau^2;
    for(j in i+1:NPred) {
      delta_pred_Sigma[i,j] = (tau^2)/2;
      delta_pred_Sigma[j,i] = delta_pred_Sigma[i,j];
    }
  }

  delta_pred = multi_normal_rng(delta_pred_mean, delta_pred_Sigma);


  for(i in 1:NPred) 
    Pred_probs[i] = inv_logit(mu + delta_pred[i]);
  
}

