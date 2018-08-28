/* A Beta-binomial model for
   pairwise meta-analysis
   Kuss (2014)
*/


data {
  int<lower=1> N;                           // num studies
  int<lower=0> rctrl[N];                    // num events, control
  int<lower=1> nctrl[N];                    // num patients, control
  int<lower=0> rtrt[N];                     // num events, treatment
  int<lower=1> ntrt[N];                     // num patients, treatment
}

parameters {
  real<lower=0,upper=1> pctrl[N];
  real<lower=0,upper=1> ptrt[N];
  real<lower=0,upper=1> rho;
  real<lower=0,upper=1> muctrl;
  real<lower=0,upper=1> mutrt;
}

transformed parameters {
  real<lower=0> alphactrl;
  real<lower=0> betactrl;
  real<lower=0> alphatrt;
  real<lower=0> betatrt;

  alphactrl = muctrl * (1 - rho) / rho;
  alphatrt  = mutrt * (1 - rho) / rho;
  betactrl = (1 - muctrl) * (1 - rho) / rho;
  betatrt = (1 - mutrt) * (1 - rho) / rho;


}

model {
  pctrl ~ beta(alphactrl, betactrl);
  ptrt ~ beta(alphatrt, betatrt);

  // likelihood
  rctrl ~ binomial(nctrl, pctrl);
  rtrt ~ binomial(ntrt, ptrt);
}

generated quantities {
  real b_0;
  real theta;                             //relative treatment effect (log odds ratio)

  b_0 = logit(muctrl);
  theta = logit(mutrt) - logit(muctrl);
}

