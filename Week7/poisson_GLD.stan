#include quantile_functions.stan
data {
  int<lower = 0> N;     // number of observations
  int<lower = 0> y[N];  // count outcome

  real<lower = 0> m;    // prior median
  real<lower = 0> r;    // prior IQR
  real<lower = -1, upper = 1> asymmetry;
  real<lower =  0, upper = 1> steepness;
}
parameters {
  real<lower = 0, upper = 1> p;
}
transformed parameters {
  real mu = GLD_icdf(p, m, r, asymmetry, steepness);
} // implicit: p has a standard uniform prior
model {
  target += poisson_lpmf(y | mu);           // log likelihood
}
