data {
  int<lower = 0> N;     // number of observations
  int<lower = 0> y[N];  // count outcome
  
  real<lower = 0> a;    // shape of gamma prior
  real<lower = 0> b;    // rate of gamma prior
}
parameters {
  real<lower = 0> mu;   // mean of DGP
}
model {
  target += gamma_lpdf(mu | a, b); // log prior PDF
  target += poisson_lpmf(y | mu);  // log likelihood
}
