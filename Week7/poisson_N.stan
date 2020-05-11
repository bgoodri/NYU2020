data {
  int<lower = 0> N;     // number of observations
  int<lower = 0> y[N];  // count outcome
  
  real<lower = 0> loc;  // location of normal prior
  real<lower = 0> scal; // scale of normal prior
}
parameters {
  real eta; // log of mean of DGP
}

transformed parameters {
  /* could do this here
  real mu = exp(eta);
  */
}

model {
  target += normal_lpdf(eta | loc, scal); // log prior PDF
  target += poisson_log_lpmf(y | eta);    // log likelihood
}
generated quantities {
  real mu = exp(eta);
}
