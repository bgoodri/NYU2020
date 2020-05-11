data {
  int<lower = 1> N; // number of studies
  vector<lower = 0>[N] se; // std. errors
}
transformed data { // for SBC
  vector[N] se2 = square(se);
  real mu_ = normal_rng(0, 1); // truth
  real tau_ = exponential_rng(1);
  vector[N] y; // estimates of mu
  for (n in 1:N) { // prior predictions
    real epsilon = normal_rng(0, 1);
    real delta = mu_ + tau_ * epsilon;
    y[n] = normal_rng(delta, se[n]);
  }
}
parameters {
  real mu;
  real<lower = 0> tau;
}
model {
  target += normal_lpdf(mu | 0, 1);
  target += exponential_lpdf(tau | 1);
  target += normal_lpdf(y | mu, sqrt(square(tau) + se2));
}
generated quantities {
  vector[N] y_ = y;  // copy of prior predictions
  vector[N] log_lik; // for loo()
  vector[2] pars_;   // copy of prior realizations
  int ranks_[2] = {mu > mu_, tau > tau_}; // for plot
  pars_[1] = mu_;    
  pars_[2] = tau_;
  for (n in 1:N) {
    real s = sqrt(square(tau) + se2[n]);
    log_lik[n] = normal_lpdf(y[n] | mu, s);
  }
}
