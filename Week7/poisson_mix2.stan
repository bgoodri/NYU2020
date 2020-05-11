data {
  int<lower = 0> N;     // number of observations
  int<lower = 0> y[N];  // count outcomes
  vector<lower = 0>[2] loc;  // location of normal prior
  vector<lower = 0>[2] scal; // scal of normal prior
}
parameters {
  ordered[2] eta;                // log of means
  real<lower = 0, upper = 1> pi; // mixture probability
}
model {
  target += normal_lpdf(eta | loc, scal);
  target += normal_lccdf(eta[1] | loc[2], scal[2]); // truncation of PDF for eta[2]
  target += log_mix(pi, poisson_log_lpmf(y | eta[1]), 
                        poisson_log_lpmf(y | eta[2]));
}
generated quantities {
  vector[2] mu = exp(eta);
}
