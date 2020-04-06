#include quantile_functions.stan
data { /* these are known and passed as a named list from R */
  int<lower = 0> n;                          // number of tests
  int<lower = 0, upper = n> y;               // number of positives
  real m;
  real<lower = 0> IQR;
  real<lower = -1, upper = 1> asymmetry;
  real<lower =  0, upper = 1> steepness;
}
parameters { /* these are unknowns whose posterior distribution is sought */
  real<lower = 0, upper = 1> p;             // CDF of positive test probability
}
transformed parameters { /* deterministic unknowns that get stored in RAM */
  real theta = GLD_icdf(p, m, IQR, asymmetry, steepness); // positive probability
} // this function ^^^ is defined in the quantile_functions.stan file
model { /* log-kernel of Bayes' Rule that essentially returns "target" */
  target += binomial_lpmf(y | n, theta); // log-likelihood (as a function of theta)
} // implicit: p ~ uniform(0, 1) <=> theta ~ GLD(m, IQR, asymmetry, steepness)
generated quantities { /* other unknowns that get stored but are not needed */
  int y_ = binomial_rng(n, theta);       // posterior predictive realization
}
