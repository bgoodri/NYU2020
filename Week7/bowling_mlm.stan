#include bowling_kernel.stan
data {
  int<lower = 0> J;                           // number of bowlers
  int<lower = 0, upper = 10> x1_x2[J, 10, 2]; // results of each bowler's frames
  vector<lower = 0>[11] a;                    // shapes for Dirichlet prior on mu
  real<lower = 0> s;                          // scale factor on top of theta
}
parameters {
  simplex[11] mu;        // overall probability of knocking down 0:10 pins
  real<lower = 0> theta; // concentration parameter across bowlers
  simplex[11] pi[J];     // bowler's probability of knocking down 0:10 pins
}
model { // target becomes the log-numerator of Bayes Rule
  vector[11] mu_theta = mu * theta * s;                  // not saved in results
  for (j in 1:J) // bowling_kernel() is defined in the functions block
    target += bowling_kernel(pi[j], mu_theta, x1_x2[j]); // note indexing
  target += dirichlet_lpdf(mu | a);                      // prior on mu
  target += exponential_lpdf(theta | 1);                 // prior on theta
}
