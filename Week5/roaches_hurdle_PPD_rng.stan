functions {
  matrix roaches_hurdle_PPD_rng(int S, vector log_roach1, vector treatment,
                                vector senior, vector offset) {
    int N = rows(log_roach1);
    matrix[S, N] PPD;
    for (s in 1:S) {
      // hurdle parameters
      real gamma = normal_rng(0, 2);
      real lambda[3] = normal_rng([0,0,0], 1);
      
      // negative binomial parameters
      real alpha = normal_rng(0, 5);
      real beta[3] = normal_rng([0,0,0], 2);
      real phi = exponential_rng(1);
      
      for (n in 1:N) {
        real log_odds = gamma + lambda[1] * (log_roach1[n] == 0) 
                      + lambda[2] * treatment[n] + lambda[3] * senior[n];
        int hurdle = bernoulli_logit_rng(log_odds); // same as bernoulli_rng(inv_logit(log_odds));
        if (hurdle == 1) {
          real eta = alpha + offset[n] + beta[1] * log_roach1[n] 
                   + beta[2] * treatment[n] + beta[3] * senior[n];
          int y_ = neg_binomial_2_log_rng(eta, phi); // may be zero
          while(y_ == 0) y_ = neg_binomial_2_log_rng(eta, phi);
          PPD[s, n] = y_;
        } else PPD[s, n] = 0;
      }
    }
    return PPD;
  }
}
