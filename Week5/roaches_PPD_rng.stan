functions {
  matrix roaches_PPD_rng(int S, vector log_roach1, vector treatment,
                         vector senior, vector offset) {
    int N = rows(log_roach1);
    vector[N] log_roach1_ = log_roach1 - mean(log_roach1);
    vector[N] treatment_ = treatment - mean(treatment);
    vector[N] senior_ = senior - mean(senior);
    vector[N] offset_ = offset - mean(offset);
    matrix[S, N] PPD;
    for (s in 1:S) {
      real alpha = normal_rng(4, 5);
      real beta[3] = normal_rng([0,0,0], 2);
      real phi[N] = rep_array(exponential_rng(1), N);
      real lambda[N] = gamma_rng(phi, phi);
      vector[N] log_lambda = to_vector(log(lambda));
      vector[N] eta = alpha + offset + beta[1] * log_roach1 + 
        beta[2] * treatment + beta[3] * senior + log_lambda;
      vector[N] mu = exp(eta);
      PPD[s, ] = to_row_vector(poisson_rng(mu));
    }
    return PPD;
  }
}
