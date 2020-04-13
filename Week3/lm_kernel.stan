functions { /* saved as lm_kernel.stan*/
  real lm_kernel(real alpha, real beta, real tau,
                 vector y, vector x) {
    int N = rows(x);
    vector[N] x_ = x - mean(x);
    vector[N] mu = alpha + beta * x_;
  
  
  
    // alpha and beta have improper priors ...
    // ... so they add nothing to the log-kernel
    //       vvv inv_sqrt(tau) = 1 / sqrt(tau)
    real sigma = inv_sqrt(tau);
    return -log(tau) // Jeffreys prior on tau
           + normal_lpdf(y | mu, sigma);
           // ^^^ log-likelihood of parameters
  }
}
