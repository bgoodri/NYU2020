functions { /* saved as OLS_rng.stan*/
  matrix OLS_rng(int S, real alpha, real beta, 
                 real sigma, vector x) {
    matrix[S, 3] out; int N = rows(x); 
    vector[N] x_ = x - mean(x);
    vector[N] mu = alpha + beta * x_; 
    real SSX = sum(square(x_)); int Nm2 = N - 2; 
    for (s in 1:S) {
      vector[N] y = to_vector(normal_rng(mu, sigma));
      real a = mean(y);
      real b = sum(y .* x_) / SSX;
      vector[N] residuals = y - (a + b * x_);
      real s2_hat = sum(square(residuals)) / Nm2;
      out[s, ] = [a, b, s2_hat];
    }
    return out;
  }
}
