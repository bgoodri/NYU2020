functions { /* saved as binormal_rng.stan in R's working directory */
  matrix binormal_rng(int S, real mu_X, real mu_Y, real sigma_X, real sigma_Y, real rho) {
    matrix[S, 2] draws; real beta = rho * sigma_Y / sigma_X; // calculate constants once ...
    real sigma = sigma_Y * sqrt(1 - square(rho));            // ... before the loop begins
    for (s in 1:S) {
      real x_ = normal_rng(mu_X, sigma_X);
      real y_ = normal_rng(mu_Y + beta * (x_ - mu_X), sigma);
      draws[s, ] = [x_, y_]; // a row_vector
    }
    return draws;
  }
}
