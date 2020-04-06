functions { /* saved as Metropolis_rng.stan in R's working directory */

  real binormal_lpdf(row_vector xy, real mu_X, real mu_Y, real sigma_X, real sigma_Y, real rho) {
    real beta = rho * sigma_Y / sigma_X; real sigma = sigma_Y * sqrt(1 - square(rho));
    if (is_inf(xy[1]) || is_inf(xy[2])) return negative_infinity();
    return normal_lpdf(xy[1] | mu_X, sigma_X) + // normal_lpdf is the logarithm of the normal PDF
           normal_lpdf(xy[2] | mu_Y + beta * (xy[1] - mu_X), sigma);
  }
  
  matrix Metropolis_rng(int S, real half_width, real mu_X, real mu_Y, real sigma_X, real sigma_Y, real rho) {
    matrix[S, 2] draws; real x = 0; real y = 0; // must initialize these before the loop so they persist
    for (s in 1:S) {
      real x_ = uniform_rng(x - half_width, x + half_width); 
      real y_ = uniform_rng(y - half_width, y + half_width); // vvv can call previously-declared functions
      real alpha_star = exp(binormal_lpdf([x_, y_] | mu_X, mu_Y, sigma_X, sigma_Y, rho) - 
                            binormal_lpdf([x , y ] | mu_X, mu_Y, sigma_X, sigma_Y, rho));
      if (alpha_star > uniform_rng(0, 1)) { // Q([x, y]) / Q[x_, y_] = 1 in this case
        x = x_; y = y_;
      } // otherwise leave x and y the same as they were on iteration s - 1
      draws[s, ] = [x, y];
    } // x_, y_, and alpha_star all get deleted here but x and y do not
    return draws;
  }
}
