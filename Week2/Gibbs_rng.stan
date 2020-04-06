functions { /* saved as Gibbs_rng.stan in R's working directory */
  matrix Gibbs_rng(int S, real mu_X, real mu_Y, real sigma_X, real sigma_Y, real rho) {
    matrix[S, 2] draws; real x = 0; // must initialize before loop so that it persists
    real beta = rho * sigma_Y / sigma_X;
    real lambda = rho * sigma_X / sigma_Y;
    real sqrt1mrho2 = sqrt(1 - square(rho));
    real sigma_YX = sigma_Y * sqrt1mrho2;
    real sigma_XY = sigma_X * sqrt1mrho2; // this is smaller than in binormal_rng.stan !
    for (s in 1:S) {
      real y = normal_rng(mu_Y + beta * (x - mu_X), sigma_YX); // y needs a persistent x
      x = normal_rng(mu_X + lambda * (y - mu_Y), sigma_XY); // overwritten not redeclared
      draws[s, ] = [x, y];
    } // y gets deleted here but x does not
    return draws;
  }
}
