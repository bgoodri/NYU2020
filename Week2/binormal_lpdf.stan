functions {
  real binormal_lpdf(row_vector xy, real mu_X, real mu_Y, real sigma_X, real sigma_Y, real rho) {
    real beta = rho * sigma_Y / sigma_X; 
    real sigma = sigma_Y * sqrt(1 - square(rho));
    if (is_inf(xy[1]) || is_inf(xy[2])) return negative_infinity();
    return normal_lpdf(xy[1] | mu_X, sigma_X) + 
           normal_lpdf(xy[2] | mu_Y + beta * (xy[1] - mu_X), sigma);
  }
}
