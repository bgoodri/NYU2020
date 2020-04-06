// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
#include <Rcpp.h>
#include <stan/math.hpp>  // pulls in everything; could be more specific with included headers

// [[Rcpp::export]]
std::vector<double> g(double x, double y, 
                      double mu_X, double mu_Y, double sigma_X, double sigma_Y, double rho) {
  auto x_var = stan::math::to_var(x); auto y_var = stan::math::to_var(y);
  std::vector<stan::math::var> theta; theta.push_back(x_var); theta.push_back(y_var);
  stan::math::var lp = stan::math::normal_lpdf(x_var, mu_X, sigma_X) + 
                       stan::math::normal_lpdf(y_var, mu_Y + rho * sigma_Y / sigma_X * (x_var - mu_X),
                                               sigma_Y * sqrt(1 - stan::math::square(rho)));
  std::vector<double> grad; lp.grad(theta, grad); return grad;
}
