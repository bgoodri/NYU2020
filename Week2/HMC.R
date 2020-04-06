library(rstan)
library(rgl)
expose_stan_functions("binormal_lpdf.stan")
Rcpp::sourceCpp("gradient.cpp")

# bivariate normal PDF in log form and negated
dbvn <- function(x, y, mu_X = 0, mu_Y = 0, sigma_X = 1, sigma_Y = 1, rho = 0.75) {
  return(-apply(cbind(x, y), MARGIN = 1, FUN = binormal_lpdf, mu_X = mu_X,
                mu_Y = mu_Y, sigma_X = sigma_X, sigma_Y = sigma_Y, rho = rho))
}

# 3D plot of dbvn. Use mouse to rotate and right-click to zoom in
persp3d(dbvn, xlim = c(-2,2), ylim = c(-2,2), alpha = 0.5, 
        xlab = "x", ylab = "y", zlab = "neg-log-density")

# same as dbvn but without vectorization and also returns gradient wrt x and y
dbvn2 <- function(initial, grad = TRUE, mu_X = 0, mu_Y = 0, sigma_X = 1, sigma_Y = 1, rho = 0.75) {
  x <- initial[1]; y <- initial[2]
  out <- binormal_lpdf(c(x, y), mu_X, mu_Y, sigma_X, sigma_Y, rho)
  if (grad) attributes(out)$grad <- g(x, y, mu_X, mu_Y, sigma_X, sigma_Y, rho)
  return(out)
}

# source some of Radford Neal's functions ( http://www.cs.utoronto.ca/~radford/GRIMS.html )
results <- sapply(c("utilities.r", "mcmc.r", "basic_hmc.r"), FUN = function(x)
  source(paste0("http://www.cs.toronto.edu/~radford/ftp/GRIMS-2012-06-07/", x)))

set.seed(12345)
HMC <- basic_hmc(dbvn2, initial = c(x = 0.9, y = 0.2), nsteps = 700, step = .65, return.traj = TRUE)
pos <- HMC$traj.q
# starting point
ID <- points3d(x = pos[1,1], y = pos[1,2], z = dbvn(pos[1,1], pos[1,2]), col = "green", size = 7)
for (i in 2:nrow(pos)) { # Hamiltonian flow
  x <- pos[i,1]; y <- pos[i,2]; rgl.pop(id = c(ID))
  ID <- points3d(x, y, z = dbvn(x, y), col = "green", size = 7)
}
