library(rstanarm)
roaches <- roaches[roaches$roach1 > 0, ]; str(roaches)

roaches$log_roach1 <- log(roaches$roach1)
roaches_ <- t(replicate(1000, {
  alpha_ <- rnorm(1, mean = log(42), sd = 3)
  beta_lag <- rnorm(1, mean = 1, sd = 1)
  beta_trt <- rnorm(1, mean = 0, sd = 0.5)
  beta_snr <- rnorm(1, mean = 0, sd = 0.5)
  
  eta_ <- with(roaches, alpha_ + log(exposure2) - mean(log(exposure2)) + 
                 beta_lag * (log_roach1 - mean(log_roach1)) + 
                 beta_trt * (treatment - mean(treatment)) + 
                 beta_snr * (senior - mean(senior)))
  mu_ <- exp(eta_) # too big!
  roaches_ <- rpois(n = length(mu_), mu_)
  roaches_
}))
