library(rstanarm)
str(wells)

switch_ <- t(replicate(1000, {
  alpha_ <- rnorm(1, mean = 0, sd = 1.5)
  beta_a <- rnorm(1, mean = 0, sd = 1)
  beta_d <- rnorm(1, mean = 0, sd = 1)

  eta_ <- with(wells, alpha_ +  
                 beta_a * (arsenic - mean(arsenic)) + 
                 beta_d * (dist - mean(dist)))
  epsilon_ <- rlogis(n = length(eta_)) # rnorm(n = length(eta_)) -> "probit" model
  utility_ <- eta_ + epsilon_
  switch_ <- utility_ > 0
  switch_
}))

prop.table(table(c(switch_)))

# this draws from an equivalent prior predictive distribution

switch_ <- t(replicate(1000, {
  alpha_ <- rnorm(1, mean = 0, sd = 1.5)
  beta_a <- rnorm(1, mean = 0, sd = 1)
  beta_d <- rnorm(1, mean = 0, sd = 1)
  
  eta_ <- with(wells, alpha_ +  
                 beta_a * (arsenic - mean(arsenic)) + 
                 beta_d * (dist - mean(dist)))
  mu_ <- plogis(eta_)
  switch_ <- rbinom(n = length(mu_), size = 1, prob = mu_)
  switch_
}))

prop.table(table(c(switch_)))
