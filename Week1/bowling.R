# computes the x-th Fibonacci number without recursion and with vectorization
F <- function(x) {
  stopifnot(is.numeric(x), all(x == as.integer(x)))
  sqrt_5 <- sqrt(5)
  golden_ratio <- 0.5 * (1 + sqrt_5)
  return(round(golden_ratio ^ (x + 1) / sqrt_5))
}

# probability of knocking down x out of n pins
Pr <- function(x, n = 10) return(ifelse(x > n, 0, F(x)) / (-1 + F(n + 2)))

Omega <- 0:10 # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
names(Omega) <- as.character(Omega)

# probability of knocking down x1 pins on the first roll and x2 on the second
joint_Pr <- matrix(0, nrow = 11, ncol = 11)
rownames(joint_Pr) <- colnames(joint_Pr) <- as.character(Omega)
for (x1 in Omega) {
  Pr_x1 <- Pr(x1)
  for (x2 in 0:(10 - x1))
    joint_Pr[x1 + 1, x2 + 1] <- Pr_x1 * Pr(x2, 10 - x1)
}

total <- row(joint_Pr) - 1L + col(joint_Pr) - 1L
