# test data (True gaussian mixture)
library(mvtnorm)
set.seed(110104)
nclust <- 4
d <- 2
mu <- list(c(.25, .25), c(.75, .75), c(.25, .75), c(.75, .25))
S <- list(matrix(.1^2 * c(1, .5, .5, 1), 2, 2), .05^2 * diag(2),
          matrix(.1^2 * c(1, -.5, -.5, 1), 2, 2), .075^2 * diag(2))
n <- 100
k <- sample(1:2, n, TRUE)
y <- matrix(0, n, d)
for (i in 1:nclust) {
  idx <- (k == i)
  print(mu[[i]])
  y[idx, ] <- rmvnorm(sum(idx), mu[[i]], S[[i]])
}
summary(y)
plot(y, xlim = c(0, 1), ylim = c(0, 1), bg = k, pch = 21,
     main = "Simulated normal mixture data")

res <- dp_gaussian_mixture(y[1:100, ], 
                           alpha = 1, 
                           lambda = runif(2), 
                           Sigma = 1 * diag(2), 
                           kappa = 10, 
                           nu = 30,
                           Omega = diag(2))



