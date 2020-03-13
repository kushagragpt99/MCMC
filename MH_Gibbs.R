#######################################################
## Aim : MH and Gibbs for bivariate normal target
## Comparing the performance of the two samplers
#######################################################

set.seed(10)
library(mvtnorm)

Ftarget <- function(X, rho)  {
  sigma <- rbind(c(1, rho), c(rho, 1))
  return(dmvnorm(X, mean = c(2, -2), sigma = sigma))
}

# Target is bivariate normal with 0 off diagonal terms
proposal <- function(x,h) {
  return(rmvnorm(1, mean = x, sigma = diag(h,2,2)))
}

# Function to calculate MH output
MH <- function(n, init, h, rho)  {
  output <- matrix(0, nrow = n, ncol = 2)
  output[1,] = init
  accept_prob <- 0
  for(i in 2:n) {
    proposed = proposal(output[i-1,], h)                                              # Working with log for numerical stability
    alpha = exp( log(Ftarget(proposed, rho)) - log(Ftarget(output[i-1,], rho)) )
    if(runif(1) < alpha)  {                                                           # Acceptance step
      output[i,] = proposed
      accept_prob = accept_prob + 1
    } else  {                                                                         # Rejection step
      output[i,] = output[i-1,]
    }
  }
  accept_prob = accept_prob/n
  print(paste0("acceptance probaility is :", accept_prob))
  return(output)
}

# Function to calculate Gibbs output
Gibbs <- function(n = 1e3, init = c(0,0), rho) {
  output <- matrix(0, nrow = n, ncol = 2)
  output[1,] = init                                                                   # Deterministic Scan Gibbs sampler
  output[1,2] = rnorm(1, mean = -2 + rho*(output[1,1] - 2), sd = sqrt(1 - rho^2))
  for(i in 2:n) {
    output[i,1] = rnorm(1, mean = 2 + rho*(output[i-1,2] + 2), sd = sqrt(1 - rho^2))
    output[i,2] = rnorm(1, mean = -2 + rho*(output[i,1] - 2), sd = sqrt(1 - rho^2))
  }
  return(output)
}

init <- numeric(length = 2)
rho_vec <- c(0, 0.5, 0.99)
h_vec <- c(3, 2.3, 0.1)
n <- 1e3
chain.mh <- matrix(nrow = n, ncol = 2)
chain.gibbs <- matrix(nrow = n, ncol = 2)
#par(mfrow=c(1,1,1,1))
t <- rnorm(n, mean = 2)


# True plots for marginals
plot(density(t), main = "True marginal for X1")
plot(density(rnorm(n, mean = -2)), main = "True marginal for X2")

# Comparsision of the two methods for different values of correlation
for(i in 1:3) {
    par(mfrow=c(2,1))
    chain.mh = MH(n, init, h_vec[i], rho_vec[i])
    chain.gibbs = Gibbs(n,init, rho_vec[i])
    print(paste0("rho is :", rho_vec[i]))
    plot(density(chain.mh[,1]), main = "marginal for X1 in chain.mh")
    plot(density(chain.mh[,2]), main = "marginal for X2 in chain.mh")
    plot.ts(chain.mh)
    plot(acf(chain.mh))
    plot(density(chain.gibbs[,1]), main = "marginal for X1 in chain.gibbs")
    plot(density(chain.gibbs[,2]), main = "marginal for X2 in chain.gibbs")
    plot.ts(chain.gibbs)
    plot(acf(chain.gibbs))
}
