###################################################
## Aim : MH for Bayes Logit Regression and 
## Gibbs sampler using data augmentation as proposed in https://arxiv.org/pdf/1205.0310.pdf
##################################################
set.seed(10)
library(mcmc)
library(mvtnorm)
library(BayesLogit)
data(logit)

# Log unnormalized posterior with multivariate normal prior
lupost <-function(m, X, y, beta)  {
  val <- 0
  for(i in 1:m) {
    val <- val + (X[i,]%*%beta)*y[i] - log(1 + exp(X[i,]%*%beta) )
  }
  val <- val - (t(beta)%*%(beta))/200
  return(val)
}

# Proposal is multivariate normal
proposal <- function(x, h, p)  {
  return(rmvnorm(1, mean = x, sigma = h*diag(1,p,p)))
}

# MH for Bayesian logistics regression
MH <- function(X, y, h, n, init)  {
  acceptance_prob <- 0
  p <- dim(X)[2]
  m <- dim(X)[1]
  output <- matrix(0, nrow = n, ncol = p)
  output[1,] = init
  for(j in 2:n) {
    proposed <- proposal(output[j-1,], h, p)
    alpha <- exp(lupost(m, X, y, t(proposed)) - lupost(m, X, y, output[j-1,]))    # Working with log for numerical stability
    if(runif(1) < alpha)  {                                                       # Acceptance step
      output[j,] = proposed
      acceptance_prob = acceptance_prob + 1
    } else  {                                                                     # Rejection step
      output[j,] = output[j-1,]
    }
  }
  print(acceptance_prob/n)
  return(output)
}

# Gibbs with data augmentation using polya-gamma latent variable, working details can be found in the link on the header
polya_gamma <- function(X, y, init, n){
  m <- dim(X)[1]
  p <- dim(X)[2]
  output <- matrix(0, nrow = n, ncol = p)
  output[1,] = init
  K <- y - 1/2
  W <- diag(0, m, m)
  for(j in 2:n) {
    for(i in 1:m) {
      W[i,i] = rpg(1,1, (X[i,]%*%output[j-1,]))
    }
    Vw <- solve(t(X)%*%W%*%X + 0.01*diag(1,p,p))
    mw <- Vw%*%(t(X)%*%K)
    output[j,] <- rmvnorm(1, mean = mw, sigma = Vw)
  }
  
  return(output)
}

# Working with logit data
data(logit)
X <- as.matrix(logit[,-1])
y <- as.matrix(logit[,1])
p <- dim(X)[2]
n <- 1e3
init <- rep(0, p)
h <- 0.2
MH_output <- MH(X, y, h, n, init)
PG_output <- polya_gamma(X, y, init, n)

# Plotting output for comparision of MH and Gibbs with augmentation 
par(mfrow = c(2,4))
acf(MH_output[,1],main="MH Comp 1")
acf(MH_output[,2],main="MH Comp 2")
acf(MH_output[,3],main="MH Comp 3")
acf(MH_output[,4],main="MH Comp 4")
acf(PG_output[,1],main="PG Comp 1")
acf(PG_output[,2],main="PG Comp 2")
acf(PG_output[,3],main="PG Comp 3")
acf(PG_output[,4],main="PG Comp 4")
