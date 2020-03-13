###################################################
## Aim: MCMC samples from a 100 dimensional normal 
## Varying the spread of the proposal
####################################################

set.seed(10)

if(!library(mvtnorm)) {
  install.packages('mvtnorm')
}
library(mvtnorm)

Ftarget <- function(X)  {
  return(dmvnorm(X))
}

#Proposal is multivariate normal
proposal <- function(x,h) {
  return(rmvnorm(1, mean = x, sigma = diag(h,100,100)))
}

# Functino for generating MCMC output using MH
MH <- function(n, init, h)  {
  acceptance_prob <- 0
  output <- matrix(nrow = n, ncol = 100)
  output[1,] = init
  for(j in 2:n) {
    proposed = proposal(output[j-1,], h)
    alpha = Ftarget(proposed)/Ftarget(output[j-1,])
    if(runif(1) < alpha)  {                             # Acceptance step
      output[j,] = proposed
      acceptance_prob = acceptance_prob+1
    } else  {                                           # Reection step
      output[j,] = output[j-1,]
    }
  }
  print(acceptance_prob/n)
  return(output)
}


init <- numeric(length = 100)
n <- 1e3
hval <- c(.001,.0005,0.0001)                          # proposal spread needs to be really small to ensure a
                                                      # reasonable acceptance rate
par(mfrow = c(3,1))
for(i in 1:3) {
  output <- MH(n, init, hval[i])
}

#print(acceptance_prob/n)

#Plotting the output

#When 100 parameters are being updated together, it just requires one of the 100 proposed values to be bad for the full vector
#to be rejected, therefore we need really small updates for them to be accepted
plot.ts(output[,1])
plot.ts(output[,2])
plot.ts(output[,3])

par(mfrow = c(3,1))
acf(output[,1])
acf(output[,2])
acf(output[,3])

par(mfrow = c(3,1))
plot(density(output[,1]))
plot(density(output[,2]))
plot(density(output[,3]))

lupost <- function(X) {
  #var <- runif(100)
  return(log(dmvnorm(X)))
}

library(mcmc)
tenet <- runif(100)
scale <- runif(100, min = 0, max = .5)
init = numeric(length = 100)
out <- metrop(lupost, initial = init, nbatch = 1e3, scale = .25)
out$accept
out$batch[900:1e3, 1:10]
