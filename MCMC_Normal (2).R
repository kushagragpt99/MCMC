#####################################
## Aim: MCMC samples from a N(0,1)
## Varying the spread of the proposal
#####################################

set.seed(10)

# Normal target distribution
Ftarget <- function(X)  {
  return(dnorm(X))
}

# Proposal is Normal(x,h)
proposal <- function(x,h) {
  return(rnorm(1, mean = x, sd = sqrt(h)))
}


init <- 0                                       # Good starting value
n <- 1e3                                        # number of samples
output <- matrix(nrow = n, ncol = 3)            # Matrix to store MCMC output
hval <- c(1,5,10)                               # Varying the spread of proposal
acceptance_prob <- numeric(length = 3)
output[1,] = init
alpha <- 0
proposed <- 0


# Running the MH algorithm for proposals N(x,1), N(x,5), N(x,10)
for(i in 1:3) {
  h <- hval[i]
  for(j in 2:n) {
    proposed = proposal(output[j-1, i], h)
    alpha = Ftarget(proposed)/Ftarget(output[j-1, i])
    if(runif(1) < alpha)  {                                     # Acceptance step
      output[j,i] = proposed
      acceptance_prob[i] = acceptance_prob[i]+1
    } else  {                                                   # Rejection step
      output[j,i] = output[j-1, i]
    }
  }
}

# Plotting MCMC output
print(acceptance_prob/n)
par(mfrow = c(3,1))
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