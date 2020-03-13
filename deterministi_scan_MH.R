##########################################################
## Aim : Deterministic Scan MH for 100 dimensional normal
#########################################################

set.seed(10)

MH <- function(n, init, h)  {
  output <- matrix(0, nrow = n, ncol = 100)                                   # matrix to store output
  output[1,] = init
  accept_prob <- numeric(length = 100)
  for(i in 2:n) {
    
    for(j in 1:100) {                                                       # Looping over every component for the update
      proposed = rnorm(1,output[i-1,j], sd = sqrt(h))                       # Working with log for numerical stability
      alpha = exp( log(dnorm(proposed)) - log(dnorm(output[i-1, j])) )
      
      if(runif(1) < alpha)  {                                               # Acceptance step
        output[i,j] = proposed
        accept_prob[j] = accept_prob[j] + 1 
      } else  {                                                             # Rejection step
        output[i,j] = output[i-1,j]
      }
    }
  }
  accept_prob = accept_prob/n
  return(accept_prob)
}

init <- numeric(length = 100)
n <- 1e3
accept.vec <- MH(n, init, 1)
summary(accept.vec)
