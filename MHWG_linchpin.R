#################################################################
## Aim : Metropolis Hasting within Gibbs (MHWG) and Linchpin variable sampler
################################################################

set.seed(10)

# Log unnormalized posterior for MHWG with gamma priors
lupost_MHWG <- function(m, a1, b1, beta, T, lambda) {
  val <- (m+a1-1)*log(beta) + (beta-1)*sum(log(T)) -lambda*sum(T^beta) - b1*beta
  return(val)
}

# Log unnormalized posterior for linchpin with gamma priors
lupost_linchpin <- function(m, a1, a0, b0, b1, beta, T)  {
  val <- (m+a1-1)*log(beta) + (beta-1)*sum(log(T)) - b1*beta - (m+a0)*log(b0+sum(T^beta))
  return(val)
}

# Function to generate output from MHWG
MHWG <- function(T, a0, a1, b0, b1, init, n, h)  {
  m = length(T)
  acceptance_prob <- 0
  output <- matrix(0, nrow = n, ncol = 2)
  output[1,] = init
  output[1,1] = rgamma(1, shape = m+a0, rate = b0 + sum(T^output[1,2]))       
  
  for(j in 2:n) {
    output[j,1] = rgamma(1, shape = m+a0, rate = b0 + sum(T^output[j-1,2]))     # Gibbs step for the first variable
    
    proposed <- rgamma(1, output[j-1,2]^2/h, output[j-1,2]/h)                   # MH for updating the second variable
    if(proposed<0)  {
      proposed = output[j-1,2]
    }
    alpha <- exp( lupost_MHWG(m, a1, b1, proposed, T, output[j,1]) + ((proposed^2)/h-1)*log(proposed) - (proposed/h)*proposed - 
                    lupost_MHWG(m, a1, b1, output[j-1,2], T, output[j,1]) - ((output[j-1,2]^2)/h-1)*log(output[j-1,2]) + (output[j-1,2]/h)*output[j-1,2])
    if(runif(1) < alpha)  {
      output[j,2] = proposed
      acceptance_prob = acceptance_prob+1
    } else  {
      output[j,2] = output[j-1,2]
    }
    
  }
  print(acceptance_prob/n)
  return(output)
}

# Function to generate output from linchpin
linchpin <- function(T, a0, a1, b0, b1, init, n, h) {
  m = length(T)
  acceptance_prob <- 0
  output <- matrix(0, nrow = n, ncol = 2)
  output[1,] = init
  for(j in 2:n) {
    proposed <- rgamma(1, output[j-1,2]^2/h, output[j-1,2]/h)                  # MH for the non-linchpin variable
    if(proposed<0)  {
      proposed = output[j-1,2]
    }
    alpha <- exp( lupost_linchpin(m, a1, a0, b0, b1, proposed, T) + ((proposed^2)/h-1)*log(proposed) - (proposed/h)*proposed - 
                  lupost_linchpin(m, a1, a0, b0, b1, output[j-1,2], T) - ((output[j-1,2]^2)/h-1)*log(output[j-1,2]) + (output[j-1,2]/h)*output[j-1,2])
    if(runif(1) < alpha)  {
      output[j,2] = proposed
      acceptance_prob = acceptance_prob+1
    } else  {
      output[j,2] = output[j-1,2]
    }
    
    output[j,1] = rgamma(1, shape = m+a0, rate = b0 + sum(T^output[j,2]))       # First variable is the linchpin variable
  }
  print(acceptance_prob/n)
  return(output)
}


# Input data
T <- c(387,182,244,600,627, 332,418,300,798,584,660,39,274,174,50, 
       34,1895,158,974,345,1755,1752,473,81,954,1407,230,464,380,131,1205)

# Hyperparameters for the gamma priors
a0= 2.5
b0= 2350
a1= 1
b1= 1
init <- c(1,1)
n=1e3

h = 0.02                                                      # h value for the proposal for reasonable acceptance for linchpin 
linchpin_output <- linchpin(T, a0, a1, b0, b1, init, n, h)

h = 0.009                                                     # h value for the proposal for reasonable acceptance for MHWG
MHWG_output <- MHWG(T, a0, a1, b0, b1, init, n, h)

# Plotting the outputs
par(mfrow = c(2,2))
plot(density(linchpin_output[,1]), xlim=c(0,0.02))
plot(density(linchpin_output[,2]))
plot(density(MHWG_output[,1]), xlim=c(0,0.02))
plot(density(MHWG_output[,2]))

plot.ts(linchpin_output)
plot.ts(MHWG_output)

par(mfrow = c(2,2))
acf(linchpin_output[,1])
acf(linchpin_output[,2])
acf(MHWG_output[,1])
acf(MHWG_output[,2])
