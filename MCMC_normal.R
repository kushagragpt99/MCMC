###########################################
## MCMC samples from a N(0,1)
## Starting from good starting value
###########################################
set.seed(1)
foo <- seq(-4, 4, length = 1e4)

# We make a function that runs the Markov chain
normal_MC <- function(T, start, h = 1)
{
  MCsamp <- numeric(length = T)
  MCsamp[1] <- start    # get starting value
  
  acc <- 0
  for(t in 2:T)
  {
    proposal <- MCsamp[t-1] + runif(1, min = -h, max = h)   #draw unif proposal
    
    alpha <- min(1, exp( MCsamp[t-1]^2/2 - proposal^2/2) )   #calculate MH ratio
    if(runif(1) < alpha)
    {
      MCsamp[t] <- proposal    #accept
      acc <- acc + 1
    } else{
      MCsamp[t] <- MCsamp[t-1]   #reject
    }
  }
  print(paste("Acceptance prob = ", acc/T))
  return(MCsamp)
}

T <- 1e4

## Running with h = 1
out1 <- normal_MC(T = T, start = 0, h = 1)   # accept prob is high

plot(density(out1))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out1)  #close to 0

## Let's see what happens when we decrease the size of the box
## Running with h = .2
out2 <- normal_MC(T = T, start = 0, h = .2)   # accept prob is high
plot(density(out2))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out2) # a bit far

## We will run again to decrease acceptance
## Running with h = 3.5
out3 <- normal_MC(T = T, start = 0, h = 3.5)   # accept prob is high
plot(density(out3))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out3) # better

## What happens when we increase the jumps even more
## Running with h = 20
out4 <- normal_MC(T = T, start = 0, h = 20)   # accept prob is high
plot(density(out4))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out4)

par(mfrow = c(2,2))
acf(out1, main = "Autocorrelation h = 1")
acf(out2, main = "Autocorrelation h =.2")
acf(out3, main = "Autocorrelation h = 3.5")   # Best ACF
acf(out4, main = "Autocorrelation h = 20")

index <- 1:1000
plot.ts(out1[index], main = "Trace Plot h = 1")
plot.ts(out2[index], main = "Trace Plot h = .2")
plot.ts(out3[index],  main = "Trace Plot h = 3.5")   #Best looking trace plot
plot.ts(out4[index],  main = "Trace Plot h = 20")

plot(density(out1))
plot(density(out2))
plot(density(out3))
plot(density(out4))


######
# In general, for univariate target distributions aim for 44% acceptance rate
# For small dimensional target distribution (2-5) aim for between 40% and 25%
# For large dimensional target (>5) aim for 23% acceptance
######







###########################################
## MCMC samples from a N(0,1)
## Starting from a bad starting point
###########################################
par(mfrow = c(1,1))

## Running with h = 1
out1 <- normal_MC(T = T, start = 10, h = 1)   # accept prob is high

plot(density(out1))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out1)  #close to 0

## Let's see what happens when we decrease the size of the box
## Running with h = .2
out2 <- normal_MC(T = T, start = 10, h = .2)   # accept prob is high
plot(density(out2))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out2) # a bit far

## We will run again to decrease acceptance
## Running with h = 3.5
out3 <- normal_MC(T = T, start = 10, h = 3.5)   # accept prob is high
plot(density(out3))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out3) # better

## What happens when we increase the jumps even more
## Running with h = 20
out4 <- normal_MC(T = T, start = 10, h = 20)   # accept prob is high
plot(density(out4))   # density from the Markov chain draws
lines(foo, dnorm(foo), col = "red") # true density

mean(out4)

par(mfrow = c(2,2))
acf(out1, main = "Autocorrelation h = 1")
acf(out2, main = "Autocorrelation h =.2")
acf(out3, main = "Autocorrelation h = 3.5")   # Best ACF
acf(out4, main = "Autocorrelation h = 20")

index <- 1:1000
plot.ts(out1[index], main = "Trace Plot h = 1")
plot.ts(out2[index], main = "Trace Plot h = .2")
plot.ts(out3[index],  main = "Trace Plot h = 3.5")   #Best looking trace plot
plot.ts(out4[index],  main = "Trace Plot h = 20")


######
# In general, for univariate target distributions aim for 44% acceptance rate
# For small dimensional target distribution (2-5) aim for between 40% and 25%
# For large dimensional target (>5) aim for 23% acceptance
######
