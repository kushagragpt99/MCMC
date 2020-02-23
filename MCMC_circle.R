###########################################
## MCMC samples from a uniform distribution
## within a circle
###########################################
set.seed(1)
unif_circle <- function(T, start, h)
{
  xsamp <- numeric(length = T)
  ysamp <- numeric(length = T)

  xsamp[1] <- start[1]   # starting values
  ysamp[1] <- start[2]
  
  acc <- 0
  for(i in 2:T)
  {
    propx <- xsamp[i-1] + runif(1, -h, h)   #uniform square proposal
    propy <- ysamp[i-1] + runif(1, -h, h)
    
    if((propx^2 + propy^2) < 1)
    {
      xsamp[i] <- propx   #both are accepted together
      ysamp[i] <- propy
      acc <- acc + 1
    } else
    {
      xsamp[i] <- xsamp[i-1]
      ysamp[i] <- ysamp[i-1]
    }
  }  
  print(paste("Acceptance prob = ", acc/T))
  return(cbind(xsamp, ysamp))
}

par(mfrow = c(2,2))
T <- 2e4
out1 <- unif_circle(T = T, start = c(.2, .6), h = .2)  # too high acceptance
plot(out1, col = rgb(0,0,0, alpha = .03), asp = 1)  # loos like more or less uniform except the corners
acf(out1[,1])

## again with a larger box to reduce acceptance
out2 <- unif_circle(T = T, start = c(.2, .6), h = 1.5)  # should be between .25 and .44
plot(out2, col = rgb(0,0,0, alpha = .03), asp = 1)  # loos like more or less uniform. Since more rejections, the points are darker
acf(out2[,1])


par(mfrow = c(1,2))
# Seeing marginal distribution distribution. We know that
# the marginal pdf of x is f(x) = 2sqrt(1-x^2)/pi
x <- seq(-1, 1, length = 1e3)
plot(density(out1[,1]))
lines(x, 2*sqrt(1-x^2)/pi, col = "red")

plot(density(out2[,1]))
lines(x, 2*sqrt(1-x^2)/pi, col = "red")

#########
# One thing you can see is that it is difficult to verify and visualize
# a distribution in 2 or more dimensions. What if we didn't know the marginal of x
# then we wouldn't know how well we are doing.
