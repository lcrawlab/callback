library(lamW)

# https://math.stackexchange.com/questions/2761563/maximum-likelihood-estimation-for-zero-inflated-poisson-distribution

estimate_zi_poisson <- function(data) {
  num.zeros <- sum(data == 0)
  r0 <- 1 /length(data) * num.zeros
  
  x.bar = mean(data)
  
  gamma <- x.bar / (1 - r0)
  
  lambda.hat <- lambertW0(-gamma * exp(-gamma)) + gamma
  
  pi.hat <- (r0 - exp(-lambda.hat)) / (1 - exp(-lambda.hat))
  
  return.list <- list("lambda.hat" = lambda.hat, "pi.hat" = pi.hat)
  return(return.list)
}

rzipoisson <- function(n, lambda, prop.zero) {
  data <- c()
  for (i in 1:n) {
    if (runif(1) < prop.zero) {
      data[i] <- 0
    }
    else {
      data[i] <- rpois(1, lambda)
    }
  }
  return(data)
} 

