
# https://en.wikipedia.org/wiki/Zero-inflated_model#Estimators_of_ZIP_parameters
# https://math.stackexchange.com/questions/2761563/maximum-likelihood-estimation-for-zero-inflated-poisson-distribution
# https://ieeexplore.ieee.org/document/9032203
estimate_zi_poisson <- function(data) {
  num.zeros <- sum(data == 0)
  r0 <- 1 / length(data) * num.zeros
  
  x.bar = mean(data)
  
  gamma <- x.bar / (1 - r0)
  
  lambda.hat <- lamW::lambertW0(-gamma * exp(-gamma)) + gamma
  
  pi.hat <- 1 - x.bar / lambda.hat


  return.list <- list("lambda.hat" = lambda.hat, "pi.hat" = pi.hat)
  return(return.list)
}

rzipoisson <- function(n, lambda, prop.zero) {
  data <- c()


  for (i in 1:n) {
    if (stats::runif(1) < prop.zero) {
      data[i] <- 0
    }
    else {
      data[i] <- stats::rpois(1, lambda)
    }
  }
  return(data)
} 

