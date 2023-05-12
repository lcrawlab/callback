# This function is based off of the paper 
# Familywise error rate control via knockoffs
# by Lucas Janson and Weijie Su

knockoff.kfwer.threshold <- function(W, k, alpha) {
    
    # determined by equation 3.1 (in Theorem 3.1)
    # NB(r,p) -> p(k) = (k + r - 1 choose k) * (1-p)^k * p^r
    # NB(r,0.5) -> p(k) = (k + r - 1 choose k) * (1/2)^(k + r)

    for (v in 0:k) {
      negative_binomial_cdf <- pnbinom(k, v, 0.5, lower.tail = FALSE)
      
      if (negative_binomial_cdf > alpha) break; # this v is too high and doesn't control the k-FWER at level alpha
    }

    if (v == 0) {
        print("Cannot control the k-FWER for this combination of k and alpha. Returning threshold of infinity.")
        return(Inf)
    }

    v <- v - 1 # go back down to the last v that controlled the k-FWER at level alpha

    sorted_W <- sort(W)

    T_v = abs(sorted_W[v])

    return(T_v)
}