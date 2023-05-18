


#' @title Returns a threshold that controls the k-FWER at level alpha for the given distribution of knockoff statistics.
#'
#' @description Given a vector of knockoff statistics, returns a threshold that controls the k-FWER at level alpha.
#'
#' @details Given \eqn{k} and \eqn{\alpha}, calculates a parameter \eqn{v} and then returns as a threshold the \eqn{v}th
#' smallest (i.e. most negative) knockoff statistic. This assumes that knockoff statistics are symmetrically distributed
#' about 0 (for example when \eqn{W = U_o - U_k}). The parameter \eqn{v} is calculated under a negative binomial model
#' for the signs of the null knockoff statistics.
#' This function is based off of the paper Familywise error rate control via knockoffs by Lucas Janson and Weijie Su
#'
#' @param W A vector of Knockoff statistics
#' @param k Defines k-FWER
#' @param alpha
#' @returns A threshold that controls the k-FWER
#' @examples
#' knockoff.fwer.threshold(W, 4, 0.05)
#' @name knockoff.kfwer.threshold
#' @export
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