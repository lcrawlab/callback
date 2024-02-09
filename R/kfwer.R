


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
#' @param alpha The desired control for the FWER.
#' @returns A threshold that controls the k-FWER
#' @name knockoff.kfwer.threshold
knockoff.kfwer.threshold <- function(W, k, alpha) {
    
    # determined by equation 3.1 (in Theorem 3.1)
    # NB(r,p) -> p(k) = (k + r - 1 choose k) * (1-p)^k * p^r
    # NB(r,0.5) -> p(k) = (k + r - 1 choose k) * (1/2)^(k + r)

    for (v in 0:k) {
      negative_binomial_cdf <- stats::pnbinom(k, v, 0.5, lower.tail = FALSE)
      
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



#' @title Returns the maximum of the threshold that controls the FDR at level q and the threshold that controls the k-FWER at level alpha.
#'
#' @description Given a vector of knockoff statistics, returns the maximum of the threshold that controls the FDR at 
#' level q and the threshold that controls the k-FWER at level alpha.
#'
#' @details Given \eqn{q}, \eqn{k}, and \eqn{\alpha}, calculates the threshold that controls the FDR at level \eqn{q} and the
#' threshold that controls the \eqn{k}-FWER at level \eqn{\alpha} and returns the maximum of the two.
#'
#' @param W A vector of Knockoff statistics
#' @param q The desired control for the FDR
#' @param k Defines k-FWER
#' @param alpha The desired control for the FWER.
#' @returns The maximum of the threshold that controls the FDR and the threshold that controls the k-FWER
#' @name knockoff.heuristic.threshold
knockoff.heuristic.threshold <- function(W, q, k, alpha) {
    kfwer_threshold <- knockoff.kfwer.threshold(W, k, alpha)

    fdr_threshold <- knockoff::knockoff.threshold(W, fdr=q, offset=1)

    heuristic_threshold <- max(fdr_threshold, kfwer_threshold)

    return(heuristic_threshold)
}