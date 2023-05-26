

# https://en.wikipedia.org/wiki/Conjugate_prior#When_the_likelihood_function_is_a_discrete_distribution
estimate_gamma_poisson <- function(data, alpha, beta) {
    alpha_post <- alpha + sum(data)
    beta_post <- beta + length(data)

    return(list("alpha_post" = alpha_post, "beta_post" = beta_post))
}

# https://en.wikipedia.org/wiki/Conjugate_prior#When_the_likelihood_function_is_a_discrete_distribution
rgamma_poisson <- function(n, alpha, beta) {
    return(rnbinom(n, size=alpha, prob = (beta / (1+beta))))

}