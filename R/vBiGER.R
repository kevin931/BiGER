#' Variational BiGER (vBiGER) for efficient rank aggregation.
#' 
#' This is the vBiGER algorithm implemented using Mean-Field Variational Inference.
#' For the original BiGER algorithm based on MCMC, see \link[BiGER]{BiGER()} or
#' \link[BiGER]{bBiGER()}.
#'
#' @param r A integer rank matrix with rows as genes and columns as studies. All
#' ties are parametrized as ties with the `min` method. See [preprocess_genelist()]
#' for constructing this matrix.
#' @param n_r A vector containing the number of ranked items in each gene list.
#' @param n_u A vector containing the number of unranked items in each gene list.
#' @param mu The initialization for the global importance vector. Optional.
#' @param sigma2_inv The initialization for study precision vector (inverse of
#' study variance). Optional.
#' @param alpha The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
#' distribution. Defaults to `1`.
#' @param beta The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
#' distribution. Defaults to `1`.
#' @param max_iter The maximum number of VI iterations to run. Defaults to `100`.
#' @param delta The numerical tolerance as the stopping criteria. This is calculated
#' by the average change (delta) of the approximated mean parameters \eqn{\mu_g}.
#' Defaults to `0.0001`.
#' @param verbose The number of iteration to print progress. No printouts are
#' shown when set to `-1` as the default.
#' @return A list containing VI-estimated latent importance \eqn{\mu},
#' study variance \eqn{\sigma^2}, and the convergence monitoring of `delta`.
#' @seealso \link[BiGER]{BiGER()} and \link[BiGER]{bBiGER()} MCMC-based algorithms.
#' @useDynLib BiGER
#' @export

vBiGER <- function(r,
                   n_r,
                   n_u,
                   mu=NULL,
                   sigma2_inv=NULL,
                   alpha=1.0,
                   beta=1.0,
                   max_iter=100,
                   delta=0.0001,
                   verbose=-1) {
  
  # Inits
  G <- nrow(r)
  J <- ncol(r)
  
  if (is.null(mu)) {
    mu <- rnorm(G)
  }
  
  if (is.null(sigma2)) {
    sigma2_inv <- rep(1, J)
  }
  
  
  return(cpp_vBiGER(r = r,
                    n_r = n_r,
                    n_u = n_u,
                    mu = mu,
                    sigma2_inv = sigma2_inv,
                    alpha = alpha,
                    beta = beta,
                    max_iter = max_iter,
                    delta = delta,
                    verbose = verbose))
}
