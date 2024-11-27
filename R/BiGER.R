#' BiGER algorithm for gene list aggregation.
#' 
#' This is the BiGER algorithm using a Gibbs sampler. For boundary approximation,
#' see \link[BiGER]{bBiGER()}; for efficiency, use \link[BiGER]{vBiGER()} for a fast variational
#' inference-based method.
#'
#' @param r A integer rank matrix with rows as genes and columns as studies. All
#' ties are parametrized as ties with the `min` method. See [preprocess_genelist()]
#' for constructing this matrix.
#' @param n_r A vector containing the number of ranked items in each gene list.
#' @param n_u A vector containing the number of unranked items in each gene list.
#' @param mu The initialization for the global importance vector. Optional.
#' @param sigma2 The initialization for study variance vector. Optional.
#' @param w The initialization for the latent weight matrix. Optional.
#' @param alpha The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
#' distribution. Defaults to `1`.
#' @param beta The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
#' distribution. Defaults to `1`.
#' @param chains The number of chains to run. Defaults to 1.
#' @param save_chains Whether to return the full chains of the MCMC run. If `FALSE`, only
#' posterior means of parameters are returned. For large studies, saving the full
#' chains can have a considerable RAM cost. Defaults to `FALSE`.
#' @param save_burnin Whether to save the burn-in period of the MCMC run, given
#' that `save_chains` is set to `TRUE`. This is used mainly for diagnostic
#' purposes. Defaults to `FALSE`.
#' @param iter The number of MCMC iterations to run. Defaults to `5000`.
#' @param burnin The number of iterations to discard as burn-in. Defaults to `2500`.
#' @param verbose The number of iteration to print progress. No printouts are
#' shown when set to `-1` as the default.
#' @return A nested list containing chains of a MCMC bBiGER run. Each chain
#' itself is a list with with posterior mean \eqn{\mu}, posterior mean of \eqn{\sigma^2},
#' and optionally the full chains of \eqn{\mu} and \eqn{\sigma^2}. All genes and
#' studies are reported in the order given in `r`. The full chains are reported
#' in matrices with rows as iterations and columns as parameters.
#' @seealso \link[BiGER]{bBiGER()} for bBiGER and \link[BiGER]{vBiGER()} VI-based vBiGER.
#' @useDynLib BiGER
#' @export

BiGER <- function(r,
                  n_r,
                  n_u,
                  mu=NULL,
                  sigma2=NULL,
                  w=NULL,
                  alpha=1.0,
                  beta=1.0,
                  save_chains=F,
                  save_burnin=F,
                  chains=1,
                  iter=5000,
                  burnin=2500,
                  verbose=-1) {
  
  # Inits
  G <- nrow(r)
  J <- ncol(r)
  
  if (is.null(mu)) {
    mu <- rnorm(G)
  }
  
  if (is.null(sigma2)) {
    sigma2 <- rep(1, J)
  }
  
  if (is.null(w)) {
    w <- init_w_r(G, J, r)
  }
  
  results <- list()
  for (i in 1:chains) {
    results[[paste0("Chain_", i)]] <- cpp_BiGER(r = r,
                                                n_r = n_r,
                                                n_u = n_u,
                                                W = w,
                                                mu = mu,
                                                sigma2 = sigma2,
                                                alpha = alpha,
                                                beta = beta,
                                                save_chains = save_chains,
                                                save_burnin = save_burnin,
                                                iter = iter,
                                                burnin = burnin,
                                                verbose = verbose)
  }
  
  return(results)
}
