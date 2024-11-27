#' Bounded BiGER (bBiGER) algorithm with fixed boundaries. 
#' 
#' This is the bounded BiGER (bBiGER) algorithm using a Gibbs sampler. This is
#' different from BiGER by fixing the boundaries of \eqn{\omega_{gj}} a priori.
#'
#' @param r A integer rank matrix with rows as genes and columns as studies. All
#' ties are parametrized as ties with the `min` method. See [preprocess_genelist()]
#' for constructing this matrix.
#' @param n_r A vector containing the number of ranked items in each gene list.
#' @param n_u A vector containing the number of unranked items in each gene list.
#' @param mu The initialization for the global importance vector. Optional.
#' @param sigma2 The initialization for study variance vector. Optional.
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
#' @seealso \link[BiGER]{BiGER()} for the vanilla BiGER with Gibbs sampling and
#' \link[BiGER]{vBiGER()} VI-based vBiGER.
#' @useDynLib BiGER
#' @export

bBiGER <- function(r,
                   n_r,
                   n_u,
                   mu=NULL,
                   sigma2=NULL,
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
  
  results <- list()
  for (i in 1:chains) {
    results[[paste0("Chain_", i)]] <- cpp_bBiGER(r = r,
                                                 n_r = n_r,
                                                 n_u = n_u,
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
  
  # if (save_chains & save_warmup) {
  #   n <- iter
  # } else if (save_chains) {
  #   n <- iter-burnin
  # } else {
  #   n <- 1
  # }
  # 
  # # Sequential chains
  # results <- array(dim = c(n, chains, G+J),
  #                  dimnames = list(1:n,
  #                                  1:chains,
  #                                  c(paste0("Mu_", 1:G), paste0("Sigma2_", 1:J)))
  #                  )
  # for (i in 1:chains) {
  #   temp <- cpp_bBiGER(r = r, n_r = n_r, n_u = n_u,
  #                      mu = mu, sigma2 = sigma2,
  #                      alpha = alpha, beta = beta,
  #                      save_chains = save_chains, save_warmup = save_warmup,
  #                      iter = iter, burnin = burnin,
  #                      verbose = verbose)
  #   
  #   # Build Array
  #   if (save_chains) {
  #     results[,i,] <- cbind(temp$mu, temp$sigma2)
  #   } else {
  #     results[1,i,] <- c(temp$mu_post, temp$sigma2_post)
  #   }
  # }
  # 
  # results <- posterior::as_draws(results)
  # return(results)
}
