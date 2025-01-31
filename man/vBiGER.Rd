% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vBiGER.R
\name{vBiGER}
\alias{vBiGER}
\title{Variational BiGER (vBiGER) for efficient rank aggregation.}
\usage{
vBiGER(
  r,
  n_r,
  n_u,
  mu = NULL,
  sigma2_inv = NULL,
  alpha = 1,
  beta = 1,
  method_bound = "normal",
  max_iter = 100,
  delta = 1e-04,
  verbose = -1
)
}
\arguments{
\item{r}{A integer rank matrix with rows as genes and columns as studies. All
ties are parametrized as ties with the `min` method. See [preprocess_genelist()]
for constructing this matrix.}

\item{n_r}{A vector containing the number of ranked items in each gene list.}

\item{n_u}{A vector containing the number of unranked items in each gene list.}

\item{mu}{The initialization for the global importance vector. Optional.}

\item{sigma2_inv}{The initialization for study precision vector (inverse of
study variance). Optional.}

\item{alpha}{The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
distribution. Defaults to `1`.}

\item{beta}{The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
distribution. Defaults to `1`.}

\item{method_bound}{A string as either 'normal' or 'uniform' to indicate the method for
finding boundaries. If the input string does not exactly match the options, am error
is thrown.}

\item{max_iter}{The maximum number of VI iterations to run. Defaults to `100`.}

\item{delta}{The numerical tolerance as the stopping criteria. This is calculated
by the average change (delta) of the approximated mean parameters \eqn{\mu_g}.
Defaults to `0.0001`.}

\item{verbose}{The number of iteration to print progress. No printouts are
shown when set to `-1` as the default.}
}
\value{
A list containing VI-estimated latent importance \eqn{\mu},
study variance \eqn{\sigma^2}, and the convergence monitoring of `delta`.
}
\description{
This is the vBiGER algorithm implemented using Mean-Field Variational Inference.
For the original BiGER algorithm based on MCMC, see \link[BiGER]{BiGER()} or
\link[BiGER]{bBiGER()}.
}
\seealso{
\link[BiGER]{BiGER()} and \link[BiGER]{bBiGER()} MCMC-based algorithms.
}
