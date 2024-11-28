#include <cmath>
#include <stdexcept>
#include <Rcpp.h>
#include <RcppTN.h>
// [[Rcpp::depends(RcppTN)]]

using namespace Rcpp;
using namespace RcppTN;


NumericVector matrix_vector_index(NumericMatrix mat, IntegerVector x, int y) {
  NumericVector vec = mat(_,y);
  return vec[x];
}


NumericVector matrix_vector_index(NumericMatrix mat, int x, IntegerVector y) {
  NumericVector vec = mat(x,_);
  return vec[y];
}


LogicalVector na_rm(LogicalVector a) {
  a[is_na(a)] = 0;
  return a;
}


List get_ranked_index(const NumericMatrix &r,
					  const NumericVector &n_r,
					  const bool na = false) {

	List L;
	for (int s=0; s < r.cols(); s++) {
		IntegerVector study;
		if (n_r[s] > 0) {
			for (int i=1; i<=n_r[s]; i++) {
				for (int j=0; j<r.rows(); j++) {
					if (r(j, s) == i) {
						study.push_back(j);
						break;
					}
				}
			}
		}
		L.push_back(study);
	}
	return L;
}


List get_unranked_index(const NumericMatrix &r,
					const NumericVector &n_r,
					const NumericVector &n_u,
					const bool na = false) {

	List L;
	for (int s=0; s < r.cols(); s++) {
		LogicalVector unranked = r(_,s) == (n_r[s] + 1);
		if (na) {
			unranked = na_rm(unranked);
		}

		IntegerVector study;
		if (n_u[s] > 0) {
			for (int i=0; i<r.rows(); i++) {
				if (unranked[i]) {
					study.push_back(i);
				}
			}
		}
		L.push_back(study);
	}
	return L;
}


List get_bottom_index(const NumericMatrix &r,
				  const NumericVector &n_r,
				  const NumericVector &n_u,
				  const bool na = false) {

	List L;
	for (int s=0; s < r.cols(); s++) {
		LogicalVector bottom = r(_,s) == (n_r[s] + n_u[s] + 1);
		if (na) {
			bottom = na_rm(bottom);
		}

		IntegerVector study;
		if (sum(bottom) > 0) {
			for (int i=0; i<r.rows(); i++) {
				if (bottom[i]) {
					study.push_back(i);
				}
			}
		}
		L.push_back(study);
	}
	return L;
}


List find_boundaries(const NumericMatrix &r,
					 const NumericVector &n_r,
					 const NumericVector &n_u,
					 const LogicalMatrix &ranked) {
	
	NumericMatrix L(r.rows(), r.cols());
	NumericMatrix U(r.rows(), r.cols());
	NumericVector boundaries(r.rows() + 1);
	boundaries = seq_along(boundaries);
	boundaries = (boundaries - 1)/(r.rows()); 

	for (int s=0; s<r.cols(); s++) {
		for (int g=0; g<r.rows(); g++) {
			// Top Ranked
			if (r(g,s) <= n_r[s]) {
				L(g, s) = boundaries(r.rows()-r(g,s));
				U(g, s) = boundaries(r.rows()-r(g,s)+1);
			} else if (r(g,s) == n_r[s] + 1 && n_u[s] > 0) {
				// Unranked
				L(g, s) = boundaries(r.rows()-n_r[s]-n_u[s]);
				U(g, s) = boundaries(r.rows()-n_r[s]);
			} else if (r(g,s) == n_r[s] + n_u[s] + 1) {
				// Bottom
				L(g, s) = boundaries(r.rows()-sum(ranked(_,s)));
				U(g, s) = boundaries(r.rows()-n_r[s]-n_u[s]);
			}
		}
	}

	List out;
	out.push_back(L);
	out.push_back(U);
	return out;
}



List find_boundaries_norm(const NumericMatrix &r,
					 	  const NumericVector &n_r,
					 	  const NumericVector &n_u,
					 	  const LogicalMatrix &ranked) {
	
	NumericMatrix L(r.rows(), r.cols());
	NumericMatrix U(r.rows(), r.cols());
	NumericVector boundaries(r.rows() + 1);
	boundaries = seq_along(boundaries);
	boundaries = (boundaries - 1)/(r.rows()); 
	double upper;
	double lower;

	for (int s=0; s<r.cols(); s++) {
		for (int g=0; g<r.rows(); g++) {
			// Top Ranked
			if (r(g,s) <= n_r[s]) {
				lower = boundaries(r.rows()-r(g,s));
				upper = boundaries(r.rows()-r(g,s)+1);
				L(g, s) = R::qnorm(lower, 0.0, 1.0, true, false);
				U(g, s) = R::qnorm(upper, 0.0, 1.0, true, false);
			} else if (r(g,s) == n_r[s] + 1 && n_u[s] > 0) {
				// Unranked
				lower = boundaries(r.rows()-n_r[s]-n_u[s]);
				upper = boundaries(r.rows()-n_r[s]);
				L(g, s) = R::qnorm(lower, 0.0, 1.0, true, false);
				U(g, s) = R::qnorm(upper, 0.0, 1.0, true, false);
			} else if (r(g,s) == n_r[s] + n_u[s] + 1) {
				// Bottom
				lower = boundaries(r.rows()-sum(ranked(_,s)));
				upper = boundaries(r.rows()-n_r[s]-n_u[s]);
				L(g, s) = R::qnorm(lower, 0.0, 1.0, true, false);
				U(g, s) = R::qnorm(upper, 0.0, 1.0, true, false);
			}
		}
	}

	List out;
	out.push_back(L);
	out.push_back(U);
	return out;
}

// [[Rcpp::export]]
NumericMatrix init_w_r(int num_genes,int num_studies,const NumericMatrix &r) {

	NumericMatrix W(r.rows(), r.cols());

	for (int i=0; i<r.cols(); i++) {
		NumericVector temp = rnorm(num_genes, 0, 1);
		std::sort(temp.begin(), temp.end(), std::greater<>());
		NumericVector order = r(_,i) - 1;
		temp = temp[order];
		W(_, i) = temp;
	}

	return W;
}


NumericMatrix init_w(const NumericMatrix &lower,
					 const NumericMatrix &upper,
					 const NumericVector &mu,
					 const NumericVector &sigma) {

	NumericMatrix e_w(lower.rows(), lower.cols());
	for (int s=0; s < lower.cols(); s++) {
		for (int g = 0; g < lower.rows(); g++) {
			e_w(g, s) =  RcppTN::etn1(mu(g), sqrt(1/sigma(s)), lower(g, s), upper(g, s));
		}
	}
	return e_w;
}


NumericMatrix init_w2(const NumericMatrix &lower,
					  const NumericMatrix &upper,
					  const NumericVector &mu,
					  const NumericVector &sigma,
					  const NumericMatrix &w) {
	NumericMatrix e_w2(lower.rows(), lower.cols());
	for (int s=0; s < lower.cols(); s++) {
		for (int g = 0; g < lower.rows(); g++) {
			e_w2(g, s) =  pow(w(g, s), 2) + RcppTN::vtn1(mu(g), sqrt(1/sigma(s)), lower(g, s), upper(g, s));	
		}
	}
	return e_w2;
}



NumericMatrix init_w_gibbs(const NumericMatrix &lower,
					 	   const NumericMatrix &upper,
					 	   const NumericVector &mu,
					 	   const NumericVector &sigma,
						   const LogicalMatrix &ranked) {

	NumericMatrix w(lower.rows(), lower.cols());
	for (int s=0; s < lower.cols(); s++) {
		for (int g = 0; g < lower.rows(); g++) {
			if (ranked(g, s)) {
				w(g, s) =  RcppTN::rtn1(mu(g), sqrt(sigma(s)), lower(g, s), upper(g, s));
			}
		}
	}
	return w;
}


// BiGER algorithm for gene list aggregation.
// 
// This is the internal C++ implementation of the BiGER algorithm
// using a Gibbs sampler. The officially R implementation adds a few quality-of-life
// improvements, such as making defaults of `mu`, `sigma2`, and `w` optional as well as removing
// empty outputs when necessary. This version may have a marginally small performance
// advantage.
//
// @param r A integer rank matrix with rows as genes and columns as studies. All
// ties are parametrized as ties with the `min` method. See \link[BiGER]{preprocess_genelist()}
// for constructing this matrix.
// @param n_r A vector containing the number of ranked items in each gene list.
// @param n_u A vector containing the number of unranked items in each gene list.
// @param mu The initialization for the global importance vector. Optional.
// @param sigma2 The initialization for study variance vector. Optional.
// @param w The initialization for the latent weight matrix. Optional.
// @param alpha The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param beta The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param chains The number of chains to run. Defaults to 1.
// @param save_chains Whether to return the full chains of the MCMC run. If `FALSE`, only
// posterior means of parameters are returned. For large studies, saving the full
// chains can have a considerable RAM cost. Defaults to `FALSE`.
// @param save_burnin Whether to save the burn-in period of the MCMC run, given
// that `save_chains` is set to `TRUE`. This is used mainly for diagnostic
// purposes. Defaults to `FALSE`.
// @param iter The number of MCMC iterations to run. Defaults to `5000`.
// @param burnin The number of iterations to discard as burn-in. Defaults to `2500`.
// @param verbose The number of iteration to print progress. No printouts are
// shown when set to `-1` as the default.
// @return A nested list containing chains of a MCMC bBiGER run. Each chain
// itself is a list with with posterior mean \eqn{\mu}, posterior mean of \eqn{\sigma^2},
// and optionally the full chains of \eqn{\mu} and \eqn{\sigma^2}. All genes and
// studies are reported in the order given in `r`. The full chains are reported
// in matrices with rows as iterations and columns as parameters.
// @seealso [bBiGER()] for bBiGER and [vBiGER()] VI-based vBiGER.
// [[Rcpp::export]]
List cpp_BiGER(const NumericMatrix &r,
           	   const NumericVector &n_r,
           	   const NumericVector &n_u,
           	   NumericMatrix &W,
           	   const NumericVector &mu,
           	   const NumericVector &sigma2,
           	   const double alpha=1.0,
		   	   const double beta=1.0,
			   const bool save_chains=false,
			   const bool save_burnin=false,
           	   const int iter=20000,
           	   const int burnin=10000,
           	   const int verbose=-1) {

	int num_genes = r.rows();
	int num_studies = r.cols();
	int iter_kept = iter - burnin;

	// Init Mu
	NumericVector Mu(num_genes);
	Mu = mu;
	NumericVector Mu_post(num_genes);
	Mu_post.fill(0);

	// Init Sigma_S2
	NumericVector Sigma_s2(num_studies);
	Sigma_s2 = sigma2;
	NumericVector Sigma_s2_post(num_studies);
	Sigma_s2_post.fill(0);

	// Chains
	int chain_dim1;
	int chain_dim2;
	int chain_dim3;

	if (save_chains) {
		if (save_burnin) {
			chain_dim1 = iter;
		} else {
			chain_dim1 = iter_kept;
		}
		chain_dim2 = num_genes;
		chain_dim3 = num_studies;
	} else {
		chain_dim1 = 0;
		chain_dim2 = 0;
		chain_dim3 = 0;
	}

	NumericMatrix Mu_keep(chain_dim1, chain_dim2);
	NumericMatrix Sigma_s2_keep(chain_dim1, chain_dim3);

	// Ranked
	LogicalMatrix ranked(num_genes, num_studies);
	for (int i=0; i<num_studies; i++) {
		ranked(_, i) = !is_na(r(_,i));
	}

	List ranked_index = get_ranked_index(r, n_r, true);
	List unranked_index = get_unranked_index(r, n_r, n_u, true);
	List bottom_index = get_bottom_index(r, n_r, n_u, true);

	for (int m=1; m<iter; m++) {

		if (verbose != -1 && m%verbose == 0) {
			Rcout << m << "\n";
		}

		// Update Mu for each gene
		for (int g=0; g<num_genes; g++) {
			NumericVector temp = W(g, _)/Sigma_s2;
			temp = temp[ranked(g,_)];
			double c = sum(temp);

			temp = 1/Sigma_s2;
			temp = temp[ranked(g,_)];
			double d =  sum(temp) + 1;
			Mu(g) = rnorm(1, c/d, sqrt(1/d))(0);

			if (save_chains & save_burnin) {
				Mu_keep(m, _) = Mu;
			}
			
			if (m >= burnin) {
				Mu_post(g) += Mu(g)/iter_kept;
				if (save_chains & !save_burnin) {
					Mu_keep(m-burnin, _) = Mu;
				}
			}
		}

		// Update sigma_s2
		for (int s=0; s<num_studies; s++) {
			NumericVector temp = pow(W(_,s) - Mu, 2);
			temp = temp[ranked(_,s)];
			double rate = sum(temp)/2+beta;
			double shape = sum(ranked(_,s))/2+alpha;
			Sigma_s2(s) = 1/rgamma(1, shape, 1/rate)(0);

			if (save_chains & save_burnin) {
				Sigma_s2_keep(m, _) = Sigma_s2;
			}

			if (m >= burnin) {
				Sigma_s2_post(s) += Sigma_s2(s)/iter_kept;
				if (save_chains & !save_burnin) {
					Sigma_s2_keep(m-burnin, _) = Sigma_s2;
				}
			}
		}

		// Update W
		for (int s=0; s<num_studies; s++) {
			double lower;
			double upper;
			int num_non_na_genes = sum(ranked(_,s));

			IntegerVector study_ranked_index = ranked_index[s];
			IntegerVector study_bottom_index = bottom_index[s];
			IntegerVector study_unranked_index = unranked_index[s];

			// Top Ranked List
			if (n_r[s] > 0) {
				NumericVector genes(n_r[s]);
				genes = seq_along(genes); 
				NumericVector index = sample(genes, n_r[s]);

				for (int g=0; g<index.length(); g++) {
					// The first Gene
					if (index[g] == 1) {
						lower = W(study_ranked_index[1], s);
						upper = R_PosInf;
					} else if (index[g] == n_r[s]) {
						// The Last gene
						if (n_r[s] < num_non_na_genes) {
							if (n_u[s] > 0) {
								lower = max(matrix_vector_index(W, study_unranked_index, s));
							} else {
								lower = max(matrix_vector_index(W, study_bottom_index, s));
							}
						} else {
							lower = R_NegInf;
						}
						upper = W(study_ranked_index[n_r[s]-2], s);

					} else {
						lower = W(study_ranked_index[index[g]], s);
						upper = W(study_ranked_index[index[g]-2], s);
					}

					W(study_ranked_index[index[g]-1], s) = RcppTN::rtn1(Mu(study_ranked_index[index[g]-1]), sqrt(Sigma_s2(s)), lower, upper);
				}
			}

			//unranked List
			if (n_u[s] > 0) {
				if (n_r[s] + n_u[s] < num_non_na_genes) {
					lower = max(matrix_vector_index(W, study_bottom_index, s));
				} else {
					lower = R_NegInf;
				}

				if (n_r[s] == 0) {
					upper = R_PosInf;
				} else {
					upper = W(study_ranked_index[n_r[s]-1], s);
				}

				for (int i=0; i<study_unranked_index.length(); i++) {
					W(study_unranked_index[i], s) = RcppTN::rtn1(Mu(study_unranked_index[i]), sqrt(Sigma_s2(s)), lower, upper);
				}
			}

			//Bottom Ties
			if (n_r[s] + n_u[s] < num_non_na_genes) {
				lower = R_NegInf;
				// Whether there are unranked items
				if (n_u[s] > 0) {
					upper = min(matrix_vector_index(W, study_unranked_index, s));
				} else {
					upper = W(study_ranked_index[n_r[s]-1], s);
				}

				for (int i=0; i<study_bottom_index.length(); i++) {
					W(study_bottom_index[i], s) = RcppTN::rtn1(Mu(study_bottom_index[i]), sqrt(Sigma_s2(s)), lower, upper);
				}
			}
		}
	}


	List L = List::create(Named("mu_post") = Mu_post,
						  _["sigma2_post"]  = Sigma_s2_post,
						  _["mu"]  = Mu_keep,
						  _["sigma2"]  = Sigma_s2_keep);

	return L;
}


// Bounded BiGER (bBiGER) algorithm with fixed boundaries. 
// 
// This is the internal C++ implementation of bounded BiGER (bBiGER) algorithm
// using a Gibbs sampler. The official R implementation adds a few quality-of-life
// improvements, such as making defaults of `mu` and `sigma2` optional as well as removing
// empty outputs when necessary. This version may have a marginally small performance
// advantage.
//
// @param r A integer rank matrix with rows as genes and columns as studies. All
// ties are parametrized as ties with the `min` method. See \link[BiGER]{preprocess_genelist()}
// for constructing this matrix.
// @param n_r A vector containing the number of ranked items in each gene list.
// @param n_u A vector containing the number of unranked items in each gene list.
// @param mu The initialization for the global importance vector. Optional.
// @param sigma2 The initialization for study variance vector. Optional.
// @param alpha The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param beta The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param method_bound A string as either 'normal' or 'uniform' to indicate the method for
// finding boundaries. If the input string does not exactly match the options, am error
// is thrown.
// @param save_chains Whether to return the full chains of the MCMC run. If `FALSE`, only
// posterior means of parameters are returned. For large studies, saving the full
// chains can have a considerable RAM cost. Defaults to `FALSE`.
// @param save_burnin Whether to save the burn-in period of the MCMC run, given
// that `save_chains` is set to `TRUE`. This is used mainly for diagnostic
// purposes. Defaults to `FALSE`.
// @param iter The number of MCMC iterations to run. Defaults to `5000`.
// @param burnin The number of iterations to discard as burn-in. Defaults to `2500`.
// @param verbose The number of iteration to print progress. No printouts are
// shown when set to `-1` as the default.
// @return A nested list with posterior mean \eqn{\mu}, posterior mean of \eqn{\sigma^2},
// and optionally the full chains of \eqn{\mu} and \eqn{\sigma^2}. All genes and
// studies are reported in the order given in `r`. The full chains are reported
// in matrices with rows as iterations and columns as parameters.
// [[Rcpp::export]]
List cpp_bBiGER(const NumericMatrix &r,
                const NumericVector &n_r,
                const NumericVector &n_u,
                const NumericVector &mu,
                const NumericVector &sigma2,
                const double alpha=1.0,
				const double beta=1.0,
				const std::string method_bound="normal",
				const bool save_chains=false,
				const bool save_burnin=false,
                const int iter=5000,
                const int burnin=2500,
                const int verbose=-1) {

	int num_genes = r.rows();
	int num_studies = r.cols();
	int iter_kept = iter - burnin;

	// Init Mu
	NumericVector Mu(num_genes);
	Mu = mu;
	NumericVector Mu_post(num_genes);
	Mu_post.fill(0);

	// Init Sigma_S2
	NumericVector Sigma_s2(num_studies);
	Sigma_s2 = sigma2;
	NumericVector Sigma_s2_post(num_studies);
	Sigma_s2_post.fill(0);

	// Chains
	int chain_dim1;
	int chain_dim2;
	int chain_dim3;

	if (save_chains) {
		if (save_burnin) {
			chain_dim1 = iter;
		} else {
			chain_dim1 = iter_kept;
		}
		chain_dim2 = num_genes;
		chain_dim3 = num_studies;
	} else {
		chain_dim1 = 0;
		chain_dim2 = 0;
		chain_dim3 = 0;
	}

	NumericMatrix Mu_keep(chain_dim1, chain_dim2);
	NumericMatrix Sigma_s2_keep(chain_dim1, chain_dim3);

	// Ranked
	LogicalMatrix ranked(num_genes, num_studies);
	for (int i=0; i<num_studies; i++) {
		ranked(_, i) = !is_na(r(_,i));
	}

	// Bounds
	List bounds;
	
	if (method_bound == "normal") {
	  bounds = find_boundaries_norm(r, n_r, n_u, ranked);
	} else if (method_bound == "uniform") {
	  bounds = find_boundaries(r, n_r, n_u, ranked);
	} else {
	  throw std::invalid_argument("The boundary method must be either 'normal' or 'uniform'.");
	}
	
	
	NumericMatrix lower = bounds[0];
	NumericMatrix upper = bounds[1];

	NumericMatrix W = init_w_gibbs(lower, upper, mu, Sigma_s2, ranked);

	for (int m=1; m<iter; m++) {

		if (verbose != -1 && m%verbose == 0) {
			Rcout << m << "\n";
		}

		// Update Mu for each gene
		for (int g=0; g<num_genes; g++) {
			NumericVector temp = W(g, _)/Sigma_s2;
			temp = temp[ranked(g,_)];
			double c = sum(temp);

			temp = 1/Sigma_s2;
			temp = temp[ranked(g,_)];
			double d =  sum(temp) + 1;
			Mu(g) = rnorm(1, c/d, sqrt(1/d))(0);

			if (save_chains & save_burnin) {
				Mu_keep(m, _) = Mu;
			}
			
			if (m >= burnin) {
				Mu_post(g) += Mu(g)/iter_kept;
				if (save_chains & !save_burnin) {
					Mu_keep(m-burnin, _) = Mu;
				}
			}
		}

		// Update sigma_s2
		for (int s=0; s<num_studies; s++) {
			NumericVector temp = pow(W(_,s) - Mu, 2);
			temp = temp[ranked(_,s)];
			double rate = sum(temp)/2+beta;
			double shape = sum(ranked(_,s))/2+alpha;
			Sigma_s2(s) = 1/rgamma(1, shape, 1/rate)(0);

			if (save_chains & save_burnin) {
				Sigma_s2_keep(m, _) = Sigma_s2;
			}

			if (m >= burnin) {
				Sigma_s2_post(s) += Sigma_s2(s)/iter_kept;
				if (save_chains & !save_burnin) {
					Sigma_s2_keep(m-burnin, _) = Sigma_s2;
				}
			}
		}

		// Update W
		for (int s=0; s<num_studies; s++) {
			for (int g=0; g<num_genes; g++) {
				if (ranked(g, s) == 0) {
					continue;
				}
				W(g, s) = RcppTN::rtn1(Mu(g), sqrt(Sigma_s2(s)), lower(g, s), upper(g, s));
			}
		}
	}


	List L = List::create(Named("mu_post") = Mu_post,
						  _["sigma2_post"]  = Sigma_s2_post,
						  _["mu"]  = Mu_keep,
						  _["sigma2"]  = Sigma_s2_keep);

	return L;
}


// Variational BiGER (vBiGER) for efficient rank aggregation.
// 
// This is the internal C++ implementation of Variational BiGER (vBiGER) algorithm
// using Mean-Field Variational Inference. The official R implementation adds a
// few quality-of-life improvements, such as making defaults of `mu` and `sigma2_inv`
// optional. This version may have a marginally small performance advantage.
//
// @param r A integer rank matrix with rows as genes and columns as studies. All
// ties are parametrized as ties with the `min` method. See \link[BiGER]{preprocess_genelist()}
// for constructing this matrix.
// @param n_r A vector containing the number of ranked items in each gene list.
// @param n_u A vector containing the number of unranked items in each gene list.
// @param mu The initialization for the global importance vector. Optional.
// @param sigma2_inv The initialization for study precision vector (inverse of study
// variance). Optional.
// @param alpha The \eqn{\alpha} (shape) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param beta The \eqn{\beta} (rate) parameter for the Inverse Gamma prior
// distribution. Defaults to `1`.
// @param method_bound A string as either 'normal' or 'uniform' to indicate the method for
// finding boundaries. If the input string does not exactly match the options, am error
// is thrown.
// @param max_iter The maximum number of VI iterations to run. Defaults to `100`.
// @param delta The numerical tolerance as the stopping criteria. This is calculated
// by the average change (delta) of the approximated mean parameters \eqn{\mu_g}.
// Defaults to `0.0001`.
// @param verbose The number of iteration to print progress. No printouts are
// shown when set to `-1` as the default.
// @return A list containing VI-estimated latent importance \eqn{\mu},
// study variance \eqn{\sigma^2}, and the convergence monitoring of `delta`.
// @seealso [BiGER()] and [bBiGER()] MCMC-based algorithms.
// [[Rcpp::export]]
List cpp_vBiGER(const NumericMatrix &r,
              const NumericVector &n_r,
              const NumericVector &n_u,
			  const NumericVector &mu,
			  const NumericVector &sigma2_inv,
              double alpha=1.0,
			  double beta=1.0,
			  const std::string method_bound="normal",
              int max_iter=100,
			  double delta=0.0001,
              int verbose=-1) {

	int num_genes = r.rows();
	int num_studies = r.cols();

	// Inits
	NumericVector m_mu(num_genes);
	NumericVector m_mu1(num_genes);
	NumericVector s2_mu(num_genes);
	NumericVector rate(num_studies);
	NumericVector shape(num_studies);
	NumericVector e_sigma2_inv(num_studies);
	NumericVector convergence;
	e_sigma2_inv = clone(sigma2_inv);

	// Ranked
	LogicalMatrix ranked(num_genes, num_studies);
	for (int i=0; i<num_studies; i++) {
		ranked(_, i) = !is_na(r(_,i));
	}

	// Bounds
	List bounds;
	if (method_bound == "normal") {
	  bounds = find_boundaries_norm(r, n_r, n_u, ranked);
	} else if (method_bound == "uniform") {
	  bounds = find_boundaries(r, n_r, n_u, ranked);
	} else {
	  throw std::invalid_argument("The boundary method must be either 'normal' or 'uniform'.");
	}
	
	NumericMatrix lower = bounds[0];
	NumericMatrix upper = bounds[1];

	NumericMatrix e_w = init_w(lower, upper, mu, e_sigma2_inv);
	NumericMatrix s2_w(lower.rows(), lower.cols());
	NumericMatrix e_w2 = init_w2(lower, upper, mu, e_sigma2_inv, e_w);

	for (int m=1; m<max_iter; m++) {

		if (verbose != -1 && m%verbose == 0) {
			Rcout << m << "\n";
		}

		// Update Mu for each gene
		for (int g=0; g<num_genes; g++) {
			NumericVector temp = e_w(g, _)*e_sigma2_inv;
			temp = temp[ranked(g,_)];
			double c = sum(temp);

			temp = e_sigma2_inv[ranked(g,_)];
			double d =  sum(temp) + 1;
			m_mu1(g) = c/d;
			s2_mu(g) = 1/d;
		}

		convergence.push_back(mean(m_mu1 - m_mu));
		m_mu = clone(m_mu1);

		// Update sigma_s2
		for (int s=0; s<num_studies; s++) {
			NumericVector temp = e_w2(_,s)/2 - e_w(_,s)*m_mu + (s2_mu + m_mu*m_mu)/2;
			temp = temp[ranked(_,s)];
			rate[s] = sum(temp)+beta;
			shape[s] = sum(ranked(_,s))/2+alpha;
			e_sigma2_inv(s) = shape[s]/rate[s];
		}

		// Update W
		for (int s=0; s<num_studies; s++) {
			for (int g=0; g<num_genes; g++) {
				if (ranked(g, s) == 0) {
					continue;
				}
				e_w(g, s) = RcppTN::etn1(m_mu(g), sqrt(1/e_sigma2_inv(s)), lower(g, s), upper(g, s));
				s2_w(g, s) = RcppTN::vtn1(m_mu(g), sqrt(1/e_sigma2_inv(s)), lower(g, s), upper(g, s));
				e_w2(g, s) = pow(e_w(g, s), 2) + s2_w(g, s);
			}
		}

		if (abs(convergence(m-1)) <= delta) {
			break;
		}
	}

	List L = List::create(Named("mu") = m_mu,
		   				  _["sigma2_inv"]  = e_sigma2_inv,
		   				  _["convergence"] = convergence);

	return L;
}
