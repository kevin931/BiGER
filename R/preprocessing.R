#' Preprocess gene lists  
#' 
#' This function is used to generate rank matrices used by the BiGER algorithms.
#' The input format is a list of ordered or unordered genes, each representing
#' results from a single study. 
#' 
#' @param genelist A list of of studies containing genes. For ranked studies,
#                 the genes are in the order shown. For unranked studies,
#                 the genes can be in any order.
#' @param ranking_status A boolean vector indicating whether each study
#' is ranked or unranked. This is unnecessary if `n_r` or `n_u` are provided,
#' or if the list has mixed gene lists.
#' @param n_r: The number of ranked items in each gene list. This is necessary
#' if there are mixed gene lists.
#' @param n_u: The number of unranked items in each gene list. This is necessary
#' if there are mixed gene lists.
#' @param na_as_ties: Whether missing genes from each study is considered as
#' bottom ties. If FALSE, they are considered as missing genes (NAs).
#' @param tie_method: The method used to rank tied genes. Acceptable methods
#' are: "min", "mean", "median", and "max".
#' @param standardize: Whether to perform min-max standardization to the rank
#'                     matrix.
#' @returns A list containing a numeric rank matrix, the number of ranked items
#' in each gene list, and the number of unranked items in each gene list.
#' @export
preprocess_genelist <- function(genelist,
                                ranking_status=NULL,
                                n_r=NULL,
                                n_u=NULL,
                                na_as_ties = FALSE,
                                tie_method = "min",
                                standardize = FALSE) {
  
  
  # Preprocess to pad with NAs
  genelist_len <- unlist(lapply(genelist, length))
  if (length(unique(genelist_len)) > 0) {
    max_len <- max(genelist_len)
    for (i in 1L:length(genelist)) {
      if (length(genelist[[i]]) != max_len) {
        genelist[[i]] <- c(genelist[[i]], rep(NA, max_len - length(genelist[[i]])))
      }
    }
    rm(max_len, i, genelist_len)
  }
  
  prev_length <- length(genelist[[1]])
  max_length <- length(genelist[[1]])
  fix <- FALSE
  
  
  # Unique Genes
  genelist_mat <- as.matrix(as.data.frame(genelist))
  all_genes <- unique(c(genelist_mat))
  if (any(is.na(all_genes))) {
    all_genes <- all_genes[!is.na(all_genes)]
  }
  
  # Sanity Check
  if (is.null(ranking_status)) {
    if (is.null(n_r) | is.null(n_u)) {
      stop("Either `ranking_status` or `n_r` and `n_u` have to be specified." )
    }
  }
  
  # Ranking Information
  if (is.null(n_r) | is.null(n_u)) {
    n_ranked <- numeric(ncol(genelist_mat))
    n_r <- numeric(ncol(genelist_mat))
    n_u <- numeric(ncol(genelist_mat))
    
    for (i in 1:ncol(genelist_mat)) {
      n_ranked[i] <- sum(!is.na(genelist_mat[,i]))
      if (ranking_status[i] == 1) {
        n_r[i] <- sum(!is.na(genelist[[i]]))
        n_u[i] <- 0
      } else {
        n_r[i] <- 0
        n_u[i] <- sum(!is.na(genelist[[i]]))
      }
    }
  }
  
  # Transformation
  r <- matrix(NA, nrow=length(all_genes), ncol=ncol(genelist_mat))
  rownames(r) <- all_genes
  
  for (i in 1:ncol(genelist_mat)) {
    if (n_r[i] > 0) {
      r[genelist_mat[,i][1:n_r[i]], i] <- 1:n_r[i]
    }
    if (n_u[i] > 0) {
      if (tie_method == "min") {
        r[genelist_mat[,i][(n_r[i]+1):(n_r[i]+n_u[i])], i] <- n_r[i]+1
      } else if (tie_method == "mean") {
        r[genelist_mat[,i][(n_r[i]+1):(n_r[i]+n_u[i])], i] <- mean((n_r[i]+1):(n_r[i]+n_u[i]))
      } else if (tie_method == "median") {
        r[genelist_mat[,i][(n_r[i]+1):(n_r[i]+n_u[i])], i] <- median((n_r[i]+1):(n_r[i]+n_u[i]))
      } else if (tie_method == "max") {
        r[genelist_mat[,i][(n_r[i]+1):(n_r[i]+n_u[i])], i] <- n_r[i]+n_u[i]
      }
    }
  }
  
  ## Bottom Ties if any
  if (na_as_ties) {
    for (i in 1:ncol(genelist_mat)) {
      if (tie_method == "min") {
        r[is.na(r[,i]), i] <- n_r[i] + n_u[i] + 1
      } else if (tie_method == "mean") {
        r[is.na(r[,i]), i] <- mean((n_r[i] + n_u[i] + 1):(length(all_genes)))
      } else if (tie_method == "median") {
        r[is.na(r[,i]), i] <- median((n_r[i] + n_u[i] + 1):(length(all_genes)))
      } else if (tie_method == "max") {
        r[is.na(r[,i]), i] <- length(all_genes)
      }
    }
  }
  
  if (standardize) {
    r <- apply(r, 2, function(x) x/max(x))
  }
  
  return(list(r = r, n_r = n_r, n_u = n_u))
}


# rank_to_genelist <- function(r, genes, n_r, n_u, mixed_method="discard") {
#   # This function converts a rank matrix to a genelist, which is the reverse
#   # process to the `preprocess_genelist` function.
#   
#   
#   # Parameters:
#   # (matrix) r: A rank matrix (num_genes x num_studies) with ranks for each gene
#   #             in eact study.
#   # (vector) genes: A list of all the genes present. The order must be the same
#   #                 as the rows in `r`.
#   # (vector) n_r: The number of ranked items in each study.
#   # (vector) n_u: The number of unranked items in each study.
#   # (bool) na_as_ties: Whether missing genes from each study is considered as
#   #                    bottom ties. If FALSE, they are considered as NAs.
#   # (str) tie_method: The method used to rank tied genes. Acceptable methods
#   #                   are: "min", "mean", "median", and "max".
#   # (bool) standardize: Whether to perform min-max standardization to the rank
#   #                     matrix.
#   
#   # (list) Returns:
#   # (matrix) r: A numeric rank matrix.
#   # (vec) n_r: The number of ranked items in each gene list.
#   # (vec) n_u: The number of unranked items in each gene list.
#   
#   na_rm <- function(vec) {
#     vec[is.na(vec)] <- F
#     return(vec)
#   }
#   
#   # Sanity Check
#   if (!mixed_method %in% c("discard", "new", "append")) {
#     stop("`mixed_method` must be of value `discard` or `new`.")
#   }
#   
#   # Gene List
#   genelist <- list()
#   i <- 1
#   
#   for (s in 1:ncol(r)) {
#     # Ranked Study
#     if (n_r[s] > 0 & n_u[s] == 0) {
#       genelist[[i]] <- vector(mode = "character", length = n_r[s])
#       for (g in 1:n_r[s]) {
#         genelist[[i]][g] <- genes[na_rm(r[,s] == g)]
#       }
#       i <- i + 1
#     }
#     
#     # Unranked Study
#     if (n_r[s] == 0 & n_u[s] > 0) {
#       genelist[[i]] <- genes[na_rm(r[,s] == min(r[,s], na.rm = T))]
#       i <- i + 1
#     }
#     
#     # Mixed Study
#     if (n_r[s] > 0 & n_u[s] > 0) {
#       # Ranked items
#       genelist[[i]] <- vector(mode = "character", length = n_r[s])
#       for (g in 1:n_r[s]) {
#         genelist[[i]][g] <- genes[na_rm(r[,s] == g)]
#       }
#       i <- i + 1
#       
#       # Unranked items
#       unranked_genes <- genes[na_rm(r[,s] > n_r[s] & r[,s] <= (n_r[s] + n_u[s]))]
#       if (mixed_method == "new") {
#         genelist[[i]] <- unranked_genes
#         i <- i + i
#       } else if (mixed_method == "append") {
#         genelist[[i-1]] <- append(genelist[[i-1]], unranked_genes)
#       }
#     }
#   }
#   return(genelist)
# }
# 
# 
# mixed_to_top <- function(r, n_r, n_u, unranked_study = "remove") {
#   # This function converts a rank matrix with mixed lists to a top-ranked-only
#   # matrix, which is supported by the original BiG.
#   
#   
#   # Parameters:
#   # (matrix) r: A rank matrix (num_genes x num_studies) with ranks for each gene
#   #             in eact study.
#   # (vector) n_r: The number of ranked items in each study.
#   # (vector) n_u: The number of unranked items in each study.
#   # (str) unranked_study: The method to preprocess unranked studies. `remove`
#   #                       removes the studies; `random` randomly assigns top
#   #                       ranks to the unranked items.
#   
#   # (list) Returns:
#   # (matrix) r: A numeric rank matrix.
#   # (vec) n_r: The number of ranked items in each gene list.
#   # (vec) n_u: The number of unranked items in each gene list.
#   
#   
#   
#   remove <- c()
#   for (s in 1:ncol(r)) {
#     # Top Ranked
#     if (n_r[s] > 0 & n_u[s] == 0) {
#       next
#     }
#     
#     # Top Unranked
#     if (n_r[s] == 0 & n_u[s] > 0) {
#       if (unranked_study == "remove") {
#         remove <- append(remove, s)
#       } else if (unranked_study == "random") {
#         r[r[,s]==1,s] <- sample(1:n_u[s])
#         n_r[s] <- n_u[s]
#       }
#     }
#     
#     # Mixed Study
#     if (n_r[s] >0 & n_u[s] > 0) {
#       r[r[,s]==n_r[s]+n_u[s]+1,s] <- n_r[s]+1
#     }
#   }
#   
#   if (length(remove)  == ncol(r)) {
#     stop("All studies are unranked! BiG does not support this configuration! Try
#          setting `unranked_study` to 'random'.")
#   }
#   
#   if (length(remove) > 0 & length(remove) < ncol(r)) {
#     r <- r[,-remove]
#     n_r <- n_r[-remove]
#     n_u <- n_u[-remove]
#   }
#   
#   return(list(r = r, n_r = n_r, n_u = rep(0, ncol(r))))
# }
# 
# init_fixed <- function(r, n_r, n_u) {
#   L <- matrix(nrow = nrow(r), ncol = ncol(r))
#   U <- matrix(nrow = nrow(r), ncol = ncol(r))
#   w <- matrix(nrow = nrow(r), ncol = ncol(r))
#   
#   boundaries <- seq(0, 1, length.out = nrow(r) + 1)
#   
#   for (s in 1:ncol(r)) {
#     for (g in 1:nrow(r)) {
#       if (r[g, s] <= n_r) {
#         U[g, s] <- boundaries[nrow(r)-r[g,s]+2]
#         L[g, s] <- boundaries[nrow(r)-r[g,s]+1]
#         w[g, s] <- runif(1, min = L[g,s], max = U[g,s])
#       } else if (r[g, s] == n_r[s] + 1) {
#         U[g, s] <- boundaries[nrow(r)-r[g,s]+2]
#         L[g, s] <- boundaries[nrow(r)-r[g,s]+1]
#         w[g, s] <- runif(1, min = L[g,s], max = U[g,s])
#       } else if (r[g, s] ==  n_r[s] + n_u[s] + 1) {
#         U[g, s] <- boundaries[nrow(r)-r[g,s]+2]
#         L[g, s] <- boundaries[nrow(r)-r[g,s]+1]
#         w[g, s] <- runif(1, min = L[g,s], max = U[g,s])
#       } else {
#         U[g, s] <- boundaries[nrow(r)-r[g,s]+2]
#         L[g, s] <- boundaries[nrow(r)-r[g,s]+1]
#         w[g, s] <- runif(1, min = L[g,s], max = U[g,s])
#       }
#       
#       
#     }
#   }
#   
#   return(list(U = U, L=U, w=w))
# }
# 
# 
# 
# genelist_to_maic <- function(genelist, ranking_status, n_genes) {
#   
#   
#   n_study <- length(genelist)
#   df <- matrix(data = "", nrow = n_study, ncol = n_genes + 4)
#   df[,1] <- paste0("C", 1:n_study)
#   df[,2] <- paste0("Sim", 1:n_study)
#   df[,3] <- ifelse(ranking_status == 1, "RANKED", "NOT_RANKED")
#   df[,4] <- "NAMED_GENES"
#   
#   for (i in 1:n_study) {
#     df[i, (1:length(genelist[[i]]))+4] <- genelist[[i]]
#   }
#   
#   df <- as.data.frame(df)
#   return(df)
# }



