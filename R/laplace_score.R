library(vegan)
library(FastKNN)
library(phyloseq)
library(kernlab)
library(stringr)
library(magrittr)
library(Matrix)

data("GlobalPatterns")

test <- matrix(rnorm(10000), nrow = 1000, ncol = 100)
distance <- dist(test)
# this function constructs a knn graph from a distance matrix
construct_knn_graph <- function(distance_matrix, k) {
  dist <- as.matrix(distance)
  adjacency <- matrix(0, nrow = nrow(dist), ncol = ncol(dist))
  for (i in 1:nrow(dist)) {
    idx <- FastKNN::k.nearest.neighbors(i = i, distance_matrix = dist, k = k)
    adjacency[i,idx] <- 1 
  }
  return(adjacency)
}

# adding weights based on kernel 
add_weights <- function(data, adjacency_mat, kern_name = "rbfdot", ...) {
  supported_kernels <- ls('package:kernlab')[ls('package:kernlab') %>% stringr::str_ends('dot')]
  if (!kern_name %in% supported_kernels){
    rlang::abort("Not a supported kernel in the kernlab package")
  }
  kern_init <- rlang::exec(kern_name, ...)
  kern_mat <- kernelMatrix(kernel = kern_init, data)
  w_adj <- adjacency_mat * kern_mat 
  return(w_adj)
}

adj <- construct_knn_graph(distance, 5)
weights <- add_weights(adj, data = test)

# https://arxiv.org/pdf/1811.07939.pdf
# https://papers.nips.cc/paper/2909-laplacian-score-for-feature-selection.pdf

laplace_score <- function(feature, weighted_adj) {
  n_samp <- length(feature)
  # constructing the degree matrix and the laplacian matrix
  D <- matrix(0, nrow = n_samp, ncol = n_samp)
  diag(D) <- weighted_adj %*% rep(1, n_samp) # degree matrix
  L <- D - weighted_adj # laplacian matrix
  
  num_tilde <- t(feature) %*% D %*% rep(1, n_samp) # numerator 
  den_tilde <- t(rep(1, n_samp)) %*% D %*% rep(1, n_samp) # denominator
  f_tilde <- feature - c(num_tilde/den_tilde) * rep(1, n_samp)
  
  score <- (t(f_tilde) %*% L %*% f_tilde)/(t(f_tilde) %*% D %*% f_tilde)
  return(as.numeric(score))
}
