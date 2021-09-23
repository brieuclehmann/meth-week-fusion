

consensus_combine_wrapper <- function(shard_samples) {
  
  if (is.vector(shard_samples[[1]])) 
    shard_samples <- lapply(shard_samples, function(x) matrix(x, ncol = 1))
  
  # number of posterior samples
  M <- nrow(shard_samples[[1]])
  # dimension of posterior
  p <- ncol(shard_samples[[1]])
  
  n_shard <- length(shard_samples)
  precision <- function(x) solve(var(x))
  precision_list <- lapply(shard_samples, precision)
  shard_precisions <- array(NA, c(p, p, n_shard))
  for (shard in seq_along(precision_list)) {
    shard_precisions[ , ,shard] <- precision_list[[shard]]
  }
  sum_prec <- apply(shard_precisions, c(1, 2), sum)
  
  consensus_samples <- matrix(NA, M, p) 
  for (i in seq(M)) {
    this_shard_samples <- matrix(sapply(shard_samples, function(x) x[i,]), nrow = p)
    
    weighted_samples <- array(NA, dim(this_shard_samples))
    for (s in seq(n_shard)) {
      weighted_samples[ ,s] <- shard_precisions[ , ,s] %*% this_shard_samples[ ,s]
    }
    
    sum_weighted_samples <- rowSums(weighted_samples)
    consensus_samples[i,] <- solve(sum_prec) %*% sum_weighted_samples
  }
  
  consensus_samples
}
