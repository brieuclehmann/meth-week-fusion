


normal_mean_posterior <- function(obs_data, obs_cov, prior_mean, prior_cov) {
  
  if (all(prior_cov == 0)) {
    
    post_cov  <- prior_cov
    post_mean <- prior_mean
    
  } else {
    
    n <- dim(obs_data)[1]
    
    post_cov  <- solve(solve(prior_cov) + n * solve(obs_cov))
    post_mean <- post_cov %*% ((solve(prior_cov) %*% prior_mean) +
                                 (n * (solve(obs_cov) %*% colMeans(obs_data))))
    
  }
  
  list(mean = post_mean, cov = post_cov)
}


#' @param shard_samples A p x s matrix of samples, one column per shard.
#' @param shard_vars An p x p x s array of shard covariance matrices.

consensus_combine <- function(shard_samples, shard_vars) {
  
  total_var <- apply(shard_vars, c(1, 2), sum)
  
  
}