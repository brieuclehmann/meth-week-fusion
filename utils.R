


normal_mean_posterior <- function(obs_data, obs_cov, prior_mean, prior_cov) {
  
  if (all(prior_cov == 0)) {
    
    post_cov  <- prior_cov
    post_mean <- prior_mean
    
  } else {
    
    if (is.matrix(obs_data)) {
      
      n <- dim(obs_data)[1]
      
      post_cov  <- solve(solve(prior_cov) + n * solve(obs_cov))
      post_mean <- post_cov %*% ((solve(prior_cov) %*% prior_mean) +
                                   (n * (solve(obs_cov) %*% colMeans(obs_data))))
      
    } else {
      
      n = length(obs_data)
      
      post_mean = (prior_cov/(obs_cov/n + prior_cov))*mean(obs_data) + (obs_cov/(obs_cov/n + prior_cov))*prior_mean 
      post_cov  = 1/(1/(prior_cov) + n/(obs_cov))
      
    }

  }
  
  list(mean = post_mean, cov = post_cov)
}


#' @param shard_samples A p x s matrix of samples, one column per shard.
#' @param shard_vars An p x p x s array of shard precision matrices.

consensus_combine <- function(shard_samples, shard_vars) {
  
  n_shards <- ncol(shard_samples)
  total_var <- apply(shard_vars, c(1, 2), sum)
  
  weighted_samples <- array(NA, dim(shard_samples))
  for (s in seq(n_shards)) {
    weighted_samples[ ,s] <- shard_vars[ , ,s] %*% shard_samples[ ,s]
  }
  
  sum_weighted_samples <- rowSums(weighted_samples)
  solve(total_var) %*% sum_weighted_samples
}


accuracy <- function(approx_samples, target_samples) {
  range.x <- range(c(approx_samples, target_samples))
  h_approx <- KernSmooth::dpik(approx_samples)
  approx_density <- KernSmooth::bkde(approx_samples, 
                                     bandwidth = h_approx,
                                     range.x = range.x)
  
  h_target <- KernSmooth::dpik(target_samples)
  target_density <- KernSmooth::bkde(target_samples,
                                     bandwidth = h_target,
                                     range.x = range.x)
  
  f <- approxfun(target_density$x, abs(target_density$y - approx_density$y))
  1 - 0.5 * integrate(f, range.x[1], range.x[2])$value
  
}

