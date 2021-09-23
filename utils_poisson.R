gamma_posterior_parameters <- function(obs_data, prior_shape, prior_scale) {
      
  n <- length(obs_data)
  
  # scale the prior to ensure the right combination
  post_shape <- sum(obs_data) + prior_shape
  post_rate  <- n + prior_scale
    
  list(shape = post_shape, rate = post_rate)
}

# Support function to draw posterior samples using lapply
poisson_draw_posterior_sample = function(obs_data, prior_shape, prior_rate, posterior_sample_size){
  
  posterior_param = gamma_posterior_parameters(obs_data, prior_shape, prior_rate)
  theta_Y    = rgamma(posterior_sample_size, shape = posterior_param$shape, rate = posterior_param$rate)
  
  theta_Y 
}

# Support function to draw posterior parameters using lapply
poisson_posterior_parameters = function(obs_data, prior_shape, prior_rate){
  
  posterior_param = gamma_posterior_parameters(obs_data, prior_shape, prior_rate)

  c(posterior_param$shape, posterior_param$rate)
}

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
