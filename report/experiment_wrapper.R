
#' @param y A n x p matrix of observations
#' @param n_shard Number of shards
#' @param sharding_func Function to split observations into n_shard shards
#' @param prior_params Vector of prior parameters
#' @param posterior_update Function taking as input prior_params, (subset of) y 
#'   and n_shards and returning posterior_params
#' @param M Number of posterior samples
#' @param posterior_sampler Function taking as input posterior_params and
#'   returning M posterior samples
#' @param sub_combiner Function taking as input subposterior samples and
#'   returning (approximate) posterior samples

run_experiment <- function(y, 
                           n_shard,
                           sharding_func,
                           prior_params, 
                           posterior_update,
                           posterior_sampler,
                           sub_combiner,
                           M = 500) {
  
  # Coerce univariate observations into matrix
  if (is.vector(y)) 
    y <- matrix(y, ncol = 1)
  n <- nrow(y)
  
  # Split observations into shards
  y_sharded <- sharding_func(y, n_shard, prior_params, 
                             posterior_update, posterior_sampler, M)
  
  
  # Draw (sub)posterior samples for each shard
  posterior_params <- lapply(y_sharded, 
                             function(x) posterior_update(x, prior_params, n_shard))
  shard_samples <- list()
  for (i in seq_along(y_sharded)) {
    shard_samples[[i]] <- posterior_sampler(M, posterior_params[[i]])
  }
  names(shard_samples) <- seq_along(shard_samples)
  approx_samples <- sub_combiner(shard_samples)
  
  list(approx_samples = approx_samples, shard_samples = shard_samples)
}


#######################################
# LOGIT
########################################

run_experiment_logit <- function(y, 
                           n_shard,
                           sharding_func,
                           rho_initial, 
                           iterations, 
                           initial_d,
                           sub_combiner,
                           M = 500) {
  
  # Coerce univariate observations into matrix
  if (is.vector(y)) 
    y <- matrix(y, ncol = 1)
  n <- nrow(y)
  
  # Split observations into shards
  y_sharded <- sharding_func(y, n_shard, 
                             rho_initial, 
                             iterations, 
                             initial_d,
                             M)
  
  
  # Draw (sub)posterior samples for each shard
  shard_samples <- list()
  for (i in seq_along(y_sharded)) {
    nc = ncol(y_sharded[[i]])
    yy = y_sharded[[i]][,nc]
    XX = y_sharded[[i]][,1:(nc-1)]
    
    GLM<-glm(yy ~ . -1 ,family=binomial, data=as.data.frame(XX) )
    
    initial_beta = GLM$coefficients
    
    shard_samples[[i]] <- draw_posterior_sample_logit(y_sharded[[i]], 
                                                      initial_beta, 
                                                      rho_initial, 
                                                      iterations, 
                                                      initial_d, )[seq((iterations-M+1),iterations,1),]
  }
  names(shard_samples) <- seq_along(shard_samples)
  approx_samples <- sub_combiner(shard_samples)
  
  list(approx_samples = approx_samples, shard_samples = shard_samples)
}
