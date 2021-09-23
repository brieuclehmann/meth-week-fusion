N = 512
K = 4

prior_params <- list(mean = 0, cov = 1)
obs_cov <- 1

mu   <- MASS::mvrnorm(1, prior_params$mean, prior_params$cov)
X    <- MASS::mvrnorm(N, mu, obs_cov)

post_update <- function(y, prior_params, n_shard) {
 normal_mean_posterior(y, obs_cov, prior_params$mean, prior_params$cov * n_shard)
}

post_sampler <- function(M, params) {
  MASS::mvrnorm(M, params$mean, params$cov)
}

out_random <- run_experiment(X, K, random_shards, prior_params, 
                             post_update, post_sampler, consensus_combine_wrapper)
                      
out_balanced <- run_experiment(X, K, balanced_shards, prior_params, 
                               post_update, post_sampler, consensus_combine_wrapper)

out_clustered <- run_experiment(X, K, clustered_shards, prior_params, 
                                post_update, post_sampler, consensus_combine_wrapper)

shard_df <- plyr::ldply(out_random$shard_samples, as_tibble, .id = "shard")

