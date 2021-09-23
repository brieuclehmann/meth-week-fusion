source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

set.seed(1)

N = 512
K = 4

prior_params <- list(mean = c(0, 0), cov = diag(2))
obs_cov <- diag(2)

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


### EXTRACT OUTPUT ###

true_posterior_params <- post_update(X, prior_params, 1)
true_post_samples <- post_sampler(500, true_posterior_params)
truth_df <- as_tibble(true_post_samples) %>%
  mutate(type = "full", shard = "1", strategy = "truth")

extract_output <- function(x, name) {
  p <- ncol(x$approx_samples)
  approx_df <- as_tibble(x$approx_samples) %>%
    mutate(type = "full", shard = "1", strategy = name)
  shard_df <- plyr::ldply(x$shard_samples, as_tibble) %>%
    mutate(type = "sub_posterior", strategy = name) %>%
    rename(shard = .id)
  
  bind_rows(approx_df, shard_df)
}

out_df <- bind_rows(extract_output(out_random, "random"),
                    extract_output(out_balanced, "balanced"), 
                    extract_output(out_clustered, "clustered"),
                    truth_df)

ggplot(out_df, aes(V1, V2, color = strategy, 
                   group = interaction(strategy, type, shard),
                   linetype = type)) +
  geom_density2d()
