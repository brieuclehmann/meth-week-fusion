source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

library(dplyr)
library(ggplot2)

set.seed(11)

N = 512
K = 4

prior_params <- list(mean = c(0, 0), cov = 2*diag(2))
obs_cov <- 0.5*diag(2)
mixture_coefficient = 0.25

mu_12   <- MASS::mvrnorm(2, prior_params$mean, prior_params$cov)

X1    <- MASS::mvrnorm(N, mu_12[1,], obs_cov)
X2    <- MASS::mvrnorm(N, mu_12[2,], obs_cov)
mask  = matrix(rbinom(2*N, 1, mixture_coefficient), nrow =N)

X = mask*X1 + (1-mask)*X2

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
true_post_samples <- post_sampler(100000, true_posterior_params)
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

dev.new()
out_df <- bind_rows(extract_output(out_random, "random"),
                    truth_df)

ggplot(out_df, aes(V1, V2, color = strategy, 
                   group = interaction(strategy, type, shard),
                   linetype = type)) +
  geom_density2d()

dev.new()
out_df <- bind_rows(extract_output(out_balanced, "balanced"),
                    truth_df)

ggplot(out_df, aes(V1, V2, color = strategy, 
                   group = interaction(strategy, type, shard),
                   linetype = type)) +
  geom_density2d()

dev.new()
out_df <- bind_rows(extract_output(out_clustered, "clustered"),
                    truth_df)

ggplot(out_df, aes(V1, V2, color = strategy, 
                   group = interaction(strategy, type, shard),
                   linetype = type)) +
  geom_density2d()