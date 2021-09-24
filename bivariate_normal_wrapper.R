source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

library(ggplot2)
library(dplyr)

extract_output <- function(x, name) {
  M <- nrow(x$approx_samples)
  p <- ncol(x$approx_samples)
  n_shard <- length(x$shard_samples)
  approx_df <- as_tibble(x$approx_samples) %>%
    mutate(type = "full", shard = "1", strategy = name, sample = seq(M))
  shard_df <- plyr::ldply(x$shard_samples, as_tibble) %>%
    mutate(type = "sub_posterior", strategy = name, 
           sample = rep(seq(M), n_shard)) %>%
    rename(shard = .id)
  
  bind_rows(approx_df, shard_df)
}

set.seed(1)

N = 512
K = 4

prior_params <- list(mean = c(0, 0), cov = diag(2))
obs_cov <- diag(2)

post_update <- function(y, prior_params, n_shard) {
  normal_mean_posterior(y, obs_cov, prior_params$mean, prior_params$cov * n_shard)
}

post_sampler <- function(M, params) {
  MASS::mvrnorm(M, params$mean, params$cov)
}

n_sim <- 10
M <- 500
results_df <- tibble()
for (sim in 1:n_sim) {
  set.seed(sim)
  
  mu   <- MASS::mvrnorm(1, prior_params$mean, prior_params$cov)
  X    <- MASS::mvrnorm(N, mu, obs_cov)
  
  out_random <- run_experiment(X, K, random_shards, prior_params, 
                               post_update, post_sampler, 
                               consensus_combine_wrapper, M)
  
  out_balanced <- run_experiment(X, K, balanced_shards, prior_params, 
                                 post_update, post_sampler,
                                 consensus_combine_wrapper, M)
  
  out_clustered <- run_experiment(X, K, clustered_shards, prior_params, 
                                  post_update, post_sampler,
                                  consensus_combine_wrapper, M)
  
  true_posterior_params <- post_update(X, prior_params, 1)
  true_post_samples <- post_sampler(M, true_posterior_params)
  truth_df <- as_tibble(true_post_samples) %>%
    mutate(type = "full", shard = "1", strategy = "truth", sample = seq(M))
  
  out_df <- bind_rows(extract_output(out_random, "random"),
                      extract_output(out_balanced, "balanced"), 
                      extract_output(out_clustered, "clustered"),
                      truth_df) %>%
    mutate(sim = sim)
  
  results_df <- bind_rows(results_df, out_df)
}

true_df <- results_df %>%
  filter(strategy == "truth") %>%
  select(true_V1 = V1, true_V2 = V2, sim, sample)
accuracy_df <- results_df %>%
  filter(type == "full" & strategy != "truth") %>%
  left_join(true_df, by = c("sim", "sample")) %>%
  group_by(strategy, sim) %>%
  summarise(V1_accuracy = accuracy(V1, true_V1),
            V2_accuracy = accuracy(V2, true_V2), .groups = "drop")

accuracy_df %>%
  pivot_longer(c(V1_accuracy, V2_accuracy), names_to = "parameter") %>%
  group_by(parameter, strategy) %>%
  summarise(min = min(value), mean = mean(value), max = max(value))

### PLOT OUTPUT ###

results_df %>% 
  filter(sim == 1) %>%
  ggplot(aes(V1, V2, color = strategy, 
             group = interaction(strategy, type, shard),
             linetype = type)) +
  geom_density2d()
