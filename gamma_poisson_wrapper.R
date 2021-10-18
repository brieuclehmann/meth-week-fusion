source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

library(dplyr)
library(xtable)
library(ggplot2)

extract_output <- function(x, name) {
  M <- nrow(x$approx_samples)
  p <- ncol(x$approx_samples)
  n_shard <- length(x$shard_samples)
  approx_df <- as_tibble(x$approx_samples) %>%
    mutate(type = "full", shard = "1", strategy = name, sample = seq(M))
  shard_df <- plyr::ldply(x$shard_samples, as_tibble) %>%
    mutate(type = "sub_posterior", strategy = name, 
           sample = rep(seq(M), n_shard)) %>%
    rename(shard = .id, V1 = value)
  
  bind_rows(approx_df, shard_df)
}

set.seed(1)

gamma_post_update <- function(y, prior_params, prior_scale) {
  
  n <- length(y)
  # scale the prior to ensure the right combination
  post_shape <- sum(y) + prior_params['shape']
  post_rate  <- n + prior_params['rate']
  
  c(post_shape, post_rate)
}

# Support function to draw posterior samples using lapply
gamma_post_sampler = function(M, params) {
  rgamma(M, shape = params['shape'], rate = params['rate'])
}

# create the data
prior_params <- c(shape = 10, rate  = 1)
N <- 512
n_sim <- 20
M <- 500
K <- 4
results_df <- tibble()
for (sim in 1:n_sim) {
  set.seed(sim)
  
  lambda <- rgamma(1, shape = prior_params['shape'], rate = prior_params['rate'])
  y      <- rpois(N, lambda)
  
  out_random <- run_experiment(y, K, random_shards, prior_params, 
                               gamma_post_update, gamma_post_sampler, 
                               consensus_combine_wrapper, M)
  
  out_balanced <- run_experiment(y, K, balanced_shards, prior_params, 
                                 gamma_post_update, gamma_post_sampler,
                                 consensus_combine_wrapper, M)
  
  out_clustered <- run_experiment(y, K, clustered_shards, prior_params, 
                                  gamma_post_update, gamma_post_sampler,
                                  consensus_combine_wrapper, M)
  
  true_posterior_params <- gamma_post_update(y, prior_params, 1)
  true_post_samples <- gamma_post_sampler(M, true_posterior_params)
  truth_df <- as_tibble(true_post_samples) %>%
    mutate(type = "full", shard = "1", strategy = "truth", sample = seq(M)) %>%
    rename(V1 = 1)
  
  out_df <- bind_rows(extract_output(out_random, "random"),
                      extract_output(out_balanced, "HeMP"), 
                      extract_output(out_clustered, "clustered"),
                      truth_df) %>%
    mutate(sim = sim)
  
  results_df <- bind_rows(results_df, out_df)
}


true_df <- results_df %>%
  filter(strategy == "truth") %>%
  select(true_V1 = V1, sim, sample)
accuracy_df <- results_df %>%
  filter(type == "full" & !strategy %in% c("truth", "clustered")) %>%
  left_join(true_df, by = c("sim", "sample")) %>%
  group_by(strategy, sim) %>%
  summarise(V1_accuracy = accuracy(V1, true_V1), .groups = "drop")

table_df <- accuracy_df %>%
  pivot_longer(V1_accuracy, names_to = "parameter") %>%
  group_by(parameter, strategy) %>%
  summarise(min = min(value), mean = mean(value), max = max(value))

print(xtable(table_df))

results_df %>% 
  filter(sim == 1 & strategy != "clustered") %>%
  ggplot(aes(V1, color = strategy, 
             group = interaction(strategy, type, shard),
             linetype = type)) +
  geom_density() +
  theme_minimal() +
  xlab("theta")
