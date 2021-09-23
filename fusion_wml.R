# Methodology week
rm(list=ls())
set.seed(124)

N = 512
K = 32

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
K_random <- 4
n_grps = K
grps_init <- grps <- split(permutation, rep_len(1:n_grps, length(permutation)))
grps_random <- lapply(seq(K_random), function(x) unlist(grps_init[seq(K) %% K_random == x - 1]))

# create the data
prior_mean_true = 0
prior_cov_true  = 1

obs_cov_true = 1

mu   = rnorm(1, prior_mean_true, prior_cov_true)
Y    = rnorm(N, mu, obs_cov_true)

#########################################################
# Deprecated
#normal_mean_posterior_uni = function(obs_data, obs_cov, prior_mean, prior_cov ){
#  n = length(obs_data)
#  post_mean = (prior_cov/(obs_cov/n + prior_cov))*mean(obs_data) + (obs_cov/(obs_cov/n + prior_cov))*prior_mean 
#  post_cov  = 1/(1/(prior_cov) + n/(obs_cov))
  
#  list(mean = post_mean, cov = post_cov)
#  }

#posterior_param = normal_mean_posterior_uni(Y[grps$`1`], 1, 0, 1)
#########################################################

source("utils.R")
#posterior_param = normal_mean_posterior(Y[grps$`1`], 1, 0, 1)
#posterior_sample_size = 100

# Support function to draw posterior samples using lapply
draw_posterior_sample = function(Y_sub, obs_cov, prior_mean, prior_cov, posterior_sample_size){
  
  posterior_param = normal_mean_posterior(Y_sub, obs_cov, prior_mean, prior_cov)
  mu_Y    = rnorm(posterior_sample_size, posterior_param$mean, sqrt(posterior_param$cov))
  
  mu_Y 
}

# Support function to draw posterior parameters using lapply
draw_posterior_parameters_mean = function(Y_sub, obs_cov, prior_mean, prior_cov){
  
  posterior_param = normal_mean_posterior(Y_sub, obs_cov, prior_mean, prior_cov)
  
  posterior_param$mean
}

# Create the groups data according to indices created before
select_indices = function(grps, data){ data[grps]}
Y_grouped = lapply(grps, select_indices, data = Y)

# mu_Y_list = lapply(Y_grouped, draw_posterior_sample)

source("pairing.R")

Levels = 3
p = 100
for (l in 1:Levels) {
  # Draw the parameters according to the groups
  n_shards <- length(Y_grouped)
  mu_Y_param_list = lapply(Y_grouped, draw_posterior_parameters_mean, 
                           obs_cov = obs_cov_true, 
                           prior_mean = prior_mean_true, 
                           prior_cov = prior_cov_true * n_shards)

  # Create the distance matriY
  D = pairwise_square_distances(unlist(mu_Y_param_list))
  # Choose the pairs to merge
  indices_pairs_list = greedy_pairs(D)

  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
          (c(data[[indices[1]]], data[[indices[2]]]))
  }

  # Merge the data
  Y_grouped <- lapply(indices_pairs_list, merge, data = Y_grouped)
  
}

mu_Y_sample_list = lapply(Y_grouped, draw_posterior_sample, 
                          obs_cov = obs_cov_true, 
                          prior_mean = prior_mean_true, 
                          prior_cov = prior_cov_true * n_shards, 
                          posterior_sample_size = p)

precision = function(x){solve(var(x))}
precision_list = lapply(mu_Y_sample_list, precision)
shard_precisions <- array(NA, c(1,1,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

consensus_samples <- double(p) 
for (i in seq(p)) {
  shard_samples <- matrix(sapply(mu_Y_sample_list, function(x) x[i]), nrow = 1)
  consensus_samples[i] <- consensus_combine(shard_samples, shard_precisions)
}




true_posterior <- normal_mean_posterior(Y, obs_cov_true, 
                                        prior_mean_true, prior_cov_true)
xgrid <- seq(min(consensus_samples), max(consensus_samples), length.out = 1000)
truth <- dnorm(xgrid, true_posterior$mean, true_posterior$cov)
plot(xgrid, truth, col = "red", type = "l")
lines(density(consensus_samples))

##############################################
# Random
N = 128
K = 4

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
n_grps = K
grps = split(permutation, rep_len(1:n_grps, length(permutation)))

# Create the groups data according to indices created before
select_indices = function(grps){ Y[grps]}
Y_grouped_random = lapply(grps_random, select_indices)

lapply(Y_grouped_random, mean)

mu_Y_random_sample_list = lapply(Y_grouped_random, draw_posterior_sample, 
                                 obs_cov = obs_cov_true, 
                                 prior_mean = prior_mean_true, 
                                 prior_cov = prior_cov_true *n_shards, 
                                 posterior_sample_size = p)

precision = function(x){solve(var(x))}
precision_random_list = lapply(mu_Y_random_sample_list, precision)
shard_random_precisions <- array(NA, c(1,1,length(precision_random_list)))
for (shard in seq_along(precision_random_list)) {
  shard_random_precisions[ , ,shard] <- precision_random_list[[shard]]
}

consensus_random_samples <- double(p) 
for (i in seq(p)) {
  shard_random_samples <- matrix(sapply(mu_Y_random_sample_list, function(x) x[i]), nrow = 1)
  consensus_random_samples[i] <- consensus_combine(shard_random_samples, shard_random_precisions)
}

lines(density(consensus_random_samples), col = "blue")

true_samples <- rnorm(p, true_posterior$mean, sqrt(true_posterior$cov))

names(mu_Y_sample_list) <- seq_along(mu_Y_sample_list)
names(mu_Y_random_sample_list) <- seq_along(mu_Y_random_sample_list)
truth_df <- tibble(x = true_samples), 
                   shard = '1', strategy = 'truth', type = 'full') %>%
  bind_rows(tibble(x = consensus_samples, shard = '1', 
                   strategy = 'merged', type = 'full')) %>%
  bind_rows(tibble(x = consensus_random_samples, shard = '1', 
                   strategy = 'random', type = 'full')) 

merged_df <- as_tibble(mu_Y_sample_list) %>%
  pivot_longer(seq(length(mu_Y_sample_list)), names_to = "shard", values_to = "x") %>%
  mutate(strategy = 'merged', type = 'subposterior')

random_df <- as_tibble(mu_Y_random_sample_list) %>%
  pivot_longer(seq(length(mu_Y_random_sample_list)), names_to = "shard", values_to = "x") %>%
  mutate(strategy = 'random', type = 'subposterior')


out_df <- bind_rows(merged_df, random_df, truth_df)
ggplot(out_df, aes(x, group = interaction(shard, type, strategy), linetype = type, color = strategy)) +
  geom_density() +
  theme_minimal()
