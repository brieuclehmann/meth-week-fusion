source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

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

N = 256
dimension = 1

# create the data
prior_params <- list(shape = 2, rate  = 1)

gamma_post_update <- function(y, prior_params, prior_scale) {
  
  n <- length(y)
  # scale the prior to ensure the right combination
  post_shape <- sum(obs_data) + prior_params$shape
  post_rate  <- n + prior_params$scale
  
  list(shape = post_shape, rate = post_rate)
}

# Support function to draw posterior samples using lapply
gamma_post_sampler = function(M, params) {
  rgamma(M, shape = params$shape, rate = params$rate)
}

# Support function to draw posterior parameters using lapply
poisson_posterior_parameters = function(obs_data, prior_shape, prior_rate){
  
  posterior_param = gamma_posterior_parameters(obs_data, prior_shape, prior_rate)
  
  c(posterior_param$shape, posterior_param$rate)
}


lambda = rgamma(1, shape = prior_shape_true, rate = prior_rate_true)
Y      = rpois(N, lambda)

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
n_shard = K
shard = split(permutation, rep_len(1:n_shard, length(permutation)), drop =T)


# Create the groups data according to indices created before
select_indices = function(shard, data){ data[shard]}
Y_grouped = lapply(shard, select_indices, data = Y)

source("utils_poisson.R")

# draw posterior samples from the whole thing
post_sample_size = 1000
lambda_Y_sample_full = poisson_draw_posterior_sample(obs_data = Y, prior_shape = prior_shape_true, prior_rate = prior_rate_true, post_sample_size)


source("pairing.r")
Levels = 3

for (l in 1:Levels) {
  # Draw the parameters according to the groups
  theta_Y_param_list = lapply(Y_grouped, poisson_posterior_parameters, prior_shape = prior_shape_true, prior_rate = prior_rate_true)
  
  # Create the distance matriY
  D = pairwise_square_distances(theta_Y_param_list)
  # Choose the pairs to merge
  indices_pairs_list = greedy_pairs(D)
  
  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
    rbind(data[[indices[1]]], data[[indices[2]]])
  }
  
  # Merge the data
  Y_grouped <- lapply(indices_pairs_list, merge, data = Y_grouped)
  
}

post_sample_size = 100
theta_Y_sample_list = lapply(Y_grouped, poisson_draw_posterior_sample, prior_shape=prior_shape_true, prior_rate= prior_rate_true, posterior_sample_size= post_sample_size)

names(theta_Y_sample_list) = sapply(1:length(theta_Y_sample_list), as.character)

library(plyr)
out_df_theta_shard = ldply(theta_Y_sample_list, data.frame)
colnames(out_df_theta_shard) = c("shard", "lambda")
out_df_theta_shard$type <- "balanced shards"

library(ggplot2)
plot_balanced = ggplot(out_df_theta_shard, aes(lambda, group = shard, color = type))+
      geom_density() +
      theme_minimal()+
  xlim(0.4, 1.2)+
  ylim(0, 8) 


precision = function(x){solve(var(x))}
precision_list = lapply(theta_Y_sample_list, precision)

shard_samples = matrix(unlist(theta_Y_sample_list), nrow =post_sample_size, byrow=FALSE)

shard_precisions <- array(NA, c(dimension,dimension,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

consensus_samples <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(theta_Y_sample_list, function(x) x[i]), nrow = dimension)
  consensus_samples[i,] <- consensus_combine(shard_samples, shard_precisions)
}

out_df_consensus = as.data.frame(consensus_samples)
colnames(out_df_consensus) = c("lambda")
out_df_consensus$type <- "balanced shards consensus"

plot_balanced =  plot_balanced+
              geom_density(data=out_df_consensus, aes(lambda, color = type)) +
              theme_minimal() 


out_df_truth = as.data.frame(lambda_Y_sample_full)
colnames(out_df_truth) = c("lambda")
out_df_truth$type <- "truth"

plot_balanced =  plot_balanced+
  geom_density(data=out_df_truth, aes(lambda, color = type)) +
  theme_minimal() 


##############################################
# Random
K_random = 4

# create a permutation over the indices
permutation = sample(N,N)

batch_size_random=N/K_random

# divide in K groups and add indices per each group
n_shard_random = K_random
shard_random = split(permutation, rep_len(1:n_shard_random, length(permutation)))

# Create the groups data according to indices created before
Y_grouped_random = lapply(shard_random, select_indices, data = Y)

# Plot
theta_Y_sample_list_random = lapply(Y_grouped_random, poisson_draw_posterior_sample, prior_shape=prior_shape_true, prior_rate= prior_rate_true, posterior_sample_size= post_sample_size)

names(theta_Y_sample_list_random) = sapply(1:length(theta_Y_sample_list_random), as.character)

library(plyr)
out_df_theta_random = ldply(theta_Y_sample_list_random, data.frame)
colnames(out_df_theta_random) = c("random", "lambda")
out_df_theta_random$type <- "random shards"

library(ggplot2)
plot_random = ggplot(out_df_theta_random, aes(lambda, group = random, color = type))+
  geom_density() +
  theme_minimal()+
  xlim(0.4, 1.2)+
  ylim(0, 8)


precision = function(x){solve(var(x))}
precision_list = lapply(theta_Y_sample_list, precision)

shard_samples = matrix(unlist(theta_Y_sample_list_random), nrow =post_sample_size, byrow=FALSE)

shard_precisions <- array(NA, c(dimension,dimension,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

consensus_samples <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(theta_Y_sample_list_random, function(x) x[i]), nrow = dimension)
  consensus_samples[i,] <- consensus_combine(shard_samples, shard_precisions)
}

out_df_consensus_random = as.data.frame(consensus_samples)
colnames(out_df_consensus_random) = c("lambda")
out_df_consensus_random$type <- "random shards consensus"

plot_random =  plot_random+
  geom_density(data=out_df_consensus_random, 
               aes(lambda, color = type), inherit.aes = FALSE,) +
  theme_minimal() 


out_df_truth = as.data.frame(lambda_Y_sample_full)
colnames(out_df_truth) = c("lambda")
out_df_truth$type <- "truth"

plot_random =  plot_random+
  geom_density(data=out_df_truth, 
               aes(lambda, color = type), inherit.aes = FALSE,) +
  theme_minimal()

plot_random
dev.new()
plot_balanced
