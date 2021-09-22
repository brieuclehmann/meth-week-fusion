# Methodology week
rm(list=ls())
set.seed(123)

N = 128
K = 32

dimension = 2

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
n_grps = K
grps = split(permutation, rep_len(1:n_grps, length(permutation)))


# create the data
prior_mean_true = matrix(0, 2, 1)
prior_cov_true  = matrix(c(1,0,0,1), 2, 2)

obs_cov_true = matrix(c(2,0.5,0.5,1), 2, 2)

library(MASS)
mu   = mvrnorm(1, prior_mean_true, prior_cov_true)
Y    = mvrnorm(N, mu, obs_cov_true)

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
  mu_Y    = mvrnorm(posterior_sample_size, posterior_param$mean, posterior_param$cov)
  
  mu_Y 
}

# Support function to draw posterior parameters using lapply
draw_posterior_parameters_mean = function(Y_sub, obs_cov, prior_mean, prior_cov){
  
  posterior_param = normal_mean_posterior(Y_sub, obs_cov, prior_mean, prior_cov)
  
  posterior_param$mean
}

# Create the groups data according to indices created before
select_indices = function(grps, data){ data[grps,]}
Y_grouped = lapply(grps, select_indices, data = Y)

# mu_Y_list = lapply(Y_grouped, draw_posterior_sample)

source("pairing.R")

Levels = 3

for (l in 1:Levels) {
  # Draw the parameters according to the groups
  mu_Y_param_list = lapply(Y_grouped, draw_posterior_parameters_mean, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true)

  # Create the distance matriY
  D = pairwise_square_distances(mu_Y_param_list)
  # Choose the pairs to merge
  indices_pairs_list = greedy_pairs(D)

  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
          (rbind(data[[indices[1]]], data[[indices[2]]]))
  }

  # Merge the data
  Y_grouped <- lapply(indices_pairs_list, merge, data = Y_grouped)
  
}

post_sample_size = 100
mu_Y_sample_list = lapply(Y_grouped, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = post_sample_size)

precision = function(x){solve(var(x))}
precision_list = lapply(mu_Y_sample_list, precision)

shard_samples = matrix(unlist(mu_Y_sample_list), nrow =post_sample_size, byrow=FALSE)

shard_precisions <- array(NA, c(dimension,dimension,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

consensus_samples <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(mu_Y_sample_list, function(x) x[i,]), nrow = dimension)
  consensus_samples[i,] <- consensus_combine(shard_samples, shard_precisions)
}


##############################################
# Random
N = 128
K_random = 4

# create a permutation over the indices
permutation = sample(N,N)

batch_size_random=N/K_random

# divide in K groups and add indices per each group
n_grps_random = K_random
grps_random = split(permutation, rep_len(1:n_grps_random, length(permutation)))

# Create the groups data according to indices created before
Y_grouped_random = lapply(grps_random, select_indices, data = Y)

# Plot
mu_Y_sample_list_random = lapply(Y_grouped_random, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = post_sample_size)

library(plyr)
out_df = ldply(mu_Y_sample_list_random)
colnames(out_df) = c("shard", "dim_1", "dim_2")

library(ggplot2)
ggplot(data = out_df)
ggplot(data=out_df, aes(x=dim_1, y=dim_2, group = shard)) + geom_density()


precision_list_random = lapply(mu_Y_sample_list_random, precision)

shard_samples_random = matrix(unlist(mu_Y_sample_list_random), nrow =post_sample_size, byrow=FALSE)

shard_precisions_random <- array(NA, c(dimension,dimension,length(precision_list_random)))
for (shard in seq_along(precision_list_random)) {
  shard_precisions_random[ , ,shard] <- precision_list_random[[shard]]
}

consensus_samples_random <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(mu_Y_sample_list_random, function(x) x[i,]), nrow = dimension)
  consensus_samples_random[i,] <- consensus_combine(shard_samples, shard_precisions_random)
}

out_df = as.data.frame(consensus_samples_random)
ggplot2(out_df, aes(x,y, group = shard)) + geom_density()


