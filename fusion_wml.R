# Methodology week
rm(list=ls())
set.seed(124)

N = 512
K = 32

dimension = 2

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
n_shard = K
shard = split(permutation, rep_len(1:n_shard, length(permutation)), drop =T)

# create the data
prior_mean_true = matrix(0, 2, 1)
prior_cov_true  = matrix(c(1,0,0,1), 2, 2)

obs_cov_true = matrix(c(1,0.5,0.5,2), 2, 2)

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

#posterior_param = normal_mean_posterior_uni(Y[shard$`1`], 1, 0, 1)
#########################################################

source("utils.R")
#posterior_param = normal_mean_posterior(Y[shard$`1`], 1, 0, 1)
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

# draw posterior samples from the whole thing
post_sample_size = 500
mu_Y_sample_full = draw_posterior_sample(Y, obs_cov_true, prior_mean_true, prior_cov_true, post_sample_size)
  


# Create the groups data according to indices created before
select_indices = function(shard, data){ data[shard,]}
Y_grouped = lapply(shard, select_indices, data = Y)

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
  D = pairwise_square_distances(mu_Y_param_list)
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
mu_Y_sample_list = lapply(Y_grouped, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = post_sample_size)

names(mu_Y_sample_list) = sapply(1:length(mu_Y_sample_list), as.character)

library(plyr)
out_df = ldply(mu_Y_sample_list)
colnames(out_df) = c("shard", "dim_1", "dim_2")
out_df$type <- "balanced shards"

library(ggplot2)
plot_homog = ggplot(data=out_df, aes(x=dim_1, y=dim_2, group = shard)) + 
  geom_density2d(size=0.5, bins =5, aes(color = type))+
  xlim(-1, 2)+
  ylim(-2.5, 3)

precision = function(x){solve(var(x))}
precision_list = lapply(mu_Y_sample_list, precision)
shard_precisions <- array(NA, c(1,1,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

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

out_df_consensus = as.data.frame(consensus_samples)
colnames(out_df_consensus) = c("dim_1", "dim_2")
out_df_consensus$type <- "balanced shards consensus"
plot_homog = plot_homog +
  geom_density2d(data=out_df_consensus, 
                 aes(dim_1,dim_2, color = type),
                 size=1, bins =5)

out_df_truth = as.data.frame(mu_Y_sample_full)
colnames(out_df_truth) = c("dim_1", "dim_2")
out_df_truth$type <- "truth"

plot_homog = plot_homog +
  geom_density2d(data=out_df_truth, 
                 aes(dim_1,dim_2, color = type),
                 size=1, bins =5)


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
mu_Y_sample_list_random = lapply(Y_grouped_random, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = post_sample_size)

library(plyr)
out_df = ldply(mu_Y_sample_list_random)
colnames(out_df) = c("random", "dim_1", "dim_2")
out_df$type <- "random shards"

library(ggplot2)
plot_random = ggplot(data=out_df, 
                     aes(x=dim_1, y=dim_2, group = random)) + 
  geom_density2d(size=0.5, bins =5, aes(color = type))+
  xlim(-1, 2)+
  ylim(-2.5, 3)


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

out_df_consensus_random = as.data.frame(consensus_samples_random)
out_df_consensus_random$random <- NA 
colnames(out_df_consensus_random) = c("dim_1", "dim_2")
out_df_consensus_random$type = "random consensus"
plot_random = plot_random + 
  geom_density2d(data=out_df_consensus_random, 
                 aes(dim_1,dim_2, color = type), inherit.aes = FALSE,
                 size=1, bins =5)

plot_random = plot_random +
  geom_density2d(data=out_df_truth, 
                 aes(dim_1,dim_2, color = type), inherit.aes = FALSE,
                 size=1, bins =5)

#dev.new()
#plot_homog
#dev.new()
#plot_random


##############################################
# Clustering
K_cluster = 4
df = as.data.frame(Y)
km <- kmeans(df, centers = K_cluster, nstart = 25)

which_shard = function(i, kmeans){which(kmeans$cluster==i)}
shard_cluster = lapply(1:K_cluster, which_shard, km)

# Create the groups data according to indices created before
Y_grouped_cluster = lapply(shard_cluster, select_indices, data = Y)

# Plot
mu_Y_sample_list_cluster = lapply(Y_grouped_cluster, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = post_sample_size)

names(mu_Y_sample_list_cluster) = sapply(1:length(mu_Y_sample_list_cluster), as.character)


library(plyr)
out_df = ldply(mu_Y_sample_list_cluster)
colnames(out_df) = c("cluster", "dim_1", "dim_2")
out_df$type <- "clustered shards"


library(ggplot2)
plot_cluster = ggplot(data=out_df, aes(x=dim_1, y=dim_2, group = cluster, col = cluster)) + 
  geom_density2d(size=0.5, aes(colour = type)) +
  xlim(-1, 2)+
  ylim(-2.5, 3)

precision_list_cluster = lapply(mu_Y_sample_list_cluster, precision)

shard_samples_cluster = matrix(unlist(mu_Y_sample_list_cluster), nrow =post_sample_size, byrow=FALSE)

shard_precisions_cluster <- array(NA, c(dimension,dimension,length(precision_list_cluster)))
for (shard in seq_along(precision_list_cluster)) {
  shard_precisions_cluster[ , ,shard] <- precision_list_cluster[[shard]]
}

consensus_samples_cluster <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(mu_Y_sample_list_cluster, function(x) x[i,]), nrow = dimension)
  consensus_samples_cluster[i,] <- consensus_combine(shard_samples, shard_precisions_cluster)
}

out_df_consensus_cluster = as.data.frame(consensus_samples_cluster)
colnames(out_df_consensus_cluster) = c("dim_1", "dim_2")
out_df_consensus_cluster$type <- "clustered consensus"
plot_cluster = plot_cluster + 
  geom_density2d(data=out_df_consensus_cluster, 
                 aes(dim_1,dim_2, colour = type), inherit.aes = FALSE,
                 size=1, bins =5)

out_df_truth = as.data.frame(mu_Y_sample_full)
colnames(out_df_truth) = c("dim_1", "dim_2")
out_df_truth$type <- "truth"
plot_cluster = plot_cluster +
  geom_density2d(data=out_df_truth, 
                 aes(dim_1,dim_2, colour = type), inherit.aes = FALSE,
                 size=1, bins =5)


dev.new()
plot_homog
dev.new()
plot_random
dev.new()
plot_cluster
