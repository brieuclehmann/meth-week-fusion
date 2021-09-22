# Methodology week
rm(list=ls())
set.seed(123)

N = 128
K = 32

# create a permutation over the indices
permutation = sample(N,N)

batch_size=N/K

# divide in K groups and add indices per each group
n_grps = K
grps = split(permutation, rep_len(1:n_grps, length(permutation)))


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
  mu_Y    = rnorm(posterior_sample_size, posterior_param$mean, posterior_param$cov)
  
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

for (l in 1:Levels) {
  # Draw the parameters according to the groups
  mu_Y_param_list = lapply(Y_grouped, draw_posterior_parameters_mean, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true)

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

p = 100
mu_Y_sample_list = lapply(Y_grouped, draw_posterior_sample, obs_cov = obs_cov_true, prior_mean = prior_mean_true, prior_cov = prior_cov_true, posterior_sample_size = p)

precision = function(x){solve(var(x))}
precision_list = lapply(mu_Y_sample_list, precision)

shard_samples = matrix(unlist(mu_Y_sample_list), nrow =p, byrow=FALSE)

lapply(Y_grouped, mean)






























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
Y_grouped_random = lapply(grps, select_indices)

lapply(Y_grouped_random, mean)
