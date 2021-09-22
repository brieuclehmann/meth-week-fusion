# Methodology week
set.seed(123)

N = 128
K = 16

permutation = sample(N,N)

batch_size=N/K

n_grps = K
grps = split(permutation, rep_len(1:n_grps, length(permutation)))


prior_mean = 0
prior_cov  = 1

mu   = rnorm(1, prior_mean, prior_cov)
X    = rnorm(N, mu, 1)

#normal_mean_posterior_uni = function(obs_data, obs_cov, prior_mean, prior_cov ){
#  n = length(obs_data)
#  post_mean = (prior_cov/(obs_cov/n + prior_cov))*mean(obs_data) + (obs_cov/(obs_cov/n + prior_cov))*prior_mean 
#  post_cov  = 1/(1/(prior_cov) + n/(obs_cov))
  
#  list(mean = post_mean, cov = post_cov)
#  }

#posterior_param = normal_mean_posterior_uni(X[grps$`1`], 1, 0, 1)

source("utils.R")
posterior_param = normal_mean_posterior(X[grps$`1`], 1, 0, 1)
posterior_sample_size = 100

draw_posterior_sample = function(X_sub, posterior_sample_size=100){
  posterior_param = normal_mean_posterior(X_sub, 1, 0, 1)
  mu_X    = rnorm(posterior_sample_size, posterior_param$mean, posterior_param$cov)
  
  mu_X 
}

draw_posterior_parameters_mean = function(X_sub){
  posterior_param = normal_mean_posterior(X_sub, 1, 0, 1)
  
  posterior_param$mean
}

select_indices = function(grps){ X[grps]}

X_grouped = lapply(grps, select_indices)

mu_X_list = lapply(X_grouped, draw_posterior_sample)

mu_X_param_list = lapply(X_grouped, draw_posterior_parameters_mean)

source("pairing.R")

M = pairwise_square_distances(unlist(mu_X_param_list))
