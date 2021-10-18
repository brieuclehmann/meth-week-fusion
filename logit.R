rm(list=ls())

source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")
library(ggplot2)
library(dplyr)

library(MASS)

set.seed(10)

dimension =2

N = 512
K = 32

d=2

X=matrix(rnorm(N*d),N,d) #continuous covariates
X[,1]=rep(1,N)

beta = seq(-1, 1, length.out = d)
y = rbinom(N,1, 1/(1+exp(-X %*% beta)) )

GLM<-glm(y ~ . -1 ,family=binomial, data=as.data.frame(X) )

#to get reference posterior
rwLogit=function(y,X,xi0,sigma,run,A, S=1){
  ptm=proc.time()
  xi=xi0
  d=length(xi0)
  U=function(q){as.numeric(sum(log(1+exp(X %*% q)))- t(y) %*% (X %*% q) + (t(q) %*% q)/(2*S))}
  tAch=t(chol(A))
  out=matrix(0,run,d)
  count=0
  for(n in 1:run){
    xinew <- xi + sigma * tAch %*% rnorm(d)
    if(runif(1)<exp(U(xi)-U(xinew))){
      xi=xinew
      count=count+1
    }
    if(n%%(run/10)==0) cat("Iteration n.",n,"\n") #track progress1
    out[n,]=xi
  }
  time=proc.time()-ptm
  cat('Run Time',round(time[1]/60,digits=2),'Min','\n')
  cat("Acc Rate",count/run,"\n")
  return(list(out,time[1], count/run))
}

draw_posterior_sample = function(Xy, initial_beta, rho_initial, iterations, initial_d, S=1){
  nc = ncol(Xy)
  y = Xy[,nc]
  X = Xy[,1:(nc-1)]
  acc = 0
  AA  = initial_d
  beta = initial_beta
  rho = rho_initial
  
  while(acc<0.2 || acc>0.3){
    resRW=rwLogit(y, X, initial_beta, rho, iterations,  AA)
    acc =resRW[[3]]
    AA=cov(resRW[[1]])
    
    if (acc>0.3){rho = 1.2*rho}
    else{rho = 0.8*rho}
  }
  
  resRW[[1]]
}

initial_d_value=diag(d)
iterations_value = 50000
rho_initial_value=1
initial_beta_value = GLM$coefficients

Xy = cbind(X, y)

beta_sample_full = draw_posterior_sample(Xy, initial_beta_value, rho_initial_value, iterations_value, initial_d_value)

# divide in K groups and add indices per each group
n_shard = K
# create a permutation over the indices
permutation = sample(N,N)
shard = split(permutation, rep_len(1:n_shard, length(permutation)), drop =T)

# Create the groups data according to indices created before
select_indices = function(shard, data){ 
  
    if(is.numeric(dim(data))){data[shard,]}
    else {data[shard]}
}
  
Xy_grouped = lapply(shard, select_indices, data = Xy)

source("pairing.R")

Levels = 3
initial_d_value=diag(d)
iterations_value = 10000
rho_initial_value=1.5
initial_beta_value = GLM$coefficients

for (l in 1:Levels) {
  # Draw the parameters according to the groups
  n_shards <- length(Xy_grouped)
  
  mu_Y_sample_list = lapply(Xy_grouped, draw_posterior_sample, 
                           initial_beta = initial_beta_value, 
                           rho_initial = rho_initial_value, 
                           iterations = iterations_value, 
                           initial_d = initial_d_value)
  
  mu_Y_param_list = lapply(mu_Y_sample_list, colMeans)
  
  # Create the distance matriY
  D = pairwise_square_distances(mu_Y_param_list)
  # Choose the pairs to merge
  indices_pairs_list = greedy_pairs(D)
  
  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
    rbind(data[[indices[1]]], data[[indices[2]]])
  }
  
  # Merge the data
  Xy_grouped <- lapply(indices_pairs_list, merge, data = Xy_grouped)
  
} 


post_sample_size = 100
beta_sample_list = lapply(Xy_grouped, draw_posterior_sample, 
                          initial_beta = initial_beta_value, 
                          rho_initial = rho_initial_value, 
                          iterations = iterations_value, 
                          initial_d = initial_d_value, S= 4)

names(beta_sample_list) = sapply(1:length(beta_sample_list), as.character)

library(plyr)
out_df = ldply(beta_sample_list)
colnames(out_df) = c("shard", "dim_1", "dim_2")
out_df$type <- "balanced shards"

library(ggplot2)
plot_homog = ggplot(data=out_df, aes(x=dim_1, y=dim_2, group = shard)) + 
  geom_density2d(size=0.5, bins =5, color = "green2")


precision = function(x){solve(var(x))}
precision_list = lapply(beta_sample_list, precision)
shard_precisions <- array(NA, c(2,2,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

shard_samples = matrix(unlist(beta_sample_list), nrow =post_sample_size, byrow=FALSE)

shard_precisions <- array(NA, c(2,2,length(precision_list)))
for (shard in seq_along(precision_list)) {
  shard_precisions[ , ,shard] <- precision_list[[shard]]
}

consensus_samples <- matrix(NA, post_sample_size, 2) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(beta_sample_list, function(x) x[i,]), nrow = 2)
  consensus_samples[i,] <- consensus_combine(shard_samples, shard_precisions)
}

out_df_consensus = as.data.frame(consensus_samples)
colnames(out_df_consensus) = c("dim_1", "dim_2")
out_df_consensus$type <- "balanced shards consensus"
plot_homog = plot_homog +
  geom_density2d(data=out_df_consensus, 
                 aes(dim_1,dim_2, color = type),
                 size=1, bins =5)

out_df_truth = as.data.frame(beta_sample_full)
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
Xy_grouped_random = lapply(shard_random, select_indices, data = Xy)

# Plot
beta_sample_list_random = lapply(Xy_grouped_random, draw_posterior_sample, 
                                 initial_beta = initial_beta_value, 
                                 rho_initial = rho_initial_value, 
                                 iterations = iterations_value, 
                                 initial_d = initial_d_value, S= 4)
library(plyr)
out_df = ldply(beta_sample_list_random)
colnames(out_df) = c("random", "dim_1", "dim_2")
out_df$type <- "random shards"

library(ggplot2)
plot_random = ggplot(data=out_df, 
                     aes(x=dim_1, y=dim_2, group = random)) + 
  geom_density2d(size=0.5, bins =5, aes(color = type, alpha =0.01))


precision_list_random = lapply(beta_sample_list_random, precision)

shard_samples_random = matrix(unlist(beta_sample_list_random), nrow =post_sample_size, byrow=FALSE)

shard_precisions_random <- array(NA, c(dimension,dimension,length(precision_list_random)))
for (shard in seq_along(precision_list_random)) {
  shard_precisions_random[ , ,shard] <- precision_list_random[[shard]]
}

consensus_samples_random <- matrix(NA, post_sample_size, dimension) 
for (i in seq(post_sample_size)) {
  shard_samples <- matrix(sapply(beta_sample_list_random, function(x) x[i,]), nrow = dimension)
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

dev.new()
is plot_homog
#dev.new()
#plot_random

  
