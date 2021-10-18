##########################
### Sharding functions ###
##########################

##################################################
### Partition observations uniformly at random ###
##################################################
random_shards <- function(y, n_shard, ...) {
  
  # Number of observations
  n <- nrow(y)
  
  # Split observations into shards
  perm <- sample(n, n)
  shard <- split(perm, rep_len(1:n_shard, length(perm)), drop = TRUE)
  
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}
  
  lapply(shard, select_indices, data = y)
  
}

#############################
### Balanced partitioning ###
#############################

balanced_shards <- function(y, n_shard, prior_params, posterior_update,
                            posterior_sampler, M, n_level = 3) {
  
  # Initial sharding
  n_shard_init <- n_shard * (2 ^ n_level)
  n <- nrow(y)
  
  # Split observations into shards
  perm <- sample(n, n)
  shard <- split(perm, rep_len(1:n_shard_init, length(perm)), drop = TRUE)
  
  
  ### Utility functions ###
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}

  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
    rbind(data[[indices[1]]], data[[indices[2]]])
  }
  
  for (l in seq(n_level)) {
    y_sharded <- lapply(shard, select_indices, data = y)
    
    # Draw the parameters according to the groups
    n_shard <- length(shard)
    posterior_params <- lapply(y_sharded, 
                               function(x) posterior_update(x, prior_params, n_shard))
    
    shard_samples <- list()
    for (i in seq_along(y_sharded)) {
      shard_samples[[i]] <- posterior_sampler(M, posterior_params[[i]])
    }
    
    if (is.vector(shard_samples[[1]])) {
      posterior_means <- lapply(shard_samples, mean)
    } else {
      posterior_means <- lapply(shard_samples, rowMeans)
    }
    
    # Create the distance matrix
    D <- pairwise_square_distances(posterior_means)
    # Choose the pairs to merge
    indices_pairs_list = greedy_pairs(D)
    
    # Merge the data
    shard <- lapply(indices_pairs_list, merge, data = shard)
  }
  
  lapply(shard, select_indices, data = y)
}

##############################
### Clustered partitioning ###
##############################

clustered_shards <- function(y, n_shard, ...) {
  
  df = as.data.frame(y)
  km <- kmeans(df, centers = n_shard, nstart = 25)
  
  which_shard = function(i, kmeans){which(kmeans$cluster==i)}
  shard_cluster = lapply(1:n_shard, which_shard, km)
  
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}
  lapply(shard_cluster, select_indices, data = y)
}



##########################
### Sharding functions for LOGIT ###
##########################

##################################################
### Partition observations uniformly at random ###
##################################################
random_shards_logit <- function(y, n_shard, ...) {
  
  # Number of observations
  n <- nrow(y)
  
  # Split observations into shards
  perm <- sample(n, n)
  shard <- split(perm, rep_len(1:n_shard, length(perm)), drop = TRUE)
  
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}
  
  lapply(shard, select_indices, data = y)
  
}

#############################
### Balanced partitioning ###
#############################

balanced_shards_logit <- function(y, n_shard, 
                                  rho_initial, 
                                  iterations, 
                                  initial_d,
                                  M, n_level = 3) {
  
  # Initial sharding
  n_shard_init <- n_shard * (2 ^ n_level)
  n <- nrow(y)
  
  # Split observations into shards
  perm <- sample(n, n)
  shard <- split(perm, rep_len(1:n_shard_init, length(perm)), drop = TRUE)
  
  
  ### Utility functions ###
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}
  
  # Create a merge function that can merge data using lapply
  merge <- function(indices, data){
    rbind(data[[indices[1]]], data[[indices[2]]])
  }
  
  for (l in seq(n_level)) {
    y_sharded <- lapply(shard, select_indices, data = y)
    
    # Draw the parameters according to the groups
    n_shard <- length(shard)

    shard_samples <- list()
    for (i in seq_along(y_sharded)) {
      nc = ncol(y_sharded[[i]])
      yy = y_sharded[[i]][,nc]
      XX = y_sharded[[i]][,1:(nc-1)]
      
      GLM<-glm(yy ~ . -1 ,family=binomial, data=as.data.frame(XX) )
      
      initial_beta = GLM$coefficients
      
      shard_samples[[i]] <- draw_posterior_sample_logit(y_sharded[[i]],
                                                        initial_beta, 
                                                        rho_initial, 
                                                        iterations, 
                                                        initial_d,)
    }
    
    if (is.vector(shard_samples[[1]])) {
      posterior_means <- lapply(shard_samples, mean)
    } else {
      posterior_means <- lapply(shard_samples, rowMeans)
    }
    
    # Create the distance matrix
    D <- pairwise_square_distances(posterior_means)
    # Choose the pairs to merge
    indices_pairs_list = greedy_pairs(D)
    
    # Merge the data
    shard <- lapply(indices_pairs_list, merge, data = shard)
  }
  
  lapply(shard, select_indices, data = y)
}

##############################
### Clustered partitioning ###
##############################

clustered_shards_logit <- function(y, n_shard, ...) {
  
  df = as.data.frame(y)
  km <- kmeans(df, centers = n_shard, nstart = 25)
  
  which_shard = function(i, kmeans){which(kmeans$cluster==i)}
  shard_cluster = lapply(1:n_shard, which_shard, km)
  
  # Create the groups data according to indices created before
  select_indices = function(shard, data) { data[shard,]}
  lapply(shard_cluster, select_indices, data = y)
}