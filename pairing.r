points = c(1, 3, 6, 10, 15) # example 'data set'

pairwise_square_distances <- function(v){ # example 'distance matrix'
  n <- length(v)
  M <- matrix(-1, nrow = n, ncol = n) 
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      M[i, j] <- (v[i] - v[j])^2
    }
  }
  return(M)
}

greedy_pairs <- function(D){ # takes distance matrix, outputs greedy 'anti-clustering'
  # check D is square
  n = dim(D)[1]
  m = floor(n/2)
  # initialise list
  # distant_pairs = matrix(NA, m, 2)
  distant_pairs = list()
  for (k in 1:m){
    # find indices corresponding to largest remaining distances
    indices = which(D == max(D), arr.ind = TRUE)
    chosen_pair = indices[1,]
    # append indices to list
    # distant_pairs[k,] <- chosen_pair
    distant_pairs[[k]] <- chosen_pair
    # hide distances involving chosen elements
    D[indices[1,1],] = -1
    D[indices[1,2],] = -1
    D[,indices[1,1]] = -1
    D[,indices[1,2]] = -1
  }
  # return list of indices
  if (2*m != n){
    leftover = which(D == max(D), arr.ind = TRUE)
    return(append(distant_pairs, leftover)) # should have some warning?
  } else{
    return(distant_pairs)
  }
}

# to work out: is this representation of the pairs appropriate?
