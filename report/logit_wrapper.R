source("experiment_wrapper.R")
source("sharding_funcs.R")
source("combine_funcs.R")
source("utils.R")
source("pairing.r")

library(tidyr)
library(dplyr)
library(readr)
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
    rename(shard = .id)
  
  bind_rows(approx_df, shard_df)
}

set.seed(11)

N = 512
K = 4

d=2

initial_d=diag(d)
iterations = 20000
rho_initial=1

n_sim <- 10
M <- 10000
results_df <- tibble()

for (sim in 1:n_sim) {
  # set.seed(sim)
  X=matrix(rnorm(N*d),N,d) #continuous covariates
  X[,1]=rep(1,N)
  
  beta = seq(-1, 1, length.out = d)
  y = rbinom(N,1, 1/(1+exp(-X %*% beta)) )
  
  
  Xy = cbind(X, y)
  
  out_random <- run_experiment_logit(Xy, K, 
                                     random_shards_logit, 
                                     rho_initial, 
                                     iterations, 
                                     initial_d,
                                     consensus_combine_wrapper, M)
  
  out_balanced <- run_experiment_logit(Xy, K, 
                                       balanced_shards_logit,
                                       rho_initial, 
                                       iterations, 
                                       initial_d,
                                       consensus_combine_wrapper, M)
  
  nc = ncol(Xy)
  yy = Xy[,nc]
  XX = Xy[,1:(nc-1)]
  
  GLM<-glm(yy ~ . -1 ,family=binomial, data=as.data.frame(XX) )
  
  initial_beta = GLM$coefficients
  
  true_post_samples <- draw_posterior_sample_logit(Xy, 
                                                    initial_beta, 
                                                    rho_initial, 
                                                    iterations, 
                                                    initial_d, )[seq((iterations-M+1),iterations,1),]
  
  truth_df <- as_tibble(true_post_samples) %>%
    mutate(type = "full", shard = "1", strategy = "truth", sample = seq(M))
  
  out_df <- bind_rows(extract_output(out_random, "random"),
                      extract_output(out_balanced, "HeMP"), 
                      truth_df) %>%
    mutate(sim = sim)
  
  results_df <- bind_rows(results_df, out_df)
}


true_df <- results_df %>%
  filter(strategy == "truth") %>%
  select(true_V1 = V1, true_V2 = V2, sim, sample)
accuracy_df <- results_df %>%
  filter(sim != 3) %>%
  filter(type == "full" & strategy != "truth") %>%
  left_join(true_df, by = c("sim", "sample")) %>%
  group_by(strategy, sim) %>%
  summarise(V1_accuracy = accuracy(V1, true_V1),
            V2_accuracy = accuracy(V2, true_V2), .groups = "drop")

table_df <- accuracy_df %>%
  pivot_longer(c(V1_accuracy, V2_accuracy), names_to = "parameter") %>%
  group_by(parameter, strategy) %>%
  summarise(min = min(value), mean = mean(value), max = max(value))

print(xtable(table_df))
### PLOT OUTPUT ###

results_df %>% 
  filter(sim == 1 & strategy != "clustered") %>%
  ggplot(aes(V1, V2, color = strategy, 
             group = interaction(strategy, type, shard),
             linetype = type)) +
  geom_density2d() +
  theme_minimal() +
  xlab("beta1") + ylab("beta2")
