# Bayesian-Hierarchical-Spatiotemporal-
Developing a Bayesian Method for Multi-pathogen outbreak

library(rstan)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

# Use all cores

# ---------- Simulation (or load your simulated arrays) ----------
set.seed(42)

# small test size — increase to full size later
N <- 20        # regions
Tt <- 100      # time
P <- 2         # pathogens
K <- 3         # covariates

Define the Population, enviromental, mobility matrix, and Lambda latent

# Populations
population <- sample(30000:70000, N, replace = TRUE)

# environmental covariates: array N x T x K
E_array <- array(rnorm(N * Tt * K), dim = c(N, Tt, K))

# Mobility matrix (stable gravity) - helper function
compute_stable_mobility <- function(pop, coords, min_dist = 10) {
    Nloc <- length(pop)
    dist_mat <- as.matrix(dist(coords))
    dist_mat[dist_mat < min_dist] <- min_dist
    log_pop <- log(pop)
    log_M <- matrix(-Inf, Nloc, Nloc)
    for (i in 1:Nloc) {
        for (j in 1:Nloc) {
            if (i != j) log_M[i,j] <- log_pop[i] + log_pop[j] - 2*log(dist_mat[i,j])
        }
    }
    log_M <- log_M - max(log_M[is.finite(log_M)])
    M <- exp(log_M)
    row_sums <- rowSums(M)
    row_sums[row_sums == 0] <- 1
    M / row_sums
}
