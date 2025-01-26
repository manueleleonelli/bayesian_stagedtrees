sample_partition <- function(c_ij, kappa_iv, csi_iv, lamb_iv, iterations = 100) {
  # c_ij: Initial cluster assignments
  # kappa_iv: Weight for creating a new cluster
  # csi_iv: Scaling parameter for distances
  # lamb_iv: Vector of weights for distance covariates
  # dist_H_iv: Pairwise distance matrix
  # dist_nums_iv: List of additional distance matrices
  # iterations: Number of iterations to run (default = 100)
  
  N_iv <- length(c_ij)  # Number of observations
  K_iv <- max(c_ij)     # Initial number of clusters
  nj_iv <- table(c_ij)  # Cluster sizes
  
  # Convert to numeric for easier manipulation
  nj_iv <- as.numeric(nj_iv)
  
  # Repeat the Gibbs sampling process for the specified number of iterations
  for (iter in 1:iterations) {
    for (i in 1:N_iv) {
      # Allocation of stage i
      aux_i <- c_ij[i]
      
      # Update number of tables serving dish k
      nj_iv[aux_i] <- nj_iv[aux_i] - 1
      
      # It could have been the only table serving dish k
      alone_i <- (nj_iv[aux_i] == 0)
      
      # Initialize probabilities for each cluster
      f_k <- rep(1, K_iv + 1)
      dist_H_j <- rep(0, K_iv)
      dist_nums_j <- vector("list", length = length(dist_nums_iv))
      
      # Calculate distances for existing clusters
      for (k in 1:K_iv) {
        ind_k <- setdiff(which(c_ij == k), i)
        dist_H_j[k] <- sum(dist_H_iv[i, ind_k])
        for (ind in 1:length(dist_nums_iv)) {
          dist_nums_j[[ind]][k] <- sum(dist_nums_iv[[ind]][i, ind_k])
        }
      }
      
      # Compute weights for clusters
      w <- c(kappa_iv, nj_iv * exp(-csi_iv * dist_H_j - Reduce(`+`, Map(function(x, y) x * y, dist_nums_j, lamb_iv))))
      
      # If observation i is alone in its cluster, remove the option to reassign
      if (alone_i) {
        w[aux_i + 1] <- 0
      }
      
      # Normalize weights to form probabilities
      f_k <- w / sum(w)
      
      # Sample a new cluster assignment
      hh <- sample.int(K_iv + 1, 1, prob = f_k)
      
      if (hh == 1) {
        if (alone_i) {
          # Keep the same number of clusters
          nj_iv[aux_i] <- 1
        } else {
          # Create a new singleton cluster
          c_ij[i] <- K_iv + 1
          K_iv <- K_iv + 1
          nj_iv <- c(nj_iv, 1)
        }
      } else {
        # Assign to an existing cluster
        h_k <- hh - 1
        c_ij[i] <- h_k
        nj_iv[h_k] <- nj_iv[h_k] + 1
        if (alone_i) {
          # Remove empty cluster and adjust labels
          K_iv <- K_iv - 1
          c_ij[c_ij > aux_i] <- c_ij[c_ij > aux_i] - 1
          nj_iv <- nj_iv[-aux_i]
        }
      }
    }
  }
  
  return(c_ij)
}
