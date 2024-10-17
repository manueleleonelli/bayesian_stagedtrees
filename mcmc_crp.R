### CRP Function
# alpha: concentration parameter
# assignments: current assignment of vertices to stages
# vertex: the index of the vertex to be reassigned
# crp_assignment <- function(assignments, alpha) {
#   unique_stages <- unique(assignments)
#   n_stages <- length(unique_stages)
#   counts <- table(assignments)
#   
#   # CRP probabilities for existing stages
#   crp_probs <- counts / (sum(counts) + alpha)
#   
#   # Probability of creating a new stage
#   new_stage_prob <- alpha / (sum(counts) + alpha)
#   
#   # Combine probabilities
#   probs <- c(crp_probs, new_stage_prob)
#   
#   # Sample a new stage assignment
#   new_stage_idx <- sample(1:(n_stages + 1), 1, prob = probs)
#   
#   if (new_stage_idx == n_stages + 1) {
#     # Create a new stage
#     new_stage <- paste0("stage_", n_stages + 1)
#   } else {
#     # Assign to an existing stage
#     new_stage <- unique_stages[new_stage_idx]
#   }
#   
#   return(new_stage)
# }

### Function to compute the prior probability of edge probabilities using imaginary sample size
comp_prior <- function(tree, a = 1){
  priors <- tree$stages
  nlev_full <- unlist(lapply(tree$tree, length)) # Number of levels of each variable
  K <-  prod(nlev_full) # Number of root to leaf paths
  paths <- rev(cumprod(rev(nlev_full)))[-1] # Number of root to leaf paths for nodes at different levels
  nlev <- nlev_full[-1] # We drop the first variable since we never use it
  depth <- length(nlev) # Number of layers where we perform model search
  for(i in 1:length(tree$stages)) priors[[i]] <- (a/K)*(paths[i])*sapply(tree$stages[[i]], function(j) sum(tree$stages[[i]] == j))
  priors_filtered_repeated_divided <- mapply(function(prior, n) {
    unique_names <- unique(names(prior))
    filtered_prior <- prior[unique_names]
    repeated_divided <- rep(filtered_prior / n, each = n)
    return(repeated_divided)
  }, priors, nlev, SIMPLIFY = FALSE)
  
  return(priors_filtered_repeated_divided)
}

### MCMC Algorithm with CRP
mcmc_crp <- function(tree, n_save, n_burn = 0, thin = 1, a = 1, kappa = 1, prior = "Friedman", scope = NULL, beta = NULL, tau = NULL, update_SM = TRUE, update_CRP= TRUE){
  
  n_tot <- n_burn + thin * n_save
  
  # Check the inputs (use check_mcmc function as before, omitted here for brevity)
  # priors <- comp_prior(stndnaming(tree), a)
  if(is.null(scope)){scope <- sevt_varnames(tree)[-1]}
  pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
  
  #Objects for CRP update
  n_v <- length(scope)
  alloc_v <- vector("list", length = n_v)
  priors_v <- vector("list", length = n_v)
  n_out_v <- rep(NA, n_v)
  N_v <- sapply(tree$stages[scope], length)
  for(iv in 1:n_v){
    v <- scope[iv]
    n_out_v[iv] <- length(tree$prob[[v]][[1]])
    alloc_v[[iv]] <- tree$stages[[v]]
    # priors_v[[iv]] <- rep(1/2^(N_v[iv]-1)/n_out_v[iv],n_out_v[iv]) 
    priors_v[[iv]] <- rep(1/2,n_out_v[iv]) 
  }
  K_v <- sapply(alloc_v, function(x){length(unique(x))})
  nj_v <- lapply(alloc_v, table)
  
  SM_accept <- 0
  SM_count <- 0
  Merge_count <-  0
  Split_count <-  0
  
  #Output arrays and lists
  alloc_v_out <- vector("list", length = n_save)
  nj_v_out <- vector("list", length = n_save)
  K_v_out <- matrix(NA, n_save, n_v) 
  chain_out <- vector("list", length = n_save)
  
  ## Cycle over iteration
  for(it in 1:n_tot){
    
    # ## Compute prior for current tree
    # priors <- comp_prior(tree, a)
    
    
    
    ## Cycle over variables that are estimated
    for(iv in 1:n_v){
      v <- scope[iv]
      
      #Update allocation labels for each node
      N_iv <- N_v[iv]
      c_ij <- as.numeric(alloc_v[[iv]])
      nj_iv <- nj_v[[iv]]
      K_iv <- K_v[iv]
      
      priors_iv <- priors_v[[iv]]
      
      if(update_CRP){
        for(i in 1:N_iv){
          
          #Allocation of stage i
          aux_i <- c_ij[i]
          
          #Update number of tables serving dish k
          nj_iv[aux_i] <- nj_iv[aux_i] - 1
          #It could have been the only table serving dish k
          alone_i <- (nj_iv[aux_i] == 0)
          
          f_k <- rep(-Inf,K_iv+1)
          
          #Prior predictive for new cluster of one element
          #This is the marginal likelihood from the Dirichlet-Multinomial model
          tree_new <- tree
          tree_new$stages[[v]][i] <- K_iv+1
          
          # ## Counts for new stage
          # priors_new <- comp_prior(tree_new, a)
          
          ix <- (tree_new$stages[[v]] == (K_iv+1))
          # pr_new <- unlist(priors_new[[v]][names(priors_new[[v]]) == (K_iv+1)])
          pr_new <- priors_iv
          tt_new <- tree_new$ctables[[v]][ix, ]
          
          f_k[1] <- lgamma(sum(pr_new)) - sum(lgamma(pr_new)) + sum(lgamma(tt_new + pr_new)) - lgamma(sum(tt_new + pr_new)) 
          
          for(k in 1:K_iv){
            
            #Probability of choosing an existing cluster
            #This is the marginal likelihood from the Dirichlet-Multinomial model
            tree_k <- tree
            tree_k$stages[[v]][i] <- k
            
            # ## Counts for k-th stage
            # priors_k <- comp_prior(tree_k, a)
            
            ix <- (tree_k$stages[[v]] == k)
            # pr_k <- unlist(priors_k[[v]][names(priors_k[[v]]) == k])
            pr_k <- priors_iv
            if (sum(ix) > 1) {
              tt_k <- apply(tree_k$ctables[[v]][ix, ], MARGIN = 2, sum)
            } else {
              tt_k <- tree_k$ctables[[v]][ix, ]
            }
            
            #Now the same but excluding the current stage (priors are unchanged)
            tt_k_minus_i <- tt_k - tree_k$ctables[[v]][i, ]
            
            f_k[k + 1] <- lgamma(sum(pr_k + tt_k_minus_i)) - sum(lgamma(pr_k + tt_k_minus_i)) + sum(lgamma(tt_k + pr_k)) - lgamma(sum(tt_k + pr_k))
          }
          
          w <- c(kappa, nj_iv)
          
          #The nj term could be < 0!
          if(alone_i){
            w[aux_i + 1] = 0
          }
          
          f_k = exp(f_k-max(f_k)) * w
          f_k = f_k/sum(f_k)
          
          hh <- sample.int(K_iv + 1, 1, prob = f_k)
          
          if(hh == 1){
            if(alone_i){
              #Same number of clusters
              nj_iv[aux_i] <- 1
            }else{
              #New singleton cluster
              c_ij[i] <- K_iv + 1
              K_iv <- K_iv + 1
              nj_iv <- c(nj_iv, 1)
            }
          }else{
            #Old cluster
            h_k <- hh - 1
            c_ij[i] <- h_k
            nj_iv[h_k] <- nj_iv[h_k] + 1
            if(alone_i){
              #Remove empty cluster and adjust labels
              K_iv <- K_iv - 1
              c_ij[c_ij > aux_i] <- c_ij[c_ij > aux_i] - 1
              nj_iv <- nj_iv[-aux_i]
            }
          }
          
          tree$stages[[v]] <- as.character(c_ij)
        }
      }
      
      
      if(update_SM){
        ## Sampling of two vertices
        vertices <- sample(length(tree$stages[[v]]), 2, replace = FALSE)
        stage_1 <- as.numeric(tree$stages[[v]][vertices[1]])
        stage_2 <- as.numeric(tree$stages[[v]][vertices[2]])
        
        # priors <-  comp_prior(tree, a)
        
        tree_new <- tree
        log_ar_SM <- 0
        if(stage_1 != stage_2){ # MERGE MOVE
          
          Merge_count <- Merge_count + 1
          K_new <- K_iv - 1
          
          # Get stages for the selected vertices
          merge_stage <- min(stage_1,stage_2)
          other_stage <- max(stage_1,stage_2)
          
          #Merge
          tree_new$stages[[v]][which(tree_new$stages[[v]] == other_stage)] <- merge_stage
          c_new <- as.numeric(tree_new$stages[[v]])
          c_new[c_new > other_stage] <- c_new[c_new > other_stage] - 1 
          tree_new$stages[[v]] <- as.character(c_new)
          nj_new <- table(c_new)
          
          ## Fit new tree
          tree_new <- sevt_fit(tree_new, scope = v)
          # ## Compute edges prior for new tree
          # priors_new <- comp_prior(tree_new, a)
          
          ## Computing counts for the first stage
          ixa <- tree$stages[[v]] == stage_1
          # pr_a <- priors[[v]][names(priors[[v]]) == stage_1]
          pr_a <- priors_iv
          if (sum(ixa) > 1) {
            tt_a <- apply(tree$ctables[[v]][ixa, ], MARGIN = 2, sum) 
          } else {
            tt_a <- tree$ctables[[v]][ixa, ]
          }
          
          ## Computing counts for the second stage
          ixb <- tree$stages[[v]] == stage_2
          # pr_b <- priors[[v]][names(priors[[v]]) == stage_2]
          pr_b <- priors_iv
          if (sum(ixb) > 1) {
            tt_b <- apply(tree$ctables[[v]][ixb, ], MARGIN = 2, sum)
          } else {
            tt_b <- tree$ctables[[v]][ixb, ]
          }
          
          
          ## Computing counts for the second stage
          ix_SM <- tree_new$stages[[v]] == merge_stage
          # pr_SM <- priors_new[[v]][names(priors_new[[v]]) == merge_stage]
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree_new$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree_new$ctables[[v]][ix_SM, ]
          }
          
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(tt_SM+pr_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM + r_12
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM + log(kappa) + lgamma(nj_new[merge_stage]) - lgamma(nj_iv[as.numeric(stage_1)]) - lgamma(nj_iv[as.numeric(stage_2)])
          
          #if(prior == "Heckerman"){lp <- -log(beta[[v]])} else if(prior == "Friedman") {lp <- log(length(tree$stages[[v]]) - length(unique(tree$stages[[v]]))) -log(length(unique(tree$stages[[v]])))} else if(prior == "Pensar"){lp <- -tau*log(nrow(tree$data_raw))*(1+tau)^(length(unique(tree$stages[[v]]))-1)} else {lp <- 1}
          
          ## Compute transition probability ratio
          log_ar_SM <- log_ar_SM + (nj_new[merge_stage]-2) * log(0.5)
          
        } else { ## SPLIT MOVE
          
          Split_count <- Split_count + 1
          K_new <- K_iv + 1
          
          ## Get the name of the current stage
          current_stage <- tree$stages[[v]][vertices[1]]
          new_stage <- K_new
          
          
          ## Get all vertices in the same stage minus those sampled
          target_indices <- setdiff(which(tree$stages[[v]] == current_stage), vertices)
          ## Create new stage
          tree_new$stages[[v]][vertices[2]] <- new_stage
          ## Randomly assign vertices in chosen stage
          assignments <- sample(c(current_stage,new_stage), length(target_indices), replace = TRUE)
          #Split
          tree_new$stages[[v]][target_indices] <- assignments
          
          ## Fit new tree
          tree_new <- sevt_fit(tree_new, scope = v)
          # ## Compute edges prior for new tree
          # priors_new <- comp_prior(tree_new, a)
          c_new <- tree_new$stages[[v]]
          nj_new <- table(c_new)
          
          
          ## Computing counts for the first stage
          ixa <- tree_new$stages[[v]] == new_stage
          # pr_a <- priors_new[[v]][names(priors_new[[v]]) == new_stage]
          pr_a <- priors_iv
          if (sum(ixa) > 1) {
            tt_a <- apply(tree_new$ctables[[v]][ixa, ], MARGIN = 2, sum) 
          } else {
            tt_a <- tree_new$ctables[[v]][ixa, ]
          }
          
          ## Computing counts for the second stage
          ixb <- tree_new$stages[[v]] == current_stage
          # pr_b <- priors_new[[v]][names(priors_new[[v]]) == current_stage]
          pr_b <- priors_iv
          if (sum(ixb) > 1) {
            tt_b <- apply(tree_new$ctables[[v]][ixb, ], MARGIN = 2, sum)
          } else {
            tt_b <- tree_new$ctables[[v]][ixb, ]
          }
          
          ## Computing counts for the second stage
          ix_SM <- tree$stages[[v]] == current_stage
          # pr_SM <- priors[[v]][names(priors[[v]]) == current_stage]
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree$ctables[[v]][ix_SM, ]
          }
          
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(pr_SM+tt_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM - r_12
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM - log(kappa) + lgamma(nj_new[current_stage]) - lgamma(nj_iv[current_stage]) - lgamma(nj_iv[new_stage])
          
          #if(prior == "Heckerman"){lp <- -log(beta[[v]])} else if(prior == "Friedman") {lp <- log(length(tree$stages[[v]]) - length(unique(tree$stages[[v]]))) -log(length(unique(tree$stages[[v]])))} else if(prior == "Pensar"){lp <- -tau*log(nrow(tree$data_raw))*(1+tau)^(length(unique(tree$stages[[v]]))-1)} else {lp <- 1}
          
          ## Compute transition probability ratio
          log_ar_SM <- log_ar_SM - (nj_iv[current_stage]-2) * log(0.5)
        }
        
        accept_SM <- 1
        if( is.finite(log_ar_SM) ){
          if(log_ar_SM < 0){
            accept_SM <- exp(log_ar_SM)
          }
        }else{
          accept_SM <- 0
        }
        
        SM_accept <- SM_accept + accept_SM
        SM_count <- SM_count + 1
        
        ## Decide on MCMC move
        if( runif(1) < accept_SM ){
          K_iv <- K_new
          nj_iv <- nj_new
          c_ij <- c_new
          tree <- tree_new
          # priors <- priors_new
        }
      }
      
      alloc_v[[iv]] <- c_ij
      nj_v[[iv]] <- nj_iv
      K_v[iv] <- K_iv
    }
    
    if(it%%10 == 0){
      print(K_v)  
    }
    
    if((it > n_burn) & ((it - n_burn)%%thin == 0)){
      iter <- (it - n_burn)%/%thin
      
      alloc_v_out[[iter]] <- alloc_v
      nj_v_out[[iter]] <- nj_v[[iv]]
      K_v_out[iter,] <- K_v 
      chain_out[[iter]] <- stndnaming(sevt_fit(tree, scope = v)) 
    }
    
    setTxtProgressBar(pb, it)
  }
  
  close(pb)
  
  if(update_SM){
    print("merge ar")
    print(Merge_count / n_tot)
    print("split ar")
    print(Split_count / n_tot)
  }
  
  OUTPUT_MCMC <- list("alloc_v_out" = alloc_v_out, "nj_v_out" = nj_v_out, "K_v_out" = K_v_out, "chain_out" = chain_out)
  return(OUTPUT_MCMC)
}

# ### Example Usage
# # Assuming you have a tree object named 'my_tree' and the appropriate setup.
# set.seed(123)  # For reproducibility
# result_chain <- mcmc_crp(tree = tree, it = 50000, a = 0.1, alpha = 100000, scope = scope)
# 
# # Print the first result of the MCMC chain
# print(result_chain[[1]])
