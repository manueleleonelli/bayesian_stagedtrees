
predloglik <- function(alpha0,beta0,kappa0, mean, var, size){
  sum(alpha0*log(beta0) + 0.5*log(kappa0) - lgamma(alpha0) - 0.5*size*log(2*pi) + lgamma(alpha0 + size/2) - 0.5*log(kappa0+size) - (alpha0 + size/2)*log(beta0+0.5*var+ (kappa0*size*(mean-mu0)^2)/(2*(kappa0+size))))
}

### MCMC Algorithm with CRP
mcmc_crp_ppmx <- function(tree, data, n_save, n_burn = 0, thin = 1, a = 1, kappa = 1, csi = 0.25, scope = NULL, update_SM = TRUE, update_CRP= TRUE, alpha0 = 1, beta0 = 1, mu0 = 0, kappa0 = 1){
  
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
    priors_v[[iv]] <- rep(a/n_out_v[iv],n_out_v[iv]) 
  }
  K_v <- sapply(alloc_v, function(x){length(unique(x))})
  nj_v <- lapply(alloc_v, table)
  
  #Compute Hamming distances between stages paths
  #these are fixed because the tree is fixed
  n_vars <- length(tree$tree)
  n_cat <- apply(matrix(1:n_vars,ncol = 1), 1, function(i){length(tree$prob[[i]][[1]])})
  n_paths <- prod(n_cat[-n_vars])
  tree_paths <- matrix(0, n_paths, n_vars)
  for(iv in 1:(n_vars-1)){
    n_cat_iv <- n_cat[iv]
    if(iv < n_vars - 1){
      exp_prod <- prod(n_cat[(iv+1):(n_vars-1)])
    }else{
      exp_prod <- 1
    }
    rep_cat <- c(apply(matrix(1:n_cat_iv, ncol = 1), 1, function(i){rep(i,exp_prod)}))
    
    #Use recycling
    tree_paths[,iv+1] <- rep_cat
  }
  
  #Now remove duplicates
  tree_paths_list <- vector("list",length = n_vars)
  for(iv in 1:n_vars){
    tree_paths_list[[iv]] <- unique(tree_paths[,1:iv])
  }
  
  
  dist_H <- vector("list", length = n_v)
  for(iv in 1:n_v){
    iv_pos <- match(scope[iv], names(tree$tree))
    N_v_iv <- N_v[iv]
    mat <- matrix(0, N_v_iv, N_v_iv)
    for(i1 in 1:(N_v_iv-1)){
      for(i2 in (i1+1):N_v_iv){
        p_v1 <- tree_paths_list[[iv_pos]][i1,]
        p_v2 <- tree_paths_list[[iv_pos]][i2,]
        mat[i1,i2] <- hamming.distance(p_v1, p_v2)
        mat[i2,i1] <- mat[i1,i2]
      }
    }
    #Use normalised distance: divide by the number of variables (taking values in (0,1))
    dist_H[[iv]] <- mat/(length(p_v1)-1)
  }
  
  # Initialize result objects as nested lists
  means <- list()
  vars <- list()
  sizes <- list()
  
  # Iterate through each target variable in the scope
  for (target_var in scope) {
    # Find the variables preceding the current target variable
    vars_order <- names(tree$tree)
    preceding_vars <- vars_order[1:(which(vars_order == target_var) - 1)]
    
    # Get all possible combinations of levels for the preceding variables
    levels_list <- tree$tree[rev(preceding_vars)]
    combinations <- expand.grid(levels_list)
    
    # Process all numerical variables
    numeric_vars <- names(data)[sapply(data, is.numeric)]
    
    # Initialize nested lists for this target variable
    means[[target_var]] <- list()
    vars[[target_var]] <- list()
    sizes[[target_var]] <- list()
    
    # Iterate over continuous (numeric) variables
    for (num_var in numeric_vars) {
      # Initialize storage for this numeric variable
      means[[target_var]][[num_var]] <- c()
      vars[[target_var]][[num_var]] <- c()
      sizes[[target_var]][[num_var]] <- c()
      
      # Iterate through combinations of preceding variables
      for (i in 1:nrow(combinations)) {
        # Subset the data based on the combination
        subset_condition <- TRUE
        for (var in preceding_vars) {
          subset_condition <- subset_condition & (data[[var]] == combinations[[var]][i])
        }
        subset_data <- data[subset_condition, ]
        
        # Compute statistics
        mean_value <- mean(subset_data[[num_var]], na.rm = TRUE)
        sum_squared_diff <- sum((subset_data[[num_var]] - mean_value)^2, na.rm = TRUE)
        sample_size <- nrow(subset_data)
        
        # Store results for this numeric variable and combination
        means[[target_var]][[num_var]] <- c(means[[target_var]][[num_var]], mean_value)
        vars[[target_var]][[num_var]] <- c(vars[[target_var]][[num_var]], sum_squared_diff)
        sizes[[target_var]][[num_var]] <- c(sizes[[target_var]][[num_var]], sample_size)
      }
    }
  }
  nums <- sapply(data, function(col) is.numeric(col) | is.integer(col))
  n_nums <- sum(nums)
  
  
  margs <- vector("list", length = n_v)
  
  for(iv in 1:n_v){
    margs[[iv]] <- rep(-Inf,length(alloc_v[[iv]]))
    sizes_iv <- sizes[[iv]]
    means_iv <- means[[iv]]
    vars_iv <- vars[[iv]]
    for(k in 1:length(margs[[iv]])){
      sizes_ivi <- sapply(sizes_iv, function(x) x[k])
      means_ivi <- sapply(means_iv, function(x) x[k])
      vars_ivi <- sapply(vars_iv, function(x) x[k])
      margs[[iv]][k] <- exp(predloglik(alpha0,beta0,kappa0,means_ivi,vars_ivi,sizes_ivi))   
      }
  }
    
    
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
    
    ## Cycle over variables that are estimated
    for(iv in 1:n_v){
      v <- scope[iv]
      
      #Update allocation labels for each node
      N_iv <- N_v[iv]
      c_ij <- as.numeric(alloc_v[[iv]])
      nj_iv <- nj_v[[iv]]
      K_iv <- K_v[iv]
      
      priors_iv <- priors_v[[iv]]
      
      dist_H_iv <- dist_H[[iv]]
      
      means_iv <- means[[iv]]
      vars_iv <- vars[[iv]]
      sizes_iv <- sizes[[iv]]
      
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
          
          ix <- (tree_new$stages[[v]] == (K_iv+1))
          pr_new <- priors_iv
          tt_new <- tree_new$ctables[[v]][ix, ]
          
          f_k[1] <- lgamma(sum(pr_new)) - sum(lgamma(pr_new)) + sum(lgamma(tt_new + pr_new)) - lgamma(sum(tt_new + pr_new)) 
          
          g_k <- rep(-Inf,K_iv+1)
          g_k[1] <- margs[[iv]][i]
          
          for(k in 1:K_iv){
            
            #Probability of choosing an existing cluster
            #This is the marginal likelihood from the Dirichlet-Multinomial model
            tree_k <- tree
            tree_k$stages[[v]][i] <- k
            
            ix <- (tree_k$stages[[v]] == k)
            pr_k <- priors_iv
            if (sum(ix) > 1) {
              tt_k <- apply(tree_k$ctables[[v]][ix, ], MARGIN = 2, sum)
            } else {
              tt_k <- tree_k$ctables[[v]][ix, ]
            }
            
            #Now the same but excluding the current stage (priors are unchanged)
            tt_k_minus_i <- tt_k - tree_k$ctables[[v]][i, ]
            
            f_k[k + 1] <- lgamma(sum(pr_k + tt_k_minus_i)) - sum(lgamma(pr_k + tt_k_minus_i)) + sum(lgamma(tt_k + pr_k)) - lgamma(sum(tt_k + pr_k))
            
            if(alone_i & aux_i == k){
              g_k[k+1] <- g_k[1]
            } else{
              sizes_new <- colSums(sapply(sizes_iv, function(x) x[ix]))
              means_new <- colSums(sapply(sizes_iv, function(x) x[ix]) * sapply(means_iv, function(x) x[ix])) / sizes_new
              vars_new <- colSums(sapply(vars_iv, function(x) x[ix])) + colSums(sapply(sizes_iv, function(x) x[ix]) * (sapply(means_iv, function(x) x[ix]) - means_new)^2)
              ixi <- ix
              ixi[i] <- FALSE
              sizes_raw <- sapply(sizes_iv, function(x) x[ixi])
              means_raw <- sapply(means_iv, function(x) x[ixi])
              vars_raw <- sapply(vars_iv, function(x) x[ixi])
              
              # Handle single-row case by converting to a matrix
              if (is.null(dim(sizes_raw))) {
                sizes_raw <- matrix(sizes_raw, nrow = 1)
              }
              if (is.null(dim(means_raw))) {
                means_raw <- matrix(means_raw, nrow = 1)
              }
              if (is.null(dim(vars_raw))) {
                vars_raw <- matrix(vars_raw, nrow = 1)
              }
              
              # Compute total sizes
              sizes_old <- colSums(sizes_raw)
              
              # Compute combined means
              means_old <- colSums(sizes_raw * means_raw) / sizes_old
              
              # Compute combined variances
              vars_old <- colSums(vars_raw) + colSums(sizes_raw * (means_raw - means_old)^2)
              
              g_k[k+1] <- exp(predloglik(alpha0,beta0,kappa0,means_new,vars_new,sizes_new) - predloglik(alpha0,beta0,kappa0,means_old,vars_old,sizes_old))
            }
          
         }
          
          dist_H_j <- rep(0,K_iv)
          for(k in 1:K_iv){
            ind_k <- setdiff(which(c_ij == k),i)
            dist_H_j[k] <- sum(dist_H_iv[i,ind_k])
          }
          
          w <- c(kappa, nj_iv * exp(-csi * dist_H_j))
          
            #The nj term could be < 0!
          if(alone_i){
            w[aux_i + 1] = 0
          }
          
          f_k = exp(f_k-max(f_k)) * w*g_k
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
        
          ## Computing counts for the first stage
          ixa <- tree$stages[[v]] == stage_1
          # pr_a <- priors[[v]][names(priors[[v]]) == stage_1]
          pr_a <- priors_iv
          if (sum(ixa) > 1) {
            tt_a <- apply(tree$ctables[[v]][ixa, ], MARGIN = 2, sum) 
          } else {
            tt_a <- tree$ctables[[v]][ixa, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_a <- sum(dist_H_iv[ixa,ixa])/2
          
          ## Computing counts for the second stage
          ixb <- tree$stages[[v]] == stage_2
          # pr_b <- priors[[v]][names(priors[[v]]) == stage_2]
          pr_b <- priors_iv
          if (sum(ixb) > 1) {
            tt_b <- apply(tree$ctables[[v]][ixb, ], MARGIN = 2, sum)
          } else {
            tt_b <- tree$ctables[[v]][ixb, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_b <- sum(dist_H_iv[ixb,ixb])/2
          
          ## Computing counts for the second stage
          ix_SM <- tree_new$stages[[v]] == merge_stage
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree_new$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree_new$ctables[[v]][ix_SM, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_SM <- sum(dist_H_iv[ix_SM,ix_SM])/2
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(tt_SM+pr_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM + r_12 - r_1 - r_2 ## fixed 
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM + log(kappa) + lgamma(nj_new[merge_stage]) - lgamma(nj_iv[as.numeric(stage_1)]) - lgamma(nj_iv[as.numeric(stage_2)])
          #Add Hamming distance part
          log_ar_SM <- log_ar_SM  - csi*(dist_SM - dist_a - dist_b)
          
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
          
          c_new <- tree_new$stages[[v]]
          nj_new <- table(c_new)
          
          
          ## Computing counts for the first stage
          ixa <- tree_new$stages[[v]] == new_stage
          pr_a <- priors_iv
          if (sum(ixa) > 1) {
            tt_a <- apply(tree_new$ctables[[v]][ixa, ], MARGIN = 2, sum) 
          } else {
            tt_a <- tree_new$ctables[[v]][ixa, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_a <- sum(dist_H_iv[ixa,ixa])/2
          
          ## Computing counts for the second stage
          ixb <- tree_new$stages[[v]] == current_stage
          pr_b <- priors_iv
          if (sum(ixb) > 1) {
            tt_b <- apply(tree_new$ctables[[v]][ixb, ], MARGIN = 2, sum)
          } else {
            tt_b <- tree_new$ctables[[v]][ixb, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_b <- sum(dist_H_iv[ixb,ixb])/2
          
          ## Computing counts for the second stage
          ix_SM <- (tree$stages[[v]] == current_stage)
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree$ctables[[v]][ix_SM, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_SM <- sum(dist_H_iv[ix_SM,ix_SM])/2
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(pr_SM+tt_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM - r_12 + r_1 + r_2 ## fixed
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM - log(kappa) + lgamma(nj_new[current_stage]) - lgamma(nj_iv[current_stage]) - lgamma(nj_iv[new_stage])
          #Add Hamming distance part
          log_ar_SM <- log_ar_SM  - csi*(- dist_SM + dist_a + dist_b)
          
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
      chain_out[[iter]] <- stndnaming(sevt_fit(tree)) 
    }
    
    setTxtProgressBar(pb, it)
  }
  
  close(pb)
  
  if(update_SM){
    print("merge ar")
    print(Merge_count / n_tot / n_v)
    print("split ar")
    print(Split_count / n_tot / n_v)
  }
  
  OUTPUT_MCMC <- list("alloc_v_out" = alloc_v_out, "nj_v_out" = nj_v_out, "K_v_out" = K_v_out, "chain_out" = chain_out)
  return(OUTPUT_MCMC)
}
