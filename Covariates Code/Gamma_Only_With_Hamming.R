

mcmc_crp_gamma_hamming <- function(tree, data, n_save, n_burn = 0, thin = 1, a = 1, kappa = 1, csi = 0.25 , scope = NULL, update_SM = TRUE, update_CRP= TRUE){
  
  
  n_tot <- n_burn + thin * n_save
  
  # Check the inputs (use check_mcmc function as before, omitted here for brevity)
  # priors <- comp_prior(stndnaming(tree), a)
  if(is.null(scope)){scope <- sevt_varnames(tree)[-1]}
  pb <- txtProgressBar(min = 0, max = n_tot, style = 3)
  
  #Objects for CRP update
  nums <- sapply(data, function(col) is.numeric(col) | is.integer(col))
  n_nums <- sum(nums)
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
  
  
  #WE NEED to extend the above for non-binary trees
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
  
  dist_nums <- vector("list", length = n_v)
  for(iv in 1:n_v){
    dist_nums[[iv]] <- vector("list", length = n_nums)
    iv_pos <- match(scope[iv], names(tree$tree))
    columns_to_filter <- names(tree$tree)[1:(iv_pos-1)]
    N_v_iv <- N_v[iv]
    means <- vector("list", length=N_v_iv)
    for(ind in 1:N_v_iv){
      p_v <- tree_paths_list[[iv_pos]][ind,]
      filter_conditions <- sapply(1:length(columns_to_filter), function(i) {
        as.numeric(data[[columns_to_filter[i]]]) == p_v[i + 1]
      })
      subset_data <- data[apply(filter_conditions, 1, all), ]
      means[[ind]] <- apply(subset_data[nums],2,mean, na.rm = TRUE)
    }
    for(ind in 1:n_nums){
      mat <- matrix(0, N_v_iv, N_v_iv)
      for(i1 in 1:(N_v_iv-1)){
        for(i2 in (i1+1):N_v_iv){
          mat[i1,i2] <- abs(means[[i1]][ind]-means[[i2]][ind])
          mat[i2,i1] <- mat[i1,i2]
        }
      }
      dist_nums[[iv]][[ind]] <- mat
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
  lambda_out <- vector("list", length = n_save)
  prop_var_out <- vector("list", length = n_save)
  
  ### 
 lambda_current <- lapply(1:n_v, function(x) rgamma(n_nums, 0.4^2/0.01, scale = 0.01/0.4))
  
  prop_var <- replicate(n_v, rep(0.04, n_nums), simplify = FALSE)
  count_var <- replicate(n_v, rep(0, n_nums), simplify = FALSE)
  count_gamma <- replicate(n_v, rep(0, n_nums), simplify = FALSE)
  num_checks <- replicate(n_v, rep(0, n_nums), simplify = FALSE)
  adapt_factor <- replicate(n_v, rep(0, n_nums), simplify = FALSE)
  
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
      csi_iv <- csi[iv]
      kappa_iv <- kappa[iv]
      
      dist_H_iv <- dist_H[[iv]]
      dist_nums_iv <- dist_nums[[iv]]
      
      lambda_current_iv <- lambda_current[[iv]]
      lamb_iv <- lambda_current_iv
      
      bb_iv <- bb[[iv]]
      cc_iv <- cc[[iv]]
      ff_iv <- ff[[iv]]
      gg_iv <- gg[[iv]]
      
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
          
          dist_H_j <- rep(0,K_iv)
          dist_nums_j <- vector("list",length= n_nums)
          
          for(k in 1:K_iv){
            ind_k <- setdiff(which(c_ij == k),i)
            dist_H_j[k] <- sum(dist_H_iv[i,ind_k])
            for(ind in 1:n_nums){
              dist_nums_j[[ind]][k] <- sum(dist_nums_iv[[ind]][i,ind_k])
            }
          }
          
          
          w <- c(kappa_iv, nj_iv * exp(-csi_iv * dist_H_j - Reduce(`+`, Map(function(x, y) x * y, dist_nums_j, lamb_iv)) ))
          
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
          # tree_new <- sevt_fit(tree_new)
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
          #Sum of pairwise Hamming distances between elements of cluster
          dist_a <- sum(dist_H_iv[ixa,ixa])/2
          dist_a_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_a_cont[ind] <- sum(dist_nums_iv[[ind]][ixa,ixa])/2
          
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
          dist_b_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_b_cont[ind] <- sum(dist_nums_iv[[ind]][ixb,ixb])/2
          
          
          ## Computing counts for the second stage
          ix_SM <- tree_new$stages[[v]] == merge_stage
          # pr_SM <- priors_new[[v]][names(priors_new[[v]]) == merge_stage]
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree_new$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree_new$ctables[[v]][ix_SM, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_SM <- sum(dist_H_iv[ix_SM,ix_SM])/2
          dist_SM_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_SM_cont[ind] <- sum(dist_nums_iv[[ind]][ix_SM,ix_SM])/2
          
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(tt_SM+pr_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM + r_12 - r_1 - r_2 ## fixed 
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM + log(kappa_iv) + lgamma(nj_new[merge_stage]) - lgamma(nj_iv[as.numeric(stage_1)]) - lgamma(nj_iv[as.numeric(stage_2)])
          #Add Hamming distance part
          log_ar_SM <- log_ar_SM  - csi_iv*(dist_SM - dist_a - dist_b) - sum(lamb_iv*(dist_SM_cont-dist_a_cont-dist_b_cont))
          
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
          # tree_new <- sevt_fit(tree_new)
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
          #Sum of pairwise Hamming distances between elements of cluster
          dist_a <- sum(dist_H_iv[ixa,ixa])/2
          dist_a_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_a_cont[ind] <- sum(dist_nums_iv[[ind]][ixa,ixa])/2
          
          ## Computing counts for the second stage
          ixb <- tree_new$stages[[v]] == current_stage
          # pr_b <- priors_new[[v]][names(priors_new[[v]]) == current_stage]
          pr_b <- priors_iv
          if (sum(ixb) > 1) {
            tt_b <- apply(tree_new$ctables[[v]][ixb, ], MARGIN = 2, sum)
          } else {
            tt_b <- tree_new$ctables[[v]][ixb, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_b <- sum(dist_H_iv[ixb,ixb])/2
          dist_b_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_b_cont[ind] <- sum(dist_nums_iv[[ind]][ixb,ixb])/2
          
          ## Computing counts for the second stage
          ix_SM <- (tree$stages[[v]] == current_stage)
          # pr_SM <- priors[[v]][names(priors[[v]]) == current_stage]
          pr_SM <- priors_iv
          if (sum(ix_SM) > 1) {
            tt_SM <- apply(tree$ctables[[v]][ix_SM, ], MARGIN = 2, sum)
          } else {
            tt_SM <- tree$ctables[[v]][ix_SM, ]
          }
          #Sum of pairwise Hamming distances between elements of cluster
          dist_SM <- sum(dist_H_iv[ix_SM,ix_SM])/2
          dist_SM_cont <- rep(0,n_nums)
          for(ind in 1:n_nums) dist_SM_cont[ind] <- sum(dist_nums_iv[[ind]][ix_SM,ix_SM])/2
          
          ## Compute terms for marginal likelihood ratio
          r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
          r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
          r_12 <- +lgamma(sum(pr_SM)) - lgamma(sum(pr_SM+tt_SM)) - sum(lgamma(pr_SM)) + sum(lgamma(tt_SM+pr_SM))
          log_ar_SM <- log_ar_SM - r_12 + r_1 + r_2 ## fixed
          
          ## Compute Tree prior ratio (Dirichlet process)
          log_ar_SM <- log_ar_SM - log(kappa_iv) + lgamma(nj_new[current_stage]) - lgamma(nj_iv[current_stage]) - lgamma(nj_iv[new_stage])
          #Add Hamming distance part
          log_ar_SM <- log_ar_SM  - csi_iv*(- dist_SM + dist_a + dist_b) - sum(lamb_iv*(- dist_SM_cont + dist_a_cont + dist_b_cont))
          
          
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
      
      
      if(update_CONT){
        
        ## CREO LE VARIABILI TEMPORANEE
        
        lambda_proposal_iv <- lambda_current_iv
        
        ## CICLO SULLE COVARIATE CONTINUE
        
        for(i in 1:n_nums){
          
          # PROPONGO GAMMA E LAMBDA
          
          
            
            lambda_proposal_iv[i] <- rgamma(1, shape = lambda_current_iv[i]^2/prop_var[[iv]][i], scale= prop_var[[iv]][i]/lambda_current_iv[i])
            while(lambda_proposal_iv[i] <= 0.0000000000000001) {lambda_proposal_iv[i] <- rgamma(1, shape = lambda_current_iv[i]^2/prop_var[[iv]][i], scale= prop_var[[iv]][i]/lambda_current_iv[i])}
        
          
          # MINIBATCH PARTIZIONE AUSILIARIA
          rho_prime <- sample_partition(c_ij, kappa_iv, csi_iv, lambda_proposal_iv[i], iterations = 100) 
          
          ## Quantita MH prior lambda
          num_lambda <- dgamma(lambda_proposal_iv[i], shape = 0.4^2/0.01, scale = 0.01/0.4,log =T)
          #(ff_iv[i]-1)*log(lambda_proposal_iv[i]) - gg_iv[i]*lambda_proposal_iv[i]
          den_lambda <- dgamma(lambda_current_iv[i], shape = 0.4^2/0.01, scale = 0.01/0.4,log =T)
          #(ff_iv[i]-1)*log(lambda_current_iv[i]) - gg_iv[i]*lambda_current_iv[i]
          
          # Quantita MH dist partizione current
          num_part1 <- length(unique(c_ij))
          num_part2 <- length(unique(rho_prime))
          
          num_dist1 <- 0
          den_dist1 <- 0
          for(k in 1:num_part1){
            ind_1 <- which(c_ij == k)
            mat1 <- dist_nums_iv[[i]][ind_1,ind_1]
            num_dist1 <- num_dist1 - lambda_proposal_iv[i]*sum(mat1)/2
            den_dist1 <- den_dist1 - lambda_current_iv[i]*sum(mat1)/2
          }
          
          # Quantita MH dist partizione ausiliaria
          num_dist2 <- 0
          den_dist2 <- 0
          for(k in 1:num_part2){
            ind_2 <- which(rho_prime == k )
            mat_2 <- dist_nums_iv[[i]][ind_2,ind_2]
            num_dist2 <- num_dist2 -lambda_current_iv[i]*sum(mat_2)/2
            den_dist2 <- den_dist2 - lambda_proposal_iv[i]*sum(mat_2)/2
          }
          
          # Transition probabilities
          
          num_trans <- dgamma(lambda_current_iv[i], shape = lambda_proposal_iv[i]^2/prop_var[[iv]][i], scale = prop_var[[iv]][i]/lambda_proposal_iv[i], log = T)
          den_trans <- dgamma(lambda_proposal_iv[i], shape = lambda_current_iv[i]^2/prop_var[[iv]][i], scale = prop_var[[iv]][i]/lambda_current_iv[i], log = T)
          
          prob <- exp(num_lambda+num_trans+num_dist1+num_dist2- den_lambda-den_trans-den_dist1-den_dist2)
          if( runif(1) < prob ){
            
            lambda_current_iv <- lambda_proposal_iv
            count_var[[iv]][i] <- count_var[[iv]][i]+1
          } else{
           
            lambda_proposal_iv <- lambda_current_iv
          }
          
        }
        lambda_current[[iv]] <- lambda_current_iv
        
        
         
        ## ADAPTIVE VARIANCE
        
        
        if (it%% 50 ==0) {
          for (idx in 1:n_nums) {
            if ( 0.01 <= 1/sqrt(num_checks[[iv]][idx]+1)) {adapt_factor[[iv]][idx] = 0.01;} else{adapt_factor[[iv]][idx] = 1/sqrt(num_checks[[iv]][idx]+1);}
            if(count_var[[iv]][idx] <= 22){prop_var[[iv]][idx] <- exp(log(prop_var[[iv]][idx]) - adapt_factor[[iv]][idx]); } else {prop_var[[iv]][idx] <- exp(log(prop_var[[iv]][idx]) + adapt_factor[[iv]][idx] );}
            num_checks[[iv]][idx] <- num_checks[[iv]][idx] + 1
            count_var[[iv]][idx] <- 0
          }
        }
        
      }
      
    }
    
    if(it%%10 == 0){
      print(K_v)  
    }
    
    if((it > n_burn) & ((it - n_burn)%%thin == 0)){
      iter <- (it - n_burn)%/%thin
      
      prop_var_out[[iter]] <- prop_var
      alloc_v_out[[iter]] <- alloc_v
      nj_v_out[[iter]] <- nj_v[[iv]]
      K_v_out[iter,] <- K_v 
      chain_out[[iter]] <- stndnaming(sevt_fit(tree)) 
      lambda_out[[iter]] <- lambda_current
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
  
  OUTPUT_MCMC <- list("alloc_v_out" = alloc_v_out, "nj_v_out" = nj_v_out, "K_v_out" = K_v_out, "chain_out" = chain_out, "lambda_out" = lambda_out, "prop_var_out" = prop_var_out)
  return(OUTPUT_MCMC)
}

