### CRP Function
# alpha: concentration parameter
# assignments: current assignment of vertices to stages
# vertex: the index of the vertex to be reassigned
crp_assignment <- function(assignments, alpha) {
  unique_stages <- unique(assignments)
  n_stages <- length(unique_stages)
  counts <- table(assignments)
  
  # CRP probabilities for existing stages
  crp_probs <- counts / (sum(counts) + alpha)
  
  # Probability of creating a new stage
  new_stage_prob <- alpha / (sum(counts) + alpha)
  
  # Combine probabilities
  probs <- c(crp_probs, new_stage_prob)
  
  # Sample a new stage assignment
  new_stage_idx <- sample(1:(n_stages + 1), 1, prob = probs)
  
  if (new_stage_idx == n_stages + 1) {
    # Create a new stage
    new_stage <- paste0("stage_", n_stages + 1)
  } else {
    # Assign to an existing stage
    new_stage <- unique_stages[new_stage_idx]
  }
  
  return(new_stage)
}

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
mcmc_crp <- function(tree, it, a = 1, alpha = 1, prior = "Friedman", burn = 0, thin = 0, scope = NULL, beta = NULL, tau = NULL){
  
  # Check the inputs (use check_mcmc function as before, omitted here for brevity)
  
  acc <- rep(0, length(scope))
  chain <- vector("list", it+1) # Create list where results will be stored
  chain[[1]] <- stndnaming(tree)
  priors <- comp_prior(chain[[1]], a)
  if(is.null(scope)){scope <- sevt_varnames(tree)[-1]}
  pb <- txtProgressBar(min = 0, max = it, style = 3)
  
  ## Cycle over iteration
  for(i in 1:(it)){
    chain[[i+1]] <- chain[[i]]
    ## Cycle over variables that are estimated
    for(v in scope){
      tree <- tree1 <- chain[i+1][[1]]
      
      ### Weight of each vertex
      counts <- table(tree$stages[[v]])[tree$stages[[v]]]
      ## Sampling of two vertices
      vertices <- sample(length(tree$stages[[v]]), 2, replace = FALSE, prob=counts)
      
      if(tree$stages[[v]][vertices[1]] != tree$stages[[v]][vertices[2]]){ # MERGE MOVE
        
        # Get stages for the selected vertices
        stage_1 <- tree$stages[[v]][vertices[1]]
        stage_2 <- tree$stages[[v]][vertices[2]]
        
        ## Merge the stages
        merged_stage <- crp_assignment(tree$stages[[v]], alpha)
        
        tree1$stages[[v]][which(tree1$stages[[v]] == stage_2)] <- merged_stage
        chain[[i+1]] <- stndnaming(sevt_fit(tree1, scope = v)) 
        acc[match(v, scope)] <- acc[match(v, scope)] + 1
        priors <- comp_prior(chain[[i+1]], a)
        
      } else { ## SPLIT MOVE
        
        ## Get the name of the current stage
        current_stage <- tree$stages[[v]][vertices[1]]
        
        ## Randomly assign some vertices to a new stage using CRP
        new_stage <- crp_assignment(tree$stages[[v]], alpha)
        
        tree1$stages[[v]][vertices[2]] <- new_stage
        ## Fit the new tree
        tree1 <- sevt_fit(tree1, scope = v)
        ## Compute prior for new tree
        priors <- comp_prior(tree1, a)
        
        chain[[i+1]] <- stndnaming(tree1)
        acc[match(v, scope)] <- acc[match(v, scope)] + 1
      }
    }
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(chain)
}

### Example Usage
# Assuming you have a tree object named 'my_tree' and the appropriate setup.
set.seed(123)  # For reproducibility
result_chain <- mcmc_crp(tree = tree, it = 50000, a = 0.1, alpha = 100000, scope = scope)

# Print the first result of the MCMC chain
print(result_chain[[1]])
