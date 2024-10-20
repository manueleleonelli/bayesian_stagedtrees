
### Function to compute the prior probability of edge probabilities using imaginary sample size
comp_prior <- function(tree, a = 16){
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
  # Remove names from both levels (first-level list and inner elements)
  priors_transformed_no_names <- lapply(priors_filtered_repeated_divided, function(vec) {
    split_list <- split(vec, names(vec))
    #    split_list_no_names <- lapply(split_list, unname)
    return(split_list)  # Remove the names from the first level
  })
  #  priors_transformed_no_names <- unname(priors_transformed_no_names)
  return(priors_transformed_no_names)
}


### Check that all inputs of the MCMC function are correct
check_mcmc <- function(tree, it, a, prior, burn, thin, scope, beta, tau) {
  
  # Check if 'prior' is valid
  if (!prior %in% c("Friedman", "Uniform", "Heckerman", "Pensar")) {
    stop("Invalid prior specified. Must be one of 'Friedman', 'Uniform', 'Heckerman', or 'Pensar'.")
  }
  
  # Check if 'tree' is of class 'sevt'
  if (!inherits(tree, "sevt")) {
    stop("'tree' must be an object of class 'sevt'.")
  }
  
  # Check if 'a' is a number greater than 0
  if (!is.numeric(a) || a <= 0) {
    stop("'a' must be a number greater than 0.")
  }
  
  # Check if 'it' is a positive integer
  if (!is.numeric(it) || it <= 0 || floor(it) != it) {
    stop("'it' must be a positive integer.")
  }
  
  # Check if 'burn' is a non-negative integer
  if (!is.numeric(burn) || burn < 0 || floor(burn) != burn) {
    stop("'burn' must be a non-negative integer.")
  }
  
  # Check if 'thin' is a non-negative integer
  if (!is.numeric(thin) || thin < 0 || floor(thin) != thin) {
    stop("'thin' must be a non-negative integer.")
  }
  
  # Check if 'scope' contains only elements from sevt_varnames(tree)
  valid_scope_vars <- sevt_varnames(tree)  # Assuming sevt_varnames is a valid function
  if (!all(scope %in% valid_scope_vars)) {
    stop("All elements of 'scope' must be present in the variable names of 'tree'.")
  }
  
  # Check additional inputs for 'Heckerman' prior
  if (prior == "Heckerman") {
    if (is.null(beta)) {
      stop("For 'Heckerman' prior, 'beta' must be provided.")
    }
    if (length(beta) != length(scope)) {
      stop("Length of 'beta' must be the same as the length of 'scope'.")
    }
    if (!all(beta > 0 & beta <= 1)) {
      stop("Each element of 'beta' must be in the range (0, 1].")
    }
  }
  
  # Check additional inputs for 'Pensar' prior
  if (prior == "Pensar") {
    if (is.null(tau)) {
      stop("For 'Pensar' prior, 'tau' must be provided.")
    }
    if (length(tau) != length(scope)) {
      stop("Length of 'tau' must be the same as the length of 'scope'.")
    }
    if (!all(tau >= 0)) {
      stop("Each element of 'tau' must be in the range [0, infinity).")
    }
  }
  
  message("All input checks passed.")
}


### MCMC Algorithm
mcmc <- function(tree, it, a = 1, prior = "Friedman", burn = 0, thin = 0, scope = NULL, beta = NULL, tau = NULL){
  check_mcmc(tree, it, a, prior, burn, thin, scope, beta, tau) 
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
      ## Sampling of two vertices
      vertices <- sample(length(tree$stages[[v]]), 2, replace = F) 
      if(tree$stages[[v]][vertices[1]] != tree$stages[[v]][vertices[2]]){# MERGE MOVE
        
        ## Deriving the staging indexes
        vertices <- tree$stages[[v]][vertices]
        
        ## Computing counts for the first stage
        ix1 <- tree$stages[[v]] == vertices[1]
        pr_a <- priors[[v]][[vertices[1]]]
        if (sum(ix1) > 1) {
          tt_a <- apply(tree$ctables[[v]][ix1, ], MARGIN = 2, sum) 
        } else {
          tt_a <- tree$ctables[[v]][ix1, ]
        }
        
        ## Computing counts for the second stage
        ix2 <- tree$stages[[v]] == vertices[2]
        pr_b <- priors[[v]][[vertices[2]]]
        if (sum(ix2) > 1) {
          tt_b <- apply(tree$ctables[[v]][ix2, ], MARGIN = 2, sum)
        } else {
          tt_b <- tree$ctables[[v]][ix2, ]
        }
        
        ## Compute terms for marginal likelihood ratio
        r_1 <- -lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) 
        r_2 <-  - lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a))
        r_12 <- +lgamma(sum(pr_b+pr_a)) - lgamma(sum(pr_b+pr_a+tt_a+tt_b)) - sum(lgamma(pr_b+pr_a)) + sum(lgamma(tt_a+tt_b+pr_a+pr_b))
        
        
        ## Compute Tree prior ratio
        if(prior == "Heckerman"){lp <- -log(beta[[v]])} else if(prior == "Friedman") {lp <- log(length(tree$stages[[v]]) - length(unique(tree$stages[[v]]))) -log(length(unique(tree$stages[[v]])))} else if(prior == "Pensar"){lp <- -tau*log(nrow(tree$data_raw))*(1+tau)^(length(unique(tree$stages[[v]]))-1)} else {lp <- 1}
        lp <- 0
        ## Compute transition probability ratio
        
        ## ChatGPT Transition  
        #trans <-   lchoose(sum(ix1+ix2),2)  + lchoose((sum(ix1+ix2)-2),(sum(ix1)-1)) +log(0.5)*(sum(ix1+ix2)-2) - log(sum(ix1)) - log(sum(ix2)) 
        
        ## Uniform transition  
        trans <-   log(0.5)*(sum(ix1+ix2)-2) -lchoose(sum(ix1+ix2),2) - log(sum(ix1+ix2)) + lchoose(length(tree$stages[[v]]),2)
       
        #trans <- 0
         
        ## Decide on MCMC move
        if(runif(1)<= r_12 + r_1 + r_2 + lp + trans){
          tree1$stages[[v]][which(tree1$stages[[v]]==vertices[2])] <- vertices[1]
          chain[[i+1]] <- stndnaming(sevt_fit(tree1, scope=v)) 
          acc[match(v, scope)] <- acc[match(v, scope)]+1
          priors <- comp_prior(chain[[i+1]], a)
        }
        
      } else { ## SPLIT MOVE
        
        ## Get name of stage
        name <- tree$stages[[v]][vertices[1]]
        ## Get all vertices in the same stage minus those sampled
        target_indices <- setdiff(which(tree$stages[[v]] == name), vertices)
        ## Create new stage
        tree1$stages[[v]][vertices[2]] <- "a"
        ## Randomly assign vertices in chosen stage
        assignments <- sample(c("a", name), length(target_indices), replace = TRUE)
        tree1$stages[[v]][target_indices] <- assignments
        ## Fit new tree
        tree1 <- sevt_fit(tree1, scope = v)
        ## Compute edges prior for new tree
        priors <- comp_prior(tree1, a)
        
        ## Counts for new stage 1
        ix1 <- tree1$stages[[v]] == name
        pr_a <- priors[[v]][[name]]
        if (sum(ix1) > 1) {
          tt_a <- apply(tree1$ctables[[v]][ix1, ], MARGIN = 2, sum) 
        } else {
          tt_a <- tree1$ctables[[v]][ix1, ]
        }
        
        ## Counts for new stage 2
        ix2 <- tree1$stages[[v]] == "a"
        pr_b <- priors[[v]][["a"]]
        if (sum(ix2) > 1) {
          tt_b <- apply(tree1$ctables[[v]][ix2, ], MARGIN = 2, sum)
        } else {
          tt_b <- tree1$ctables[[v]][ix2, ]
        }
        
        ## Compute terms for marginal likelihood ratio
        r_1 <- -(-lgamma(sum(pr_a)) + lgamma(sum(tt_a+pr_a)) + sum(lgamma(pr_a)) - sum(lgamma(tt_a+pr_a)) )
        r_2 <-  -(- lgamma(sum(pr_b)) + lgamma(sum(tt_b+pr_b)) + sum(lgamma(pr_b)) - sum(lgamma(tt_b+pr_a)))
        r_12 <- -(+lgamma(sum(pr_b+pr_a)) - lgamma(sum(pr_b+pr_a+tt_a+tt_b)) - sum(lgamma(pr_b+pr_a)) + sum(lgamma(tt_a+tt_b+pr_a+pr_b)))
        
        ## Compute terms for tree prior ratio
        if(prior == "Heckerman"){lp <- -(-log(beta[[v]]))} else if(prior == "Friedman") {lp <- -(log(length(tree$stages[[v]]) - length(unique(tree$stages[[v]]))) -log(length(unique(tree$stages[[v]]))))} else if(prior == "Pensar"){lp <- -(-tau*log(nrow(tree$data_raw))*(1+tau)^(length(unique(tree$stages[[v]]))-1))} else {lp <- 1}
        
        #lp <- 0
        #trans <- 0
        ## Compute ratio of transition probabilities
        
        ## ChatGPT Transition
        #trans <-  - (lchoose(sum(ix1+ix2),2)  + lchoose((sum(ix1+ix2)-2),(sum(ix1)-1)) +log(0.5)*(sum(ix1+ix2)-2) - log(sum(ix1)) - log(sum(ix2)) )
        
        ## Uniform Transition
         trans <- - lchoose(length(tree$stages[[v]]),2) - (log(0.5)*(sum(ix1+ix2)-2) -lchoose(sum(ix1+ix2),2) - log(sum(ix1+ix2)))
        
        if(runif(1)<= r_12 + r_1 + r_2 + lp + trans){
          chain[[i+1]] <- stndnaming(tree1) 
          acc[match(v, scope)] <- acc[match(v, scope)]+1
          priors <- comp_prior(chain[[i+1]], a)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(chain)
}
