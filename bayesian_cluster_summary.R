# FUNCTION FOR BAYESIAN CLUSTERING

## Create plots chatgpt
process_variable <- function(v, result, estimate_B, estimate_VI, lower, upper, horizontal) {
  # Create a matrix for each 'v' across all 'result'
  mat <- t(sapply(result, function(res) res$stages[[v]]))
  
  # Perform salso with Binder loss
  estimate_B$stages[[v]] <- salso(mat, loss = "binder")
  
  # Perform salso with VI loss
  estimate_VI$stages[[v]] <- salso(mat, loss = "VI")
  
  # Compute credible ball for VI distance
  ball <- credibleball(estimate_VI$stages[[v]], mat, c.dist = "VI", alpha = 0.1)
  
  # Update lower, upper, and horizontal stages
  lower$stages[[v]] <- ball$c.lowervert[1,]
  upper$stages[[v]] <- ball$c.uppervert[1,]
  horizontal$stages[[v]] <- ball$c.horiz[1,]
  
  # Return the updated objects
  list(estimate_B = estimate_B, estimate_VI = estimate_VI, lower = lower, upper = upper, horizontal = horizontal)
}