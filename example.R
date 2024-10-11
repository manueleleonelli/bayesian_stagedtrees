# List of CRAN packages
cran_packages <- c("ATbounds", "stagedtrees", "salso", "mcclust.ext")

# Install missing CRAN packages
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Install mcclust from GitHub if not installed
if (!require("mcclust", character.only = TRUE)) {
  # Install devtools if not installed
  if (!require("devtools", character.only = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("sarawade/mcclust")
}

library(ATbounds) ## Package with dataset
library(stagedtrees)
library(salso)
library(mcclust.ext)

data("EFM")

# Convert the dataset to a data.frame and drop the 'year' column
data <- as.data.frame(EFM[ , !names(EFM) %in% "year"])

# Convert each column to a factor with levels 'Yes' and 'No'
data <- lapply(data, function(x) {
  factor(x, levels = c(0, 1), labels = c("No", "Yes"))
})

# Convert the list back to a data.frame and order variables
data <- as.data.frame(data)
data <- data[,c(3,4,5,1,2)]


### Choose a starting tree for the MCMC and subset of data for uncertainty
tree <- stndnaming(stages_hc(full(data)))



tree <- stndnaming(stages_hc(indep(data)))
### tree: starting input of the MCMC
### a: imaginary sample size for edge probability priors
### prior: choice of prior distribution over staged tree space
### scope: which variables to consider
### beta: additional input for Heckerman prior

it <- 50000
a <- 1
prior <- "Heckerman"
scope <- c("monitor")
beta <- list(monitor = 0.3)

ciao <- mcmc(tree,it = it,a = a,prior = prior, scope = scope, beta=beta)

## Burn-in and thinning
indices <- seq(10000, length(ciao), by = 100)
result <- ciao[indices]

## Check some measure of convergence
x <- lapply(result,logLik)
plot(1:length(x),x,type="l")

nstages <- function(tree) sum(unlist(lapply(tree$stages,function(i) length(unique(i)))))
y <- lapply(result,nstages)
plot(1:length(y),y,type="l")

## This should be already within the MCMC but to be sure I do it here too, just renaming
result <- lapply(result,stndnaming)


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

# Initialize the objects
estimate_VI <- estimate_B <- lower <- upper <- horizontal <- result[[1]]

# Use lapply to loop over 'scope' and process each variable
for (v in scope) {
  results <- process_variable(v, result, estimate_B, estimate_VI, lower, upper, horizontal)
   estimate_B <- results$estimate_B
  estimate_VI <- results$estimate_VI
  lower <- results$lower
  upper <- results$upper
  horizontal <- results$horizontal
}

par(mfrow=c(2,2),mar=c(0,0,0,0))
plot(estimate_VI)
plot(horizontal)
plot(lower)
plot(upper)

