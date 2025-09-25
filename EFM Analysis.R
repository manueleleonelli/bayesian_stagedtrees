library(stagedtrees)
library(ATbounds)
source("effect_sevt.R")
source("ppmx_hamming.R")
source("ppmx_hamming_Z.R")
source("bayesian_cluster_summary.R")
data <- EFM[,-6]
data$cesarean <- factor(data$cesarean)
data$monitor <- factor(data$monitor)
data$arrest <- factor(data$arrest)
data$breech <- factor(data$breech)
data$nullipar <- factor(data$nullipar)

tree <- stndnaming(full(data,order = c("nullipar","breech","arrest","monitor","cesarean")))

scope <- c("monitor","cesarean")
n_v <- length(scope)
a <- 1
n_burn <- 1000
thin <- 5
n_save <- 2000
csi <- 0.25
kappa <- 1
mu0 <- 0
kappa0 <- 1
alpha0 <- 1
beta0 <- 1
update_SM = TRUE
update_CRP= TRUE

set.seed(1988)
model_hamming <- mcmc_crp_ppmx_Hamming(tree, data, n_save, n_burn = n_burn, thin = thin, a = a, 
                              kappa = kappa, csi = csi,
                              scope = scope, update_SM = update_SM, update_CRP= update_CRP)


result <- model_hamming$chain_out

## Check some measure of convergence
x <- lapply(result,logLik)
plot(1:length(x),x,type="l")

nstages <- function(tree) sum(unlist(lapply(tree$stages,function(i) length(unique(i)))))
y <- lapply(result,nstages)
plot(1:length(y),y,type="l")


## This should be already within the MCMC but to be sure I do it here too, just renaming
result <- lapply(result,stndnaming)

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

ce_diffs <- sapply(result, function(model) {
  diff(ce_randomized(sevt_fit(model, data), outcome = "cesarean", treatment = "monitor")[, 2])
})
hist(ce_diffs)


# Initialize a list to store the differences for each result
cesarean_differences_list <- list()

# Loop over all models in result
for (i in seq_along(result)) {
  
  # Fit the staged tree model for this result
  fitted_model <- sevt_fit(result[[i]], data)
  
  # Extract the stage assignments and stage probabilities for cesarean
  stages_vec <- result[[i]]$stages$cesarean
  prob_list <- fitted_model$prob$cesarean
  
  # Extract the probability of "1" for each stage in order
  prob_vec <- sapply(stages_vec, function(stage) prob_list[[as.character(stage)]]["1"])
  
  # Compute the differences between second and first, fourth and third, etc.
  diffs <- prob_vec[seq(2, length(prob_vec), by = 2)] - prob_vec[seq(1, length(prob_vec), by = 2)]
  
  # Store in list
  cesarean_differences_list[[i]] <- diffs
}

# Optionally, convert to a matrix or data frame for easy inspection
cesarean_differences_matrix <- do.call(rbind, cesarean_differences_list)

# View result
print(cesarean_differences_matrix)

apply(cesarean_differences_matrix,2,mean)
apply(cesarean_differences_matrix,2,sd)
apply(cesarean_differences_matrix, 2, function(x) mean(x > 0))
apply(cesarean_differences_matrix, 2, function(x) mean(x == 0))
apply(cesarean_differences_matrix, 2, function(x) mean(x < 0))
