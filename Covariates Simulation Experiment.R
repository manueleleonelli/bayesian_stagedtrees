library(stagedtrees)
library("e1071")
library(salso)
library(fossil)
library(mcclust.ext)

set.seed(93)
## Definisco l'albero
tree_def <- list(X1 = c("no","yes"),
                 X2 = c("no","yes"),
                 X3 = c("no","yes"),
                 X4 = c("no","yes"))

## Lo creo come oggetto per pacchetto stagedtrees
tree_def <- sevt(tree_def, full= F)

tree_def$stages$X2 <- c("1","2")
tree_def$stages$X3 <- c("1","2","1","3")
tree_def$stages$X4 <- c("1","1","2","2","1","1","3","3")

tree_def <- stndnaming(tree_def)                        

## Assegno probabilita all'albero
## Per tutte le variabili precedenti metto 50/50 cosi da avere osservazioni in tutti i paths
tree_def$prob <- list()
tree_def$prob$X1 <- list("NA"=c("no"= 0.6,"yes" = 0.4))
tree_def$prob$X2 <- list("1" = c("no"= 0.7,"yes" = 0.3),
                         "2" = c("no"= 0.3,"yes" = 0.7))
tree_def$prob$X3 <- list("1" = c("no"= 0.4,"yes" = 0.6),
                         "2" = c("no"= 0.7,"yes" = 0.3),
                         "3" = c("no"= 0.75,"yes" = 0.25))
tree_def$prob$X4 <- list("1" = c("no" = 0.7, "yes" = 0.3),
                         "2" = c("no" = 0.3, "yes" = 0.7),
                         "3" = c("no" = 0.75, "yes" = 0.25))

sample_size <- 500
scope <- c("X3","X4")
n_v <- length(scope)
a <- 1
n_burn <- 1000
thin <- 1
n_save <- 1000
csi <- 0.25
kappa <- 1
mu0 <- 0
kappa0 <- 1
alpha0 <- 1
beta0 <- 1

set.seed(123)
data <- sample_from(tree_def, size = sample_size)
tree <- stndnaming(stages_fbhc(full(data)))

# Ensure data is valid (no empty counts for X4)
while (any(tree$ctables$X4[, 1] == 0 & tree$ctables$X4[, 2] == 0)) {
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
}

# #Non-informative covariate
# data$Z1 <- rnorm(sample_size)
#well-separated covariates within clusters
data$Z1 <- rep(NA,sample_size)
data$Z1[which(data$X2=="no")] <- rnorm(length(which(data$X2=="no")),0,sd =0.5)
data$Z1[which(data$X1=="yes" & data$X2 == "yes")] <- rnorm(length(which(data$X1=="yes" & data$X2 == "yes")),-1,sd =0.1)
data$Z1[which(data$X1=="no" & data$X2 == "yes")] <- rnorm(length(which(data$X1=="no" & data$X2 == "yes")),1,sd =0.1)
hist(data$Z1,freq=F)

data$group <- NA

data$group[data$X2 == "no"] <- "A"
data$group[data$X1 == "yes" & data$X2 == "yes"] <- "B"
data$group[data$X1 == "no" & data$X2 == "yes"] <- "C"

# Optionally convert to factor with desired order
data$group <- factor(data$group, levels = c("A", "B", "C"))

hist(data$Z1[which(data$X2=="no")], col = rgb(1,0,0,alpha = 0.5), main = "", xlab = "Z", ylab = "", xlim = c(-2,2), freq = FALSE)
hist(data$Z1[which(data$X1=="yes" & data$X2 == "yes")], col = rgb(0,1,0,alpha = 0.5), add = TRUE, breaks = 35, freq = FALSE)
hist(data$Z1[which(data$X1=="no" & data$X2 == "yes")], col = rgb(0,0,1,alpha = 0.5), add = TRUE, freq = FALSE)
legend("topright", legend = c("Cluster 1", "Cluster 2", "Cluster 3"), fill = c(2,3,4), bty = "n")


data <- data[,1:5]
#Number of continuous variables
q_Z <- sum(sapply(data, is.numeric))
#Penalty coefficient for continuous variables in PPMx
lambda_Z <- rep(0.5, q_Z)


#Algorithm specs
update_SM = TRUE
update_CRP= TRUE

ciao <- mcmc_crp_ppmx_Hamming_Z(tree, data, n_save, n_burn = n_burn, thin = thin, a = a, 
                                kappa = kappa, csi = csi, lambda_Z = lambda_Z,
                                scope = scope, update_SM = update_SM, update_CRP= update_CRP, 
                                mu0 = mu0, kappa0 = kappa0, alpha0 = alpha0, beta0 = beta0)


## Burn-in and thinning
result <- ciao$chain_out

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

est_cov <- estimate_VI





ciao <- mcmc_crp_ppmx_Hamming(tree, data, n_save, n_burn = n_burn, thin = thin, a = a, 
                              kappa = kappa, csi = csi,
                              scope = scope, update_SM = update_SM, update_CRP= update_CRP)


## Burn-in and thinning
result <- ciao$chain_out

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



png("Simulated_HammingStages.png")
plot(estimate_VI)
dev.off()

