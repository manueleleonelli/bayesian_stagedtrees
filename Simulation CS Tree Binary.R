library(stagedtrees)
library("e1071")
library(salso)
library(fossil)

## Definisco l'albero
tree_def <- list(X1 = c("no","yes"),
                 X2 = c("no","yes"),
                 X3 = c("no","yes"),
                 X4 = c("no","yes"),
                 X5 = c("no","yes"),
                 X6 = c("no","yes"))

## Lo creo come oggetto per pacchetto stagedtrees
tree_def <- sevt(tree_def, full= F)

tree_def$stages$X2 <- c("1","2")
tree_def$stages$X3 <- c("1","2","1","3")
tree_def$stages$X4 <- c("1","2","3","4","1","1","3","4")
tree_def$stages$X5 <- rep(c("1","2","3","4"),4)


## PARTIZIONE X6 DEFINITA SOLO DA CONTEXT-SPECIFIC & SYMMETRIC INDEPENDENCES
# X6 INDEP X1 | X2, X3, X4, X5
# X6 INDEP X4 | X2 = no
# X6 INDEP X3 | X5 = yes, X4 = yes
# X6 INDEP X5 | X2 = yes , X4 = no

tree_def$stages$X6 <- rep(c("1","2","1","2",
                            "3","2","3","2",
                            "9","9","11","12",
                            "13","13","15","12"),2)

tree_def <- stndnaming(tree_def)                        

## Assegno probabilita all'albero
## Per tutte le variabili precedenti metto 50/50 cosi da avere osservazioni in tutti i paths
tree_def$prob <- list()
tree_def$prob$X1 <- list("NA"=c("no"= 0.6,"yes" = 0.4))
tree_def$prob$X2 <- list("1" = c("no"= 0.5,"yes" = 0.5),
                         "2" = c("no"= 0.4,"yes" = 0.6))
tree_def$prob$X3 <- list("1" = c("no"= 0.4,"yes" = 0.6),
                         "2" = c("no"= 0.5,"yes" = 0.5),
                         "3" = c("no"= 0.7,"yes" = 0.3))
tree_def$prob$X4 <- list("1" = c("no"= 0.5,"yes" = 0.5),
                         "2" = c("no"= 0.4,"yes" = 0.6),
                         "3" = c("no"= 0.3,"yes" = 0.7),
                         "4" = c("no"= 0.7,"yes" = 0.3))
tree_def$prob$X5 <- list("1" = c("no"= 0.5,"yes" = 0.5),
                         "2" = c("no"= 0.4,"yes" = 0.6),
                         "3" = c("no"= 0.3,"yes" = 0.7),
                         "4" = c("no"= 0.7,"yes" = 0.3))
tree_def$prob$X6 <- list("1" = c("no"= 0.2,"yes" = 0.8),
                         "2" = c("no"= 0.7,"yes" = 0.3),
                         "3" = c("no"= 0.3,"yes" = 0.7),
                         "4" = c("no"= 0.5,"yes" = 0.5),
                         "5" = c("no"= 0.6,"yes" = 0.4),
                         "6" = c("no"= 0.4,"yes" = 0.6),
                         "7" = c("no"= 0.8,"yes" = 0.2),
                         "8" = c("no"= 0.9,"yes" = 0.1))

n_burn <- 1000
thin <- 2
n_save <- 1000
update_CRP <- TRUE
update_SM <- TRUE
scope <- "X6"
beta <- list(monitor = 0.3)
a <- 1
prior <- "Heckerman"
kappa <- 1
csi <- 0.4

source("mcmc_crp_distances.R")

# Define vectors of values for csi and kappa
# Define vectors of values for csi, kappa, and sample sizes
csi_values <- c(0, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1)
kappa_values <- c(0.01, 0.05, 0.1, 0.5, 1, 5)
sample_sizes <- c(500, 1000, 2500, 5000,7500, 10000)

# Initialize a list to store estimate_VI results
estimate_VI_results <- list()

# Loop over all sample sizes
for (sample_size in sample_sizes) {
  # Generate new data for each sample size
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
  
  # Check and regenerate data if necessary
  while (any(tree$ctables$X6[,1] == 0 & tree$ctables$X6[,2] == 0)) {
    data <- sample_from(tree_def, size = sample_size)
    tree <- stndnaming(stages_fbhc(full(data)))
  }
  
  # Loop over all combinations of csi and kappa for the current sample size
  for (csi_val in csi_values) {
    for (kappa_val in kappa_values) {
      
      estimate_VI_results[[paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)]] <- list()
      
      # Run MCMC for each combination of csi and kappa
      ciao <- mcmc_crp_distances(tree, n_save = n_save, n_burn = n_burn, thin = thin, 
                                 a = a, kappa = kappa_val, csi = csi_val, prior = prior, 
                                 scope = scope, beta = beta, update_SM = update_SM, 
                                 update_CRP = update_CRP)
      
      # Burn-in and thinning
      result <- ciao$chain_out
      
      # Initialize the storage for current combination
      estimate_VI <- result[[1]]
      
      # Store data object
      estimate_VI_results[[paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)]][["data"]] <- data
      
      # Store salso results for each variable in scope
      for (v in scope) {
        # Create matrix for variable v
        mat <- t(sapply(result, function(res) res$stages[[v]]))
        
        # Perform salso with VI loss for estimate_VI
        estimate_VI$stages[[v]] <- salso(mat, loss = "VI")
      }
      
      # Store the ce_randomized diff() output for each model
      ce_diffs <- sapply(result, function(model) {
        diff(ce_randomized(sevt_fit(model, data, lambda = a/2), outcome = "X6", treatment = "X5")[, 2])
      })
      
      # Add ce_diffs to the estimate_VI object
      estimate_VI_results[[paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)]][["ce_diffs"]] <- ce_diffs
      
      # Store the result in the main list
      estimate_VI_results[[paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)]][["tree"]] <- estimate_VI
    }
  }
}


library(ggplot2)
library(reshape2)
library(gridExtra)

# Extract unique stage counts for each entry
results <- data.frame(Sample_Size = integer(), Csi = numeric(), Kappa = numeric(), Unique_Stages = integer(), Hamming = integer(), Local = integer(),Parents = integer(), Rand= numeric())

for (key in names(estimate_VI_results)) {
  # Parse the sample, csi, and kappa values from the key
  parts <- unlist(strsplit(key, "_"))
  sample_size <- as.integer(parts[2])
  csi <- as.numeric(parts[4])
  kappa <- as.numeric(parts[6])
  
  # Calculate the number of unique stages for X6
  unique_stages <- length(unique(estimate_VI_results[[key]]$tree$stages$X6))
  hamming <- sum(hamming_stages(tree_def,estimate_VI_results[[key]]$tree,return_tree = T)$X6)
  locals <- length(as_parentslist(estimate_VI_results[[key]]$tree)$X6$local)
  parents <- length(as_parentslist(estimate_VI_results[[key]]$tree)$X6$parents)
  rand <- rand.index(as.numeric(tree_def$stages$X6),as.numeric(estimate_VI_results[[key]]$tree$stages$X6))
  # Add the result to the dataframe
  results <- rbind(results, c( sample_size, csi, kappa, unique_stages,hamming,locals,parents,rand))
}

results <- as.data.frame(results)
colnames(results) <- c("Sample_Size", "Csi","Kappa","Unique_Stages","Hamming","Locals","Parents","Rand")
results$Unique_Stages <- results$Unique_Stages - length(unique(tree_def$stages$X6))

heatmap_list <- list()
for (sample_size in sample_sizes) {
  # Subset data for the current sample size
  data_subset <- subset(results, Sample_Size == sample_size)
  
  # Create heatmap
  p <- ggplot(data_subset, aes(x = as.factor(Kappa), y = as.factor(Csi), fill = Hamming)) +
    geom_tile() +
    scale_fill_gradient2(name = "Rand", low = "red", mid = "white", high = "blue", midpoint = 0) +
    labs(title = paste("Sample Size:", sample_size), x = "Kappa", y = "Csi") +
    theme_minimal()
  
  # Append to list
  heatmap_list[[as.character(sample_size)]] <- p
}

# Arrange heatmaps in a grid
do.call("grid.arrange", c(heatmap_list, ncol = 2))

