library(stagedtrees)
library("e1071")

## Definisco l'albero
tree_def <- list(X1 = c("no","maybe","yes"),
             X2 = c("no","maybe","yes"),
             X3 = c("no","maybe","yes"),
             X4 = c("no","yes"),
             X5 = c("no","yes"))

## Lo creo come oggetto per pacchetto stagedtrees
tree_def <- sevt(tree_def, full= F)

tree_def$stages$X2 <- c("1","1","2")
tree_def$stages$X3 <- c("1","2","3","1","2","3","1","1","1")
tree_def$stages$X4 <- rep(c("1","2","3","1","2","3","4","5","6"),3)

## PARTIZIONE X5 DEFINITA SOLO DA CONTEXT-SPECIFIC & SYMMETRIC INDEPENDENCES
## NO PARTIAL RELATIONSHIP BETWEEN VARIABLES
# X5 INDEP X1 | X2, X3, X4
# X5 INDEP X4, X3 | X2 = no
# X5 INDEP X2 | X3 = maybe
tree_def$stages$X5 <- rep(c("1","1","1","1","1","1",
                        "7","8","1","1","11","12",
                        "13","14","1","1","17","18"),3)

tree_def <- stndnaming(tree_def)                        

## Assegno probabilita all'albero
tree_def$prob <- list()
tree_def$prob$X1 <- list("NA"=c(0.4,0.3,0.3))
tree_def$prob$X2 <- list("1" = c(0.5,0.25,0.25),
                         "2" = c(0.3,0.3,0.4))
tree_def$prob$X3 <- list("1" = c(0.2,0.2,0.6),
                         "2" = c(0.3,0.4,0.3),
                         "3" = c(0.2,0.4,0.4))
tree_def$prob$X4 <- list("1" = c(0.5,0.5),
                         "2" = c(0.4,0.6),
                         "3" = c(0.3,0.7),
                         "4" = c(0.7,0.3),
                         "5" = c(0.6,0.4),
                         "6" = c(0.55,0.45))
tree_def$prob$X5 <- list("1" = c(0.2,0.8),
                         "2" = c(0.3,0.7),
                         "3" = c(0.4,0.6),
                         "4" = c(0.5,0.5),
                         "5" = c(0.6,0.4),
                         "6" = c(0.7,0.3),
                         "7" = c(0.8,0.2),
                         "8" = c(0.9,0.1),
                         "9" = c(0.1,0.9))

## Simulo i dati dall'albero
data <- sample_from(tree_def, size = 2000, seed=2024)
tree <- stndnaming(stages_fbhc(full(data)))
while(any(tree$ctables$X5[,1]== 0 & tree$ctables$X5[,2]== 0)){
  data <- sample_from(tree_def, size = N)
  tree <- stndnaming(stages_fbhc(full(data)))
}

n_burn <- 1000
thin <- 2
n_save <- 1000
update_CRP <- TRUE
update_SM <- TRUE
scope <- "X5"
beta <- list(monitor = 0.3)
a <- 1
prior <- "Heckerman"
kappa <- 1
csi <- 0.4

source("mcmc_crp_distances.R")


csi_values <- c(0, 1/112, 1/56, 1/28, 1/14, 1/7, 1/2, 1)
kappa_values <- c(0.01, 0.05, 0.1, 0.5, 1, 5)
sample_sizes <- c(500, 1000, 2500, 5000, 10000)

# Initialize a list to store estimate_VI results
estimate_VI_results <- list()

# Loop over all sample sizes
for (sample_size in sample_sizes) {
  # Generate new data for each sample size
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
  
  # Check and regenerate data if necessary
  while (any(tree$ctables$X5[,1] == 0 & tree$ctables$X5[,2] == 0)) {
    data <- sample_from(tree_def, size = sample_size)
    tree <- stndnaming(stages_fbhc(full(data)))
  }
  
  # Loop over all combinations of csi and kappa for the current sample size
  for (csi_val in csi_values) {
    for (kappa_val in kappa_values) {
      # Run MCMC for each combination of csi and kappa
      ciao <- mcmc_crp_distances(tree, n_save = n_save, n_burn = n_burn, thin = thin, 
                                 a = a, kappa = kappa_val, csi = csi_val, prior = prior, 
                                 scope = scope, beta = beta, update_SM = update_SM, 
                                 update_CRP = update_CRP)
      
      # Burn-in and thinning
      result <- ciao$chain_out
      
      # Initialize estimate_VI for the current run
      estimate_VI <- result[[1]]
      
      # Process each variable in scope and calculate estimate_VI
      for (v in scope) {
        # Create matrix for variable v
        mat <- t(sapply(result, function(res) res$stages[[v]]))
        
        # Perform salso with VI loss for estimate_VI
        estimate_VI$stages[[v]] <- salso(mat, loss = "VI")
      }
      
      # Store estimate_VI in the results list with a unique name for each combination
      estimate_VI_results[[paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)]] <- estimate_VI
    }
  }
}

library(ggplot2)
library(gridExtra)

# Initialize a dataframe for storing the median results
results <- data.frame(Sample_Size = integer(), Csi = numeric(), Kappa = numeric(),
                      Median_Unique_Stages = numeric(), Median_Hamming = numeric(),
                      Median_Locals = numeric(), Median_Parents = numeric(), Median_Rand = numeric())

# Loop through each combination of parameters in the estimate_VI_results list
for (key in names(estimate_VI_results)) {
  # Parse the sample, csi, and kappa values from the key
  parts <- unlist(strsplit(key, "_"))
  sample_size <- as.integer(parts[2])
  csi <- as.numeric(parts[4])
  kappa <- as.numeric(parts[6])
  
  # Initialize storage for metrics across the 5 replications
  unique_stages_list <- c()
  hamming_list <- c()
  locals_list <- c()
  parents_list <- c()
  rand_list <- c()
  
  # Loop through each replication (rep_1 to rep_5)
  for (rep_key in names(estimate_VI_results[[key]])) {
    rep_result <- estimate_VI_results[[key]][[rep_key]]$VI  # Get the VI object
    
    # Compute metrics for this replication
    unique_stages_list <- c(unique_stages_list, length(unique(rep_result$stages$X6)))
    hamming_list <- c(hamming_list, sum(hamming_stages(tree_def, rep_result, return_tree = TRUE)$X6))
    locals_list <- c(locals_list, length(as_parentslist(rep_result)$X6$local))
    parents_list <- c(parents_list, length(as_parentslist(rep_result)$X6$parents))
    rand_list <- c(rand_list, rand.index(as.numeric(tree_def$stages$X6), as.numeric(rep_result$stages$X6)))
  }
  
  # Compute median values across the 5 replications
  median_unique_stages <- median(unique_stages_list) - length(unique(tree_def$stages$X6))  # Adjusted baseline
  median_hamming <- median(hamming_list)
  median_locals <- median(locals_list)
  median_parents <- median(parents_list)
  median_rand <- median(rand_list)
  
  # Append results for this parameter combination
  results <- rbind(results, data.frame(
    Sample_Size = sample_size,
    Csi = csi,
    Kappa = kappa,
    Median_Unique_Stages = median_unique_stages,
    Median_Hamming = median_hamming,
    Median_Locals = median_locals,
    Median_Parents = median_parents,
    Median_Rand = median_rand
  ))
}

# Visualize results using heatmaps
heatmap_list <- list()
for (sample_size in unique(results$Sample_Size)) {
  # Subset data for the current sample size
  data_subset <- subset(results, Sample_Size == sample_size)
  
  # Create heatmap for median Rand index
  p <- ggplot(data_subset, aes(x = as.factor(Kappa), y = as.factor(Csi), fill = Median_Rand)) +
    geom_tile() +
    scale_fill_gradient2(name = "Rand", low = "red", mid = "white", high = "blue", midpoint = 0) +
    labs(title = paste("Sample Size:", sample_size), x = "Kappa", y = "Csi") +
    theme_minimal()
  
  # Append the heatmap to the list
  heatmap_list[[as.character(sample_size)]] <- p
}

# Arrange heatmaps in a grid
do.call("grid.arrange", c(heatmap_list, ncol = 2))


