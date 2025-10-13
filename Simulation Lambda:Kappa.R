library(stagedtrees)
library("e1071")
library(salso)
library(fossil)
source("effect_sevt.R")
source("ppmx_hamming.R")
source("ppmx_hamming_Z.R")
source("bayesian_cluster_summary.R")

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


# Define vectors of values for csi, kappa, and sample sizes
csi_values <- c(0, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1)
kappa_values <- c(0.01, 0.05, 0.1, 0.5, 1, 5)
sample_sizes <- c(500, 1000, 2500, 5000, 7500, 10000)
num_reps <- 5

# Initialize a list to store estimate_VI results
estimate_VI_results <- list()

# Loop over all sample sizes
for (sample_size in sample_sizes) {
  
  # Loop over all combinations of csi and kappa
  for (csi_val in csi_values) {
    for (kappa_val in kappa_values) {
      
      # Initialize the storage for the current combination
      key <- paste0("sample_", sample_size, "_csi_", csi_val, "_kappa_", kappa_val)
      estimate_VI_results[[key]] <- list()
      
      # Loop for 25 repetitions
      for (rep in 1:num_reps) {
        # Generate new data for each repetition
        data <- sample_from(tree_def, size = sample_size)
        tree <- stndnaming(stages_fbhc(full(data)))
        
        # Ensure data is valid (no empty counts for X6)
        while (any(tree$ctables$X6[, 1] == 0 & tree$ctables$X6[, 2] == 0)) {
          data <- sample_from(tree_def, size = sample_size)
          tree <- stndnaming(stages_fbhc(full(data)))
        }
        
        # Run MCMC for the current combination of csi and kappa
        ciao <- mcmc_crp_ppmx_Hamming(tree, n_save = n_save, n_burn = n_burn, thin = thin, 
                                   a = a, kappa = kappa_val, csi = csi_val, 
                                   scope = scope,  update_SM = update_SM, 
                                   update_CRP = update_CRP)
        
        # Burn-in and thinning
        result <- ciao$chain_out
        
        # Initialize the storage for current repetition
        estimate_VI <- result[[1]]
        estimate_Binder <- result[[1]]
        
        # Store salso results for each variable in scope
        for (v in scope) {
          # Create matrix for variable v
          mat <- t(sapply(result, function(res) res$stages[[v]]))
          
          # Perform salso with VI and Binder loss
          estimate_VI$stages[[v]] <- salso(mat, loss = "VI")
          estimate_Binder$stages[[v]] <- salso(mat, loss = "binder")
        }
        
        # Store the ce_randomized diff() output for each model
        ce_diffs <- sapply(result, function(model) {
          diff(ce_randomized(sevt_fit(model, data, lambda = a/2), outcome = "X6", treatment = "X5")[, 2])
        })
        
        # Save results for current repetition
        estimate_VI_results[[key]][[paste0("rep_", rep)]] <- list(
          VI = estimate_VI,
          Binder = estimate_Binder,
          ce_diffs = ce_diffs
        )
        
        # Clear data to save memory
        rm(data)
      }
    }
  }
}

library(ggplot2)
library(reshape2)
library(gridExtra)

# Initialize a dataframe for storing median results
results <- data.frame(Sample_Size = integer(), Csi = numeric(), Kappa = numeric(), 
                      Median_Unique_Stages = numeric(), Median_Hamming = numeric(), 
                      Median_Locals = numeric(), Median_Parents = numeric(), Median_Rand = numeric())

for (key in names(estimate_VI_results)) {
  # Parse the sample, csi, and kappa values from the key
  parts <- unlist(strsplit(key, "_"))
  sample_size <- as.integer(parts[2])
  csi <- as.numeric(parts[4])
  kappa <- as.numeric(parts[6])
  
  # Collect metrics across replications
  unique_stages_list <- c()
  hamming_list <- c()
  locals_list <- c()
  parents_list <- c()
  rand_list <- c()
  
  for (rep_key in names(estimate_VI_results[[key]])) {
    # Retrieve tree and data for this replication
    rep_result <- estimate_VI_results[[key]][[rep_key]]
    
    unique_stages_list <- c(unique_stages_list, length(unique(rep_result$VI$stages$X6)))
    hamming_list <- c(hamming_list, sum(hamming_stages(tree_def, rep_result$VI, return_tree = TRUE)$X6))
    locals_list <- c(locals_list, length(as_parentslist(rep_result$VI)$X6$local))
    parents_list <- c(parents_list, length(as_parentslist(rep_result$VI)$X6$parents))
    rand_list <- c(rand_list, rand.index(as.numeric(tree_def$stages$X6), as.numeric(rep_result$VI$stages$X6)))
  }
  
  # Compute median values for this combination of sample size, csi, and kappa
  median_unique_stages <- median(unique_stages_list)- length(unique(tree_def$stages$X6))
  median_hamming <- median(hamming_list)
  median_locals <- median(locals_list)
  median_parents <- median(parents_list)
  median_rand <- median(rand_list)
  
  # Append to the results dataframe
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

# Visualize using heatmaps
heatmap_list <- list()
for (sample_size in unique(results$Sample_Size)) {
  # Subset data for the current sample size
  data_subset <- subset(results, Sample_Size == sample_size)
  
  # Create heatmap
  p <- ggplot(data_subset, aes(x = as.factor(Kappa), y = as.factor(Csi), fill = Median_Unique_Stages)) +
    geom_tile() +
    scale_fill_gradient2(name = "Number of Stages", low = "red", mid = "white", high = "blue", midpoint = 0) +
 # scale_fill_viridis(name = "Number of Stages", option = "D")+
    labs(title = paste("Sample Size:", sample_size), x = "Kappa", y = "Csi") +
    theme_minimal()
  
  # Append to list
  heatmap_list[[as.character(sample_size)]] <- p
}

# Arrange heatmaps in a grid
do.call("grid.arrange", c(heatmap_list, ncol = 3))


## HISTOGRAMS CAUSAL EFFECTS

csi_target <- 0.25
kappa_target <- 1
vertical_line <- diff(ce_randomized(tree_def, outcome = "X6", treatment = "X5")[, 2])

# Initialize a list to store histograms
histograms <- list()

# Loop over all sample sizes to generate histograms
for (sample_size in sample_sizes) {
  # Access the specific ce_diffs
  key <- paste0("sample_", sample_size, "_csi_", csi_target, "_kappa_", kappa_target)
  if (!is.null(estimate_VI_results[[key]]$rep_1$ce_diffs)) {
    ce_diffs <- estimate_VI_results[[key]]$rep_1$ce_diffs
    
    # Create a histogram using ggplot2
    p <- ggplot(data = data.frame(ce_diffs = ce_diffs), aes(x = ce_diffs)) +
      geom_histogram(bins = 20, color = "black", fill = "blue", alpha = 0.7) +
      geom_vline(xintercept = vertical_line, color = "red", linetype = "dashed", size = 1) +
      labs(
        title = paste("Sample Size:", sample_size),
        x = "Causal Effect",
        y = "Frequency"
      ) +
      theme_minimal()
    
    # Add the plot to the list
    histograms[[as.character(sample_size)]] <- p
  }
}

# Arrange the histograms in a grid
do.call("grid.arrange", c(histograms, ncol = 2))
