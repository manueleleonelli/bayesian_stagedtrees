library(stagedtrees)
library(e1071)
library(salso)
library(mcclust.ext)
library(pcalg)
library(bnlearn)
library(BiDAG)

source("effect_sevt.R")
source("ppmx_hamming.R")
source("bayesian_cluster_summary.R")


## SIMULAZIONE DATA
# Tunables
sample_size   <- 1000
q             <- 0.7
max_resamples <- 20     # if we hit this in the inner while, restart everything
max_restarts  <- 50     # safety cap; set to Inf if you truly never want to stop

restart <- 0
repeat {
  ## (Re)build everything
 # model_gt <- random_sevt(list( X1 = c("y","n"), X2 = c("y","n"),X3 = c("y","n"), X4 = c("y","n"), X5 = c("y","n"), X6 = c("y","n") ), q = q)
  model_gt <- random_sevt(random_parentslist(6,maxp=5),q=q,rfun=rexp)
  ate_true <- diff(ce_randomized(model_gt, outcome = "X6", treatment = "X5")[, 2])
  
  model_gt <- stndnaming(model_gt)
  data <- sample_from(model_gt, sample_size)
  tree <- full(data)
  
  ## Try to fix zero rows, but give up and restart after max_resamples tries
  resamples <- 0
  while (any(tree$ctables$X6[, 1] == 0 & tree$ctables$X6[, 2] == 0)) {
    data <- sample_from(model_gt, size = sample_size)
    tree <- stndnaming(stages_fbhc(full(data)))
    resamples <- resamples + 1
    if (resamples >= max_resamples) {
      restart <- restart + 1
      message(sprintf(
        "Hit %d resamples without success; restarting everything (restart %d).",
        max_resamples, restart
      ))
      break
    }
  }
  
  ## Success if we exited the while before hitting the limit
  if (resamples < max_resamples) break
  if (restart >= max_restarts) stop("Gave up after ", max_restarts, " full restarts.")
}

# At this point: model_gt, data, tree are acceptable; ate_true matches this model.

### ESTIMATION TREE BHC
tree_est <- stndnaming(stages_bhc(full(data)))
 abs(diff(ce_randomized(tree_est, outcome = "X6", treatment = "X5")[, 2])
- ate_true)

for(i in 1:4){
  tree_est$stages[[i]] <- model_gt$stages[[i]]
}

hamming_stages(tree_est,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(tree_est$stages$X6))


### ESTIMATION TREE CSBHC
tree_est <-stndnaming(stages_csbhc(full(data)))
abs(diff(ce_randomized(tree_est, outcome = "X6", treatment = "X5")[, 2])
    - ate_true)
for(i in 1:4){
  tree_est$stages[[i]] <- model_gt$stages[[i]]
}
hamming_stages(tree_est,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(tree_est$stages$X6))

### ESTIMATION TREE BAYESIAN
scope <- c("X5","X6")
n_v <- length(scope)
a <- 1
n_burn <- 500
thin <- 5
n_save <- 1000
csi <- 0.25
kappa <- 1
mu0 <- 0
kappa0 <- 1
alpha0 <- 1
beta0 <- 1
update_SM = TRUE
update_CRP= TRUE

tree_est <- mcmc_crp_ppmx_Hamming(
  tree, data, n_save,
  n_burn = n_burn, thin = thin, a = a,
  kappa = kappa, csi = csi,
  scope = scope, update_SM = update_SM, update_CRP = update_CRP
)

result <- tree_est$chain_out 
ce_diffs <- sapply(result, function(model) {
  diff(ce_randomized(sevt_fit(model, data,lambda=0.1), outcome = "X6", treatment = "X5")[, 2])
})

abs(ate_true - mean(ce_diffs))

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

for(i in 1:4){
  estimate_VI$stages[[i]] <- model_gt$stages[[i]]
}
hamming_stages(estimate_VI,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(estimate_VI$stages$X6))

## ESTIMATION TABU

dag <- tabu(data,blacklist = ordering2blacklist(colnames(data)))
amat(dag)[5,6] <- 1
tree_tabu <- sevt_fit(as_sevt(bn.fit(dag,data),order=colnames(data)),data,lambda=0.1)

abs(diff(ce_randomized(tree_tabu, outcome = "X6", treatment = "X5")[, 2])
    - ate_true)

for(i in 1:4){
  tree_tabu$stages[[i]] <- model_gt$stages[[i]]
}
hamming_stages(tree_tabu,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(tree_tabu$stages$X6))

## ESTIMATION PC

dag <- pc.stable(data,blacklist = ordering2blacklist(colnames(data)),alpha=0.01)
amat(dag)[5,6] <- 1
tree_pc <- sevt_fit(as_sevt(bn.fit(dag,data),order=colnames(data)),data,lambda=0.1)

abs(diff(ce_randomized(tree_pc, outcome = "X6", treatment = "X5")[, 2])
    - ate_true)

for(i in 1:4){
  tree_pc$stages[[i]] <- model_gt$stages[[i]]
}
hamming_stages(tree_pc,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(tree_pc$stages$X6))


## ESTIMATION DAG BAYESIAN

pairs_to_blacklist_matrix <- function(pairs, nodes, sparse = FALSE) {
  if (!all(c("from","to") %in% tolower(colnames(pairs)))) {
    colnames(pairs) <- c("from","to")  # ensure names
  }
  from <- pairs[,1]; to <- pairs[,2]
  # index lookup
  i <- match(from, nodes); j <- match(to, nodes)
  if (any(is.na(i) | is.na(j))) stop("pairs contain names not in 'nodes'.")
  
  if (sparse) {
    if (!requireNamespace("Matrix", quietly = TRUE))
      stop("Install 'Matrix' for sparse output or use sparse=FALSE.")
    Matrix::sparseMatrix(i = i, j = j, x = 1L,
                         dims = c(length(nodes), length(nodes)),
                         dimnames = list(nodes, nodes))
  } else {
    M <- matrix(0L, length(nodes), length(nodes),
                dimnames = list(nodes, nodes))
    M[cbind(i, j)] <- 1L
    M
  }
}

# --- example with your objects ---
nodes <- colnames(data)
pairs <- ordering2blacklist(nodes)   # your from/to matrix

bl <- pairs_to_blacklist_matrix(pairs, nodes, sparse = FALSE)

scores <- scoreparameters(scoretype = "bdecat", data)
sample <- sampleBN(scores, algorithm = c("partition"), blacklist=bl)

dag <- empty.graph(colnames(data))
amat(dag) <- modelp(sample,0.5)
amat(dag)[5,6] <- 1

tree_bayes <- sevt_fit(as_sevt(bn.fit(dag,data),order=colnames(data)),data,lambda=0.01)

abs(diff(ce_randomized(tree_bayes, outcome = "X6", treatment = "X5")[, 2])
    - ate_true)

for(i in 1:4){
  tree_bayes$stages[[i]] <- model_gt$stages[[i]]
}
hamming_stages(tree_bayes,model_gt)
rand.index(as.numeric(model_gt$stages$X6),as.numeric(tree_bayes$stages$X6))

library(bnlearn)
library(stagedtrees)  # for as_sevt, sevt_fit, ce_randomized
# library(Matrix)     # if you want to be explicit with sparse matrices

nodes    <- colnames(data)
treat    <- "X5"
outcome  <- "X6"
lambda   <- 0.01      # smoothing for sevt_fit
force_RX <- TRUE      # force X5 -> X6 like in your example

inc_list <- sample$traceadd$incidence  # list of adjacency matrices (sparse)

# helper: turn a (possibly sparse) adjacency into a bnlearn DAG
adj_to_dag <- function(A, nodes, force_edge = FALSE, from = treat, to = outcome) {
  A <- as.matrix(A) * 1L                    # ensure ordinary 0/1 matrix
  A <- A[nodes, nodes, drop = FALSE]        # match node order
  diag(A) <- 0L
  if (force_edge) A[from, to] <- 1L         # optionally enforce treat -> outcome
  dag <- empty.graph(nodes)
  amat(dag) <- A
  dag
}

# loop over all sampled DAGs, compute ATE for each; skip invalid ones safely
get_ate <- function(A) {
  dag <- try(adj_to_dag(A, nodes, force_edge = force_RX), silent = TRUE)
  if (inherits(dag, "try-error")) return(NA_real_)
  fit <- try(bn.fit(dag, data), silent = TRUE)   # parameter fit on your data
  if (inherits(fit, "try-error")) return(NA_real_)
  st  <- try(sevt_fit(as_sevt(fit, order = nodes), data, lambda = lambda), silent = TRUE)
  if (inherits(st, "try-error")) return(NA_real_)
  ce  <- ce_randomized(st, outcome = outcome, treatment = treat)
  diff(ce[, 2])  # same contrast you used for ATE
}

# compute ATEs with a simple progress bar
pb <- txtProgressBar(min = 0, max = length(inc_list), style = 3)
ates <- vapply(seq_along(inc_list), function(i) {
  res <- get_ate(inc_list[[i]])
  setTxtProgressBar(pb, i)
  res
}, numeric(1L))
close(pb)

# clean and summarize
ates <- ates[is.finite(ates)]
summary_ates <- c(
  mean = mean(ates),
  sd   = sd(ates),
  quantile(ates, c(0.025, 0.25, 0.5, 0.75, 0.975))
)


# if you want error vs. ground truth for each sampled DAG:
abs(ates - ate_true)


