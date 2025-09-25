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


# Compute 16 CATEs by taking, for each covariate profile,
# P(Y="yes" | R=yes) - P(Y="yes" | R=no) from consecutive X6 vertices.
cate16_pairs <- function(tree,
                         outcome = "X6",
                         treatment = "X5",
                         covariates = c("X1","X2","X3","X4"),
                         outcome_yes = "yes",
                         treated_level = "yes") {
  st <- tree$stages[[outcome]]
  pr <- tree$prob[[outcome]]
  n <- length(st)
  if (n %% 2 != 0) stop("Outcome stage vector length is not even.")
  
  # Which treatment level is first/second along the path? default assumes c("no","yes")
  tr_levels <- tryCatch(tree$ll[[treatment]], error = function(e) NULL)
  if (is.null(tr_levels)) tr_levels <- c("no","yes")
  treated_is_second <- match(treated_level, tr_levels) == 2
  
  get_yes <- function(stage_id) {
    p <- pr[[as.character(stage_id)]]
    y <- p[[outcome_yes]]
    if (is.null(y)) stop(sprintf("Outcome level '%s' not found in stage %s.", outcome_yes, stage_id))
    unname(y)
  }
  
  idx <- matrix(seq_len(n), ncol = 2, byrow = TRUE)  # (1,2), (3,4), ...
  cate <- apply(idx, 1, function(ix) {
    y1 <- get_yes(st[ix[1]])  # R = first level in tree
    y2 <- get_yes(st[ix[2]])  # R = second level in tree
    if (isTRUE(treated_is_second)) y2 - y1 else y1 - y2
  })
  
  # Optional labels for X1..X4 (fallback to c("no","yes") if levels not stored)
  levs <- lapply(covariates, function(v) {
    lv <- tryCatch(tree$ll[[v]], error = function(e) NULL)
    if (is.null(lv)) c("no","yes") else lv
  })
  cov_df <- expand.grid(levs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(cov_df) <- covariates
  
  out <- cbind(cov_df, CATE = as.numeric(cate))
  rownames(out) <- NULL
  out
}


sample_size <- 500
scope <- c("X5","X6")
n_v <- length(scope)
a <- 1
n_burn <- 1000
thin <- 2
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
while (any(tree$ctables$X6[, 1] == 0 & tree$ctables$X6[, 2] == 0)) {
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
}


run_ppmx_pipeline <- function(tree, data,
                              n_save,
                              n_burn,
                              thin,
                              a,
                              kappa,
                              csi,
                              scope,
                              update_SM = TRUE,
                              update_CRP = TRUE,
                              lambda_fit = 0.1,
                              alpha = 0.1,
                              do_plots = TRUE) {
  #--- helpers ---------------------------------------------------------------
  nstages_fun <- function(tr) sum(unlist(lapply(tr$stages, function(i) length(unique(i)))))

  
  #--- 1) Run MCMC -----------------------------------------------------------
  ciao <- mcmc_crp_ppmx_Hamming(
    tree, data, n_save,
    n_burn = n_burn, thin = thin, a = a,
    kappa = kappa, csi = csi,
    scope = scope, update_SM = update_SM, update_CRP = update_CRP
  )
  
  # Pull chain output
  result <- ciao$chain_out
  
  #--- 2) Simple diagnostics (optional plots) --------------------------------
  # log-lik over draws (coerce safely to numeric)
  x <- sapply(result, function(tr) {
    ll <- tryCatch(logLik(tr), error = function(e) NA_real_)
    as.numeric(ll)
  })
  
  # #stages over draws
  y <- sapply(result, nstages_fun)
  
  if (isTRUE(do_plots)) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op), add = TRUE)
    plot(seq_along(x), x, type = "l", xlab = "Iteration", ylab = "logLik", main = "MCMC logLik")
    plot(seq_along(y), y, type = "l", xlab = "Iteration", ylab = "# stages (total)", main = "Complexity trace")
  }
  
  #--- 3) Standardize stage naming across draws ------------------------------
  result <- lapply(result, stndnaming)
  
  #--- 4) Build estimates per variable ---------------------------------------
  # initialize holders as a copy of the first draw's structure
  estimate_VI <- estimate_B <- lower <- upper <- horizontal <- result[[1]]
  
  for (v in scope) {
    out <- process_variable(v, result, estimate_B, estimate_VI, lower, upper, horizontal, alpha = alpha)
    estimate_B  <- out$estimate_B
    estimate_VI <- out$estimate_VI
    lower       <- out$lower
    upper       <- out$upper
    horizontal  <- out$horizontal
  }
  
  #--- 5) Optional refit with stage estimates (smooth counts) ----------------
  estimate_VI_fit <- sevt_fit(estimate_VI, data, lambda = lambda_fit)
  estimate_B_fit  <- sevt_fit(estimate_B,  data, lambda = lambda_fit)
  
  #--- 6) Return everything useful ------------------------------------------
  list(
    mcmc = list(
      draws = result,
      logLik = x,
      nstages = y
    ),
    estimates = list(
      VI = estimate_VI,
      Binder = estimate_B,
      VI_fit = estimate_VI_fit,
      Binder_fit = estimate_B_fit
    ),
    credible_balls = list(
      lower = lower,
      upper = upper,
      horizontal = horizontal,
      alpha = alpha
    )
  )
}


## ---------- helpers to evaluate ATE/CATE errors ----------
eval_effect_errors <- function(reference_tree,
                               estimated_tree,
                               outcome = "X6",
                               treatment = "X5",
                               cate_fun = cate16_pairs,
                               cate_name = "CATE") {
  # ATE abs error
  ref_ate <- diff(ce_randomized(reference_tree, outcome = outcome, treatment = treatment)[, 2])
  est_ate <- diff(ce_randomized(estimated_tree,  outcome = outcome, treatment = treatment)[, 2])
  ate_abs_err <- abs(est_ate - ref_ate)
  
  # CATE abs errors: join on covariate columns to be safe
  ref_c <- as.data.frame(cate_fun(reference_tree, outcome = outcome, treatment = treatment))
  est_c <- as.data.frame(cate_fun(estimated_tree,  outcome = outcome, treatment = treatment))
  
  # keys = all columns except the CATE column
  cate_col <- which(names(ref_c) == cate_name)
  key_cols <- setdiff(names(ref_c), cate_name)
  
  merged <- merge(ref_c, est_c, by = key_cols, suffixes = c(".ref", ".est"))
  cate_abs_err <- abs(merged[[paste0(cate_name, ".est")]] - merged[[paste0(cate_name, ".ref")]])
  
  list(
    ATE_abs_error = as.numeric(ate_abs_err),
    CATE_abs_errors = as.numeric(cate_abs_err),
    CATE_median_abs = median(as.numeric(cate_abs_err), na.rm = TRUE),
    CATE_mean_abs   = mean(as.numeric(cate_abs_err), na.rm = TRUE)
  )
}

## ---------- simulate data until valid (no empty X6 contexts) ----------
simulate_valid_data <- function(tree_def, n, max_tries = 200) {
  tries <- 0
  repeat {
    tries <- tries + 1
    data <- sample_from(tree_def, size = n)
    tr   <- tryCatch(stndnaming(stages_fbhc(full(data))), error = function(e) NULL)
    if (!is.null(tr)) {
      # validity: no row in X6 contingency table has both 0s
      ok <- !any(tr$ctables$X6[, 1] == 0 & tr$ctables$X6[, 2] == 0)
      if (isTRUE(ok)) return(list(data = data, start_tree = tr))
    }
    if (tries >= max_tries) stop("simulate_valid_data: could not build a valid starting tree in max_tries.")
  }
}

## --- replicate wrapper: keep full CATEs too ---
do_one_replicate <- function(n,
                             tree_def,
                             scope = c("X5","X6"),
                             n_save = 1000, n_burn = 1000, thin = 2,
                             a = 1, kappa = 1, csi = 0.25,
                             update_SM = TRUE, update_CRP = TRUE,
                             lambda_fit = 0.1, alpha = 0.05,
                             seed = NULL, do_plots = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  sim <- simulate_valid_data(tree_def, n)
  data <- sim$data
  start_tree <- sim$start_tree
  
  out <- run_ppmx_pipeline(
    tree = start_tree,
    data = data,
    n_save = n_save,
    n_burn = n_burn,
    thin = thin,
    a = a,
    kappa = kappa,
    csi = csi,
    scope = scope,
    update_SM = update_SM,
    update_CRP = update_CRP,
    lambda_fit = lambda_fit,
    alpha = alpha,
    do_plots = do_plots
  )
  
  est_tree <- out$estimates$VI_fit
  
  errs <- eval_effect_errors(
    reference_tree = tree_def,
    estimated_tree = est_tree,
    outcome   = "X6",
    treatment = "X5"
  )
  
  # return both aggregate metrics + the full CATE vector
  list(
    ATE_abs_error   = errs$ATE_abs_error,
    CATE_median_abs = errs$CATE_median_abs,
    CATE_mean_abs   = errs$CATE_mean_abs,
    CATE_errors     = errs$CATE_abs_errors  # length 16
  )
}

## --- main simulation loop: keep all 16 CATEs ---
sizes <- c(500, 1000, 2500, 5000, 10000)
n_reps <- 25

# tidy long format with 16 CATE errors stored
all_results <- do.call(
  rbind,
  lapply(sizes, function(n) {
    do.call(rbind, lapply(seq_len(n_reps), function(r) {
      seed <- 1e6 + n*100 + r
      cat("Running n =", n, "rep =", r, "\n")
      o <- do_one_replicate(
        n = n, tree_def = tree_def, scope = c("X5","X6"),
        seed = seed, do_plots = FALSE
      )
      data.frame(
        n = n,
        rep = r,
        cate_id = seq_along(o$CATE_errors),         # 1..16
        CATE_abs_error = o$CATE_errors,
        ATE_abs_error  = o$ATE_abs_error,           # repeat per row
        CATE_median_abs = o$CATE_median_abs,        # repeat per row
        CATE_mean_abs   = o$CATE_mean_abs           # repeat per row
      )
    }))
  })
)


## Example usage:
head(all_results)


# quick glance:
# boxplot(CATE_abs_error ~ as.factor(n), data = all_cate_raw,
#         xlab = "n", ylab = "CATE absolute error", main = "CATE error by sample size")


out <- run_ppmx_pipeline(
  tree, data,
  n_save = 1000,
  n_burn = 1000,
  thin = 2,
  a = a,
  kappa = kappa,
  csi = csi,
  scope = scope,
  update_SM = TRUE,
  update_CRP = TRUE,
  lambda_fit = 0.1,
  alpha = 0.05,
  do_plots = F
)

# Access results
out$estimates$VI_fit

ATE <- diff(ce_randomized(tree_def, outcome = "X6", treatment = "X5")[, 2])
CATEs <- cate16_pairs(tree_def)
abs(diff(ce_randomized(out$estimates$VI_fit, outcome = "X6", treatment = "X5")[, 2])- ATE)
abs(cate16_pairs(out$estimates$VI_fit)[,5]-CATEs[,5])
