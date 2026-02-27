if (packageVersion("stagedtrees") <= "2.3"){
  # install devel version from github, 
  # at time of writing it is at least 2.3.0.9999 
  # otherwise some functions are not available
  remotes::install_github("stagedtrees/stagedtrees")
}
library("stagedtrees")
library(dplyr)

source("effect_sevt.R")
source("ppmx_hamming.R")
source("ppmx_hamming_Z.R")
source("bayesian_cluster_summary.R")

# create output data
dir.create("results/CTRCD", showWarnings = FALSE, recursive = TRUE)


## data available from (last accessed 27/02/2026)
## https://figshare.com/articles/dataset/BC_cardiotox_A_cardiotoxicity_dataset_for_breast_cancer_patients/22650748/4?file=41888007
data <- read.csv2("data/BC_cardiotox_clinical_variables.csv")

## ANALISI STAGED TREE HAMMING

# preprocess variables
data$past_treat <- as.numeric(data$ACprev == 1 | 
                                data$antiHER2prev == 1 |
                                data$RTprev == 1)

data$smok <- as.integer(data$smoker == 1 | data$exsmoker == 1)

data$CTRCD <- factor(data$CTRCD)
data$AC <- factor(data$AC)
data$HTA <- factor(data$HTA)
data$DL <- factor(data$DL)
data$past_treat <- factor(data$past_treat)
data$smok <- factor(data$smok)
data$BMI <- data$weight / ( (data$height / 100)^2 )

# select variables
d <- data |> dplyr::select(CTRCD,AC,HTA,DL,past_treat)

# remove missing data
d <- na.omit(d)

# create event tree and 
# fit the full probability tree with given order
tree <- stndnaming(full(d, order = c("past_treat","DL","HTA","AC","CTRCD")))

# define parameters for the MCMC
scope <- c("AC","CTRCD")
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

model_hamming <- mcmc_crp_ppmx_Hamming(tree, d, n_save, n_burn = n_burn, thin = thin, a = a, 
                              kappa = kappa, csi = csi,
                              scope = scope, update_SM = update_SM, update_CRP= update_CRP)


## Burn-in and thinning
result <- model_hamming$chain_out

## Check some measure of convergence
x <- lapply(result,logLik)
plot(1:length(x),x,type="l", main = "logLik")

nstages <- function(tree) sum(unlist(lapply(tree$stages,function(i) length(unique(i)))))
y <- lapply(result,nstages)
plot(1:length(y),y,type="l", main = "nstages")


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

par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(estimate_VI, main = "estimate_VI", var_names = FALSE) 
plot(horizontal, main = "horizontal", var_names = FALSE)
plot(lower, main = "lower", var_names = FALSE)
plot(upper, main = "upper", var_names = FALSE)

## compute ATE and ATT for each model in the output of the MCMC
# ATE
ATEs <- sapply(result, function(model) {
   ATE(model, outcome = "CTRCD", treatment = "AC")
   })


write.csv(ATEs, file = "results/CTRCD/ATE_CTRCD.csv")

# ATT
ATTs <- sapply(result, function(model) {
  ATT(model, outcome = "CTRCD", treatment = "AC")
})
write.csv(ATTs, file = "results/CTRCD/ATT_CTRCD.csv")

par(mfrow=c(1,2),mar=c(1,1,1,1))
hist(ATEs, main = "ATEs distribution")
hist(ATTs, main = "ATTs distribution")


## ANALISI COVARIATES
d <- data |> dplyr::select(CTRCD,AC,HTA,DL,past_treat,BMI,age,heart_rate,LVEF)

d <- na.omit(d)

q_Z <- sum(sapply(d, is.numeric))
#Penalty coefficient for continuous variables in PPMx
lambda_Z <- rep(1, q_Z)

d$BMI <- scale(d$BMI)
d$age <- scale(d$age)
d$heart_rate <- scale(d$heart_rate)
d$LVEF <- scale(d$LVEF)

na.omit(d)

model_covariates <- mcmc_crp_ppmx_Hamming_Z(tree, d, n_save, n_burn = n_burn, thin = thin, a = a, 
                                kappa = kappa, csi = csi, lambda_Z = lambda_Z,
                                scope = scope, update_SM = update_SM, update_CRP= update_CRP, 
                                mu0 = mu0, kappa0 = kappa0, alpha0 = alpha0, beta0 = beta0)



## Burn-in and thinning
result <- model_covariates$chain_out

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

## compute ATE and ATT for each model in the output of the MCMC
# ATE
ATEs_cov <- sapply(result, function(model) {
  ATE(model, outcome = "CTRCD", treatment = "AC")
})


write.csv(ATEs_cov, file = "results/CTRCD/ATE_CTRCD_cov.csv")

# ATT
ATTs_cov <- sapply(result, function(model) {
  ATT(model, outcome = "CTRCD", treatment = "AC")
})
write.csv(ATTs_cov, file = "results/CTRCD/ATT_CTRCD_cov.csv")



## OTHER CAUSAL QUANTITIES

# Initialize a list to store the differences for each result
differences_list <- list()

# Loop over all models in result
for (i in seq_along(result)) {
  
  # Fit the staged tree model for this result
  fitted_model <- sevt_fit(result[[i]], d, lambda =0.1)
  
  # Extract the stage assignments and stage probabilities for cesarean
  stages_vec <- result[[i]]$stages$CTRCD
  prob_list <- fitted_model$prob$CTRCD
  
  # Extract the probability of "1" for each stage in order
  prob_vec <- sapply(stages_vec, function(stage) prob_list[[as.character(stage)]]["1"])
  
  # Compute the differences between second and first, fourth and third, etc.
  diffs <- prob_vec[seq(2, length(prob_vec), by = 2)] - prob_vec[seq(1, length(prob_vec), by = 2)]
  
  # Store in list
  differences_list[[i]] <- diffs
}

# Optionally, convert to a matrix or data frame for easy inspection
differences_matrix <- do.call(rbind, differences_list)


apply(differences_matrix,2,mean)
apply(differences_matrix,2,sd)
apply(differences_matrix, 2, function(x) mean(x > 0))
apply(differences_matrix, 2, function(x) mean(x == 0))
apply(differences_matrix, 2, function(x) mean(x < 0))

estimate_VI <- sevt_fit(estimate_VI,data= d, lambda = 0.1)
estimate_VI$stages
estimate_VI$prob



## COVARIATES PROFILES

d <- data |> dplyr::select(CTRCD,AC,HTA,DL,past_treat,BMI,age,heart_rate,LVEF)

d <- na.omit(d)

cont_vars <- c("age", "BMI", "heart_rate", "LVEF")

# Step 1. Summarise by parent profile
prof_summary <- d %>%
  group_by(past_treat, DL, HTA) %>%
  summarise(across(all_of(cont_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x,   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            n = n(),
            .groups = "drop") %>%
  arrange(past_treat, DL, HTA)

# Step 2. Attach stage labels (from estimate_VI$stages$AC)
prof_summary$AC_stage <- c(1,2,3,2,1,2,2,2)  # order must match contexts

# Step 3. Collapse into stage summaries with proper weighted mean + pooled sd
stage_summary <- function(df, cont_vars){
  results <- list()
  for (v in cont_vars){
    mean_col <- paste0(v, "_mean")
    sd_col   <- paste0(v, "_sd")
    n_col    <- "n"
    
    # weighted mean
    w_mean <- weighted.mean(df[[mean_col]], w = df[[n_col]])
    
    # pooled variance (between + within)
    pooled_var <- (sum((df[[n_col]] - 1) * (df[[sd_col]]^2), na.rm = TRUE) +
                     sum(df[[n_col]] * (df[[mean_col]] - w_mean)^2, na.rm = TRUE)) /
      (sum(df[[n_col]]) - 1)
    results[[paste0(v, "_mean")]] <- w_mean
    results[[paste0(v, "_sd")]]   <- sqrt(pooled_var)
  }
  results$n <- sum(df$n)
  as.data.frame(results)
}

summary_AC <- prof_summary %>%
  group_by(AC_stage) %>%
  group_modify(~ stage_summary(.x, cont_vars)) %>%
  ungroup()

summary_AC


# Step 1. Summarise by parent profile (past_treat, DL, HTA, AC)
prof_summary_C <- d %>%
  group_by(past_treat, DL, HTA, AC) %>%
  summarise(across(all_of(cont_vars),
                   list(mean = ~mean(as.numeric(as.character(.x)), na.rm = TRUE),
                        sd   = ~sd(as.numeric(as.character(.x)),   na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            n = n(),
            .groups = "drop") %>%
  arrange(past_treat, DL, HTA, AC)

# Step 2. Attach CTRCD stage labels
prof_summary_C$C_stage <- c(1,2,3,3,4,5,3,3,2,6,3,7,5,3,3,8)

# Step 3. Collapse into stage summaries with weighted mean + pooled SD
stage_summary_covs <- function(df, cont_vars){
  results <- list()
  for (v in cont_vars){
    mean_col <- paste0(v, "_mean")
    sd_col   <- paste0(v, "_sd")
    n_col    <- "n"
    
    # weighted mean
    w_mean <- weighted.mean(df[[mean_col]], w = df[[n_col]])
    
    # pooled variance (within + between)
    pooled_var <- (sum((df[[n_col]] - 1) * (df[[sd_col]]^2), na.rm = TRUE) +
                     sum(df[[n_col]] * (df[[mean_col]] - w_mean)^2, na.rm = TRUE)) /
      (sum(df[[n_col]]) - 1)
    
    results[[paste0(v, "_mean")]] <- w_mean
    results[[paste0(v, "_sd")]]   <- sqrt(pooled_var)
  }
  results$n <- sum(df$n)
  as.data.frame(results)
}

summary_CTRCD <- prof_summary_C %>%
  group_by(C_stage) %>%
  group_modify(~ stage_summary_covs(.x, cont_vars)) %>%
  ungroup()

summary_CTRCD


