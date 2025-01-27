library(stagedtrees)
library("e1071")
library(salso)
library(fossil)

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
tree_def$prob$X2 <- list("1" = c("no"= 0.5,"yes" = 0.5),
                         "2" = c("no"= 0.4,"yes" = 0.6))
tree_def$prob$X3 <- list("1" = c("no"= 0.4,"yes" = 0.6),
                         "2" = c("no"= 0.5,"yes" = 0.5),
                         "3" = c("no"= 0.7,"yes" = 0.3))
tree_def$prob$X4 <- list("1" = c("no" = 0.7, "yes" = 0.3),
                         "2" = c("no" = 0.3, "yes" = 0.7),
                         "3" = c("no" = 0.4, "yes" = 0.6))



sample_size <- 1000
scope <- c("X3","X4")
a <- 1
n_burn <- 1000
thin <- 5
n_save <- 5000
csi <- 0.25
kappa <- 1
mu0 <- 0
kappa0 <- 1
alpha0 <- 1
beta0 <- 1

update_CRP <- T
update_SM <- T

data <- sample_from(tree_def, size = sample_size)
tree <- stndnaming(stages_fbhc(full(data)))

# Ensure data is valid (no empty counts for X6)
while (any(tree$ctables$X4[, 1] == 0 & tree$ctables$X4[, 2] == 0)) {
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
}

  
data$Z1 <- scale(rnorm(sample_size,40, sd= 10))
data$Z2 <- 0
data$Z2[which(data$X2=="no")] <- rnorm(length(which(data$X2=="no")),20,sd =3)
data$Z2[which(data$X1=="yes" & data$X2 == "yes")] <- rnorm(length(which(data$X1=="yes" & data$X2 == "yes")),60,sd =3)
data$Z2[which(data$X1=="no" & data$X2 == "yes")] <- rnorm(length(which(data$X1=="no" & data$X2 == "yes")),100,sd =3)
data$Z2 <- scale(data$Z2)


mcmc_crp_ppmx(tree, data, n_save, n_burn = n_burn, thin = thin, a = a, kappa = kappa, csi = csi, scope = scope, update_SM = TRUE, update_CRP= TRUE)

