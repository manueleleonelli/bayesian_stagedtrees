library(stagedtrees)
library("e1071")
library(salso)
library(fossil)

source("Gamma_Only_With_Hamming.R")
source("sample_partition_function.R")
source("Gamma_Only_No_Hamming.R")
source("Only_Hamming.R")

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
scope <- c("X4")
a <- 1
n_burn <- 1000
thin <- 5
n_save <- 5000
csi <- 0.2
kappa <- 1

bb <- list(c(1,1))
cc <- list(c(1,1))
ff <- list(c(0.2^2/0.2,0.2^2/0.2))
gg <- list(c(0.2/0.2,0.2/0.2))

update_CONT <-  T
update_CRP <- T
update_SM <- T

data <- sample_from(tree_def, size = sample_size)
tree <- stndnaming(stages_fbhc(full(data)))

# Ensure data is valid (no empty counts for X6)
while (any(tree$ctables$X4[, 1] == 0 & tree$ctables$X4[, 2] == 0)) {
  data <- sample_from(tree_def, size = sample_size)
  tree <- stndnaming(stages_fbhc(full(data)))
}

data$Z1 <- rnorm(sample_size,40, sd= 10) 
data$Z1 <- (data$Z1 - min(data$Z1))/(max(data$Z1)-min(data$Z1))
data$Z2 <- 0
data$Z2[which(data$X2=="no")] <- rnorm(length(which(data$X2=="no")),20,sd =3)
data$Z2[which(data$X1=="yes" & data$X2 == "yes")] <- rnorm(length(which(data$X1=="yes" & data$X2 == "yes")),60,sd =3)
data$Z2[which(data$X1=="no" & data$X2 == "no")] <- rnorm(length(which(data$X1=="no" & data$X2 == "no")),100,sd =3)
data$Z2 <- (data$Z2 - min(data$Z2))/(max(data$Z2)-min(data$Z2))


d <- data[,c("X1","X2","X3","X4")]
tree <- stndnaming(stages_hc(full(d)))

HO <- mcmc_crp_distances(tree,n_save =n_save,n_burn= n_burn,thin =thin,scope=scope)
GNH <- mcmc_crp_gamma_no_hamming(tree, data, n_save, n_burn = n_burn, thin = thin, scope=scope)
GH <- mcmc_crp_gamma_hamming(tree, data, n_save, n_burn = n_burn, thin = thin, scope=scope)


