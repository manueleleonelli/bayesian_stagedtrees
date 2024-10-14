## Definisco l'albero
tree_def <- list(X1 = c("no","yes"),
             X2 = c("no","yes"),
             X3 = c("no","yes"),
             X4 = c("no","yes"),
             X5 = c("no","yes"))

## Lo creo come oggetto per pacchetto stagedtrees
tree_def <- sevt(tree_def, full= F)

## Creo partizione per X5 (dove andiamo a fare clustering)
tree_def$stages$X5 <- c("1","2","3","4","1","2","5","6","1","2","3","4","1","2","5","6")

## Assegno probabilita all'albero
## Per tutte le variabili precedenti metto 50/50 cosi da avere osservazioni in tutti i paths
tree_def$prob <- list()
tree_def$prob$X1 <- list("NA"=c(0.5,0.5))
tree_def$prob$X2 <- list("1" = c(0.5,0.5))
tree_def$prob$X3 <- list("1" = c(0.5,0.5))
tree_def$prob$X4 <- list("1"=c(0.5,0.5))
tree_def$prob$X5 <- list("1" = c(0.2,0.8),
                     "2" = c(0.3,0.7),
                     "3" = c(0.5,0.5),
                     "4" = c(0.7,0.3),
                     "5" = c(0.8,0.2),
                     "6" = c(0.9,0.1))

## Simulo i dati dall'albero
data <- sample_from(tree_def, size = 2000, seed=2024)

## Costruisco un albero iniziale se vuoi utilizzarlo
tree <- stages_hc(full(data))
