library(stagedtrees)
library(salso)

tree_def <- list(X1 = c("no","yes"),
                 X2 = c("no","yes"),
                 X3 = c("no","yes"),
                 X4 = c("no","yes"))

tree_def <- sevt(tree_def, full= F)

tree_def$stages$X2 <- c("1","2")
tree_def$stages$X3 <- c("1","2","1","3")
tree_def$stages$X4 <- c("1","1","2","2","1","1","3","3")

tree_def <- stndnaming(tree_def)                        
tree_def$prob <- list()

tree_def$prob$X1 <- list("NA"=c("no"= 0.6,"yes" = 0.4))
tree_def$prob$X2 <- list("1" = c("no"= 0.6,"yes" = 0.4),
                         "2" = c("no"= 0.4,"yes" = 0.6))

tree_def$prob$X3 <- list("1" = c("no"= 0.4,"yes" = 0.6),
                         "2" = c("no"= 0.5,"yes" = 0.5),
                         "3" = c("no"= 0.7,"yes" = 0.3))

tree_def$prob$X4 <- list("1" = c("no"= 0.7,"yes" = 0.3),
                         "2" = c("no"= 0.3,"yes" = 0.7),
                         "3" = c("no"= 0.4,"yes" = 0.6))


