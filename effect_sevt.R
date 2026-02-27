if (packageVersion("stagedtrees") <= "2.3"){
  # install devel version from github, 
  # at time of writing it is at least 2.3.0.9999 
  # otherwise some functions are not available
  remotes::install_github("stagedtrees/stagedtrees")
}
library("stagedtrees")

# we directly use predict.sevt 
# which can handle now marginalization

#' function for ATE
#'
#' @param object a fitted staged event tree
#' @param outcome the name of the outcome variable
#' @param treatment the name of the treatment variable 
#' @returns The Average treatmen effect estimated 
#'          from the model in \code{object}.
#' @details
#' This function assume that the second level 
#' of the treament variable is the treated and that the 
#' second level of the outcome variable 
#' is the positive outcome (e.g. survived)
#' @examples
#' model <- Titanic |> full(join_unobserved = FALSE, lambda = 1) |> stages_bhc()
#' outcome <- "Survived"
#' treatment <- "Sex"
#' ATE(model, outcome, treatment)
ATE <- function(object, outcome, treatment){
  object0 <- randomize_sevt(object, treatment)
  diff(predict(object = object0, newdata = data.frame(object$tree[treatment]),
               class = outcome, prob = TRUE)[,2])
}

########## for treated sub population
## this return the probabilities of the outcomes under the
## two possible treatments for units in the treated group
ATT <- function(object, outcome, treatment){
  ## we use as x pre treatment variables affecting the treatment
  parents <- as_parentslist(object, silent = TRUE)[[treatment]]$parents
  ct <- list(object$tree[[treatment]][2]) ## second level is treated
  names(ct) <- treatment
  object0 <- randomize_sevt(object, treatment)
  if (length(parents) == 0 ){
    ## no variable affect treatment we just default to ATE
    return(ATE(object, outcome, treatment))
    } else {
    xx <- expand.grid(object$tree[parents])
    res <- lapply(seq_len(nrow(xx)), function(ii){
      x <- xx[ii,,drop=FALSE]
      predict(object0, class = outcome,
              newdata = data.frame(x, as.data.frame(object$tree[treatment]), row.names = NULL),
              prob = TRUE) * prob(object, x, conditional_on = ct)
    })
    diff(Reduce("+", x = res, accumulate = FALSE)[,2])
  }
}


