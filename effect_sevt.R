library("stagedtrees")

randomize <- function(object, var, p = NULL, ignore = object$name_unobserved){
  kk <- length(object$tree[[var]])
  if (is.null(p)){
    p <- rep.int(1/kk, kk)
  }
  names(p) <- object$tree[[var]]
  tmp <- object$stages[[var]]
  object$stages[[var]][!(tmp %in% ignore)] <- "randomized"
  object$prob[[var]] <- c(list(randomized = p), object$prob[[var]][ignore])
  object$prob[[var]] <- object$prob[[var]][!is.na(names(object$prob[[var]]))] 
  return(object)
}


ce_randomized <- function(object, outcome, treatment){
  object0 <- randomize(object, treatment)
  xx <- c(NA)
  names(xx) <- outcome
  res <- sapply(object$tree[[outcome]], function(vo){
    xx[1] <- vo
    prob(object0, xx, conditional_on = as.data.frame(object$tree[treatment]), na0 = FALSE)
  })
  dimnames(res) <- object$tree[c(treatment, outcome)]
  return(res)
}

get_stages_ <- function(x, i, context){
  tree <- x$tree
  ixv <- which(i == names(tree))
  tree1 <- tree[1:(ixv - 1)]
  tree1[names(context)] <- as.list(context)
  more <- rev(expand.grid(rev(tree1), stringsAsFactors = FALSE))
  ixs <- vapply(1:nrow(more), function(rr) {
    cpath <- unlist(more[rr, ])
    stagedtrees:::tree_idx(cpath, tree)
  }, FUN.VALUE = 1)
  return(ixs)
}

ce_propscorestrat <- function(model, outcome, treatment, ignore = model$name_unobserved){
  stgs <- unique(model$stages[[treatment]])
  stgs <- stgs[!(stgs %in% ignore)]
  strata <- sapply(stgs, function(stage){
    get_path(model, treatment, stage)
  }, simplify = FALSE)
  
  tmp <- model$stages[[outcome]]
  kk <- sum(! (tmp %in% ignore))
  model$stages[[outcome]][! (tmp %in% ignore)] <- paste0(1:kk)
  for (i in 1:length(strata)){
    nn <- names(strata)[i]
    for (j in 1:nrow(strata[[i]])){
      todo <- get_stages_(model, outcome, strata[[i]][j,])
      model$stages[[outcome]][todo] <- paste0(nn, 1:length(todo))
    }
  }
  
  model <- sevt_fit(model)
  
  probstrata <- sapply(strata, function(str) sum(prob(model, str)))
  cestrata <- sapply(names(strata), function(ns){
    str <- strata[[ns]]
    todo <- model$stages[[outcome]][get_stages_(model, outcome, str[1,])]
    A <- simplify2array(model$prob[[outcome]][todo])
    colnames(A) <- model$tree$Treatment
    A * probstrata[ns]
  }, simplify = FALSE)
  
  cestrata <- simplify2array(cestrata)
  
  wok <- apply(cestrata, 3, function(xx) !any(is.na(xx)))
  
  t(apply(cestrata[,,wok], c(1,2), sum) / sum(probstrata[wok]))
  
}

