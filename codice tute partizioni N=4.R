parts <- list(c(1,2,3,4),c(1,1,2,3),c(1,2,1,3),c(1,2,3,1),c(1,2,2,3),c(1,2,3,2),c(1,2,3,3),c(1,1,2,2),c(1,2,1,2),c(1,2,2,1),c(1,1,1,2),c(1,1,2,1),c(1,2,1,1),c(1,2,2,2),c(1,1,1,1))
x <- c(1,2,3,4)

dist <- matrix(0,4,4)
dist[1,1] <- dist[2,2] <- dist[3,3] <- dist[4,4] <- 0
dist[1,2] <- dist[2,1] <- dist[3,4] <- dist[4,3] <- dist[1,3] <- dist[3,1] <- dist[4,2] <- dist[2,4] <- 0.5
dist[1,4] <- dist[4,1] <- dist[2,3] <- dist[3,2] <- 1

prior <- function(x,kappa,xi,dist){
  M <- length(unique(x))
  nj <- table(x)
  dp <- kappa^M*prod(gamma(nj))
  expo <- 0
  for(i in 1:M){
   expo <- expo + sum(dist[x==i,x==i])
  }
  dp*exp(-xi*expo)
}

round(sapply(parts,prior, kappa=0.5,xi=0.2,dist=dist)/sum(sapply(parts,prior, kappa=1,xi=0.1,dist=dist)),3)
