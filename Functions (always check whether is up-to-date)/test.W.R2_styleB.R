# ************************************************************************************************************
test.W.R2 <- function(Y, nb, xy, MEM.autocor = c("positive", "negative"), f = NULL, ...) 
{
  mycall <- pairlist(...)   
  res <- list()  
  MEM.autocor <- match.arg(MEM.autocor)   
  
  if (!(is.null(f))) {
    nbdist <- nbdists(nb, as.matrix(xy))
    if (!(is.null(mycall))) {   
      param <- expand.grid(as.list(mycall))
      m1 <- match(names(param), names(formals(f)))
      for (i in 1:nrow(param)) {
        formals(f)[m1] <- unclass(param[i, ])
        res[[i]] <- scores.listw(nb2listw(nb, style = "B", 
                                          glist = lapply(nbdist, f), zero.policy = TRUE), MEM.autocor = MEM.autocor)
      }
    }
    else {   
      res[[1]] <- scores.listw(nb2listw(nb, style = "B", 
                                        glist = lapply(nbdist, f)), MEM.autocor = MEM.autocor)
    }
  }
  else {   
    res[[1]] <- scores.listw(nb2listw(nb, style = "B"), MEM.autocor = MEM.autocor)
  }
  
  if(length(res) > 1) {
    res2 <- lapply(res, function(z) RsquareAdj(rda(Y, z))$adj.r.squared)
    thebest <- which.max(res2)
    return(list(param = param[thebest, ], best = list(MEM = res[[thebest]])))
  }
  else return(list(best = list(MEM = res[[1]])))
}                                                                                        # End of the function
# ************************************************************************************************************