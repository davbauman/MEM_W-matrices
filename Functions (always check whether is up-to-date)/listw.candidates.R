# Trouver comment gérer la question des DBMEM et des threshold distances

# Pour les listes dans la liste (f2 et f3), la valeur du 'y' est le nom (names) de 
# l'élément liste (v. listw.candidates)

# *********************************************************************************** #
listw.candidates <- function (coord, style = "W", del = TRUE, gab = TRUE, rel = TRUE, 
                             mst = TRUE, PCNM = TRUE,  DB = TRUE, thresh = "?", 
                             binary = TRUE, flin = TRUE, fconcdown = TRUE, 
                             fconcup = TRUE, y_fconcdown = 5, y_fconcup = 0.5)
{
  
  if (any(is.na(coord))) stop("NA entries in coord")
  
  if (length(which(c(flin, fconcdown, fconcup) == TRUE)) != 0) weightfun = TRUE
  else weightfun = FALSE
  
  # Definition of the weighting functions:
  if (weightfun == TRUE) {
    f1 <- function (D, dmax)     { 1 - (D/dmax) }       # Linear function
    f2 <- function (D, dmax, y)  { 1 - (D/dmax)^y }     # Concave-down function
    f3 <- function (D, y)        { 1 / D^y }            # Concave-up function
  }
  
  xy.d1 <- dist(coord)
  
  # Total nb of W matrices to build:
  nbB <- length(which(c(del, gab, rel, mst) == TRUE))  # AJOUTER DB + nb thresh dist.
  nbw <- nbB
  control_BinLin <- FALSE
  control_f2 <- FALSE
  if (length(which(c(binary, flin) == TRUE)) != 0) {
    control_BinLin <- TRUE
    nbw <- nbB * length(which(c(binary, flin) == TRUE))
  }
  if (weightfun == TRUE) {
    if (fconcdown == TRUE) {
      control_f2 <- TRUE
      yf2 <- length(y_fconcdown)
      if (control_BinLin == TRUE) nbw <- nbw + (nbB * yf2) 
      else nbw <- nbB * yf2
    }
    if (fconcup == TRUE) {
      yf3 <- length(y_fconcup)
      if (control_BinLin == TRUE) nbw <- nbw + (nbB * yf3) 
      else if (control_f2 == TRUE) nbw <- nbw + nbB * yf3 else nbw <- nbB * yf3
    }
  }
  if (PCNM == TRUE) nbw <- nbw + 1
  
  # List of the W matrix candidates
  listwcand <- vector("list", nbw)
  count <- 0

  # Construction of the list of W matrix candidates: 
  # ************************************************
  if (PCNM == TRUE) {
    count <- count + 1
    f <- function (D, t) { 1-(D/(4*t))^2 }           # PCNM criterion
    lowlim <- give.thresh(xy.d1)
    matB <- dnearneigh(lowlim, x = as.matrix(coord), d1 = 0)
    listwcand[[count]] <- nb2listw(matB, style = style, 
                                   glist = lapply(nbdists(matB, as.matrix(coord)), f, 
                                                  t = lowlim))
    names(listwcand)[count] <- "DBMEM_PCNM"
  }
  if (del == TRUE) {
    Y.del <- tri2nb(coord)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.del, style = style)
      names(listwcand)[count] <- "Delaunay_Binary weighting"
    }
    if (weightfun == TRUE) {
      max.del <- max(unlist(nbdists(Y.del, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.del, style = style, 
                                       glist = lapply(nbdists(Y.del, as.matrix(coord)), 
                                                      f1, dmax = max.del))
        names(listwcand)[count] <- "Delaunay_Linear weighting"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcdown))
        count_list <- 1
        for (i in y_fconcdown) {
          listf[[count_list]] <- nb2listw(Y.del, style = style, 
                                       glist = lapply(nbdists(Y.del, as.matrix(coord)),
                                                      f2, y = i, dmax = max.del))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcdown
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Delaunay_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcup))
        count_list <- 1
        for (i in y_fconcup) {
          listf[[count_list]] <- nb2listw(Y.del, style = style, 
                                 glist = lapply(nbdists(Y.del, as.matrix(coord)),
                                                f3, y = i))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcup
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Delaunay_Concave up weighting"
      }
    }
  }
  if (gab == TRUE) {
    Y.gab <- graph2nb(gabrielneigh(as.matrix(coord), nnmult = 5), sym = TRUE)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.gab, style = style)
      names(listwcand)[count] <- "Gabriel_Binary weighting"
    }
    if (weightfun == TRUE) {
      max.gab <- max(unlist(nbdists(Y.gab, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.gab, style = style, 
                                       glist = lapply(nbdists(Y.gab, as.matrix(coord)), 
                                                      f1, dmax = max.gab))
        names(listwcand)[count] <- "Gabriel_Linear weighting"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcdown))
        count_list <- 1
        for (i in y_fconcdown) {
          listf[[count_list]] <- nb2listw(Y.gab, style = style, 
                                 glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                f2, y = i, dmax = max.gab))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcdown
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Gabriel_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcup))
        count_list <- 1
        for (i in y_fconcup) {
          listf[[count_list]] <- nb2listw(Y.gab, style = style, 
                                 glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                f3, y = i))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcup
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Gabriel_Concave up weighting"
      }
    }
  }
  if (rel == TRUE) {
    Y.rel <- graph2nb(relativeneigh(as.matrix(coord), nnmult = 5), sym = TRUE)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.rel, style = style)
      names(listwcand)[count] <- "Relative neighbourhood_Binary weighting"
    }
    if (weightfun == TRUE) {
      max.rel <- max(unlist(nbdists(Y.rel, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.rel, style = style, 
                                       glist = lapply(nbdists(Y.rel, as.matrix(coord)), 
                                                      f1, dmax = max.rel))
        names(listwcand)[count] <- "Relative neighbourhood_Linear weighting"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcdown))
        count_list <- 1
        for (i in y_fconcdown) {
          listf[[count_list]] <- nb2listw(Y.rel, style = style, 
                                 glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                f2, y = i, dmax = max.rel))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcdown
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Relative neighbourhood_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcup))
        count_list <- 1
        for (i in y_fconcup) {
          listf[[count_list]] <- nb2listw(Y.rel, style = style, 
                                 glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                f3, y = i))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcup
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Relative neighbourhood_Concave up weighting"
      }
    }
  }
  if (mst == TRUE) {
    Y.mst <- mst.nb(xy.d1)
    if (binary == TRUE) {
      count <- count + 1
      listwcand[[count]] <- nb2listw(Y.mst, style = style)
      names(listwcand)[count] <- "Minimum spanning tree_Binary weighting"
    }
    if (weightfun == TRUE) {
      max.mst <- max(unlist(nbdists(Y.mst, as.matrix(coord))))
      if (flin == TRUE) {
        count <- count + 1
        listwcand[[count]] <- nb2listw(Y.mst, style = style, 
                                       glist = lapply(nbdists(Y.mst, as.matrix(coord)), 
                                                      f1, dmax = max.mst))
        names(listwcand)[count] <- "Minimum spanning tree_Linear weighting"
      } 
      if (fconcdown == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcdown))
        count_list <- 1
        for (i in y_fconcdown) {
          listf[[count_list]] <- nb2listw(Y.mst, style = style, 
                                 glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                f2, y = i, dmax = max.mst))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcdown
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Minimum spanning tree_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", length(y_fconcup))
        count_list <- 1
        for (i in y_fconcup) {
          listf[[count_list]] <- nb2listw(Y.mst, style = style, 
                                 glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                f3, y = i))
          count_list <- count_list + 1
        }
        names(listf) <- y_fconcup
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Minimum spanning tree_Concave up weighting"
      }
    }
  }
  return(listwcand)
} # End of the lisw.candidates function
# *********************************************************************************** #
  