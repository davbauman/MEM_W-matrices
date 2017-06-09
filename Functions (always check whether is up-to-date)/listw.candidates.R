# Trouver comment gérer la question des DBMEM et des threshold distances

# Idée pour f2 et f3. Avons une liste d'objets de nb2listw, mais pour f2 et f3 avons
# une liste d'autant d'éléments qu'on teste de y (exp). Avons donc une liste dans la
# liste. Dans la fonction MEM.modsel, si l'objet[[i]] de la list listwcand est une liste,
# alors trouvons le meilleur 'y' sur base du R2 (v. test.W.R2) puis faisons le test
# uniquement sur cette matrice W.

# Pour les listes dans la liste (f2 et f3), la valeur du 'y' est le nom (names) de 
# l'élément liste (v. listw.candidates)

# Problème est que cette list-output doit ressembler à n'importe quelle liste que
# l'utilisateur voudrait créer lui-même pour que MEM.modsel puisse traiter la liste
# lui était fournie de manière égale.

# *********************************************************************************** #
listw.candidates <- function (coord, style = "W", del = TRUE, gab = TRUE, rel = TRUE, 
                             mst = TRUE, PCNM = TRUE,  DB = TRUE, thresh = "?", 
                             binary = TRUE, weightfun = TRUE, flin = TRUE, 
                             fconcdown = TRUE, fconcup = TRUE, ymax = 5)
{
  
  if (any(is.na(coord))) stop("NA entries in coord")
  
  # Definition of the weighting functions:
  if (weightfun == TRUE) {
    f1 <- function (D, dmax)     { 1 - (D/dmax) }       # Linear function
    f2 <- function (D, dmax, y)  { 1 - (D/dmax)^y }     # Concave-down function
    f3 <- function (D, dmax, y)  { 1 / (D/dmax)^y }     # Concave-up function
  }
  
  xy.d1 <- dist(coord)
  
  # Total nb of W matrices to build:
  nbw <- length(which(c(del, gab, rel, mst) == TRUE))  # AJOUTER DB + nb thresh dist.
  if (weightfun == TRUE) {
    if (binary == TRUE) {
      nbw <- nbw * (length(which(c(flin, fconcdown, fconcup) == TRUE))+1)
    } else nbw <- nbw * length(which(c(flin, fconcdown, fconcup) == TRUE))
  }
  if (PCNM == TRUE) nbw <- nbw + 1
  
  # List of the W matrix candidates
  listwcand <- vector("list", nbw)
  listwfeatures <- vector("list", nbw)
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
        listf <- vector("list", ymax-1)
        for (i in 1:(ymax-1)) {
          listf[[i]] <- nb2listw(Y.del, style = style, 
                                       glist = lapply(nbdists(Y.del, as.matrix(coord)),
                                                      f2, y = i+1, dmax = max.del))
        }
        names(listf) <- c(2:ymax)
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Delaunay_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", ymax)
        for (i in 1:ymax) {
          listf[[i]] <- nb2listw(Y.del, style = style, 
                                 glist = lapply(nbdists(Y.del, as.matrix(coord)),
                                                f3, y = i, dmax = max.del))
        }
        names(listf) <- c(1:ymax)
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
        listf <- vector("list", ymax-1)
        for (i in 1:(ymax-1)) {
          listf[[i]] <- nb2listw(Y.gab, style = style, 
                                 glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                f2, y = i+1, dmax = max.gab))
        }
        names(listf) <- c(2:ymax)
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Gabriel_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", ymax)
        for (i in 1:ymax) {
          listf[[i]] <- nb2listw(Y.gab, style = style, 
                                 glist = lapply(nbdists(Y.gab, as.matrix(coord)),
                                                f3, y = i, dmax = max.gab))
        }
        names(listf) <- c(1:ymax)
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
        listf <- vector("list", ymax-1)
        for (i in 1:(ymax-1)) {
          listf[[i]] <- nb2listw(Y.rel, style = style, 
                                 glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                f2, y = i+1, dmax = max.rel))
        }
        names(listf) <- c(2:ymax)
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Relative neighbourhood_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", ymax-1)
        for (i in 1:ymax) {
          listf[[i]] <- nb2listw(Y.rel, style = style, 
                                 glist = lapply(nbdists(Y.rel, as.matrix(coord)),
                                                f3, y = i, dmax = max.rel))
        }
        names(listf) <- c(1:ymax)
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
        listf <- vector("list", ymax-1)
        for (i in 1:(ymax-1)) {
          listf[[i]] <- nb2listw(Y.mst, style = style, 
                                 glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                f2, y = i+1, dmax = max.mst))
        }
        names(listf) <- c(2:ymax)
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Minimum spanning tree_Concave down weighting"
      }
      if (fconcup == TRUE) {      
        count <- count + 1
        listf <- vector("list", ymax-1)
        for (i in 1:ymax) {
          listf[[i]] <- nb2listw(Y.mst, style = style, 
                                 glist = lapply(nbdists(Y.mst, as.matrix(coord)),
                                                f3, y = i, dmax = max.mst))
        }
        names(listf) <- c(1:ymax)
        listwcand[[count]] <- listf
        names(listwcand)[count] <- "Minimum spanning tree_Concave up weighting"
      }
    }
  }
  return(listwcand)
} # End of the lisw.candidates function
# *********************************************************************************** #
  