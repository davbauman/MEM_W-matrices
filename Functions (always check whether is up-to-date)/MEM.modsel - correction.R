# Function MEM.modsel:
# ********************
# The function requires the packages vegan, adespatial and spdep to be installed. The MEM.modsel function will load them.
# INPUT:
# ******
# Arguments:
# - x is the vector or dataframe of response variable(s).
# - coord is the matrix or dataframe of X and Y spatial coordinates. It can also be a vector in the case of a transect.
# autocor corresponds to the type of spatial correlation modelled by the MEM eigenfunctions to be considered: autocor is 
# either "positive" (default argument), "negative" or "all" --> positively, negatively or both types of MEM eigenfunctions
# are considered. In the latter case, the positive and negative models are built, tested and selected separately (a global
# adjusted R2 cannot be obtained and therefore cannot but used as second stopping criterion of the forward selection for
# saturated models (n-1 variables here, if positive and negative MEM eigenfunctions considered together).
# - PCNM: whether a distance-based MEM based on the PCNM criterion (Dray et al. 2006) should be computed and tested.
# - del, gab, rel, mst: Connectivity matrices to be tested (Delaunay triangulation, Gabriel's graph, relative neighbourhood 
# graph, and minimum spanning tree, respectively). The default arguments test all connectivity matrices.
# - weightfun: Whether weighting matrices should be added to the connectivity matrices (del, gab, rel, mst). If weightfun
# = TRUE, three weighting functions can be tested: a decreasing linear (flin = TRUE), a concave-down (fconcdown = TRUE), 
# and a concave-up function (fconcup = TRUE). 
# - ymax: The two latter functions are tested by default for a 'y' exponent parameter taking the value 2 to 10 and 1 to 10, 
# respectively. The maximum value of the 'y' exponent can be set to a different value through the 'ymax' argument.
# If weightfun = FALSE, only the binary MEM based on the selected connectivity matrices are tested.
# - alpha_thresh: Significance threshold value.
# OUTPUT:
# *******
# If only a positive or negative spatial model is tested, or if both are tested but only one was significant, the function
# returns a list containing all characteristics of the best spatial model selected: the connectivity and weighting matrices 
# used for constructing the spatial weighting matrix (W), the 'y' exponent parameter (if a concave-down or -up weighting 
# function was used), the p-value of the global model (corrected by a Sidak correction depending on the number of W matrices
# tested), the adjusted R2 of the final model (after performing the forward selection with double stopping criterion), the 
# R2 of each MEM variable and the eigenvectors (i.e., the selected MEM variables). 
# The latter can directly be used as spatial explanatory variables in statistical models.
# If both positive and negative spatial models are tested and if both are significant, the function returns a list of two
# lists (MEM.pos and MEM.neg), each one containing the same information as described above. 

MEM.modsel <- function(x, coord, autocor = c("positive", "negative", "all"), PCNM = TRUE, del = TRUE, gab = TRUE, 
   rel = TRUE, mst = TRUE, weightfun = TRUE, flin = TRUE, fconcdown = TRUE, fconcup = TRUE, ymax = 10, alpha_thresh = 0.05)
{

   library(vegan)
   library(adespatial)
   library(spdep)

   x <- as.data.frame(x)
   if(any(is.na(x)) | any(is.na(coord))) 
      stop("NA entries in x or coord")
   if(nrow(x) != nrow(coord)) 
      stop("different number of rows")

   options(warn = 2)

   # *************************************************************************************************************
   MEM.test <- function(a = x, b, c = autocor, d = nbtest, alpha = alpha_thresh)
   {
      pval <- anova.cca(rda(a, b))$Pr[1]
      pval <- 1-(1-pval)^d                      # Sidak correction (nb of MEM model tested and compared) 
      if(c == "all") pval <- 1-(1-pval)^2       # Sidak correction for positive and negative MEM model computation 
      if(pval <= alpha){  
         R2adj <- RsquareAdj(rda(a, b))$adj.r.squared
         class <- class(try(fsel <- forward.sel(a, b, adjR2thresh = R2adj, nperm = 999), TRUE))
         if(class != "try-error"){ 
            sign <- sort(fsel$order)
            MEM.select <- b[, c(sign)] 
            list(MEM.select = MEM.select, NbVar = length(sign), pval = pval, 
               R2adj = RsquareAdj(rda(a, MEM.select))$adj.r.squared, AdjR2Cum = fsel$AdjR2Cum)
            } 
      } else return(NA)
   }                                                                             # End of the MEM.test() function
   # ************************************************************************************************************

   # ************************************************************************************************************
   test.W.R2 <- function (Y, nb, xy, MEM.autocor = c("positive", "negative"), f = NULL, ...) 
   {
      mycall <- pairlist(...)   
      res <- list()  # List of the MEM eigenvectors computed for each W matrix tested within test.W.R2()
      MEM.autocor <- match.arg(MEM.autocor)   

      nbdist <- nbdists(nb, as.matrix(xy))
      param <- expand.grid(as.list(mycall))
      m1 <- match(names(param), names(formals(f)))
      for (i in 1:nrow(param)) {
         formals(f)[m1] <- unclass(param[i, ])
         res[[i]] <- scores.listw(nb2listw(nb, style = "B", 
            glist = lapply(nbdist, f), zero.policy = TRUE), MEM.autocor = MEM.autocor)
      }

      res2 <- lapply(res, function(z) RsquareAdj(rda(Y, z))$adj.r.squared)
      thebest <- which.max(res2)
      return(list(param = param[thebest, ], best = list(MEM = res[[thebest]])))
   }                                                                            # End of the test.W.R2() function
   # ************************************************************************************************************

   autocor <- match.arg(autocor) 

   if(weightfun == TRUE){
      f1 <- function (D, dmax)    {1-(D/dmax)}        # Linear function
      f2 <- function(D, dmax, y) { 1 - (D/dmax)^y }   # Concave-down function
      f3 <- function (D, y) {1/(D^y)}                 # Concave-up function
   }

   xy.d1 <- dist(coord)

   # Since the loop is entered only once if autocor = "positive" or "negative" but twice if autocor = "all" (once for the 
   # positive and once for the negative MEM models), we begin by setting up the loop.

   lenlist <- c()   # Will help with the result output  
   if(autocor != "all"){ 
      k <- 1   # number of times the loop will run
      if(autocor == "positive"){ 
         cor <- c("positive", "negative") 
      }  else cor <- c("negative", "positive") 
   } else k <- 2 ; cor <- c("positive", "negative")

   # A multitest p-value correction is needed for controling the type-I error rate. We define the total nb of tests:
   nbtest <- length(which(c(del, gab, rel, mst) == TRUE))
   if(weightfun == TRUE) nbtest <- nbtest * length(which(c(flin, fconcdown, fconcup) == TRUE))
   if(PCNM == TRUE) nbtest <- nbtest + 1   

   for(h in 1:k){

      # For model comparison and selection
      results <- as.data.frame(matrix(nrow = 17, ncol = 6))
      colnames(results) <- c("Matrix B", "Matrix A", "pval_sidak", "R2adj", "NbMEM", "y")     
      results[, 1] <- c("DBMEM", rep(c("Delaunay triangulation", "Gabriel's graph", "relative neigh.", 
                        "minimum spanning tree"), each = 4))
      results[, 2] <- c("PCNM criterion", rep(c("binary", "linear", "concave-down", "concave-up"), times = 4))

      # List of the MEM.select matrices
      listMEM <- vector("list", 17)
      # and corresponding R2
      listR2 <- vector("list", 17)
      # List of 'MEM.modsel' results
      listtest <- vector("list", 17)
      # List of global W matrices
      listW <- vector("list", 17)

      if(PCNM == TRUE){
         f <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion
         lowlim <- give.thresh(xy.d1)
         matB <- dnearneigh(lowlim, x = as.matrix(coord), d1 = 0)
         W <- scores.listw(nb2listw(matB, style = "B", glist = lapply(nbdists(matB, as.matrix(coord)), f, t = lowlim)), 
            MEM.autocor = cor[h])
         listtest[[1]] <- MEM.test(x, W)
         listW[[1]] <- W
      } 
      if(del == TRUE){
         Y.del <- tri2nb(coord)
         W <- scores.listw(nb2listw(Y.del, style = "B"), MEM.autocor = cor[h])
         listtest[[2]] <- MEM.test(x, W)
         listW[[2]] <- W
         if(weightfun == TRUE){
            max.del <- max(unlist(nbdists(Y.del, as.matrix(coord))))
            if(flin == TRUE){
               class <- class(try(W <- scores.listw(nb2listw(Y.del, style = "B", 
                  glist = lapply(nbdists(Y.del, as.matrix(coord)), f1, dmax = max.del)), MEM.autocor = cor[h]), TRUE))
               if(class[1] != "try-error"){
                  listtest[[3]] <- MEM.test(x, W)
                  listW[[3]] <- W
               }  
            } 
            if(fconcdown == TRUE){      
               class <- class(try(W <- test.W.R2(Y = x, nb = Y.del, xy = as.matrix(coord), MEM.autocor = cor[h], f = f2, 
                  y = 2:ymax, dmax = max.del), TRUE))
               if(class[1] != "try-error"){
                  listtest[[4]] <- MEM.test(x, W$best$MEM)
                  listW[[4]] <- W
               }
            }
            if(fconcup == TRUE){
               class = "try-error" ; yconcup = ymax
               while(class == "try-error"){  
                  class <- class(try(W <- test.W.R2(Y = x, nb = Y.del, xy = as.matrix(coord), MEM.autocor = cor[h], f = f3, 
                     y = 1:yconcup), TRUE))
                  if(class[1] == "try-error"){ yconcup <- yconcup - 1
                  } else {
                     listtest[[5]] <- MEM.test(x, W$best$MEM)
                     listW[[5]] <- W
                  }
                  if(yconcup == 0) break
               }
            }
         }
      }
      if(gab == TRUE){
         Y.gab <- graph2nb(gabrielneigh(as.matrix(coord), nnmult = 4), sym = TRUE)
         W <- scores.listw(nb2listw(Y.gab, style = "B"), MEM.autocor = cor[h])
         listtest[[6]] <- MEM.test(x, W)
         listW[[6]] <- W
         if(weightfun == TRUE){
            max.gab <- max(unlist(nbdists(Y.gab, as.matrix(coord))))
            if(flin == TRUE){
               class <- class(try(W <- scores.listw(nb2listw(Y.gab, style = "B", 
                  glist = lapply(nbdists(Y.gab, as.matrix(coord)), f1, dmax = max.gab)), MEM.autocor = cor[h]), TRUE))
               if(class[1] != "try-error"){
                  listtest[[7]] <- MEM.test(x, W)  
                  listW[[7]] <- W
               }   
            }
            if(fconcdown == TRUE){         
               class <- class(try(W <- test.W.R2(Y = x, nb = Y.gab, xy = as.matrix(coord), MEM.autocor = cor[h], f = f2, 
                  y = 2:ymax, dmax = max.gab), TRUE))
               if(class[1] != "try-error"){
                  listtest[[8]] <- MEM.test(x, W$best$MEM)
                  listW[[8]] <- W
               }
            }
            if(fconcup == TRUE){
               class = "try-error" ; yconcup = ymax
               while(class == "try-error"){    
                  class <- class(try(W <- test.W.R2(Y = x, nb = Y.gab, xy = as.matrix(coord), MEM.autocor = cor[h], f = f3, 
                     y = 1:yconcup), TRUE))
                  if(class[1] == "try-error"){ yconcup <- yconcup - 1
                  } else {
                     listtest[[9]] <- MEM.test(x, W$best$MEM)
                     listW[[9]] <- W
                  }
                  if(yconcup == 0) break
               }
            }
         }
      }
      if(rel == TRUE){
         Y.rel <- graph2nb(relativeneigh(as.matrix(coord), nnmult = 4), sym = TRUE)
         W <- scores.listw(nb2listw(Y.rel, style = "B"), MEM.autocor = cor[h])
         listtest[[10]] <- MEM.test(x, W)
         listW[[10]] <- W
         if(weightfun == TRUE){
            max.rel <- max(unlist(nbdists(Y.rel, as.matrix(coord)))) 
            if(flin == TRUE){
               class <- class(try(W <- scores.listw(nb2listw(Y.rel, style = "B", 
                  glist = lapply(nbdists(Y.rel, as.matrix(coord)), f1, dmax = max.rel)), MEM.autocor = cor[h]), TRUE))
               if(class[1] != "try-error"){
                  listtest[[11]] <- MEM.test(x, W)  
                  listW[[11]] <- W 
               }
            }
            if(fconcdown == TRUE){           
               class <- class(try(W <- test.W.R2(Y = x, nb = Y.rel, xy = as.matrix(coord), MEM.autocor = cor[h], f = f2, 
                  y = 2:ymax, dmax = max.rel), TRUE))
               if(class[1] != "try-error"){
                  listtest[[12]] <- MEM.test(x, W$best$MEM)
                  listW[[12]] <- W
               }
            }
            if(fconcup == TRUE){ 
               class = "try-error" ; yconcup = ymax
               while(class == "try-error"){  
                  class <- class(try(W <- test.W.R2(Y = x, nb = Y.rel, xy = as.matrix(coord), MEM.autocor = cor[h], f = f3,
                     y = 1:yconcup), TRUE))
                  if(class[1] == "try-error"){ yconcup <- yconcup - 1
                  } else {
                     listtest[[13]] <- MEM.test(x, W$best$MEM)
                     listW[[13]] <- W 
                  } 
                  if(yconcup == 0) break                
               }
            }
         }
      }
      if(mst == TRUE){
         Y.mst <- mst.nb(xy.d1)
         W <- scores.listw(nb2listw(Y.mst, style = "B"), MEM.autocor = cor[h])
         listtest[[14]] <- MEM.test(x, W)
         listW[[14]] <- W         
         if(weightfun == TRUE){
            max.mst <- max(unlist(nbdists(Y.mst, as.matrix(coord)))) 
            if(flin == TRUE){
               class <- class(try(W <- scores.listw(nb2listw(Y.mst, style = "B", 
                  glist = lapply(nbdists(Y.mst, as.matrix(coord)), f1, dmax = max.mst)), MEM.autocor = cor[h]), TRUE))
               if(class[1] != "try-error"){
                  listtest[[15]] <- MEM.test(x, W)
                  listW[[15]] <- W
               }
            }
            if(fconcdown == TRUE){         
               class <- class(try(W <- test.W.R2(Y = x, nb = Y.mst, xy = as.matrix(coord), MEM.autocor = cor[h], f = f2, 
                  y = 2:ymax, dmax = max.mst), TRUE))
               if(class[1] != "try-error"){
                  listtest[[16]] <- MEM.test(x, W$best$MEM)
                  listW[[16]] <- W
               }
            }  
            if(fconcup == TRUE){    
               class = "try-error" ; yconcup = ymax
               while(class == "try-error"){       
                  class <- class(try(W <- test.W.R2(Y = x, nb = Y.mst, xy = as.matrix(coord), MEM.autocor = cor[h], f = f3, 
                     y = 1:yconcup), TRUE))
                  if(class[1] == "try-error"){ yconcup <- yconcup - 1
                  } else {
                     listtest[[17]] <- MEM.test(x, W$best$MEM)
                     listW[[17]] <- W 
                  } 
                  if(yconcup == 0) break                  
               }  
            }
         }
      }
      # Save the results in order to compare them and choose the most parcimonious model:
      for(i in 1:nrow(results)){
         if(is.list(listtest[[i]]) == TRUE){
            results[i, 3] <- listtest[[i]]$pval
            results[i, 4] <- listtest[[i]]$R2adj
            results[i, 5] <- listtest[[i]]$NbVar
            if(length(listW[[i]]) == 2) results[i, 6] <- listW[[i]]$param[1]
            listMEM[[i]] <- listtest[[i]]$MEM.select
            listR2[[i]] <- listtest[[i]]$AdjR2Cum
         }
      }
      # Selection of the best model:
      if(length(which(results[, 3] <= alpha_thresh)) > 0){
         best <- which.max(results[, 4])
         if(autocor != "all"){
            lenlist <- c(lenlist, cor[h])
            L <- list(MEM.vec = listMEM[[best]], MEM.AdjR2Cum = listR2[[best]], Connectivity_Matrix = results[best, 1],
                    Weighting_fun = results[best, 2], pval = results[best, 3], R2adj = results[best, 4], 
                    NbVar = results[best, 5], y = results[best, 6])
         } else {
            lenlist <- c(lenlist, cor[h])
            if(h == 1){
               L1 <- list(MEM.vec = listMEM[[best]], MEM.AdjR2Cum = listR2[[best]], Connectivity_Matrix = results[best, 1],
                        Weighting_fun = results[best, 2], pval = results[best, 3], R2adj = results[best, 4], 
                        NbVar = results[best, 5], y = results[best, 6])
               
            } else {
               L2 <- list(MEM.vec = listMEM[[best]], MEM.AdjR2Cum = listR2[[best]], Connectivity_Matrix = results[best, 1],
                        Weighting_fun = results[best, 2], pval = results[best, 3], R2adj = results[best, 4], 
                        NbVar = results[best, 5], y = results[best, 6])
            }
         }
      }
   }    # End of the for loop

   options(warn = 1)

   if(length(lenlist) == 2){
      cat("\n", "\n", "************************************************************************", "\n",
          "************************************************************************", "\n",
         "Significant positive (p-value = ", L1$pval, ", R2adj = ", L1$R2adj, ") and negative (p-value = ", L2$pval, 
         ", R2adj = ", L2$R2adj, ")", "\n", "MEM models were detected and selected.", "\n", 
         "The best positive and negative models were built using ", L1$Connectivity_Matrix, " and ", 
         L2$Connectivity_Matrix, "\n", "(connectivity matrix), and ", L1$Weighting_fun, " and ", L2$Weighting_fun, 
         " weighting functions, respectively.", "\n",
         "************************************************************************", "\n",
         "The output of the function is a list of two lists (MEM.pos and MEM.neg).", "\n",
         "Within each list, the MEM variables are available in $MEM.vec.", "\n", "\n", sep = "")
      list(MEM.pos = L1, MEM.neg = L2)
   } else if(length(lenlist) == 1){
             if(autocor == "all") { if(lenlist == "positive") L <- L1 else L <- L2 }
             cat("\n", "\n", "************************************************************************", "\n",
                 "************************************************************************", "\n",
             "A significant ", lenlist, " MEM model was selected (p-value = ", L$pval, ", R2adj = ", L$R2adj, ").", "\n", 
             "The best model was built using ", L$Connectivity_Matrix, " (connectivity matrix) and a ", L$Weighting_fun, 
             " weighting function.", "\n",
             "************************************************************************", "\n",
             "The output of the function is a list containing the MEM variables ($MEM.vec) and other", "\n", 
             "characteristics and results of the model ($pval, $R2adj etc.)", "\n", "\n", sep = "")
             list(MEM.vec = L$MEM.vec, MEM.AdjR2Cum = L$MEM.AdjR2Cum, Connectivity_Matrix = L$Connectivity_Matrix,
                Weighting_fun = L$Weighting_fun, pval = L$pval, R2adj = L$R2adj, NbVar = L$NbVar, y = L$y)
          } else cat("\n", "\n", "****************************************************************", "\n",
                     "****************************************************************", "\n",
                     "No significant spatial structure was detected in the data.", "\n", "\n", sep = "")   
}

# Written by:
# David Bauman, in Bauman et al. (2017)
