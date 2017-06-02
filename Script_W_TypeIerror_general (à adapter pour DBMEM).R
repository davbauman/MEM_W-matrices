
POUR LES DBMEM --> METTRE DANS UN MEME SCRIPT LE TEST ERREUR TYPE I DES DBMEM SANS CORRECTION (exactement comme est fait ici donc) ET
LE TEST AVEC CORRECTION DE SIDAK METTANT EN EVIDENCE QUE PROBLEME REGLE QUAND CORRIGEONS LA PVALUE



################################################################################################################
# Code for computing the Type I Error Rate of severall spatial weighting matrices (W) obtained by              #
# combination of different connectivity and weighting matrices, using uni- and multivariate response variables #
# in a regular and an irregular sampling design.                                                               #
################################################################################################################

# Useful packages and functions:
# ******************************

library(vegan)
library(adespatial)
library(spdep)

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
}                                                                            # End of the test.W.R2() function
# ************************************************************************************************************

# Construction of a result matrix:                               
# ********************************

# For each B matrix, no A matrix and three A matrices tested; Two additional lines for
# the original PCNM method and for the PCNM method computed through MEM (DBEM, see further).
# For the columns: lines 3 = type I error; line 4 = mean R2adj; line 5 = sd of the R2adj;
# 1000 permutations, so that lines 6 to 1005 contain p-values, and lines 1006 to 2005 
# contain R2adj.

results <- as.data.frame(matrix(nrow = 21, ncol = 2005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
   paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), sep = ""))
results[,1] <- c(rep(c("del", "gab", "rel", "mst", "dnn"), each = 4), "DBMEM (PCNM)")
results[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), times = 5), "1-(D/4t)^2")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"
#MEM_model = "negative"   ; autocor <- "neg" 

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

# Generation of the 117 quadrats:
#################################

if(design == "regular"){

   C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m

   # We define the quadrat coordinates
   X1 <- c()
   Y1 <- c()
   for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
   for(i in 1:50){ Y1 <- c(Y1, c(1:25))}

   C[, 1] <- X1
   C[, 2] <- Y1

   # We choose the 117 regularly spaced quadrats of the grid

   tx <- seq(from = 1, to = 50, by = 4)
   ty <- seq(from = 1, to = 25, by = 3)
   C <- C[C[, 1] %in% tx, ]
   C <- C[C[, 2] %in% ty, ]

} else {

   # We choose 117 randomly spaced quadrats inside the grid 

   C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m

   set.seed(123)
   C[, 1] <- runif(117, min = 1, max = 50)
   C[, 2] <- runif(117, min = 1, max = 25)
}

xy.d1 <- dist(C)

# II. Connectivity matrices (B):
################################

Y.del <- tri2nb(C) 
Y.gab <- graph2nb(gabrielneigh(as.matrix(C), nnmult = 4), sym = T)
Y.rel <- graph2nb(relativeneigh(as.matrix(C), nnmult = 4), sym = T)
Y.mst <- mst.nb(xy.d1)

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
uplim <- max(xy.d1)
thresh10 <- seq(lowlim, uplim, le=10)   # 10 tested threshold distances
Y.list10dnn <- lapply(thresh10, dnearneigh, x=as.matrix(C), d1=0)

# III. Generation of MEM eigenvectors based on these B matrices: with a binary weighting matrix:
################################################################################################

Y.del.MEM <- test.W.R2(nb = Y.del, xy = C, MEM.autocor = MEM_model)
Y.gab.MEM <- test.W.R2(nb = Y.gab, xy = C, MEM.autocor = MEM_model)
Y.rel.MEM <- test.W.R2(nb = Y.rel, xy = C, MEM.autocor = MEM_model)
Y.mst.MEM <- test.W.R2(nb = Y.mst, xy = C, MEM.autocor = MEM_model)

# Compute the 10 sets of DBMEM variables (one for each threshold distance tested, with a binary weighting).
# In the simulation loop, one of the 10 sets will be selected for each Y based on the highest R2adj.

Y.list10dnn.MEM.list <- lapply(Y.list10dnn, function(x) test.W.R2(nb = x, xy = C, MEM.autocor = MEM_model))

# Weighting functions and fixed parametres:
# *****************************************

f1 <- function (D, dmax)    {1-(D/dmax)}        # Linear function
f2 <- function (D, dmax, y) {1-(D/dmax)^y}      # Concave-down function
f3 <- function (D, y)       {1/(D^y)}           # Concave-up function
f4 <- function (D, t)       {1-(D/(4*t))^2}     # PCNM criterion

max.del <- max(unlist(nbdists(Y.del, as.matrix(C)))) 
max.gab <- max(unlist(nbdists(Y.gab, as.matrix(C))))
max.rel <- max(unlist(nbdists(Y.rel, as.matrix(C)))) 
max.mst <- max(unlist(nbdists(Y.mst, as.matrix(C))))
 
nbdist <- lapply(Y.list10dnn, coords = as.matrix(C), nbdists)
 unlist <- lapply(nbdist, unlist)
  max.list10dnn <- lapply(unlist, max)

# Generation of MEM eigenvectors using the linear weighting function, and DBMEM with the PCNM criterion.
# These A matrices are separated from the others because they can be computed independently of a response
# variable, and therefore do not need to be computed inside the simulation loop.
# *******************************************************************************************************

   # DBMEM (PCNM criteria)
   # *********************
matB <- dnearneigh(lowlim, x = as.matrix(C), d1 = 0)
Y.list10dnn.MEM.w.f4 <- test.W.R2(nb = matB, xy = C, f = f4, t = lowlim, MEM.autocor = MEM_model)

   # Linear weighting function
   # *************************
Y.del.MEM.w.f1 <- test.W.R2(nb = Y.del, xy = C, f = f1, dmax = max.del, MEM.autocor = MEM_model)
Y.gab.MEM.w.f1 <- test.W.R2(nb = Y.gab, xy = C, f = f1, dmax = max.gab, MEM.autocor = MEM_model)
Y.rel.MEM.w.f1 <- test.W.R2(nb = Y.rel, xy = C, f = f1, dmax = max.rel, MEM.autocor = MEM_model)
Y.mst.MEM.w.f1 <- test.W.R2(nb = Y.mst, xy = C, f = f1, dmax = max.mst, MEM.autocor = MEM_model)

# Compute the 10 sets of DBMEM variables (one for each threshold distance tested, with the linear weighting f1).
# In the simulation loop, one of the 10 sets will be selected for each Y based on the highest R2adj.

Y.list10dnn.MEM.w.f1.list <- lapply(Y.list10dnn, function(x) test.W.R2(nb = x, xy = C, f = f1, 
   dmax = max(unlist(nbdists(x, as.matrix(C)))), MEM.autocor = MEM_model))

   # ***********************************************************************************
   # The simulation begins here 
   # The MEM variables weighted by a dist. function will be generated after generating Y
   # (because 'test.W.R2' requires a response vector/matrix)
   # ***********************************************************************************

# Simulation procedure:
#######################

for(i in 1:nperm){   

   set.seed(i)

if(framework == "univariate"){
   
   Y <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"             # Random (uniform)
#   Y <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" # Random (normal)
#   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                        # Exponential (1)
#   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                    # Exponential cubed

} else {

   Y <- matrix(ncol = 5, nrow = nrow(C))
   for(b in 1:ncol(Y)){

     Y[, b] <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"    
#     Y[, b] <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" 
#     Y[, b] <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
#     Y[, b] <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3" 
   }                  
}

# Construction of MEM variables using an A matrix:
# ************************************************

# del
Y.del.MEM.w.f2 <- test.W.R2(Y = Y, nb = Y.del, xy = C, MEM.autocor = MEM_model, f = f2, dmax = max.del, y = 2:10)
Y.del.MEM.w.f3 <- test.W.R2(Y = Y, nb = Y.del, xy = C, MEM.autocor = MEM_model, f = f3, y = 2:10)

# gab
Y.gab.MEM.w.f2 <- test.W.R2(Y = Y, nb = Y.gab, xy = C, MEM.autocor = MEM_model, f = f2, dmax = max.gab, y = 2:10)
Y.gab.MEM.w.f3 <- test.W.R2(Y = Y, nb = Y.gab, xy = C, MEM.autocor = MEM_model, f = f3, y = 2:10)

# rel
Y.rel.MEM.w.f2 <- test.W.R2(Y = Y, nb = Y.rel, xy = C, MEM.autocor = MEM_model, f = f2, dmax = max.rel, y = 2:10)
Y.rel.MEM.w.f3 <- test.W.R2(Y = Y, nb = Y.rel, xy = C, MEM.autocor = MEM_model, f = f3, y = 2:10)

# mst
Y.mst.MEM.w.f2 <- test.W.R2(Y = Y, nb = Y.mst, xy = C, MEM.autocor = MEM_model, f = f2, dmax = max.mst, y = 2:10)
Y.mst.MEM.w.f3 <- test.W.R2(Y = Y, nb = Y.mst, xy = C, MEM.autocor = MEM_model, f = f3, y = 2:10)

# dnn
R2.list <- lapply(Y.list10dnn.MEM.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  Y.list10dnn.MEM <- Y.list10dnn.MEM.list[[thebest]]

R2.list <- lapply(Y.list10dnn.MEM.w.f1.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  Y.list10dnn.MEM.w.f1 <- Y.list10dnn.MEM.w.f1.list[[thebest]]

Y.list10dnn.MEM.w.f2.list <- lapply(Y.list10dnn, function(x) test.W.R2(x, Y = Y, xy = C, f = f2, y = 2:10,
   dmax = max(unlist(nbdists(x, C))), MEM.autocor = MEM_model))
R2.list <- lapply(Y.list10dnn.MEM.w.f2.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  Y.list10dnn.MEM.w.f2 <- Y.list10dnn.MEM.w.f2.list[[thebest]]

Y.list10dnn.MEM.w.f3.list <- lapply(Y.list10dnn, function(x) test.W.R2(x, Y = Y, xy = C, f = f3, y = 2:10,
   MEM.autocor = MEM_model))
R2.list <- lapply(Y.list10dnn.MEM.w.f3.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  Y.list10dnn.MEM.w.f3 <- Y.list10dnn.MEM.w.f3.list[[thebest]]

# Significance test and MEM variable selection (forward selection with double stopping criterion):
# ************************************************************************************************

# del
R2adj <- RsquareAdj(rda(Y, Y.del.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.del.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.del.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.del.MEM$best$MEM[, c(sign)]
   results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[1, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[1, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[1, i+1005] <- 0
     }  
R2adj <- RsquareAdj(rda(Y, Y.del.MEM.w.f1$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.del.MEM.w.f1$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.w.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.del.MEM.w.f1$best$MEM[, c(sign)]
   results[2, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[2, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[2, i+5] <- 1  
          results[2, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.del.MEM.w.f2$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.del.MEM.w.f2$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.w.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.del.MEM.w.f2$best$MEM[, c(sign)]
   results[3, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[3, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[3, i+5] <- 1   
          results[3, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.del.MEM.w.f3$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.del.MEM.w.f3$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.w.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.del.MEM.w.f3$best$MEM[, c(sign)]
   results[4, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[4, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[4, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[4, i+1005] <- 0
     }

# gab
R2adj <- RsquareAdj(rda(Y, Y.gab.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.gab.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.gab.MEM$best$MEM[, c(sign)]
   results[5, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[5, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[5, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[5, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.w.f1$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.gab.MEM.w.f1$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.w.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.gab.MEM.w.f1$best$MEM[, c(sign)]
   results[6, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[6, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[6, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[6, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.w.f2$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.gab.MEM.w.f2$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.w.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.gab.MEM.w.f2$best$MEM[, c(sign)]
   results[7, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[7, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[7, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[7, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.w.f3$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.gab.MEM.w.f3$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.w.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.gab.MEM.w.f3$best$MEM[, c(sign)]
   results[8, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[8, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[8, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[8, i+1005] <- 0
     }

# rel
R2adj <- RsquareAdj(rda(Y, Y.rel.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.rel.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.rel.MEM$best$MEM[, c(sign)]
   results[9, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[9, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[9, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[9, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.w.f1$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.rel.MEM.w.f1$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.w.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.rel.MEM.w.f1$best$MEM[, c(sign)]
   results[10, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[10, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[10, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[10, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.w.f2$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.rel.MEM.w.f2$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.w.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.rel.MEM.w.f2$best$MEM[, c(sign)]
   results[11, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[11, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[11, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[11, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.w.f3$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.rel.MEM.w.f3$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.w.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.rel.MEM.w.f3$best$MEM[, c(sign)]
   results[12, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[12, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[12, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[12, i+1005] <- 0
     }

# mst
R2adj <- RsquareAdj(rda(Y, Y.mst.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.mst.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.mst.MEM$best$MEM[, c(sign)]
   results[13, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[13, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[13, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[13, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.w.f1$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.mst.MEM.w.f1$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.w.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.mst.MEM.w.f1$best$MEM[, c(sign)]
   results[14, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[14, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[14, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[14, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.w.f2$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.mst.MEM.w.f2$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.w.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.mst.MEM.w.f2$best$MEM[, c(sign)]
   results[15, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[15, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[15, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[15, i+1005] <- 0
     }
R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.w.f3$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.mst.MEM.w.f3$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.w.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- Y.mst.MEM.w.f3$best$MEM[, c(sign)]
   results[16, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[16, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[16, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
          results[16, i+1005] <- 0
     }

# dnn
R2adj <- RsquareAdj(rda(Y, Y.list10dnn.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.list10dnn.MEM$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.list10dnn.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
       sign <- sort(fsel$order)
       MEM.FwdSel <- Y.list10dnn.MEM$best$MEM[, c(sign)]
       results[17, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
       results[17, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    }
 } else { 
    results[17, i+5] <- 1   # p-val made equal to 1 
    results[17, i+1005] <- 0
   }

R2adj <- RsquareAdj(rda(Y, Y.list10dnn.MEM.w.f1$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.list10dnn.MEM.w.f1$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.list10dnn.MEM.w.f1$best$MEM, adjR2thresh = R2adj, nperm = 999),
       TRUE))
    if(class != "try-error"){
       sign <- sort(fsel$order)
       MEM.FwdSel <- Y.list10dnn.MEM.w.f1$best$MEM[, c(sign)]
       results[18, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
       results[18, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    }
 } else { 
    results[18, i+5] <- 1   # p-val made equal to 1 
    results[18, i+1005] <- 0
   }

R2adj <- RsquareAdj(rda(Y, Y.list10dnn.MEM.w.f2$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.list10dnn.MEM.w.f2$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.list10dnn.MEM.w.f2$best$MEM, adjR2thresh = R2adj, nperm = 999),
       TRUE))
    if(class != "try-error"){
       sign <- sort(fsel$order)
       MEM.FwdSel <- Y.list10dnn.MEM.w.f2$best$MEM[, c(sign)]
       results[19, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
       results[19, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    }
 } else { 
    results[19, i+5] <- 1   # p-val made equal to 1 
    results[19, i+1005] <- 0
   }

R2adj <- RsquareAdj(rda(Y, Y.list10dnn.MEM.w.f3$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.list10dnn.MEM.w.f3$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.list10dnn.MEM.w.f3$best$MEM, adjR2thresh = R2adj, nperm = 999),
       TRUE))
    if(class != "try-error"){
       sign <- sort(fsel$order)
       MEM.FwdSel <- Y.list10dnn.MEM.w.f3$best$MEM[, c(sign)]
       results[20, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
       results[20, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    }
 } else { 
    results[20, i+5] <- 1   # p-val made equal to 1 
    results[20, i+1005] <- 0
   }

# DBMEM (with PCNM criteria: B = give.thresh(xy.d1) and A = f4)
R2adj <- RsquareAdj(rda(Y, Y.list10dnn.MEM.w.f4$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.list10dnn.MEM.w.f4$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.list10dnn.MEM.w.f4$best$MEM, adjR2thresh = R2adj, nperm = 999),
       TRUE))
    if(class != "try-error"){
       sign <- sort(fsel$order)
       MEM.FwdSel <- Y.list10dnn.MEM.w.f4$best$MEM[, c(sign)]
       results[21, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
       results[21, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } 
 } else { 
    results[21, i+5] <- 1   # p-val made equal to 1 
    results[21, i+1005] <- 0
   }

} # End of the 'for' simulation loop

# Type I error, median and sd of R2adj:
#######################################

for(i in 1:nrow(results)){
   results[i, 3] <- length(which(results[i, c(6:(nperm + 5))] <= 0.05)) / nperm
   results[i, 4] <- median(as.numeric(results[i, c(1006:(nperm + 1005))]))
   results[i, 5] <- sd(as.numeric(results[i, c(1006:(nperm + 1005))]))
}
