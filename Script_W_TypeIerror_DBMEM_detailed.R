
#############################################################################################
#############################################################################################
### We test the type I error rate of the DBMEM process by applying the Sidak correction
### only the W matrix to which corresponds the distance threshold displaying the highest R2.
#############################################################################################
#############################################################################################


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
# For DBMEM (10 different dmax thresholds), a binary and a linear weighting function:

results <- as.data.frame(matrix(nrow = 2, ncol = 2005))

colnames(results) <- c("Matrix A", "/", "type I error", "mean R2adj", "sd R2adj",
   paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), sep = ""))
results[,1] <- c("bin", "linear")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"
#MEM_model = "negative"   ; autocor <- "neg" 

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 100

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

# II. Distance-based B (radius around points):
# ********************************************
lowlim <- give.thresh(xy.d1)
uplim <- max(xy.d1)
thresh10 <- seq(lowlim, uplim, le=10)   # 10 tested distances

Y.list10dnn <- lapply(thresh10, dnearneigh, x=as.matrix(C), d1=0)

# III. Generation of MEM eigenvectors based on these B matrices: with a binary weighting matrix:
################################################################################################

Y.list10dnn.MEM.list <- lapply(Y.list10dnn, function(x) test.W.R2(nb = x, xy = C, MEM.autocor = MEM_model))

# Weighting functions and fixed parametres:
# *****************************************

f1 <- function (D, dmax)    {1-(D/dmax)}        # Linear function
f2 <- function(D, dmax, y) { 1 - (D/dmax)^y }   # Concave-down function
f3 <- function (D, y) {1/(D^y)}                 # Concave-up function

nbdist <- lapply(Y.list10dnn, coords = as.matrix(C), nbdists)
 unlist <- lapply(nbdist, unlist)
  max.list10dnn <- lapply(unlist, max)

# Generation of MEM eigenvectors using the linear weighting function, and DBMEM with the PCNM criterion.
# These A matrices are separated from the others because they can be computed independently of a response
# variable, and therefore do not need to be computed inside the simulation loop.
# *******************************************************************************************************

   # Linear weighting function
   # *************************

# Compute the 10 sets of DBMEM variables (one for each threshold distance tested, with the linear weighting f1).
# In the simulation loop, one of the 10 sets will be selected for each Y based on the highest R2adj.

Y.list10dnn.MEM.w.f1.list <- lapply(Y.list10dnn, function(x) test.W.R2(nb = x, xy = C, f = f1, 
   dmax = max(unlist(nbdists(x, as.matrix(C)))), MEM.autocor = MEM_model))

# Correction of the p-value (Sidak):
# **********************************
# Define the number of distance threshold tested for the DBMEM:

nbtest <- length(thresh10)
corrected_alpha <- 1-((1-0.05)^(1/nbtest))   # Sidak correction

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

# Construction of the DBMEM for the concdown and concup functions:
# ****************************************************************

Y.list10dnn.MEM.w.f2.list <- lapply(Y.list10dnn, function(x) test.W.R2(x, Y = Y, xy = C, f = f2, y = 2:10,
   dmax = max(unlist(nbdists(x, C))), MEM.autocor = MEM_model))

Y.list10dnn.MEM.w.f3.list <- lapply(Y.list10dnn, function(x) test.W.R2(x, Y = Y, xy = C, f = f3, y = 2:10,
   MEM.autocor = MEM_model))

# Significance test and MEM variable selection (forward selection with double stopping criterion):
# ************************************************************************************************
# dnn
     # bin
     # ***
R2.list <- lapply(Y.list10dnn.MEM.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  MEM.best <- Y.list10dnn.MEM.list[[thebest]]$best$MEM
    if(anova.cca(rda(Y, MEM.best))$Pr[1] <= corrected_alpha){
      R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
      class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, nperm = 999), TRUE))
      if(class != "try-error"){
         sign <- sort(fsel$order)
         MEM.FwdSel <- MEM.best[, c(sign)]
         results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
         results[1, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
      } 
    } else { 
         results[1, i+5] <- 1   # p-val made equal to 1 
         results[1, i+1005] <- 0
      }

     # linear
     # ******
R2.list <- lapply(Y.list10dnn.MEM.w.f1.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  MEM.best <- Y.list10dnn.MEM.w.f1.list[[thebest]]$best$MEM
    if(anova.cca(rda(Y, MEM.best))$Pr[1] <= corrected_alpha){
      R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
      class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, nperm = 999), TRUE))
      if(class != "try-error"){
         sign <- sort(fsel$order)
         MEM.FwdSel <- MEM.best[, c(sign)]
         results[2, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
         results[2, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
      } 
    } else { 
         results[2, i+5] <- 1   # p-val made equal to 1 
         results[2, i+1005] <- 0
      }

     # Concave-down
     # ************
R2.list <- lapply(Y.list10dnn.MEM.w.f2.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  MEM.best <- Y.list10dnn.MEM.w.f2.list[[thebest]]$best$MEM
    if(anova.cca(rda(Y, MEM.best))$Pr[1] <= corrected_alpha){
      R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
      class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, nperm = 999), TRUE))
      if(class != "try-error"){
         sign <- sort(fsel$order)
         MEM.FwdSel <- MEM.best[, c(sign)]
         results[3, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
         results[3, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
      } 
    } else { 
         results[3, i+5] <- 1   # p-val made equal to 1 
         results[3, i+1005] <- 0
      }

     # Concave-up
     # **********
R2.list <- lapply(Y.list10dnn.MEM.w.f3.list, function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
 thebest <- which.max(R2.list)
  MEM.best <- Y.list10dnn.MEM.w.f3.list[[thebest]]$best$MEM
    if(anova.cca(rda(Y, MEM.best))$Pr[1] <= corrected_alpha){
      R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
      class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, nperm = 999), TRUE))
      if(class != "try-error"){
         sign <- sort(fsel$order)
         MEM.FwdSel <- MEM.best[, c(sign)]
         results[4, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
         results[4, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
      } 
    } else { 
         results[4, i+5] <- 1   # p-val made equal to 1 
         results[4, i+1005] <- 0
      }

} # End of the 'for' simulation loop

# Type I error, median and sd of R2adj:
#######################################

for(i in 1:nrow(results)){
   results[i, 3] <- length(which(results[i, c(6:(nperm + 5))] <= 0.05)) / nperm
   results[i, 4] <- median(as.numeric(results[i, c(1006:(nperm + 1005))]))
   results[i, 5] <- sd(as.numeric(results[i, c(1006:(nperm + 1005))]))
}

# Output of the results:
# **********************

res_file_name <- paste("Results_typeIerror_DBMEM_SIDAK_OneStep", paste(design, ".txt", sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")
