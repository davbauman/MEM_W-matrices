#############################################################################
### ****************** Corrected version of MEM.modsel ****************** ###

library(vegan)
library(adespatial)
library(spdep)

# Functions:

source("lmp.R")
source("MEM.modsel - multiple_test_correction - adespatial - bestR2 - V3.R")

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

# Function to retrieve the pvalue of a lm() function:
# ***************************************************

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Construction of a results matrix for each scale and for both W matrix selection:
# ********************************************************************************

resultsB_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsB_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsM_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsM_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsF_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsF_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_modsel[, 1] <- c("sim", "real_R2", "real_pval")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"   # "positive", "negative", "all"

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

# I. Generation of the 117 quadrats:
####################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[,1] <- X1
  C[,2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[,1] %in% tx, ]
  C <- C[C[,2] %in% ty, ]
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(123)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# II. We generate the MEM variables with the DBMEM (PCNM) and generate nperm simulated species structured
# at 1) broad, 2) intermediate, and 3) fine scale.
#########################################################################################################

f4 <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
matB <- dnearneigh(lowlim, x = as.matrix(C), d1 = 0)
Y.DBMEM <- test.W.R2(nb = matB, xy = C, f = f4, t = lowlim, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM$best$MEM)

# MEM est la référence utilisée pour construire les variables réponses. Il sera aussi utilisé avec les
# autres matrices W en tant que matrice explicative.

spesimB <- vector("list", nperm)   # Liste de nperm espèces structur?es simul?es (Broad scale) 
spesimM <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Medium scale)
spesimF <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Fine scale)

n <- nrow(C)

# Générons nperm réalisation d'une espèce structur?e ? large, moyenne et fine ?chelle :

if(design == "regular"){
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 25] * 0.35) + (MEM[, 26] * 0.32) + (MEM[, 27] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 56] * 0.35) + (MEM[, 57] * 0.32) + (MEM[, 58] * 0.29) + rnorm(n, mean = 0, sd = 1)
  } 
} else {                  # Irregular design
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 19] * 0.35) + (MEM[, 20] * 0.32) + (MEM[, 21] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 36] * 0.35) + (MEM[, 37] * 0.32) + (MEM[, 38] * 0.29) + rnorm(n, mean = 0, sd = 1)
  }
}

# III. Test of the MEM.modsel.R function (with Sidak correction)
########################################
########################################
########################################

   # Broad scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimB[[i]]
   x <- MEM[, 1:3]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsB_modsel[2, 1005+i] <- R2adj
   resultsB_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsB_modsel[1, 5+i] <- memsel$pval
      resultsB_modsel[1, 1005+i] <- memsel$R2adj - resultsB_modsel[2, 1005+i]
   } else { resultsB_modsel[1, 5+i] <- 1 ; resultsB_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsB_modsel[1, 3] <- length(which(resultsB_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_modsel[1, 4] <- median(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   resultsB_modsel[1, 5] <- sd(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsB_modsel[2, 4] <- median(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   resultsB_modsel[2, 5] <- sd(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsB_modsel[3, 3] <- length(which(resultsB_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Broad", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsB_modsel, file = res_file_name, sep = "\t")

   # Medium scale
   ##############
   ##############

for(i in 1:nperm){

   Y <- spesimM[[i]]
   if(design == "regular") { x <- MEM[, 25:27]
   } else { x <- MEM[, 19:21] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_modsel[2, 1005+i] <- R2adj
   resultsM_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsM_modsel[1, 5+i] <- memsel$pval
      resultsM_modsel[1, 1005+i] <- memsel$R2adj - resultsM_modsel[2, 1005+i]
   } else { resultsM_modsel[1, 5+i] <- 1 ; resultsM_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsM_modsel[1, 3] <- length(which(resultsM_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_modsel[1, 4] <- median(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   resultsM_modsel[1, 5] <- sd(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsM_modsel[2, 4] <- median(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   resultsM_modsel[2, 5] <- sd(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsM_modsel[3, 3] <- length(which(resultsM_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Medium", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsM_modsel, file = res_file_name, sep = "\t")


   # Fine scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimF[[i]]
   if(design == "regular") { x <- MEM[, 56:58]
   } else { x <- MEM[, 36:38] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsF_modsel[2, 1005+i] <- R2adj
   resultsF_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsF_modsel[1, 5+i] <- memsel$pval
      resultsF_modsel[1, 1005+i] <- memsel$R2adj - resultsF_modsel[2, 1005+i]
   } else { resultsF_modsel[1, 5+i] <- 1 ; resultsF_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsF_modsel[1, 3] <- length(which(resultsF_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsF_modsel[1, 4] <- median(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   resultsF_modsel[1, 5] <- sd(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsF_modsel[2, 4] <- median(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   resultsF_modsel[2, 5] <- sd(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsF_modsel[3, 3] <- length(which(resultsF_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Fine", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsF_modsel, file = res_file_name, sep = "\t")

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

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

# Function to retrieve the pvalue of a lm() function:
# ***************************************************

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Construction of a results matrix for each scale and for both W matrix selection:
# ********************************************************************************

resultsB_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsB_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsM_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsM_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsF_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsF_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_modsel[, 1] <- c("sim", "real_R2", "real_pval")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"   # "positive", "negative", "all"

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "irregular"   # or "irregular"

nperm <- 1000

# I. Generation of the 117 quadrats:
####################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[,1] <- X1
  C[,2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[,1] %in% tx, ]
  C <- C[C[,2] %in% ty, ]
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(123)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# II. We generate the MEM variables with the DBMEM (PCNM) and generate nperm simulated species structured
# at 1) broad, 2) intermediate, and 3) fine scale.
#########################################################################################################

f4 <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
matB <- dnearneigh(lowlim, x = as.matrix(C), d1 = 0)
Y.DBMEM <- test.W.R2(nb = matB, xy = C, f = f4, t = lowlim, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM$best$MEM)

# MEM est la référence utilisée pour construire les variables réponses. Il sera aussi utilisé avec les
# autres matrices W en tant que matrice explicative.

spesimB <- vector("list", nperm)   # Liste de nperm espèces structur?es simul?es (Broad scale) 
spesimM <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Medium scale)
spesimF <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Fine scale)

n <- nrow(C)

# Générons nperm réalisation d'une espèce structur?e ? large, moyenne et fine ?chelle :

if(design == "regular"){
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 25] * 0.35) + (MEM[, 26] * 0.32) + (MEM[, 27] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 56] * 0.35) + (MEM[, 57] * 0.32) + (MEM[, 58] * 0.29) + rnorm(n, mean = 0, sd = 1)
  } 
} else {                  # Irregular design
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 19] * 0.35) + (MEM[, 20] * 0.32) + (MEM[, 21] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 36] * 0.35) + (MEM[, 37] * 0.32) + (MEM[, 38] * 0.29) + rnorm(n, mean = 0, sd = 1)
  }
}

# III. Test of the MEM.modsel.R function (with Sidak correction)
########################################
########################################
########################################

   # Broad scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimB[[i]]
   x <- MEM[, 1:3]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsB_modsel[2, 1005+i] <- R2adj
   resultsB_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsB_modsel[1, 5+i] <- memsel$pval
      resultsB_modsel[1, 1005+i] <- memsel$R2adj - resultsB_modsel[2, 1005+i]
   } else { resultsB_modsel[1, 5+i] <- 1 ; resultsB_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsB_modsel[1, 3] <- length(which(resultsB_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_modsel[1, 4] <- median(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   resultsB_modsel[1, 5] <- sd(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsB_modsel[2, 4] <- median(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   resultsB_modsel[2, 5] <- sd(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsB_modsel[3, 3] <- length(which(resultsB_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Broad", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsB_modsel, file = res_file_name, sep = "\t")

   # Medium scale
   ##############
   ##############

for(i in 1:nperm){

   Y <- spesimM[[i]]
   if(design == "regular") { x <- MEM[, 25:27]
   } else { x <- MEM[, 19:21] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_modsel[2, 1005+i] <- R2adj
   resultsM_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsM_modsel[1, 5+i] <- memsel$pval
      resultsM_modsel[1, 1005+i] <- memsel$R2adj - resultsM_modsel[2, 1005+i]
   } else { resultsM_modsel[1, 5+i] <- 1 ; resultsM_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsM_modsel[1, 3] <- length(which(resultsM_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_modsel[1, 4] <- median(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   resultsM_modsel[1, 5] <- sd(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsM_modsel[2, 4] <- median(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   resultsM_modsel[2, 5] <- sd(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsM_modsel[3, 3] <- length(which(resultsM_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Medium", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsM_modsel, file = res_file_name, sep = "\t")


   # Fine scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimF[[i]]
   if(design == "regular") { x <- MEM[, 56:58]
   } else { x <- MEM[, 36:38] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsF_modsel[2, 1005+i] <- R2adj
   resultsF_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsF_modsel[1, 5+i] <- memsel$pval
      resultsF_modsel[1, 1005+i] <- memsel$R2adj - resultsF_modsel[2, 1005+i]
   } else { resultsF_modsel[1, 5+i] <- 1 ; resultsF_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsF_modsel[1, 3] <- length(which(resultsF_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsF_modsel[1, 4] <- median(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   resultsF_modsel[1, 5] <- sd(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsF_modsel[2, 4] <- median(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   resultsF_modsel[2, 5] <- sd(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsF_modsel[3, 3] <- length(which(resultsF_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Corrected_INTERM_Fine", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsF_modsel, file = res_file_name, sep = "\t")

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
################################ UNCORRECTED VERSION OF MEM.MODSEL #####################################
### ************************************************************************************************ ###

# Functions:

source("lmp.R")
source("MEM.modsel - no_correction - adespatial - bestR2 - V3.R")

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

# Function to retrieve the pvalue of a lm() function:
# ***************************************************

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Construction of a results matrix for each scale and for both W matrix selection:
# ********************************************************************************

resultsB_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsB_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsM_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsM_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsF_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsF_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_modsel[, 1] <- c("sim", "real_R2", "real_pval")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"   # "positive", "negative", "all"

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

# I. Generation of the 117 quadrats:
####################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[,1] <- X1
  C[,2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[,1] %in% tx, ]
  C <- C[C[,2] %in% ty, ]
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(123)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# II. We generate the MEM variables with the DBMEM (PCNM) and generate nperm simulated species structured
# at 1) broad, 2) intermediate, and 3) fine scale.
#########################################################################################################

f4 <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
matB <- dnearneigh(lowlim, x = as.matrix(C), d1 = 0)
Y.DBMEM <- test.W.R2(nb = matB, xy = C, f = f4, t = lowlim, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM$best$MEM)

# MEM est la référence utilisée pour construire les variables réponses. Il sera aussi utilisé avec les
# autres matrices W en tant que matrice explicative.

spesimB <- vector("list", nperm)   # Liste de nperm espèces structur?es simul?es (Broad scale) 
spesimM <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Medium scale)
spesimF <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Fine scale)

n <- nrow(C)

# Générons nperm réalisation d'une espèce structur?e ? large, moyenne et fine ?chelle :

if(design == "regular"){
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 25] * 0.35) + (MEM[, 26] * 0.32) + (MEM[, 27] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 56] * 0.35) + (MEM[, 57] * 0.32) + (MEM[, 58] * 0.29) + rnorm(n, mean = 0, sd = 1)
  } 
} else {                  # Irregular design
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 19] * 0.35) + (MEM[, 20] * 0.32) + (MEM[, 21] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 36] * 0.35) + (MEM[, 37] * 0.32) + (MEM[, 38] * 0.29) + rnorm(n, mean = 0, sd = 1)
  }
}

# III. Test of the MEM.modsel.R function (with Sidak correction)
########################################
########################################
########################################

   # Broad scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimB[[i]]
   x <- MEM[, 1:3]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsB_modsel[2, 1005+i] <- R2adj
   resultsB_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsB_modsel[1, 5+i] <- memsel$pval
      resultsB_modsel[1, 1005+i] <- memsel$R2adj - resultsB_modsel[2, 1005+i]
   } else { resultsB_modsel[1, 5+i] <- 1 ; resultsB_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsB_modsel[1, 3] <- length(which(resultsB_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_modsel[1, 4] <- median(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   resultsB_modsel[1, 5] <- sd(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsB_modsel[2, 4] <- median(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   resultsB_modsel[2, 5] <- sd(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsB_modsel[3, 3] <- length(which(resultsB_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Broad", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsB_modsel, file = res_file_name, sep = "\t")

   # Medium scale
   ##############
   ##############

for(i in 1:nperm){

   Y <- spesimM[[i]]
   if(design == "regular") { x <- MEM[, 25:27]
   } else { x <- MEM[, 19:21] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_modsel[2, 1005+i] <- R2adj
   resultsM_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsM_modsel[1, 5+i] <- memsel$pval
      resultsM_modsel[1, 1005+i] <- memsel$R2adj - resultsM_modsel[2, 1005+i]
   } else { resultsM_modsel[1, 5+i] <- 1 ; resultsM_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsM_modsel[1, 3] <- length(which(resultsM_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_modsel[1, 4] <- median(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   resultsM_modsel[1, 5] <- sd(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsM_modsel[2, 4] <- median(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   resultsM_modsel[2, 5] <- sd(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsM_modsel[3, 3] <- length(which(resultsM_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Medium", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsM_modsel, file = res_file_name, sep = "\t")


   # Fine scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimF[[i]]
   if(design == "regular") { x <- MEM[, 56:58]
   } else { x <- MEM[, 36:38] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsF_modsel[2, 1005+i] <- R2adj
   resultsF_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsF_modsel[1, 5+i] <- memsel$pval
      resultsF_modsel[1, 1005+i] <- memsel$R2adj - resultsF_modsel[2, 1005+i]
   } else { resultsF_modsel[1, 5+i] <- 1 ; resultsF_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsF_modsel[1, 3] <- length(which(resultsF_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsF_modsel[1, 4] <- median(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   resultsF_modsel[1, 5] <- sd(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsF_modsel[2, 4] <- median(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   resultsF_modsel[2, 5] <- sd(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsF_modsel[3, 3] <- length(which(resultsF_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Fine", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsF_modsel, file = res_file_name, sep = "\t")

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

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

# Function to retrieve the pvalue of a lm() function:
# ***************************************************

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Construction of a results matrix for each scale and for both W matrix selection:
# ********************************************************************************

resultsB_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsB_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsM_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsM_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_modsel[, 1] <- c("sim", "real_R2", "real_pval")

resultsF_modsel <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsF_modsel) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_modsel[, 1] <- c("sim", "real_R2", "real_pval")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"   # "positive", "negative", "all"

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "irregular"   # or "irregular"

nperm <- 1000

# I. Generation of the 117 quadrats:
####################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[,1] <- X1
  C[,2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[,1] %in% tx, ]
  C <- C[C[,2] %in% ty, ]
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(123)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# II. We generate the MEM variables with the DBMEM (PCNM) and generate nperm simulated species structured
# at 1) broad, 2) intermediate, and 3) fine scale.
#########################################################################################################

f4 <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
matB <- dnearneigh(lowlim, x = as.matrix(C), d1 = 0)
Y.DBMEM <- test.W.R2(nb = matB, xy = C, f = f4, t = lowlim, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM$best$MEM)

# MEM est la référence utilisée pour construire les variables réponses. Il sera aussi utilisé avec les
# autres matrices W en tant que matrice explicative.

spesimB <- vector("list", nperm)   # Liste de nperm espèces structur?es simul?es (Broad scale) 
spesimM <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Medium scale)
spesimF <- vector("list", nperm)   # Liste de nperm esp?ces structur?es simul?es (Fine scale)

n <- nrow(C)

# Générons nperm réalisation d'une espèce structur?e ? large, moyenne et fine ?chelle :

if(design == "regular"){
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 25] * 0.35) + (MEM[, 26] * 0.32) + (MEM[, 27] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 56] * 0.35) + (MEM[, 57] * 0.32) + (MEM[, 58] * 0.29) + rnorm(n, mean = 0, sd = 1)
  } 
} else {                  # Irregular design
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.35) + (MEM[, 2] * 0.32) + (MEM[, 3] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 19] * 0.35) + (MEM[, 20] * 0.32) + (MEM[, 21] * 0.29) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 36] * 0.35) + (MEM[, 37] * 0.32) + (MEM[, 38] * 0.29) + rnorm(n, mean = 0, sd = 1)
  }
}

# III. Test of the MEM.modsel.R function (with Sidak correction)
########################################
########################################
########################################

   # Broad scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimB[[i]]
   x <- MEM[, 1:3]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsB_modsel[2, 1005+i] <- R2adj
   resultsB_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsB_modsel[1, 5+i] <- memsel$pval
      resultsB_modsel[1, 1005+i] <- memsel$R2adj - resultsB_modsel[2, 1005+i]
   } else { resultsB_modsel[1, 5+i] <- 1 ; resultsB_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsB_modsel[1, 3] <- length(which(resultsB_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_modsel[1, 4] <- median(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   resultsB_modsel[1, 5] <- sd(na.omit(as.numeric(resultsB_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsB_modsel[2, 4] <- median(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   resultsB_modsel[2, 5] <- sd(na.omit(as.numeric(resultsB_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsB_modsel[3, 3] <- length(which(resultsB_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Broad", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsB_modsel, file = res_file_name, sep = "\t")

   # Medium scale
   ##############
   ##############

for(i in 1:nperm){

   Y <- spesimM[[i]]
   if(design == "regular") { x <- MEM[, 25:27]
   } else { x <- MEM[, 19:21] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_modsel[2, 1005+i] <- R2adj
   resultsM_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsM_modsel[1, 5+i] <- memsel$pval
      resultsM_modsel[1, 1005+i] <- memsel$R2adj - resultsM_modsel[2, 1005+i]
   } else { resultsM_modsel[1, 5+i] <- 1 ; resultsM_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsM_modsel[1, 3] <- length(which(resultsM_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_modsel[1, 4] <- median(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   resultsM_modsel[1, 5] <- sd(na.omit(as.numeric(resultsM_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsM_modsel[2, 4] <- median(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   resultsM_modsel[2, 5] <- sd(na.omit(as.numeric(resultsM_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsM_modsel[3, 3] <- length(which(resultsM_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Medium", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsM_modsel, file = res_file_name, sep = "\t")


   # Fine scale
   #############
   #############

for(i in 1:nperm){

   Y <- spesimF[[i]]
   if(design == "regular") { x <- MEM[, 56:58]
   } else { x <- MEM[, 36:38] }
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsF_modsel[2, 1005+i] <- R2adj
   resultsF_modsel[3, 5+i] <- lmp(lm)

   memsel <- MEM.modsel(Y, C, autocor = "positive")
   if(class(memsel) != "NULL"){
      resultsF_modsel[1, 5+i] <- memsel$pval
      resultsF_modsel[1, 1005+i] <- memsel$R2adj - resultsF_modsel[2, 1005+i]
   } else { resultsF_modsel[1, 5+i] <- 1 ; resultsF_modsel[1, 1005+i] <- NA }
}

# Power, median and sd of R2adj:
################################

   resultsF_modsel[1, 3] <- length(which(resultsF_modsel[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsF_modsel[1, 4] <- median(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   resultsF_modsel[1, 5] <- sd(na.omit(as.numeric(resultsF_modsel[1, c(1006:(nperm + 1005))])))
   # Median and SD of the real R2adj
   resultsF_modsel[2, 4] <- median(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   resultsF_modsel[2, 5] <- sd(na.omit(as.numeric(resultsF_modsel[2, c(1006:(nperm + 1005))])))
   # Correct significance detection rate
   resultsF_modsel[3, 3] <- length(which(resultsF_modsel[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power_MEM.modsel_Uncorrected_INTERM_Fine", paste(design, ".txt", sep = ""), sep = "_")
write.table(resultsF_modsel, file = res_file_name, sep = "\t")
