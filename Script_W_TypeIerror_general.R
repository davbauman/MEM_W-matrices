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

source("test.W.R2.R")

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
results[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed", "DBmax"), each = 4), "DBMEM_PCNM")
results[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), times = 7), "1-(D/4t)^2")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Sampling design:

design <- "random"   # "random" or "clustered"

nperm <- 1000

# Generation of the 117 quadrats:
#################################

xy <- expand.grid(x = seq(1, 150, 1), y = seq(1, 75, 1))

set.seed(1)

if (design == "clustered") {
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  x1 <- runif (39, min = sample(c(1:15), 1), max = sample(c(36:42), 1))
  y1 <- runif (39, min = sample(c(39:51), 1), max = sample(c(66:75), 1))
  x2 <- runif (39, min = sample(c(54:63), 1), max = sample(c(81:93), 1))
  y2 <- runif (39, min = sample(c(36:49), 1), max = sample(c(66:75), 1))
  x3 <- runif (39, min = sample(c(99:114), 1), max = sample(c(135:148), 1))
  y3 <- runif (39, min = sample(c(1:15), 1), max = sample(c(30:45), 1))
  
  C[, 1] <- rbind(x1, x2, x3)
  C[, 2] <- rbind(y1, y2, y3)
  
} else {          # design = "random"
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  C[, 1] <- runif (117, min = 1, max = 148)
  C[, 2] <- runif (117, min = 1, max = 75) 
  
}

# Attribute each sampled point to one of the grid cells:
# ******************************************************

grid.size <- 1
tri <- c()
for (k in 1:nrow(C)) {
  x <- floor((C[k, 1]) / (grid.size))
  y <- floor((C[k, 2]) / (grid.size))
  N <- y * 150 + x + 1
  tri <- c(tri, N)
}

# We can only have one sampled point by grid cell:
# ************************************************

sort <- sort(tri)
control <- length(levels(as.factor(sort)))
while (control != nrow(C)) {
  cible <- c()
  for (k in 1:(length(sort)-1)) if (sort[k+1] == sort[k]) cible <- c(cible, k+1)
  for (k in cible) {
    if (length(which(seq(50, 1250, 50) == sort[k])) == 0) {
      sort[k] = sort[k] + 1
    } else {
      sort[k] = sort[k] - 1
    }
  }
  control <- length(levels(as.factor(sort)))
}
# We rearange 'C' so that all sampled point are in different grid cells ('sort'):
C <- xy[sort, ]
xy.d1 <- dist(C)

# II. Connectivity matrices (B):
################################

Y.del <- tri2nb(C) 
Y.gab <- graph2nb(gabrielneigh(as.matrix(C), nnmult = 5), sym = TRUE)
Y.rel <- graph2nb(relativeneigh(as.matrix(C), nnmult = 5), sym = TRUE)
Y.mst <- mst.nb(xy.d1)

# Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
uplim <- max(xy.d1)
thresh <- seq(lowlim, uplim, le = 3)   # 3 tested distances
Y.listDB <- lapply(thresh, dnearneigh, x = as.matrix(C), d1 = 0)

# Weighting functions and fixed parametres:
# *****************************************
f1 <- function (D, dmax)    {1-(D/dmax)}        # Linear function
f2 <- function (D, dmax, y) {1-(D/dmax)^y}      # Concave-down function
f3 <- function (D, dmax, y) {1/(D/dmax)^y}      # Concave-up function
f4 <- function (D, t)       {1-(D/(4*t))^2}     # PCNM criterion

max.del <- max(unlist(nbdists(Y.del, as.matrix(C)))) 
max.gab <- max(unlist(nbdists(Y.gab, as.matrix(C))))
max.rel <- max(unlist(nbdists(Y.rel, as.matrix(C)))) 
max.mst <- max(unlist(nbdists(Y.mst, as.matrix(C))))

nbdist <- lapply(Y.listDB, coords = as.matrix(C), nbdists)
unlist <- lapply(nbdist, unlist)
max.DB.list <- lapply(unlist, max)

# Generation of MEM variables:
# ****************************
# del:
# ****
Y.del.MEM <- test.W.R2(nb = Y.del, xy = C, style = style, MEM.autocor = MEM_model)
Y.del.MEM.f1 <- test.W.R2(nb = Y.del, xy = C, f = f1, dmax = max.del, style = style, 
                          MEM.autocor = MEM_model)
# gab:
# ****
Y.gab.MEM <- test.W.R2(nb = Y.gab, xy = C, style = style, MEM.autocor = MEM_model)
Y.gab.MEM.f1 <- test.W.R2(nb = Y.gab, xy = C, f = f1, dmax = max.gab, style = style, 
                          MEM.autocor = MEM_model)
# rel:
# ****
Y.rel.MEM <- test.W.R2(nb = Y.rel, xy = C, style = style, MEM.autocor = MEM_model)
Y.rel.MEM.f1 <- test.W.R2(nb = Y.rel, xy = C, f = f1, dmax = max.rel, style = style, 
                          MEM.autocor = MEM_model)
# mst:
# ****
Y.mst.MEM <- test.W.R2(nb = Y.mst, xy = C, style = style, MEM.autocor = MEM_model)
Y.mst.MEM.f1 <- test.W.R2(nb = Y.mst, xy = C, f = f1, dmax = max.mst, style = style, 
                          MEM.autocor = MEM_model)
# DB:
# ***
Y.DB1.MEM <- test.W.R2(nb = Y.listDB[[1]], xy = C, style = style, 
                       MEM.autocor = MEM_model)
Y.DB1.MEM.f1 <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f1, dmax = max.DB.list[[1]], 
                          style = style, MEM.autocor = MEM_model)
Y.DB2.MEM <- test.W.R2(nb = Y.listDB[[2]], xy = C, style = style, 
                       MEM.autocor = MEM_model)
Y.DB2.MEM.f1 <- test.W.R2(nb = Y.listDB[[2]], xy = C, f = f1, dmax = max.DB.list[[2]], 
                          style = style, MEM.autocor = MEM_model)
Y.DB3.MEM <- test.W.R2(nb = Y.listDB[[3]], xy = C, style = style, 
                       MEM.autocor = MEM_model)
Y.DB3.MEM.f1 <- test.W.R2(nb = Y.listDB[[3]], xy = C, f = f1, dmax = max.DB.list[[3]], 
                          style = style, MEM.autocor = MEM_model)
# DBMEM with PCNM criteria:
# *************************
Y.DB.PCNM <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f4, t = lowlim, style = style, 
                       MEM.autocor = MEM_model)

# Simulation procedure:
#######################

for (i in 1:nperm) {   
  
  set.seed(i)
  
  if (framework == "univariate") {
    
    Y <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"             # Random (uniform)
    #   Y <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" # Random (normal)
    #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                        # Exponential (1)
    #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                    # Exponential cubed
    
  } else {
    
    Y <- matrix(ncol = 5, nrow = nrow(C))
    for (b in 1:ncol(Y)) {
      Y[, b] <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"    
      #     Y[, b] <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
      #     Y[, b] <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
      #     Y[, b] <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3" 
    }                  
  }
  
  # Generation of the remaining W matrices (with concdown and concup functions):
  # ****************************************************************************
  # del:
  # ****
  Y.del.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.del, y = 2:5)
  Y.del.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.del, y = 1:5)
  # gab:
  # ****
  Y.gab.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.gab, y = 2:5)
  Y.gab.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.gab, y = 1:5)
  # rel:
  # ****
  Y.rel.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.rel, y = 2:5)
  Y.rel.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.rel, y = 1:5)
  # mst:
  # ****
  Y.mst.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.mst, y = 2:5)
  Y.mst.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.mst, y = 1:5)
  # DB:
  # ***
  Y.DB1.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, 
                            style = style, MEM.autocor = MEM_model, f = f2, 
                            dmax = max.DB.list[[1]], y = 2:5)
  Y.DB1.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.DB.list[[1]],
                            y = 1:5)
  
  Y.DB2.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[2]], 
                            y = 2:5)
  Y.DB2.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.DB.list[[2]],
                            y = 1:5)
  
  Y.DB3.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[3]], 
                            y = 2:5)
  Y.DB3.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, style = style, 
                            MEM.autocor = MEM_model, f = f3, dmax = max.DB.list[[3]],
                            y = 1:5)
  
  # Significance test and MEM variable selection (forward selection with double stopping criterion):
  # ************************************************************************************************
  
  # del
  R2adj <- RsquareAdj(rda(Y, Y.del.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.del.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.del.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM$best$MEM[, c(sign)]
      results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[1, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[1, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[1, i+1005] <- NA
    }  
  R2adj <- RsquareAdj(rda(Y, Y.del.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.del.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f1$best$MEM[, c(sign)]
      results[2, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[2, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[2, i+5] <- 1  
      results[2, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.del.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.del.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f2$best$MEM[, c(sign)]
      results[3, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[3, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[3, i+5] <- 1   
      results[3, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.del.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.del.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.del.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f3$best$MEM[, c(sign)]
      results[4, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[4, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[4, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[4, i+1005] <- NA
    }
  
  # gab
  R2adj <- RsquareAdj(rda(Y, Y.gab.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.gab.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM$best$MEM[, c(sign)]
      results[5, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[5, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[5, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[5, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.gab.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f1$best$MEM[, c(sign)]
      results[6, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[6, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[6, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[6, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.gab.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f2$best$MEM[, c(sign)]
      results[7, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[7, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[7, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[7, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.gab.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.gab.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.gab.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f3$best$MEM[, c(sign)]
      results[8, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[8, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[8, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[8, i+1005] <- NA
    }
  
  # rel
  R2adj <- RsquareAdj(rda(Y, Y.rel.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.rel.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM$best$MEM[, c(sign)]
      results[9, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[9, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[9, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[9, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.rel.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f1$best$MEM[, c(sign)]
      results[10, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[10, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[10, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[10, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.rel.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f2$best$MEM[, c(sign)]
      results[11, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[11, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[11, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[11, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.rel.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.rel.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.rel.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f3$best$MEM[, c(sign)]
      results[12, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[12, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[12, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[12, i+1005] <- NA
    }
  
  # mst
  R2adj <- RsquareAdj(rda(Y, Y.mst.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.mst.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM$best$MEM[, c(sign)]
      results[13, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[13, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[13, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[13, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.mst.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f1$best$MEM[, c(sign)]
      results[14, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[14, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[14, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[14, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.mst.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f2$best$MEM[, c(sign)]
      results[15, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[15, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[15, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[15, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.mst.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.mst.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.mst.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f3$best$MEM[, c(sign)]
      results[16, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[16, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[16, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[16, i+1005] <- NA
    }
  
  # DBmin
  R2adj <- RsquareAdj(rda(Y, Y.DB1.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB1.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB1.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM$best$MEM[, c(sign)]
      results[17, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[17, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[17, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[17, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB1.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB1.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB1.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f1$best$MEM[, c(sign)]
      results[18, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[18, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[18, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[18, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB1.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB1.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB1.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f2$best$MEM[, c(sign)]
      results[19, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[19, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[19, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[19, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB1.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB1.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB1.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f3$best$MEM[, c(sign)]
      results[20, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[20, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[20, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[20, i+1005] <- NA
    }
  
  # DBmin
  R2adj <- RsquareAdj(rda(Y, Y.DB2.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB2.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB2.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM$best$MEM[, c(sign)]
      results[21, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[21, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[21, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[21, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB2.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB2.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB2.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f1$best$MEM[, c(sign)]
      results[22, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[22, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[22, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[22, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB2.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB2.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB2.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f2$best$MEM[, c(sign)]
      results[23, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[23, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[23, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[23, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB2.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB2.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB2.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f3$best$MEM[, c(sign)]
      results[24, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[24, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[24, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[24, i+1005] <- NA
    }
  
  # DBmin
  R2adj <- RsquareAdj(rda(Y, Y.DB3.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB3.MEM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB3.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM$best$MEM[, c(sign)]
      results[25, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[25, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[25, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[25, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB3.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB3.MEM.f1$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB3.MEM.f1$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f1$best$MEM[, c(sign)]
      results[26, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[26, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[26, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[26, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB3.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB3.MEM.f2$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB3.MEM.f2$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f2$best$MEM[, c(sign)]
      results[27, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[27, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[27, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[27, i+1005] <- NA
    }
  R2adj <- RsquareAdj(rda(Y, Y.DB3.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB3.MEM.f3$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB3.MEM.f3$best$MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f3$best$MEM[, c(sign)]
      results[28, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[28, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } } else{ 
      results[28, i+5] <- 1   # p-val made equal to 1 (we did not even enter the fwd sel)
      results[28, i+1005] <- NA
    }
  
  # DBMEM (with PCNM criteria: B = give.thresh(xy.d1) and A = f4)
  R2adj <- RsquareAdj(rda(Y, Y.DB.PCNM$best$MEM))$adj.r.squared
  if (anova.cca(rda(Y, Y.DB.PCNM$best$MEM))$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(Y, Y.DB.PCNM$best$MEM, adjR2thresh = R2adj, nperm = 999),
                       TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB.PCNM$best$MEM[, c(sign)]
      results[29, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[29, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } 
  } else { 
    results[29, i+5] <- 1   # p-val made equal to 1 
    results[29, i+1005] <- 0
  }
  
} # End of the 'for' simulation loop

# Type I error, median and sd of R2adj:
#######################################

for (i in 1:nrow(results)) {
  results[i, 3] <- length(which(results[i, c(6:(nperm + 5))] <= 0.05)) / nperm
  results[i, 4] <- median(na.omit(as.numeric(results[i, c(1006:(nperm + 1005))])))
  results[i, 5] <- sd(na.omit(as.numeric(results[i, c(1006:(nperm + 1005))])))
}

# Output of the results:
# **********************

res_file_name <- paste("Results_tIerror_general", framework, ran,
                       paste(design, ".txt", sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")
