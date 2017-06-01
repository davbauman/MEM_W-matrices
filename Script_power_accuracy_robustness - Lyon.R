##### Script du test de puissance, précision et robustesse de différentes #####
### matrices W vis-à-vis du type de design d'échantillonnage, de l'échelle ####
### spatiale de la structuration, d'un point de vue de la population réelle ###
###################### et de l'échantillon réel. ##############################
# *************************************************************************** #

# rm(list=ls())

# Usefull packages and functions:
# *******************************

library(vegan)
library(adespatial)
library(spdep)
library(ade4)

source("lmp.R")

# Definition of the simulation parameters:
# ****************************************

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"   # Either "positive" or "negative"

# Sampling design:
design <- "clustered"    # Either "clustered" or "random"

nperm <- 1000

# Structuring Intensity (low or high):
a <- 0.35   # 0.35 or 0.55

style <- "B"             # Either "B" or "W"

if (style == "B") { 
    source("test.W.R2_styleB.R") 
} else source("test.W.R2_styleW.R") 

# Construction of a results matrix for each scale:
# ************************************************
# Les résultats _pop font référence à ce qui se passe en gardant 'y' fixe et
# en faisant varier le design d'échantillonnage. 
# Les résultats _sub font référence à ce qui se passe en maintenant le design
# fixe et en faisant varier y (ce qu'on faisait avant).
# Chaque matrice contient en ligne les différentes matrices W, et en colonne :
# 1000 valeurs de p-values, puis 1000 valeurs R2_pop, R2_sub, R2_W, 
# deltaR2_pop, deltaR2_sub, deltaR2_subpop et les médianes et écart-types des
# trois deltaR2.

resultsB_pop <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsB_pop) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsB_pop[, 1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "pvalReal", "column_count")
resultsB_pop[c(1:32), 2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/")
resultsB_pop[33, ] <- c(1:ncol(resultsB_pop))

resultsM_pop <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsM_pop) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_pop[, 1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "pvalReal", "column_count")
resultsM_pop[c(1:32), 2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/")
resultsM_pop[33,] <- c(1:ncol(resultsM_pop))

resultsB_sub <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsB_sub) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsB_sub[, 1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "pvalReal", "column_count")
resultsB_sub[c(1:32), 2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/")
resultsB_sub[33, ] <- c(1:ncol(resultsB_sub))

resultsM_sub <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsM_sub) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_sub[, 1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "pvalReal", "column_count")
resultsM_sub[c(1:32), 2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/")
resultsM_sub[33, ] <- c(1:ncol(resultsM_sub))

### I. The response remains unchanged and the sampling design varies:
#####################################################################
#####################################################################
#####################################################################

# The MEM are built for a full grid (50 x 25 cells):
# **************************************************

xy <- expand.grid(x = seq(1, 50, 1), y = seq(1, 25, 1))
# plot(xy, cex = 1)

nb <- cell2nb(nrow = 25, ncol = 50, "queen")
nb2 <- nb2listw(nb)
MEM <- scores.listw(nb2, MEM.autocor = MEM_model)

# xy_attr <- attr(nb, "region.id")
# XY <- matrix(as.integer(unlist(strsplit(xy_attr, ":"))), ncol = 2, byrow = TRUE)
# plot(nb, XY)

# To know from where and in which direction the cells are considered when building MEM
# ************************************************************************************

# s.label(xy, neig = nb2neig(nb), clab = 0.5)

# Creation of the response variable:
# **********************************
   # Creation of 'y_spa' and 'y_noise' and standardisation:
   # ******************************************************

set.seed(1)

y_spa_broad <- MEM[, 1] + MEM[, 3] + MEM[, 5]
y_spa_med <- MEM[, 211] + MEM[, 212] + MEM[, 215]
y_noise <- rnorm(n = nrow(MEM), mean = 0, sd = 1)

y_spa_broad_st <- scale(y_spa_broad)
y_spa_med_st <- scale(y_spa_med)
y_noise_st <- scale(y_noise)

# par(mfrow = c(1, 3))
# for(i in c(1, 3, 5)) s.value(xy, MEM[,i])
# for(i in c(211, 212, 215)) s.value(xy, MEM[,i])

   # Creation of the response variable 'y' at the whole population level (pop):
   # **************************************************************************

y_broad <- (a * y_spa_broad_st) + ((1-a) * y_noise_st)
y_med <- (a * y_spa_med_st) + ((1-a) * y_noise_st)

# s.value(xy, y_broad)

R2_pop_broad <- cor(y_broad, y_spa_broad_st)^2
R2_pop_med <- cor(y_med, y_spa_med_st)^2

resultsB_pop[30, c(1010:2009, 3010:4009)] <- R2_pop_broad
resultsM_pop[30, c(1010:2009, 3010:4009)] <- R2_pop_med

# Begining of the simulation process:
# ***********************************

# The simulation is run first at the broad scale, and then at the medium scale.
# Since the simulation at the medium scale have to be done exactly on the same
# subsampled sites, and therefore on the same 'C' matrices and same 'MEMsub', we 
# keep all the 'C' and 'MEMsub' computed at the broad scale in two lists, and we
# reuse them at the medium scale to spare time.

C.list <- vector("list", nperm)
MEMsub.list <- vector("list", nperm)

   ######################
   ######################
   ### I. Broad scale ###
   ######################
   ######################

for (i in 1:nperm) {

# Sampling scheme:
# ****************

if (design == "clustered") {
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  set.seed(i)
  x1 <- runif(39, min = 2, max = 14)
  y1 <- runif(39, min = 17, max = 25)
  x2 <- runif(39, min = 19, max = 34)
  y2 <- runif(39, min = 15, max = 25)
  x3 <- runif(39, min = 30, max = 50)
  y3 <- runif(39, min = 1, max = 10)
  
  C[,1] <- rbind(x1, x2, x3)
  C[,2] <- rbind(y1, y2, y3)

#  x1 <- rnorm(39, mean = 0.3, sd = 0.07)
#  y1 <- rnorm(39, mean = 0.7, sd = 0.07)
#  x2 <- rnorm(39, mean = 0.65, sd = 0.07)
#  y2 <- rnorm(39, mean = 0.7, sd = 0.07)
#  x3 <- rnorm(39, mean = 1.1, sd = 0.07)
#  y3 <- rnorm(39, mean = 0.4, sd = 0.07)
#  C[,1] <- rbind(x1, x2, x3)
#  C[,2] <- rbind(y1, y2, y3)
#  C[, 1] <- C[, 1] / min(C[, 1])
#  C[, 2] <- C[, 2] / min(C[, 2])
#  C[, 1] <- (C[, 1] / max(C[, 1])) * 50
#  C[, 2] <- (C[, 2] / max(C[, 2])) * 25
  
  xy.d1 <- dist(C)
  
} else {          # design = "random"
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  set.seed(i)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25) 
  
  xy.d1 <- dist(C)
}

# Attribute each sampled point to one of the grid cells:
# ******************************************************

# /25 = taille d'un côté du quadrat ; pour que la numérotation des quadrats se
# fasse de gauche à droite en partant du coin inférieur gauche de la grille :
# N <- y * (nb horizontal de quadrats) + x + 1 ; pour que la numérotation des 
# quadrats se fasse de bas en haut en partant du coin inférieur gauche :
# N <- x * (nb vertical de quadrats) + y + 1

grid.size <- 1
tri <- c()
for (k in 1:nrow(C)) {
  x <- floor(C[k, 1]/grid.size)
  y <- floor(C[k, 2]/grid.size)
  N <- y * 50 + x + 1
  tri <- c(tri, N)
}

# We can only have one sampled point by grid cell:
# ************************************************
# We sort the cells and 1) add 1 to a cell if it has the same number than the one 
# before it, given that the cell is not at the right border. Otherwise, we substract
# 1 to it. We repeat the procedure until all the sampled point are in different cells:

sort <- sort(tri)
control <- length(levels(as.factor(sort)))
while (control != nrow(C)) {
  cible <- c()
  for(k in 1:(length(sort)-1)) if (sort[k+1] == sort[k]) cible <- c(cible, k+1)
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
C.list[[i]] <- C

# We keep the lines of MEM that correspond to the sampled cells ('tri'):
# **********************************************************************
MEMsub <- MEM[sort, ]
MEMsub.list[[i]] <- MEMsub

# We sample the response variable within the sampled cells ('y_sub'):
# *******************************************************************
y_sub <- y_broad[sort, ]

# Real p-value and R2_sub:
# ************************
x <- MEMsub[, c(1, 3, 5)]
lm <- lm(y_sub ~ ., data = x)
R2_sub <- summary(lm)$adj.r.squared
resultsB_pop[31, c(2009+i, 3009+i)] <- R2_sub
resultsB_pop[32, 9+i] <- lmp(lm)

# Visualisation:
# **************
# par(mfrow = c(1, 3))
# for(k in c(1, 3, 5)) s.value(C, MEMsub[, k])
# for(k in c(211, 212, 215)) s.value(C, MEMsub[, k])

# Construction of the different W matrices:
# #########################################

   # Connectivity matrices (B):
   # **************************

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
f3 <- function (D, y)       {1/(D^y)}           # Concave-up function
f4 <- function (D, t)       {1-(D/(4*t))^2}     # PCNM criterion

max.del <- max(unlist(nbdists(Y.del, as.matrix(C)))) 
max.gab <- max(unlist(nbdists(Y.gab, as.matrix(C))))
max.rel <- max(unlist(nbdists(Y.rel, as.matrix(C)))) 
max.mst <- max(unlist(nbdists(Y.mst, as.matrix(C))))

nbdist <- lapply(Y.list10dnn, coords = as.matrix(C), nbdists)
unlist <- lapply(nbdist, unlist)
max.DB.list <- lapply(unlist, max)

# Generation of MEM variables:
# ****************************
# del:
# ****
Y.del.MEM <- test.W.R2(nb = Y.del, xy = C, MEM.autocor = MEM_model)
Y.del.MEM.f1 <- test.W.R2(nb = Y.del, xy = C, f = f1, dmax = max.del, 
                          MEM.autocor = MEM_model)
Y.del.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, MEM.autocor = MEM_model, 
                          f = f2, dmax = max.del, y = 2:10)
Y.del.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, MEM.autocor = MEM_model, 
                          f = f3, y = 2:10)
# gab:
# ****
Y.gab.MEM <- test.W.R2(nb = Y.gab, xy = C, MEM.autocor = MEM_model)
Y.gab.MEM.f1 <- test.W.R2(nb = Y.gab, xy = C, f = f1, dmax = max.gab, 
                          MEM.autocor = MEM_model)
Y.gab.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, MEM.autocor = MEM_model, 
                          f = f2, dmax = max.gab, y = 2:10)
Y.gab.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, MEM.autocor = MEM_model, 
                          f = f3, y = 2:10)
# rel:
# ****
Y.rel.MEM <- test.W.R2(nb = Y.rel, xy = C, MEM.autocor = MEM_model)
Y.rel.MEM.f1 <- test.W.R2(nb = Y.rel, xy = C, f = f1, dmax = max.rel, 
                          MEM.autocor = MEM_model)
Y.rel.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, MEM.autocor = MEM_model, 
                          f = f2, dmax = max.rel, y = 2:10)
Y.rel.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, MEM.autocor = MEM_model, 
                          f = f3, y = 2:10)
# mst:
# ****
Y.mst.MEM <- test.W.R2(nb = Y.mst, xy = C, MEM.autocor = MEM_model)
Y.mst.MEM.f1 <- test.W.R2(nb = Y.mst, xy = C, f = f1, dmax = max.mst, 
                          MEM.autocor = MEM_model)
Y.mst.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, MEM.autocor = MEM_model, 
                            f = f2, dmax = max.mst, y = 2:10)
Y.mst.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, MEM.autocor = MEM_model, 
                            f = f3, y = 2:10)
# DB:
# ***
Y.DB1.MEM <- test.W.R2(nb = Y.listDB[[1]], xy = C, MEM.autocor = MEM_model)
Y.DB1.MEM.f1 <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f1, dmax = max.DB.list[[1]], 
                          MEM.autocor = MEM_model)
Y.DB1.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[1]], 
                          y = 2:10)
Y.DB1.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, 
                          MEM.autocor = MEM_model, f = f3, y = 2:10)

Y.DB2.MEM <- test.W.R2(nb = Y.listDB[[2]], xy = C, MEM.autocor = MEM_model)
Y.DB2.MEM.f1 <- test.W.R2(nb = Y.listDB[[2]], xy = C, f = f1, dmax = max.DB.list[[2]], 
                          MEM.autocor = MEM_model)
Y.DB2.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[2]], 
                          y = 2:10)
Y.DB2.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, 
                          MEM.autocor = MEM_model, f = f3, y = 2:10)

Y.DB3.MEM <- test.W.R2(nb = Y.listDB[[3]], xy = C, MEM.autocor = MEM_model)
Y.DB3.MEM.f1 <- test.W.R2(nb = Y.listDB[[3]], xy = C, f = f1, dmax = max.DB.list[[3]], 
                          MEM.autocor = MEM_model)
Y.DB3.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[3]], 
                          y = 2:10)
Y.DB3.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, 
                          MEM.autocor = MEM_model, f = f3, y = 2:10)

# DBMEM with PCNM criteria:
# *************************
Y.DB.PCNM <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f4, t = lowlim, 
                       MEM.autocor = MEM_model)

# Significance test and MEM var. selection (fwd.sel with double stopping criterion):
# **********************************************************************************
# **********************************************************************************
# del
R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.del.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.del.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[1, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[1, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[1, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[1, i+9] <- 1     # p-val set at a value of 1 (for power computation)   
  resultsB_pop[1, i+1009] <- NA
  resultsB_pop[1, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.del.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.del.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[2, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[2, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[2, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[2, i+9] <- 1         
  resultsB_pop[2, i+1009] <- NA
  resultsB_pop[2, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.del.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.del.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[3, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[3, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[3, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[3, i+9] <- 1         
  resultsB_pop[3, i+1009] <- NA
  resultsB_pop[3, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.del.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.del.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[4, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[4, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[4, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[4, i+9] <- 1         
  resultsB_pop[4, i+1009] <- NA
  resultsB_pop[4, i+2009] <- NA
}

# gab
R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.gab.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.gab.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[5, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[5, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[5, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[5, i+9] <- 1         
  resultsB_pop[5, i+1009] <- NA
  resultsB_pop[5, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.gab.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.gab.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[6, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[6, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[6, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[6, i+9] <- 1         
  resultsB_pop[6, i+1009] <- NA
  resultsB_pop[6, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.gab.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.gab.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[7, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[7, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[7, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[7, i+9] <- 1         
  resultsB_pop[7, i+1009] <- NA
  resultsB_pop[7, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.gab.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.gab.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[8, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[8, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[8, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[8, i+9] <- 1         
  resultsB_pop[8, i+1009] <- NA
  resultsB_pop[8, i+2009] <- NA
}

# rel
R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.rel.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.rel.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[9, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[9, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[9, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[9, i+9] <- 1         
  resultsB_pop[9, i+1009] <- NA
  resultsB_pop[9, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.rel.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.rel.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[10, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[10, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[10, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[10, i+9] <- 1         
  resultsB_pop[10, i+1009] <- NA
  resultsB_pop[10, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.rel.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.rel.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[11, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[11, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[11, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[11, i+9] <- 1         
  resultsB_pop[11, i+1009] <- NA
  resultsB_pop[11, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.rel.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.rel.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[12, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[12, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[12, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[12, i+9] <- 1         
  resultsB_pop[12, i+1009] <- NA
  resultsB_pop[12, i+2009] <- NA
}

# mst
R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.mst.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.mst.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[13, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[13, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[13, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[13, i+9] <- 1         
  resultsB_pop[13, i+1009] <- NA
  resultsB_pop[13, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.mst.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.mst.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[14, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[14, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[14, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[14, i+9] <- 1         
  resultsB_pop[14, i+1009] <- NA
  resultsB_pop[14, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.mst.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.mst.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[15, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[15, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[15, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[15, i+9] <- 1         
  resultsB_pop[15, i+1009] <- NA
  resultsB_pop[15, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.mst.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.mst.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[16, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[16, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[16, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[16, i+9] <- 1         
  resultsB_pop[16, i+1009] <- NA
  resultsB_pop[16, i+2009] <- NA
}

# DB1
R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB1.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB1.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[17, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[17, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[17, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[17, i+9] <- 1         
  resultsB_pop[17, i+1009] <- NA
  resultsB_pop[17, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB1.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB1.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[18, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[18, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[18, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[18, i+9] <- 1         
  resultsB_pop[18, i+1009] <- NA
  resultsB_pop[18, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB1.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB1.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[19, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[19, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[19, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[19, i+9] <- 1         
  resultsB_pop[19, i+1009] <- NA
  resultsB_pop[19, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB1.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB1.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[20, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[20, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[20, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[20, i+9] <- 1         
  resultsB_pop[20, i+1009] <- NA
  resultsB_pop[20, i+2009] <- NA
}

# DB2
R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB2.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB2.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[21, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[21, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[21, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[21, i+9] <- 1         
  resultsB_pop[21, i+1009] <- NA
  resultsB_pop[21, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB2.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB2.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[22, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[22, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[22, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[22, i+9] <- 1         
  resultsB_pop[22, i+1009] <- NA
  resultsB_pop[22, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB2.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB2.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[23, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[23, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[23, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[23, i+9] <- 1         
  resultsB_pop[23, i+1009] <- NA
  resultsB_pop[23, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB2.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB2.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[24, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[24, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[24, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[24, i+9] <- 1         
  resultsB_pop[24, i+1009] <- NA
  resultsB_pop[24, i+2009] <- NA
}

# DB3
R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB3.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB3.MEM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[25, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[25, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[25, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[25, i+9] <- 1         
  resultsB_pop[25, i+1009] <- NA
  resultsB_pop[25, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f1$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB3.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f1$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB3.MEM.f1$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[26, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[26, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[26, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[26, i+9] <- 1      
  resultsB_pop[26, i+1009] <- NA
  resultsB_pop[26, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f2$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB3.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f2$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB3.MEM.f2$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[27, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[27, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[27, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[27, i+9] <- 1      
  resultsB_pop[27, i+1009] <- NA
  resultsB_pop[27, i+2009] <- NA
}
R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f3$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB3.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f3$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB3.MEM.f3$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[28, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[28, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[28, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[28, i+9] <- 1      
  resultsB_pop[28, i+1009] <- NA
  resultsB_pop[28, i+2009] <- NA
}

# DB_PCNM
R2adj <- RsquareAdj(rda(y_sub, Y.DB.PCNM$best$MEM))$adj.r.squared
if (anova.cca(rda(y_sub, Y.DB.PCNM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(y_sub, Y.DB.PCNM$best$MEM, 
                                         adjR2thresh = R2adj, nperm = 999), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DB.PCNM$best$MEM[, c(sign)]
    R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
    resultsB_pop[29, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                    permutations = 9999))$Pr[1]
    resultsB_pop[29, i+1009] <- R2_W - R2_pop_broad
    resultsB_pop[29, i+2009] <- R2_W - R2_sub
  } 
} else { 
  R2_W <- NA
  resultsB_pop[29, i+9] <- 1      
  resultsB_pop[29, i+1009] <- NA
  resultsB_pop[29, i+2009] <- NA
}

} # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

for (i in 1:nrow(resultsB_pop)) {
  resultsB_pop[i, 3] <- length(which(resultsB_pop[i, c(10:(nperm + 9))] 
                                     <= 0.05)) / nperm
  resultsB_pop[i, 4] <- median(na.omit(as.numeric(resultsB_pop[i, c(1010:(nperm + 
                                                                            1009))])))
  resultsB_pop[i, 5] <- sd(na.omit(as.numeric(resultsB_pop[i, c(1010:(nperm + 
                                                                        1009))])))
  resultsB_pop[i, 6] <- median(na.omit(as.numeric(resultsB_pop[i, c(2010:(nperm + 
                                                                            2009))])))
  resultsB_pop[i, 7] <- sd(na.omit(as.numeric(resultsB_pop[i, c(2010:(nperm + 
                                                                        2009))])))
  resultsB_pop[i, 8] <- median(na.omit(as.numeric(resultsB_pop[i, c(3010:(nperm + 
                                                                            3009))])))
  resultsB_pop[i, 9] <- sd(na.omit(as.numeric(resultsB_pop[i, c(3010:(nperm + 
                                                                        3009))])))
}

# Output of the results:
# **********************
res_file_name <- paste("Results_pop", "Broad", design, paste("style", style, ".txt",
                                                             sep = ""), sep = "_")
write.table(resultsB_pop, file = res_file_name, sep = "\t")

   ########################
   ########################
   ### II. Medium scale ###
   ########################
   ########################

for (i in 1:nperm) {
  
  # Sampling scheme:
  # ****************
  C <- C.list[[i]]
  
  # We keep the lines of MEM that correspond to the sampled cells ('tri'):
  # **********************************************************************
  MEMsub <- MEMsub.list[[i]]
  
  # We sample the response variable within the sampled cells ('y_sub'):
  # *******************************************************************
  y_sub <- y_med[as.numeric(row.names(C)), ]
  
  # Real p-value and R2_sub:
  # ************************
  x <- MEMsub[, c(211, 212, 215)]
  lm <- lm(y_sub ~ ., data = x)
  R2_sub <- summary(lm)$adj.r.squared
  resultsM_pop[31, c(2009+i, 3009+i)] <- R2_sub
  resultsM_pop[32, 9+i] <- lmp(lm)
  
  # Construction of the different W matrices:
  # #########################################
  
  # Connectivity matrices (B):
  # **************************
  
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
  f3 <- function (D, y)       {1/(D^y)}           # Concave-up function
  f4 <- function (D, t)       {1-(D/(4*t))^2}     # PCNM criterion
  
  max.del <- max(unlist(nbdists(Y.del, as.matrix(C)))) 
  max.gab <- max(unlist(nbdists(Y.gab, as.matrix(C))))
  max.rel <- max(unlist(nbdists(Y.rel, as.matrix(C)))) 
  max.mst <- max(unlist(nbdists(Y.mst, as.matrix(C))))
  
  nbdist <- lapply(Y.list10dnn, coords = as.matrix(C), nbdists)
  unlist <- lapply(nbdist, unlist)
  max.DB.list <- lapply(unlist, max)
  
  # Generation of MEM variables:
  # ****************************
  # del:
  # ****
  Y.del.MEM <- test.W.R2(nb = Y.del, xy = C, MEM.autocor = MEM_model)
  Y.del.MEM.f1 <- test.W.R2(nb = Y.del, xy = C, f = f1, dmax = max.del, 
                            MEM.autocor = MEM_model)
  Y.del.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, MEM.autocor = MEM_model, 
                            f = f2, dmax = max.del, y = 2:10)
  Y.del.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.del, xy = C, MEM.autocor = MEM_model, 
                            f = f3, y = 2:10)
  # gab:
  # ****
  Y.gab.MEM <- test.W.R2(nb = Y.gab, xy = C, MEM.autocor = MEM_model)
  Y.gab.MEM.f1 <- test.W.R2(nb = Y.gab, xy = C, f = f1, dmax = max.gab, 
                            MEM.autocor = MEM_model)
  Y.gab.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, MEM.autocor = MEM_model, 
                            f = f2, dmax = max.gab, y = 2:10)
  Y.gab.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.gab, xy = C, MEM.autocor = MEM_model, 
                            f = f3, y = 2:10)
  # rel:
  # ****
  Y.rel.MEM <- test.W.R2(nb = Y.rel, xy = C, MEM.autocor = MEM_model)
  Y.rel.MEM.f1 <- test.W.R2(nb = Y.rel, xy = C, f = f1, dmax = max.rel, 
                            MEM.autocor = MEM_model)
  Y.rel.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, MEM.autocor = MEM_model, 
                            f = f2, dmax = max.rel, y = 2:10)
  Y.rel.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.rel, xy = C, MEM.autocor = MEM_model, 
                            f = f3, y = 2:10)
  # mst:
  # ****
  Y.mst.MEM <- test.W.R2(nb = Y.mst, xy = C, MEM.autocor = MEM_model)
  Y.mst.MEM.f1 <- test.W.R2(nb = Y.mst, xy = C, f = f1, dmax = max.mst, 
                            MEM.autocor = MEM_model)
  Y.mst.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, MEM.autocor = MEM_model, 
                            f = f2, dmax = max.mst, y = 2:10)
  Y.mst.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.mst, xy = C, MEM.autocor = MEM_model, 
                            f = f3, y = 2:10)
  # DB:
  # ***
  Y.DB1.MEM <- test.W.R2(nb = Y.listDB[[1]], xy = C, MEM.autocor = MEM_model)
  Y.DB1.MEM.f1 <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f1, dmax = max.DB.list[[1]], 
                            MEM.autocor = MEM_model)
  Y.DB1.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[1]], 
                            y = 2:10)
  Y.DB1.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[1]], xy = C, 
                            MEM.autocor = MEM_model, f = f3, y = 2:10)
  
  Y.DB2.MEM <- test.W.R2(nb = Y.listDB[[2]], xy = C, MEM.autocor = MEM_model)
  Y.DB2.MEM.f1 <- test.W.R2(nb = Y.listDB[[2]], xy = C, f = f1, dmax = max.DB.list[[2]], 
                            MEM.autocor = MEM_model)
  Y.DB2.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[2]], 
                            y = 2:10)
  Y.DB2.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[2]], xy = C, 
                            MEM.autocor = MEM_model, f = f3, y = 2:10)
  
  Y.DB3.MEM <- test.W.R2(nb = Y.listDB[[3]], xy = C, MEM.autocor = MEM_model)
  Y.DB3.MEM.f1 <- test.W.R2(nb = Y.listDB[[3]], xy = C, f = f1, dmax = max.DB.list[[3]], 
                            MEM.autocor = MEM_model)
  Y.DB3.MEM.f2 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, 
                            MEM.autocor = MEM_model, f = f2, dmax = max.DB.list[[3]], 
                            y = 2:10)
  Y.DB3.MEM.f3 <- test.W.R2(Y = y_sub, nb = Y.listDB[[3]], xy = C, 
                            MEM.autocor = MEM_model, f = f3, y = 2:10)
  
  # DBMEM with PCNM criteria:
  # *************************
  Y.DB.PCNM <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f4, t = lowlim, 
                         MEM.autocor = MEM_model)
  
  # Significance test and MEM var. selection (fwd.sel with double stopping criterion):
  # **********************************************************************************
  # **********************************************************************************
  # del
  R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.del.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[1, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[1, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[1, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[1, i+9] <- 1     # p-val set at a value of 1 (for power computation)   
    resultsM_pop[1, i+1009] <- NA
    resultsM_pop[1, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.del.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[2, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[2, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[2, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[2, i+9] <- 1         
    resultsM_pop[2, i+1009] <- NA
    resultsM_pop[2, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.del.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[3, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[3, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[3, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[3, i+9] <- 1         
    resultsM_pop[3, i+1009] <- NA
    resultsM_pop[3, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.del.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.del.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.del.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.del.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[4, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[4, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[4, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[4, i+9] <- 1         
    resultsM_pop[4, i+1009] <- NA
    resultsM_pop[4, i+2009] <- NA
  }
  
  # gab
  R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.gab.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[5, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[5, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[5, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[5, i+9] <- 1         
    resultsM_pop[5, i+1009] <- NA
    resultsM_pop[5, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.gab.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[6, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[6, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[6, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[6, i+9] <- 1         
    resultsM_pop[6, i+1009] <- NA
    resultsM_pop[6, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.gab.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[7, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[7, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[7, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[7, i+9] <- 1         
    resultsM_pop[7, i+1009] <- NA
    resultsM_pop[7, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.gab.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.gab.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.gab.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.gab.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[8, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[8, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[8, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[8, i+9] <- 1         
    resultsM_pop[8, i+1009] <- NA
    resultsM_pop[8, i+2009] <- NA
  }
  
  # rel
  R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.rel.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[9, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                      permutations = 9999))$Pr[1]
      resultsM_pop[9, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[9, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[9, i+9] <- 1         
    resultsM_pop[9, i+1009] <- NA
    resultsM_pop[9, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.rel.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[10, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[10, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[10, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[10, i+9] <- 1         
    resultsM_pop[10, i+1009] <- NA
    resultsM_pop[10, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.rel.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[11, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[11, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[11, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[11, i+9] <- 1         
    resultsM_pop[11, i+1009] <- NA
    resultsM_pop[11, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.rel.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.rel.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.rel.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.rel.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[12, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[12, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[12, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[12, i+9] <- 1         
    resultsM_pop[12, i+1009] <- NA
    resultsM_pop[12, i+2009] <- NA
  }
  
  # mst
  R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.mst.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[13, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[13, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[13, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[13, i+9] <- 1         
    resultsM_pop[13, i+1009] <- NA
    resultsM_pop[13, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.mst.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[14, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[14, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[14, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[14, i+9] <- 1         
    resultsM_pop[14, i+1009] <- NA
    resultsM_pop[14, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.mst.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[15, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[15, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[15, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[15, i+9] <- 1         
    resultsM_pop[15, i+1009] <- NA
    resultsM_pop[15, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.mst.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.mst.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.mst.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.mst.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[16, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[16, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[16, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[16, i+9] <- 1         
    resultsM_pop[16, i+1009] <- NA
    resultsM_pop[16, i+2009] <- NA
  }
  
  # DB1
  R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB1.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[17, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[17, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[17, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[17, i+9] <- 1         
    resultsM_pop[17, i+1009] <- NA
    resultsM_pop[17, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB1.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[18, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[18, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[18, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[18, i+9] <- 1         
    resultsM_pop[18, i+1009] <- NA
    resultsM_pop[18, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB1.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[19, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[19, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[19, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[19, i+9] <- 1         
    resultsM_pop[19, i+1009] <- NA
    resultsM_pop[19, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB1.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB1.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB1.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB1.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[20, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[20, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[20, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[20, i+9] <- 1         
    resultsM_pop[20, i+1009] <- NA
    resultsM_pop[20, i+2009] <- NA
  }
  
  # DB2
  R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB2.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[21, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[21, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[21, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[21, i+9] <- 1         
    resultsM_pop[21, i+1009] <- NA
    resultsM_pop[21, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB2.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[22, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[22, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[22, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[22, i+9] <- 1         
    resultsM_pop[22, i+1009] <- NA
    resultsM_pop[22, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB2.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[23, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[23, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[23, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[23, i+9] <- 1         
    resultsM_pop[23, i+1009] <- NA
    resultsM_pop[23, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB2.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB2.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB2.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB2.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[24, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[24, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[24, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[24, i+9] <- 1         
    resultsM_pop[24, i+1009] <- NA
    resultsM_pop[24, i+2009] <- NA
  }
  
  # DB3
  R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB3.MEM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[25, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[25, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[25, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[25, i+9] <- 1         
    resultsM_pop[25, i+1009] <- NA
    resultsM_pop[25, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f1$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB3.MEM.f1$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f1$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f1$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[26, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[26, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[26, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[26, i+9] <- 1      
    resultsM_pop[26, i+1009] <- NA
    resultsM_pop[26, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f2$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB3.MEM.f2$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f2$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f2$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[27, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[27, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[27, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[27, i+9] <- 1      
    resultsM_pop[27, i+1009] <- NA
    resultsM_pop[27, i+2009] <- NA
  }
  R2adj <- RsquareAdj(rda(y_sub, Y.DB3.MEM.f3$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB3.MEM.f3$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB3.MEM.f3$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB3.MEM.f3$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[28, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[28, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[28, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[28, i+9] <- 1      
    resultsM_pop[28, i+1009] <- NA
    resultsM_pop[28, i+2009] <- NA
  }
  
  # DB_PCNM
  R2adj <- RsquareAdj(rda(y_sub, Y.DB.PCNM$best$MEM))$adj.r.squared
  if (anova.cca(rda(y_sub, Y.DB.PCNM$best$MEM), permutations = 9999)$Pr[1] <= 0.05) {
    class <- class(try(fsel <- forward.sel(y_sub, Y.DB.PCNM$best$MEM, 
                                           adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- Y.DB.PCNM$best$MEM[, c(sign)]
      R2_W <- RsquareAdj(rda(y_sub, MEM.FwdSel))$adj.r.squared
      resultsM_pop[29, i+9] <- as.data.frame(anova.cca(rda(y_sub, MEM.FwdSel),
                                                       permutations = 9999))$Pr[1]
      resultsM_pop[29, i+1009] <- R2_W - R2_pop_broad
      resultsM_pop[29, i+2009] <- R2_W - R2_sub
    } 
  } else { 
    R2_W <- NA
    resultsM_pop[29, i+9] <- 1      
    resultsM_pop[29, i+1009] <- NA
    resultsM_pop[29, i+2009] <- NA
  }
  
} # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

for (i in 1:nrow(resultsM_pop)) {
  resultsM_pop[i, 3] <- length(which(resultsM_pop[i, c(10:(nperm + 9))] 
                                     <= 0.05)) / nperm
  resultsM_pop[i, 4] <- median(na.omit(as.numeric(resultsM_pop[i, c(1010:(nperm + 
                                                                            1009))])))
  resultsM_pop[i, 5] <- sd(na.omit(as.numeric(resultsM_pop[i, c(1010:(nperm + 
                                                                        1009))])))
  resultsM_pop[i, 6] <- median(na.omit(as.numeric(resultsM_pop[i, c(2010:(nperm + 
                                                                            2009))])))
  resultsM_pop[i, 7] <- sd(na.omit(as.numeric(resultsM_pop[i, c(2010:(nperm + 
                                                                        2009))])))
  resultsM_pop[i, 8] <- median(na.omit(as.numeric(resultsM_pop[i, c(3010:(nperm + 
                                                                            3009))])))
  resultsM_pop[i, 9] <- sd(na.omit(as.numeric(resultsM_pop[i, c(3010:(nperm + 
                                                                        3009))])))
}

# Output of the results:
# **********************
res_file_name <- paste("Results_pop", "Medium", design, paste("style", style, ".txt",
                                                              sep = ""), sep = "_")
write.table(resultsM_pop, file = res_file_name, sep = "\t")


### CONTINUE NOW FOR RESULTS_SUB
