# ********************************************** #
###### PART I. Comparison of the W matrices ######
# ********************************************** #

##### Script du test de puissance, précision et robustesse de différentes #####
### matrices W vis-à-vis du type de design d'échantillonnage, de l'échelle ####
### spatiale de la structuration, d'un point de vue de la population réelle ###
###################### et de l'échantillon réel. ##############################
# *************************************************************************** #

rm(list=ls()[-match(c("xy", "nb", "nb2", "MEM", "MEM_backup"), ls())])

# Usefull packages and functions:
# *******************************

# library(ade4)
library(vegan)
library(adespatial)
library(spdep)

source("lmp.R")
source("test.W.R2.R")
# Function to compute the global test and FwdSel on the MEM variables (for power and accuracy):
MEMfwd.test <- function (y, mem) {
  pval <- as.data.frame(anova.cca(rda(y, mem), permutations = 9999))$Pr[1]
  if (pval <= 0.05) {
    R2adj <- RsquareAdj(rda(y, mem))$adj.r.squared
    class <- class(try(fsel <- forward.sel(y, mem, adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- mem[, c(sign)]
      R2_W <- RsquareAdj(rda(y, MEM.FwdSel))$adj.r.squared
      delta_pop <- R2_W - R2_pop_broad
      delta_sub <- R2_W - R2_sub
    } else { 
      pval <- 1      
      delta_pop <- NA
      delta_sub <- NA
    }
  } else { 
    pval <- 1      
    delta_pop <- NA
    delta_sub <- NA
  }
  list(pval = pval, delta_pop = delta_pop, delta_sub = delta_sub)
}

# Definition of the simulation parameters:
# ****************************************

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"   # Either "positive" or "negative"

# Sampling design:
design <- "random"    # Either "clustered" or "random"

nperm <- 1000

# Structuring Intensity (low or high):
a <- 0.35   # 0.35 or 0.55
if (a < 0.5) intensity = "Weak" else intensity = "Strong"

style <- "B"             # Either "B" or "W"

# Construction of the results matrices :
# **************************************
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

# The MEM are built for a full grid (50 x 25 cells):
# **************************************************

xy <- expand.grid(x = seq(1, 90, 1), y = seq(1, 90, 1))

nb <- cell2nb(nrow = 90, ncol = 90, "queen")
nb2 <- nb2listw(nb, style = style)
MEM <- scores.listw(nb2, MEM.autocor = MEM_model)

# xy_attr <- attr(nb, "region.id")
# XY <- matrix(as.integer(unlist(strsplit(xy_attr, ":"))), ncol = 2, byrow = TRUE)
# plot(nb, XY)

# To know from where and in which direction the cells are considered when building MEM
# ************************************************************************************
# s.label(xy, neig = nb2neig(nb), clab = 0.5)


# ***************************************************************** #
### I. The response remains unchanged and the sampling design varies:
### ==> Robustness of the W matrices and of R2_sub to sampling scheme
### varations. ######################################################
#####################################################################
#####################################################################
# ***************************************************************** #

# Creation of the response variable:
# **********************************
   # Creation of 'y_spa' and 'y_noise' and standardisation:
   # ******************************************************

set.seed(1)

y_spa_broad <- MEM[, 6] + MEM[, 7] + MEM[, 8]
y_spa_med <- MEM[, 40] + MEM[, 41] + MEM[, 42]
y_noise <- rnorm(n = nrow(MEM), mean = 0, sd = 1)

y_spa_broad_st <- scale(y_spa_broad)
y_spa_med_st <- scale(y_spa_med)
y_noise_st <- scale(y_noise)

#par(mfrow = c(1, 3), mar = c(2, 2, 1, 1))
#for (i in 1:3) image(matrix(MEM[, i+39], ncol = 90, byrow = F))
#par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))
#image(matrix(y_spa_med_st, ncol = 90, byrow = F))

   # Creation of the response variable 'y' at the whole population level (pop):
   # **************************************************************************

y_broad <- (a * y_spa_broad_st) + ((1-a) * y_noise_st)
y_med <- (a * y_spa_med_st) + ((1-a) * y_noise_st)

#par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))
#image(matrix(y_broad, ncol = 90, byrow = F))

R2_pop_broad <- cor(y_broad, y_spa_broad_st)^2
R2_pop_med <- cor(y_med, y_spa_med_st)^2

resultsB_pop[30, c(1010:2009, 3010:4009)] <- R2_pop_broad
resultsM_pop[30, c(1010:2009, 3010:4009)] <- R2_pop_med

# Begining of the simulation process:
# ***********************************

# Will save the nperm sampling designs ('C') and corresponding lists of 29 W matrices
# to reuse them at the medium scale (to spare time):
C.list <- vector("list", nperm)  
listW.list <- vector("list", nperm)

   ######################
   ######################
   ### I. Broad scale ###
   ######################
   ######################

for (i in 1:nperm) {

  # Sampling scheme:
  # ****************
  set.seed(i)
  
  if (design == "clustered") {
    zones <- matrix(c(1:9, rep(c(0, 30, 60), times = 3),rep(c(0, 30, 60), each = 3)), ncol = 3)
    sampled.zones <- sample(c(1:9), 3)
    grid <- expand.grid(c(6:25), c(6:25))
    for (w in 1:3) {
      points <- grid[sample(c(1:nrow(grid)), 40), ]
      points[, 1] <- points[, 1] + zones[sampled.zones[w], 2]
      points[, 2] <- points[, 2] + zones[sampled.zones[w], 3]
      if (w == 1) C <- points else C <- rbind(C, points)
    }
    # Attribute each sampled point to one of the 'xy' grid cells:
    grid.size <- 1
    tri <- c()
    for (k in 1:nrow(C)) {
      x <- floor((C[k, 1]) / (grid.size))
      y <- floor((C[k, 2]) / (grid.size))
      N <- y * 90 + x + 1
      tri <- c(tri, N)
    }
    rownames(C) <- tri
  } else C <- xy[sample(c(1:nrow(xy)), 120), ]   # design = "random"

xy.d1 <- dist(C)
C.list[[i]] <- C

# We keep the lines of y_spa_broad_st that correspond to the sampled cells ('tri'):
# *********************************************************************************
MEMsub <- y_spa_broad_st[as.numeric(rownames(C))]

# We sample the response variable within the sampled cells ('y_sub'):
# *******************************************************************
y_sub <- y_broad[as.numeric(rownames(C))]

# Real p-value and R2_sub:
# ************************
lm <- lm(y_sub ~ MEMsub)
resultsB_pop[32, 9+i] <- lmp(lm)
R2_sub <- cor(y_sub, MEMsub)^2                                          
resultsB_pop[31, c(2009+i, 3009+i)] <- R2_sub
resultsB_pop[c(1:29), 3009+i] <- as.numeric(R2_sub - R2_pop_broad)

# Visualisation:
# **************
# par(mfrow = c(1, 1))
# s.value(C, MEMsub)

# Construction of the different W matrices:
# #########################################

   # Connectivity matrices (B):
   # **************************

Y.del <- tri2nb(jitter(as.matrix(C), factor = 0.001)) 
Y.gab <- graph2nb(gabrielneigh(as.matrix(C), nnmult = 5), sym = TRUE)
Y.rel <- graph2nb(relativeneigh(as.matrix(C), nnmult = 5), sym = TRUE)
Y.mst <- mst.nb(xy.d1)

   # Distance-based B (radius around points):
lowlim <- give.thresh(xy.d1)
uplim <- max(xy.d1)
thresh <- seq(lowlim, uplim, le = 10)   # 3 tested distances: thresh[1, 5, 9]
Y.listDB <- lapply(thresh[c(1, 5, 9)], dnearneigh, x = as.matrix(C), d1 = 0)

# Weighting functions and fixed parametres:
# *****************************************
f1 <- function (D, dmax)    { 1 - (D/dmax) }        # Linear function
f2 <- function (D, dmax, y) { 1 - (D/dmax)^y }      # Concave-down function
f3 <- function (D, y)       { 1 / D^y }             # Concave-up function
f4 <- function (D, t)       { 1 - (D/(4*t))^2 }     # PCNM criterion

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
Y.del.MEM.f1 <- test.W.R2(nb = Y.del, xy = C, f = f1, dmax = max.del, 
                          style = style, MEM.autocor = MEM_model)
Y.del.MEM.f2 <- test.W.R2(nb = Y.del, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.del, y = 5)
Y.del.MEM.f3 <- test.W.R2(nb = Y.del, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f3, y = 0.5)
# gab:
# ****
Y.gab.MEM <- test.W.R2(nb = Y.gab, xy = C, style = style, MEM.autocor = MEM_model)
Y.gab.MEM.f1 <- test.W.R2(nb = Y.gab, xy = C, f = f1, dmax = max.gab, 
                          style = style, MEM.autocor = MEM_model)
Y.gab.MEM.f2 <- test.W.R2(nb = Y.gab, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.gab, y = 5)
Y.gab.MEM.f3 <- test.W.R2(nb = Y.gab, xy = C, style = style,  
                          MEM.autocor = MEM_model, f = f3, y = 0.5)
# rel:
# ****
Y.rel.MEM <- test.W.R2(nb = Y.rel, xy = C, style = style, MEM.autocor = MEM_model)
Y.rel.MEM.f1 <- test.W.R2(nb = Y.rel, xy = C, f = f1, dmax = max.rel, 
                          style = style, MEM.autocor = MEM_model)
Y.rel.MEM.f2 <- test.W.R2(nb = Y.rel, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.rel, y = 5)
Y.rel.MEM.f3 <- test.W.R2(nb = Y.rel, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f3, y = 0.5)
# mst:
# ****
Y.mst.MEM <- test.W.R2(nb = Y.mst, xy = C, style = style, MEM.autocor = MEM_model)
Y.mst.MEM.f1 <- test.W.R2(nb = Y.mst, xy = C, f = f1, dmax = max.mst, 
                          style = style, MEM.autocor = MEM_model)
Y.mst.MEM.f2 <- test.W.R2(nb = Y.mst, xy = C, style = style, 
                          MEM.autocor = MEM_model, f = f2, dmax = max.mst, y = 5)
Y.mst.MEM.f3 <- test.W.R2(nb = Y.mst, xy = C, style = style,  
                          MEM.autocor = MEM_model, f = f3, y = 0.5)
# DB:
# ***
Y.DB1.MEM <- test.W.R2(nb = Y.listDB[[1]], xy = C, style = style, MEM.autocor = MEM_model)
Y.DB1.MEM.f1 <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f1, dmax = max.DB.list[[1]], 
                          style = style, MEM.autocor = MEM_model)
Y.DB1.MEM.f2 <- test.W.R2(nb = Y.listDB[[1]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f2, 
                          dmax = max.DB.list[[1]], y = 5)
Y.DB1.MEM.f3 <- test.W.R2(nb = Y.listDB[[1]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f3, y = 0.5)

Y.DB2.MEM <- test.W.R2(nb = Y.listDB[[2]], xy = C, style = style, 
                       MEM.autocor = MEM_model)
Y.DB2.MEM.f1 <- test.W.R2(nb = Y.listDB[[2]], xy = C, f = f1, dmax = max.DB.list[[2]], 
                          style = style, MEM.autocor = MEM_model)
Y.DB2.MEM.f2 <- test.W.R2(nb = Y.listDB[[2]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f2, 
                          dmax = max.DB.list[[2]], y = 5)
Y.DB2.MEM.f3 <- test.W.R2(nb = Y.listDB[[2]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f3, y = 0.5)

Y.DB3.MEM <- test.W.R2(nb = Y.listDB[[3]], xy = C, style = style, 
                       MEM.autocor = MEM_model)
Y.DB3.MEM.f1 <- test.W.R2(nb = Y.listDB[[3]], xy = C, f = f1, dmax = max.DB.list[[3]], 
                          style = style, MEM.autocor = MEM_model)
Y.DB3.MEM.f2 <- test.W.R2(nb = Y.listDB[[3]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f2, 
                          dmax = max.DB.list[[3]], y = 5)
Y.DB3.MEM.f3 <- test.W.R2(nb = Y.listDB[[3]], xy = C, 
                          style = style, MEM.autocor = MEM_model, f = f3, y = 0.5)

# DBMEM with PCNM criteria:
# *************************
Y.DB.PCNM <- test.W.R2(nb = Y.listDB[[1]], xy = C, f = f4, t = lowlim, style = style, 
                       MEM.autocor = MEM_model)

listW <- list(Y.del.MEM$MEM, Y.del.MEM.f1$MEM, Y.del.MEM.f2$MEM, Y.del.MEM.f3$MEM, 
              Y.gab.MEM$MEM, Y.gab.MEM.f1$MEM, Y.gab.MEM.f2$MEM, Y.gab.MEM.f3$MEM, 
              Y.rel.MEM$MEM, Y.rel.MEM.f1$MEM, Y.rel.MEM.f2$MEM, Y.rel.MEM.f3$MEM, 
              Y.mst.MEM$MEM, Y.mst.MEM.f1$MEM, Y.mst.MEM.f2$MEM, Y.mst.MEM.f3$MEM, 
              Y.DB1.MEM$MEM, Y.DB1.MEM.f1$MEM, Y.DB1.MEM.f2$MEM, Y.DB1.MEM.f3$MEM, 
              Y.DB2.MEM$MEM, Y.DB2.MEM.f1$MEM, Y.DB2.MEM.f2$MEM, Y.DB2.MEM.f3$MEM,
              Y.DB3.MEM$MEM, Y.DB3.MEM.f1$MEM, Y.DB3.MEM.f2$MEM, Y.DB3.MEM.f3$MEM, 
              Y.DB.PCNM$MEM)
listW.list[[i]] <- listW

# Significance test and MEM var. selection (fwd.sel with double stopping criterion):
# **********************************************************************************
# **********************************************************************************
for (q in 1:length(listW)) {
  test <- MEMfwd.test(y_sub, listW[[q]])
  resultsB_pop[q, i+9] <- test$pval
  resultsB_pop[q, i+1009] <- test$delta_pop
  resultsB_pop[q, i+2009] <- test$delta_sub
}

} # End of the simulation ('for') loop

# Summary of the results:
# ***********************

for (j in 1:(nrow(resultsB_pop)-1)) {
  resultsB_pop[j, 3] <- length(which(resultsB_pop[j, c(10:(nperm + 9))] 
                                     <= 0.05)) / nperm
  resultsB_pop[j, 4] <- median(na.omit(as.numeric(resultsB_pop[j, c(1010:(nperm + 
                                                                            1009))])))
  resultsB_pop[j, 5] <- sd(na.omit(as.numeric(resultsB_pop[j, c(1010:(nperm + 
                                                                        1009))])))
  resultsB_pop[j, 6] <- median(na.omit(as.numeric(resultsB_pop[j, c(2010:(nperm + 
                                                                            2009))])))
  resultsB_pop[j, 7] <- sd(na.omit(as.numeric(resultsB_pop[j, c(2010:(nperm + 
                                                                        2009))])))
  resultsB_pop[j, 8] <- median(na.omit(as.numeric(resultsB_pop[j, c(3010:(nperm + 
                                                                            3009))])))
  resultsB_pop[j, 9] <- sd(na.omit(as.numeric(resultsB_pop[j, c(3010:(nperm + 
                                                                        3009))])))
}

# Output of the results:
# **********************
res_file_name <- paste("Results_pop", intensity, "Broad", design, 
                       paste("style", style, ".txt", sep = ""), sep = "_")
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
  listW <- listW.list[[i]]
  
  # We keep the lines of y_spa_med_st that correspond to the sampled cells ('tri'):
  # *******************************************************************************
  MEMsub <- y_spa_med_st[as.numeric(rownames(C))]
  
  # We sample the response variable within the sampled cells ('y_sub'):
  # *******************************************************************
  y_sub <- y_med[as.numeric(rownames(C))]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsM_pop[32, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, MEMsub)^2                                          
  resultsM_pop[31, c(2009+i, 3009+i)] <- R2_sub
  resultsM_pop[c(1:29), 3009+i] <- as.numeric(R2_sub - R2_pop_med)
  
  # Significance test and MEM var. selection (fwd.sel with double stopping criterion):
  # **********************************************************************************
  # **********************************************************************************
  for (q in 1:length(listW)) {
    test <- MEMfwd.test(y_sub, listW[[q]])
    resultsM_pop[q, i+9] <- test$pval
    resultsM_pop[q, i+1009] <- test$delta_pop
    resultsM_pop[q, i+2009] <- test$delta_sub
  }
} # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

for (i in 1:(nrow(resultsM_pop)-1)) {
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
res_file_name <- paste("Results_pop", intensity, "Medium", design, 
                       paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsM_pop, file = res_file_name, sep = "\t")





# ****************************************************************** #
### II. The sampling design remains unchanged and the response varies:
# ****************************************************************** #

C <- C.list[[5]]
listW <- listW.list[[5]]

######################
######################
### I. Broad scale ###
######################
######################

# We keep the rows of y_spa_broad_st that correspond to the sampled cells (rownames(C)):
# **************************************************************************************
MEMsub <- y_spa_broad_st[as.numeric(rownames(C))]

for (i in 1:nperm) { 
  
  set.seed(i)
  y_noise_st <- scale(rnorm(n = nrow(MEM), mean = 0, sd = 1))
  
  # Creation of the response variable 'y' at the whole population level (pop):
  # **************************************************************************
  
  y_broad <- (a * y_spa_broad_st) + ((1-a) * y_noise_st)
  R2_pop_broad <- cor(y_broad, y_spa_broad_st)^2
  resultsB_sub[30, c(1009+i, 3009+i)] <- R2_pop_broad
  
  y_sub <- y_broad[as.numeric(rownames(C))]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsB_sub[32, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, MEMsub)^2                                          
  resultsB_sub[31, c(2009+i, 3009+i)] <- R2_sub
  resultsB_sub[c(1:29), 3009+i] <- as.numeric(R2_sub - R2_pop_broad)
  
  # Significance test and MEM var. selection (fwd.sel with double stopping criterion):
  # **********************************************************************************
  # **********************************************************************************
  for (q in 1:length(listW)) {
    test <- MEMfwd.test(y_sub, listW[[q]])
    resultsB_sub[q, i+9] <- test$pval
    resultsB_sub[q, i+1009] <- test$delta_pop
    resultsB_sub[q, i+2009] <- test$delta_sub
  }
} # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

for (i in 1:(nrow(resultsB_sub)-1)) {
  resultsB_sub[i, 3] <- length(which(resultsB_sub[i, c(10:(nperm + 9))] 
                                     <= 0.05)) / nperm
  resultsB_sub[i, 4] <- median(na.omit(as.numeric(resultsB_sub[i, c(1010:(nperm + 
                                                                            1009))])))
  resultsB_sub[i, 5] <- sd(na.omit(as.numeric(resultsB_sub[i, c(1010:(nperm + 
                                                                        1009))])))
  resultsB_sub[i, 6] <- median(na.omit(as.numeric(resultsB_sub[i, c(2010:(nperm + 
                                                                            2009))])))
  resultsB_sub[i, 7] <- sd(na.omit(as.numeric(resultsB_sub[i, c(2010:(nperm + 
                                                                        2009))])))
  resultsB_sub[i, 8] <- median(na.omit(as.numeric(resultsB_sub[i, c(3010:(nperm + 
                                                                            3009))])))
  resultsB_sub[i, 9] <- sd(na.omit(as.numeric(resultsB_sub[i, c(3010:(nperm + 
                                                                        3009))])))
}

# Output of the results:
# **********************
res_file_name <- paste("Results_sub", intensity, "Broad", design, 
                       paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsB_sub, file = res_file_name, sep = "\t")


########################
########################
### II. Medium scale ###
########################
########################

# We keep the lines of MEM that correspond to the sampled cells ('tri'):
# **********************************************************************
MEMsub <- y_spa_med_st[as.numeric(rownames(C))]

for (i in 1:nperm) {
  
  set.seed(i)
  y_noise_st <- scale(rnorm(n = nrow(MEM), mean = 0, sd = 1))
  
  # Creation of the response variable 'y' at the whole population level (pop):
  # **************************************************************************
  
  y_med <- (a * y_spa_med_st) + ((1-a) * y_noise_st)
  R2_pop_med <- cor(y_med, y_spa_med_st)^2
  resultsM_sub[30, c(1009+i, 3009+i)] <- R2_pop_med
  
  y_sub <- y_med[as.numeric(rownames(C))]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsM_sub[32, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, MEMsub)^2                                          
  resultsM_sub[31, c(2009+i, 3009+i)] <- R2_sub
  resultsM_sub[c(1:29), 3009+i] <- as.numeric(R2_sub - R2_pop_med)
  
  # Significance test and MEM var. selection (fwd.sel with double stopping criterion):
  # **********************************************************************************
  # **********************************************************************************
  for (q in 1:length(listW)) {
    test <- MEMfwd.test(y_sub, listW[[q]])
    resultsM_sub[q, i+9] <- test$pval
    resultsM_sub[q, i+1009] <- test$delta_pop
    resultsM_sub[q, i+2009] <- test$delta_sub
  }
} # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

for (i in 1:(nrow(resultsM_sub)-1)) {
  resultsM_sub[i, 3] <- length(which(resultsM_sub[i, c(10:(nperm + 9))] 
                                     <= 0.05)) / nperm
  resultsM_sub[i, 4] <- median(na.omit(as.numeric(resultsM_sub[i, c(1010:(nperm + 
                                                                            1009))])))
  resultsM_sub[i, 5] <- sd(na.omit(as.numeric(resultsM_sub[i, c(1010:(nperm + 
                                                                        1009))])))
  resultsM_sub[i, 6] <- median(na.omit(as.numeric(resultsM_sub[i, c(2010:(nperm + 
                                                                            2009))])))
  resultsM_sub[i, 7] <- sd(na.omit(as.numeric(resultsM_sub[i, c(2010:(nperm + 
                                                                        2009))])))
  resultsM_sub[i, 8] <- median(na.omit(as.numeric(resultsM_sub[i, c(3010:(nperm + 
                                                                            3009))])))
  resultsM_sub[i, 9] <- sd(na.omit(as.numeric(resultsM_sub[i, c(3010:(nperm + 
                                                                        3009))])))
}

# Output of the results:
# **********************
res_file_name <- paste("Results_sub", intensity, "Medium", design, 
                       paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsM_sub, file = res_file_name, sep = "\t")
