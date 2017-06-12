################################################################
### ****************** Test of MEM.modsel ****************** ###

# rm(list=ls())

# Usefull packages and functions:
# *******************************

library(vegan)
library(adespatial)
library(spdep)

source("lmp.R")
source("test.W.R2.R") 

# Definition of the simulation parameters:
# ****************************************

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"   # Either "positive" or "negative"

# Sampling design:
design <- "clustered"    # Either "clustered" or "random"

nperm <- 1000

# Structuring Intensity (low or high):
a <- 0.55   # 0.35 or 0.55
if (a < 0.5) intensity = "Weak" else intensity = "Strong"

style <- "B"                  # Either "B" or "W"
correction <- "Corrected"     # Either "Corrected" or "Uncorrected"

if (correction == "Corrected") {
  source("MEM.modsel - correction.R")
} else source("MEM.modsel - no_correction.R")

# Construction of a results matrix for each scale and for both W matrix selection:
# ********************************************************************************

resultsB_pop <- as.data.frame(matrix(nrow = 5, ncol = 4009))
colnames(resultsB_pop) <- c("Matrix B", "Matrix A", "Power", 
                            "Median deltaR2_pop", "sd deltaR2_pop", 
                            "Median deltaR2_sub", "sd deltaR2_sub",
                            "Median deltaR2_subpop", "sd deltaR2_subpop",
                            paste("pval", c(1:1000), sep = ""), 
                            paste("delR2_pop", c(1:1000), sep = ""),
                            paste("delR2_sub", c(1:1000), sep = ""),
                            paste("delR2_subpop", c(1:1000), sep = ""))
resultsB_pop[, 1] <- c("sim", "R2_pop", "R2_sub", "pvalReal", "column_count")
resultsB_pop[5, ] <- c(1:ncol(resultsB_pop))

resultsM_pop <- as.data.frame(matrix(nrow = 5, ncol = 4009))
colnames(resultsM_pop) <- c("Matrix B", "Matrix A", "Power", 
                            "Median deltaR2_pop", "sd deltaR2_pop", 
                            "Median deltaR2_sub", "sd deltaR2_sub",
                            "Median deltaR2_subpop", "sd deltaR2_subpop",
                            paste("pval", c(1:1000), sep = ""), 
                            paste("delR2_pop", c(1:1000), sep = ""),
                            paste("delR2_sub", c(1:1000), sep = ""),
                            paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_pop[, 1] <- c("sim", "R2_pop", "R2_sub", "pvalReal", "column_count")
resultsM_pop[5, ] <- c(1:ncol(resultsM_pop))

resultsB_sub <- as.data.frame(matrix(nrow = 5, ncol = 4009))
colnames(resultsB_sub) <- c("Matrix B", "Matrix A", "Power", 
                            "Median deltaR2_pop", "sd deltaR2_pop", 
                            "Median deltaR2_sub", "sd deltaR2_sub",
                            "Median deltaR2_subpop", "sd deltaR2_subpop",
                            paste("pval", c(1:1000), sep = ""), 
                            paste("delR2_pop", c(1:1000), sep = ""),
                            paste("delR2_sub", c(1:1000), sep = ""),
                            paste("delR2_subpop", c(1:1000), sep = ""))
resultsB_sub[, 1] <- c("sim", "R2_pop", "R2_sub", "pvalReal", "column_count")
resultsB_sub[5, ] <- c(1:ncol(resultsB_sub))

resultsM_sub <- as.data.frame(matrix(nrow = 5, ncol = 4009))
colnames(resultsM_sub) <- c("Matrix B", "Matrix A", "Power", 
                            "Median deltaR2_pop", "sd deltaR2_pop", 
                            "Median deltaR2_sub", "sd deltaR2_sub",
                            "Median deltaR2_subpop", "sd deltaR2_subpop",
                            paste("pval", c(1:1000), sep = ""), 
                            paste("delR2_pop", c(1:1000), sep = ""),
                            paste("delR2_sub", c(1:1000), sep = ""),
                            paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_sub[, 1] <- c("sim", "R2_pop", "R2_sub", "pvalReal", "column_count")
resultsM_sub[5, ] <- c(1:ncol(resultsM_sub))

# The MEM are built for a full grid (50 x 25 cells):
# **************************************************

xy <- expand.grid(x = seq(1, 150, 1), y = seq(1, 75, 1))
# plot(xy, cex = 1)

nb <- cell2nb(nrow = 75, ncol = 150, "queen")
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
#####################################################################
#####################################################################
# ***************************************************************** #

# Creation of the response variable:
# **********************************
# Creation of 'y_spa' and 'y_noise' and standardisation:
# ******************************************************

set.seed(1)

y_spa_broad <- MEM[, 1] + MEM[, 2] + MEM[, 3]
y_spa_med <- MEM[, 306] + MEM[, 307] + MEM[, 308]
y_noise <- rnorm(n = nrow(MEM), mean = 0, sd = 1)

y_spa_broad_st <- scale(y_spa_broad)
y_spa_med_st <- scale(y_spa_med)
y_noise_st <- scale(y_noise)

# Creation of the response variable 'y' at the whole population level (pop):
# **************************************************************************

y_broad <- (a * y_spa_broad_st) + ((1-a) * y_noise_st)
y_med <- (a * y_spa_med_st) + ((1-a) * y_noise_st)

# s.value(xy, y_broad)

R2_pop_broad <- cor(y_broad, y_spa_broad_st)^2
R2_pop_med <- cor(y_med, y_spa_med_st)^2

resultsB_pop[2, c(1010:2009, 3010:4009)] <- R2_pop_broad
resultsM_pop[2, c(1010:2009, 3010:4009)] <- R2_pop_med

# Begining of the simulation process:
# ***********************************

# The simulation is run first at the broad scale, and then at the medium scale.
# Since the simulation at the medium scale have to be done exactly on the same
# subsampled sites, and therefore on the same 'C' matrices and same 'MEMsub', we 
# keep all the 'C' and 'MEMsub' computed at the broad scale in two lists, and we
# reuse them at the medium scale to spare time.

C.list <- vector("list", nperm)
candidates.list <- vector("list", nperm)

# To compute the selection percentages of each W matrix:
namesw <- c("DBMEM_PCNM", "Delaunay_Binary weighting", "Delaunay_Linear weighting",
           "Delaunay_Concave down weighting", "Delaunay_Concave up weighting",
           "Gabriel_Binary weighting", "Gabriel_Linear weighting", 
           "Gabriel_Concave down weighting", "Gabriel_Concave up weighting",
           "Relative neighbourhood_Binary weighting", 
           "Relative neighbourhood_Linear weighting", 
           "Relative neighbourhood_Concave down weighting",
           "Relative neighbourhood_Concave up weighting", 
           "Minimum spanning tree_Binary weighting",
           "Minimum spanning tree_Linear weighting",
           "Minimum spanning tree_Concave down weighting",
           "Minimum spanning tree_Concave up weighting")
bestmod_indexB <- as.data.frame(matrix(nrow = length(namesw), ncol = nperm + 2))
colnames(bestmod_indexB) <- c("W matrix", "Proportion", paste("best", c(1:nperm),
                                                              sep = ""))
bestmod_indexB[, 1] <- namesw

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
    x1 <- runif(39, min = sample(c(1:15), 1), max = sample(c(36:42), 1))
    y1 <- runif(39, min = sample(c(39:51), 1), max = sample(c(66:75), 1))
    x2 <- runif(39, min = sample(c(54:63), 1), max = sample(c(81:93), 1))
    y2 <- runif(39, min = sample(c(36:49), 1), max = sample(c(66:75), 1))
    x3 <- runif(39, min = sample(c(99:114), 1), max = sample(c(135:148), 1))
    y3 <- runif(39, min = sample(c(1:15), 1), max = sample(c(30:45), 1))
    
    C[, 1] <- rbind(x1, x2, x3)
    C[, 2] <- rbind(y1, y2, y3)
    
  } else {          # design = "random"
    
    C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
    set.seed(i)
    C[, 1] <- runif(117, min = 1, max = 148)
    C[, 2] <- runif(117, min = 1, max = 75) 
    
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
  C.list[[i]] <- C
  
  # We keep the lines of MEM that correspond to the sampled cells ('tri'):
  # **********************************************************************
  MEMsub <- y_spa_broad_st[sort]

  # We sample the response variable within the sampled cells ('y_sub'):
  # *******************************************************************
  y_sub <- y_broad[sort]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsB_pop[4, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, y_spa_broad_st[sort])^2                                          
  resultsB_pop[3, c(2009+i, 3009+i)] <- R2_sub
  resultsB_pop[1, 3009+i] <- R2_sub - R2_pop_broad

  # MEM.modsel function: Optimisation of the W matrix:
  # **************************************************
  # **************************************************
   candidates <- listw.candidates(C, style = style)
   candidates.list[[i]] <- candidates
  
   memsel <- MEM.modsel(y_sub, candidates, autocor = MEM_model)
   if (class(memsel) != "NULL") {
      resultsB_pop[1, 9+i] <- memsel$pval
      resultsB_pop[1, 1009+i] <- memsel$R2adj - R2_pop_broad
      resultsB_pop[1, 2009+i] <- memsel$R2adj - R2_sub
      bestmod_indexB[which(bestmod_indexB[, 1] == memsel$name), 2+i] <- 1
   } else resultsB_pop[1, 9+i] <- 1
}         # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

resultsB_pop[1, 3] <- length(which(resultsB_pop[1, c(10:(nperm + 9))] <= 0.05)) / nperm
resultsB_pop[1, 4] <- median(na.omit(as.numeric(resultsB_pop[1, c(1010:(nperm + 
                                                                          1009))])))
resultsB_pop[1, 5] <- sd(na.omit(as.numeric(resultsB_pop[1, c(1010:(nperm + 1009))])))
resultsB_pop[1, 6] <- median(na.omit(as.numeric(resultsB_pop[1, c(2010:(nperm + 
                                                                          2009))])))
resultsB_pop[1, 7] <- sd(na.omit(as.numeric(resultsB_pop[1, c(2010:(nperm + 2009))])))
resultsB_pop[1, 8] <- median(na.omit(as.numeric(resultsB_pop[1, c(3010:(nperm + 
                                                                          3009))])))
resultsB_pop[1, 9] <- sd(na.omit(as.numeric(resultsB_pop[1, c(3010:(nperm + 3009))])))

# Correct significance detection rate
resultsB_pop[4, 3] <- length(which(resultsB_pop[4, c(10:(nperm + 9))] <= 0.05)) / nperm

# Proportion of selection for each W matrix:
for (j in 1:nrow(bestmod_indexB)) {
  bestmod_indexB[j, 2] <- length(which(as.numeric(bestmod_indexB[j, c(3:nperm)]) 
                                       == 1)) / nperm
}

# Output of the results:
# **********************
res_file_name <- paste("Results_MEM.modsel_pop", correction, intensity, "Broad", 
                       design, paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsB_pop, file = res_file_name, sep = "\t")

bestmod_index_file_name <- paste("Bestmod-prop_MEM.modsel_pop", correction, 
                            intensity, "Broad", design, paste("style", style, ".txt", 
                                                              sep = ""), sep = "_")
write.table(bestmod_indexB, file = bestmod_index_file_name, sep = "\t")

   ########################
   ########################
   ### II. Medium scale ###
   ########################
   ########################

bestmod_indexM <- as.data.frame(matrix(nrow = length(namesw), ncol = nperm + 2))
colnames(bestmod_indexM) <- c("W matrix", "Proportion", paste("best", c(1:nperm),
                                                              sep = ""))
bestmod_indexM[, 1] <- namesw

for (i in 1:nperm) {
  
  # Sampling scheme:
  # ****************
  C <- C.list[[i]]
  
  # We keep the lines of MEM that correspond to the sampled cells ('tri'):
  # **********************************************************************
  MEMsub <- y_spa_med_st[as.numeric(row.names(C))]
  
  # We sample the response variable within the sampled cells ('y_sub'):
  # *******************************************************************
  y_sub <- y_med[as.numeric(row.names(C))]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsM_pop[4, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, y_spa_med_st[sort])^2     
  resultsM_pop[3, c(2009+i, 3009+i)] <- R2_sub
  resultsM_pop[1, 3009+i] <- R2_sub - R2_pop_med

  # MEM.modsel function: Optimisation of the W matrix:
  # **************************************************
  # **************************************************
  candidates <- candidates.list[[i]]
    
  bestmod_indexM <- as.data.frame(matrix(nrow = len, ncol = nperm + 2))
  colnames(bestmod_indexM) <- c("W matrix", "Proportion", paste("best", c(1:nperm),
                                                                sep = ""))
  bestmod_indexM[, 1] <- names(candidates)
  
  memsel <- MEM.modsel(y_sub, candidates, autocor = MEM_model)
  if (class(memsel) != "NULL") {
    resultsM_pop[1, 9+i] <- memsel$pval
    resultsM_pop[1, 1009+i] <- memsel$R2adj - R2_pop_med
    resultsM_pop[1, 2009+i] <- memsel$R2adj - R2_sub
    bestmod_indexM[which(bestmod_indexM[, 1] == memsel$name), 2+i] <- 1
  } else resultsM_pop[1, 9+i] <- 1
}         # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

resultsM_pop[1, 3] <- length(which(resultsM_pop[1, c(10:(nperm + 9))] <= 0.05)) / nperm
resultsM_pop[1, 4] <- median(na.omit(as.numeric(resultsM_pop[1, c(1010:(nperm + 
                                                                          1009))])))
resultsM_pop[1, 5] <- sd(na.omit(as.numeric(resultsM_pop[1, c(1010:(nperm + 1009))])))
resultsM_pop[1, 6] <- median(na.omit(as.numeric(resultsM_pop[1, c(2010:(nperm + 
                                                                          2009))])))
resultsM_pop[1, 7] <- sd(na.omit(as.numeric(resultsM_pop[1, c(2010:(nperm + 2009))])))
resultsM_pop[1, 8] <- median(na.omit(as.numeric(resultsM_pop[1, c(3010:(nperm + 
                                                                          3009))])))
resultsM_pop[1, 9] <- sd(na.omit(as.numeric(resultsM_pop[1, c(3010:(nperm + 3009))])))

# Correct significance detection rate
resultsM_pop[4, 3] <- length(which(resultsM_pop[4, c(10:(nperm + 9))] <= 0.05)) / nperm

# Proportion of selection for each W matrix:
for (j in 1:nrow(bestmod_indexM)) {
  bestmod_indexM[j, 2] <- length(which(as.numeric(bestmod_indexM[j, c(3:nperm)]) 
                                       == 1)) / nperm
}

# Output of the results:
# **********************
res_file_name <- paste("Results_MEM.modsel_pop", correction, intensity, "Medium", 
                       design, paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsM_pop, file = res_file_name, sep = "\t")

bestmod_index_file_name <- paste("Bestmod-prop_MEM.modsel_pop", correction, intensity, 
                                 "Medium", design, paste("style", style, ".txt", 
                                                         sep = ""), sep = "_")
write.table(bestmod_indexM, file = bestmod_index_file_name, sep = "\t")



# ***************************************************************** #
### II. The sampling design remains unchanged and the response varies:
#####################################################################
#####################################################################
# ***************************************************************** #

# Sampling scheme:
# ****************

set.seed(1)

# Sampling scheme:
# ****************

if (design == "clustered") {
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  x1 <- runif(39, min = sample(c(1:15), 1), max = sample(c(36:42), 1))
  y1 <- runif(39, min = sample(c(39:51), 1), max = sample(c(66:75), 1))
  x2 <- runif(39, min = sample(c(54:63), 1), max = sample(c(81:93), 1))
  y2 <- runif(39, min = sample(c(36:49), 1), max = sample(c(66:75), 1))
  x3 <- runif(39, min = sample(c(99:114), 1), max = sample(c(135:148), 1))
  y3 <- runif(39, min = sample(c(1:15), 1), max = sample(c(30:45), 1))
  
  C[, 1] <- rbind(x1, x2, x3)
  C[, 2] <- rbind(y1, y2, y3)
  
} else {          # design = "random"
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
  C[, 1] <- runif(117, min = 1, max = 148)
  C[, 2] <- runif(117, min = 1, max = 75) 
  
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
candidates <- listw.candidates(C, style = style)

   ######################
   ######################
   ### I. Broad scale ###
   ######################
   ######################

# We keep the lines of MEM that correspond to the sampled cells ('tri'):
# **********************************************************************
MEMsub <- y_spa_broad_st[sort]

# To compute the selection percentages of each W matrix:
bestmod_indexB <- as.data.frame(matrix(nrow = length(namesw), ncol = nperm + 2))
colnames(bestmod_indexB) <- c("W matrix", "Proportion", paste("best", c(1:nperm),
                                                              sep = ""))
bestmod_indexB[, 1] <- namesw

for (i in 1:nperm) { 
  
  set.seed(i)
  
  y_noise <- rnorm(n = nrow(MEM), mean = 0, sd = 1)
  y_noise_st <- scale(y_noise)
  
  # Creation of the response variable 'y' at the whole population level (pop):
  # **************************************************************************
  
  y_broad <- (a * y_spa_broad_st) + ((1-a) * y_noise_st)
  R2_pop_broad <- cor(y_broad, y_spa_broad_st)^2
  resultsB_sub[2, c(1009+i, 3009+i)] <- R2_pop_broad
  
  y_sub <- y_broad[sort]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsB_sub[4, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, y_spa_med_st[sort])^2                                          
  resultsB_sub[3, c(2009+i, 3009+i)] <- R2_sub
  resultsB_sub[1, 3009+i] <- R2_sub - R2_pop_broad

  # MEM.modsel function: Optimisation of the W matrix:
  # **************************************************
  # **************************************************
  
  memsel <- MEM.modsel(y_sub, candidates, autocor = MEM_model)
  if (class(memsel) != "NULL") {
    resultsB_sub[1, 9+i] <- memsel$pval
    resultsB_sub[1, 1009+i] <- memsel$R2adj - R2_pop_broad
    resultsB_sub[1, 2009+i] <- memsel$R2adj - R2_sub
    bestmod_indexB[which(bestmod_indexB[, 1] == memsel$name), 2+i] <- 1
  } else resultsB_sub[1, 9+i] <- 1
}         # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

resultsB_sub[1, 3] <- length(which(resultsB_sub[1, c(10:(nperm + 9))] <= 0.05)) / nperm
resultsB_sub[1, 4] <- median(na.omit(as.numeric(resultsB_sub[1, c(1010:(nperm + 
                                                                          1009))])))
resultsB_sub[1, 5] <- sd(na.omit(as.numeric(resultsB_sub[1, c(1010:(nperm + 1009))])))
resultsB_sub[1, 6] <- median(na.omit(as.numeric(resultsB_sub[1, c(2010:(nperm + 
                                                                          2009))])))
resultsB_sub[1, 7] <- sd(na.omit(as.numeric(resultsB_sub[1, c(2010:(nperm + 2009))])))
resultsB_sub[1, 8] <- median(na.omit(as.numeric(resultsB_sub[1, c(3010:(nperm + 
                                                                          3009))])))
resultsB_sub[1, 9] <- sd(na.omit(as.numeric(resultsB_sub[1, c(3010:(nperm + 3009))])))

# Correct significance detection rate
resultsB_sub[4, 3] <- length(which(resultsB_sub[4, c(10:(nperm + 9))] <= 0.05)) / nperm

# Proportion of selection for each W matrix:
for (j in 1:nrow(bestmod_indexB)) {
  bestmod_indexB[j, 2] <- length(which(as.numeric(bestmod_indexB[j, c(3:nperm)]) 
                                       == 1)) / nperm
}

# Output of the results:
# **********************
res_file_name <- paste("Results_MEM.modsel_sub", correction, intensity, "Broad", 
                       design, paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsB_sub, file = res_file_name, sep = "\t")

bestmod_index_file_name <- paste("Bestmod-prop_MEM.modsel_sub", correction, 
                                 intensity, "Broad", design, paste("style", style, ".txt", 
                                                                   sep = ""), sep = "_")
write.table(bestmod_indexB, file = bestmod_index_file_name, sep = "\t")

   ########################
   ########################
   ### II. Medium scale ###
   ########################
   ########################

# We keep the lines of MEM that correspond to the sampled cells ('tri'):
# **********************************************************************
MEMsub <- y_spa_med_st[sort]

# To compute the selection percentages of each W matrix:
bestmod_indexM <- as.data.frame(matrix(nrow = length(namesw), ncol = nperm + 2))
colnames(bestmod_indexM) <- c("W matrix", "Proportion", paste("best", c(1:nperm),
                                                              sep = ""))
bestmod_indexM[, 1] <- namesw

for (i in 1:nperm) { 
  
  set.seed(i)
  
  y_noise <- rnorm(n = nrow(MEM), mean = 0, sd = 1)
  y_noise_st <- scale(y_noise)
  
  # Creation of the response variable 'y' at the whole population level (pop):
  # **************************************************************************
  
  y_med <- (a * y_spa_med_st) + ((1-a) * y_noise_st)
  R2_pop_med <- cor(y_med, y_spa_med_st)^2
  resultsM_sub[2, c(1009+i, 3009+i)] <- R2_pop_med
  
  y_sub <- y_med[sort]
  
  # Real p-value and R2_sub:
  # ************************
  lm <- lm(y_sub ~ MEMsub)
  resultsM_sub[4, 9+i] <- lmp(lm)
  R2_sub <- cor(y_sub, y_spa_med_st[sort])^2     
  resultsM_sub[3, c(2009+i, 3009+i)] <- R2_sub
  resultsM_sub[1, 3009+i] <- R2_sub - R2_pop_med
  
  # MEM.modsel function: Optimisation of the W matrix:
  # **************************************************
  # **************************************************
  
  memsel <- MEM.modsel(y_sub, candidates, autocor = MEM_model)
  if (class(memsel) != "NULL") {
    resultsM_sub[1, 9+i] <- memsel$pval
    resultsM_sub[1, 1009+i] <- memsel$R2adj - R2_pop_broad
    resultsM_sub[1, 2009+i] <- memsel$R2adj - R2_sub
    bestmod_indexM[which(bestmod_indexM[, 1] == memsel$name), 2+i] <- 1
  } else resultsM_sub[1, 9+i] <- 1
}         # End of the simulation ('for') loop

# Median and standard deviation of the deltaR2:
# *********************************************

resultsM_sub[1, 3] <- length(which(resultsM_sub[1, c(10:(nperm + 9))] <= 0.05)) / nperm
resultsM_sub[1, 4] <- median(na.omit(as.numeric(resultsM_sub[1, c(1010:(nperm + 
                                                                          1009))])))
resultsM_sub[1, 5] <- sd(na.omit(as.numeric(resultsM_sub[1, c(1010:(nperm + 1009))])))
resultsM_sub[1, 6] <- median(na.omit(as.numeric(resultsM_sub[1, c(2010:(nperm + 
                                                                          2009))])))
resultsM_sub[1, 7] <- sd(na.omit(as.numeric(resultsM_sub[1, c(2010:(nperm + 2009))])))
resultsM_sub[1, 8] <- median(na.omit(as.numeric(resultsM_sub[1, c(3010:(nperm + 
                                                                          3009))])))
resultsM_sub[1, 9] <- sd(na.omit(as.numeric(resultsM_sub[1, c(3010:(nperm + 3009))])))

# Correct significance detection rate
resultsM_sub[4, 3] <- length(which(resultsM_sub[4, c(10:(nperm + 9))] <= 0.05)) / nperm

# Proportion of selection for each W matrix:
for (j in 1:nrow(bestmod_indexM)) {
  bestmod_indexM[j, 2] <- length(which(as.numeric(bestmod_indexM[j, c(3:nperm)]) 
                                       == 1)) / nperm
}

# Output of the results:
# **********************
res_file_name <- paste("Results_MEM.modsel_sub", correction, intensity, "Medium", 
                       design, paste("style", style, ".txt", sep = ""), sep = "_")
write.table(resultsM_sub, file = res_file_name, sep = "\t")

bestmod_index_file_name <- paste("Bestmod-prop_MEM.modsel_sub", correction, 
                                 intensity, "Medium", design, paste("style", style, ".txt", 
                                                                   sep = ""), sep = "_")
write.table(bestmod_indexB, file = bestmod_index_file_name, sep = "\t")

