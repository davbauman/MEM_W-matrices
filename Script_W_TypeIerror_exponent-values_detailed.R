#######################################################################################
#######################################################################################
### We test the type I error rate of the DBMEM process for one given distance threshold
### (fixed connectivity matrix) for which we test several exponent values in the 
### weighting matrix.
#######################################################################################
#######################################################################################


# Useful packages and functions:
# ******************************

library(vegan)
library(adespatial)
library(spdep)

source("test.W.R2.R")

# Construction of a result matrix:                               
# ********************************
# For DBMEM (10 different 'y' exponent value), in a concave-down and a concave-up
# weighting function:

results <- as.data.frame(matrix(nrow = 2, ncol = 2005))

colnames(results) <- c("Matrix A", "/", "type I error", "mean R2adj", "sd R2adj",
                       paste("p-val", c(1:1000), sep = ""), 
                       paste("R2_", c(1:1000), sep = ""))
results[, 1] <- c("conc-down", "conc-up")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"
#MEM_model = "negative"   ; autocor <- "neg" 

correction = FALSE    # TRUE (Sidak correction) or FALSE (no correction: p-val = 0.05)
if (correction == FALSE) corr <- "Uncorrected" else corr <- "Corrected"

design <- "clustered"   # or "clustered"

y_exp_f2 <- c(2:5)      # Exponent values to be compared
y_exp_f3 <- c(0.4, 1)   # Exponent values to be compared

nperm <- 1000

f2test <- TRUE # Si voulons tester f2 aussi. Si FALSE, ne testons que f3.

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

# II. Distance-based B (radius around points) using the largest edge of the MST:
# ******************************************************************************
lowlim <- give.thresh(xy.d1)
dmax <- max(xy.d1)

Y.DBMEM <- dnearneigh(x = as.matrix(C), d2 = lowlim, d1 = 0)

# Weighting functions and fixed parametres:
# *****************************************

f2 <- function (D, dmax, y) {1 - (D/dmax)^y}      # Concave-down function
#f3 <- function (D, dmax, y) {1 / (D/dmax)^y}      # Concave-up function
f3 <- function (D, y) {1 / D^y}

# III. Generation of MEM eigenvectors:
######################################

Y.DBMEM.f2.list <- lapply(y_exp_f2, function(x) test.W.R2(nb= Y.DBMEM, xy = C,
                                                           f = f2, y = x, dmax = dmax,
                                                           MEM.autocor = MEM_model))
Y.DBMEM.f3.list <- lapply(y_exp_f3, function(x) test.W.R2(nb= Y.DBMEM, xy = C,
                                                           f = f3, y = x, 
                                                           MEM.autocor = MEM_model))

# Correction of the p-value (Sidak):
# **********************************
# Define the number of distance threshold tested for the DBMEM:

nbtestf2 <- length(y_exp_f2)
nbtestf3 <- length(y_exp_f3)
alphaf2 <- 0.05
alphaf3 <- 0.05
if (correction == TRUE) {
  alphaf2 <- 1-((1-alphaf2)^(1/nbtestf2))   # Sidak correction
  alphaf3 <- 1-((1-alphaf3)^(1/nbtestf3))   # Sidak correction
} 

# Simulation procedure:
#######################

plot(c(1:1000), c(1:1000), type = "n")

for (i in 1:nperm) {   
  points(x = i, y = 1)
  set.seed(i)
   
    Y <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"             
    #   Y <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
    #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                      
    #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                   
  
  # Significance test and MEM variable selection:
  # *********************************************
  # f2:
  # ***
  if(f2test == TRUE) {
  R2.list <- lapply(Y.DBMEM.f2.list, 
                    function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
  thebest <- which.max(R2.list)
  MEM.best <- Y.DBMEM.f2.list[[thebest]]$best$MEM
  if (anova.cca(rda(Y, MEM.best), permutations = 9999)$Pr[1] <= alphaf2) {
    R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
    class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, 
                                           nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- MEM.best[, c(sign)]
      results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      results[1, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } 
  } else { 
    results[1, i+5] <- 1   # p-val made equal to 1 
    results[1, i+1005] <- NA
  }
  }
  # f3:
  # ***
  R2.list <- lapply(Y.DBMEM.f3.list, 
                    function(z) RsquareAdj(rda(Y, z$best$MEM))$adj.r.squared)
  thebest <- which.max(R2.list)
  MEM.best <- Y.DBMEM.f3.list[[thebest]]$best$MEM
  if (anova.cca(rda(Y, MEM.best), permutations = 9999)$Pr[1] <= alphaf3) {
    R2adj <- RsquareAdj(rda(Y, MEM.best))$adj.r.squared
    class <- class(try(fsel <- forward.sel(Y, MEM.best, adjR2thresh = R2adj, 
                                           nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- MEM.best[, c(sign)]
      results[2, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.best)))$Pr[1]
      results[2, i+1005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
    } 
  } else { 
    results[2, i+5] <- 1   # p-val made equal to 1 
    results[2, i+1005] <- NA
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

res_file_name <- paste("Results_tIerror_Exponent-selection", corr,
                       paste(design, ".txt", sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")
