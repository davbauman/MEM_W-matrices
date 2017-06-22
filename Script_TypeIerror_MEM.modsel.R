######################################################################################
#*****# Test and comparison of different W matrices (MEM) - A simulation study #*****#
######################################################################################

# Type I error rate of the MEM.modsel function with and without p-value correction:
# #################################################################################

rm(list=ls()[-match(c("xy", "nb", "nb2", "nb2_backup", "MEM", "MEM_backup"), ls())])

# Usefull packages and functions:
# *******************************

library(vegan)
library(adespatial)
library(spdep)

source("lmp.R")
source("listw.candidates.R")
source("MEM.modsel.R")

# Definition of the simulation parameters:
# ****************************************

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"          # Either "positive" or "negative"

# Framework:
framework <- "univariate"       # Either "univariate" or "multivariate"

# Sampling design:
design <- "clustered"           # Either "clustered" or "random"

nperm <- 1000

style <- "B"                  # Either "B" or "W"
correction <- "Corrected"     # Either "Corrected" or "Uncorrected"

if (correction == "Corrected") {
  corr <- TRUE
} else corr <- FALSE

# Construction of a result matrix:
# ********************************

# For each B matrix, no A matrix and three A matrices tested; Two additional lines for
# the original PCNM method and for the PCNM method computed through MEM (DBEM, see further).
# For the columns: lines 3 = type I error; line 4 = mean R2adj; line 5 = sd of the R2adj;
# 1000 permutations, so that lines 6 to 1005 contain p-values, and lines 1006 to 2005 
# contain R2adj.

# For DBMEM (10 different dmax thresholds), a binary and a linear weighting function:

results <- as.data.frame(matrix(nrow = 1, ncol = 2006))

colnames(results) <- c("MEM.modsel", "/", "type I error", "mean R2adj", "sd R2adj", "max R2adj",
   paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), sep = ""))

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
# We build the nb2listW list:
candidates <- listw.candidates(C, style = style)

# Here begins the simulation process:
# ***********************************

# Simulation procedure:
#######################

for (i in 1:nperm) {   
   
     set.seed(i)
     
     if (framework == "univariate") {
       
       Y <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"  
       #   Y <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
       #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
       #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"            
       
     } else {
       
       Y <- matrix(ncol = 5, nrow = nrow(C))
       for (b in 1:ncol(Y)) {
         Y[, b] <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"    
         #     Y[, b] <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
         #     Y[, b] <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
         #     Y[, b] <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3" 
       }                  
     }
     
   best.spa <- MEM.modsel(Y, candidates, autocor = MEM_model, correction = corr)

   if (is.list(best.spa) == TRUE) {
      results[1, i+6] <- best.spa$pval
      results[1, i+1006] <- best.spa$R2adj
   } else { 
      results[1, i+6] <- 1
      results[1, i+1006] <- NA
   }

}

# Summary of the results:
#########################

   results[1, 3] <- length(which(results[i, c(7:(nperm + 6))] <= 0.05)) / nperm
   results[1, 4] <- mean(as.numeric(results[i, c(1007:(nperm + 1006))]))
   results[1, 5] <- sd(as.numeric(results[i, c(1007:(nperm + 1006))]))
   results[1, 6] <- max(as.numeric(results[i, c(1007:(nperm + 1006))]))

# Output of the results:
# **********************

res_file_name <- paste("Results_tIerror_MEM.modsel", correction, framework, ran, 
                       paste(design, ".txt", sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")
