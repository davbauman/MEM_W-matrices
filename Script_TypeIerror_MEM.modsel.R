######################################################################################
#*****# Test and comparison of different W matrices (MEM) - A simulation study #*****#
######################################################################################

# Type I error rate of the MEM.modsel function with and without p-value correction:
# #################################################################################

rm(list=ls()[-match(c("xy", "nb", "nb2", "MEM", "MEM_backup"), ls())])

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

xy <- expand.grid(x = seq(1, 90, 1), y = seq(1, 90, 1))

set.seed(5)

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

# We build the nb2listW list:
candidates <- listw.candidates(C, style = style)

# Here begins the simulation process:
# ***********************************

# Simulation procedure:
#######################

for (i in 1:nperm) {   
   
     set.seed(i)

       Y <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"  
       #   Y <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
       #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
       #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"            
       
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

   results[1, 3] <- length(which(results[1, c(7:(nperm + 6))] <= 0.05)) / nperm
   results[1, 4] <- mean(as.numeric(results[1, c(1007:(nperm + 1006))]))
   results[1, 5] <- sd(as.numeric(results[1, c(1007:(nperm + 1006))]))
   results[1, 6] <- max(as.numeric(results[1, c(1007:(nperm + 1006))]))

# Output of the results:
# **********************
res_file_name <- paste("Results_tIerror_MEM.modsel", correction, ran, 
                       paste(design, ".txt", sep = ""), sep = "_")
write.table(results, file = res_file_name, sep = "\t")
