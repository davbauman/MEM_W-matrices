######################################################################################
# Code for computing the Type I Error Rate of severall spatial weighting matrices (W) 
# obtained by combination of different connectivity and weighting matrices, using uni- 
# and multivariate response variables in a regular and an irregular sampling design.
######################################################################################

# Useful packages and functions:
# ******************************

library(vegan)
library(adespatial)
library(spdep)

source("test.W.R2.R")
# Function to compute the global test and FwdSel on the MEM variables (for typeIerror):
MEMfwd.test <- function (y, mem) {
  pval <- as.data.frame(anova.cca(rda(y, mem), permutations = 9999))$Pr[1]
  if (pval <= 0.05) {
    R2adj <- RsquareAdj(rda(y, mem))$adj.r.squared
    class <- class(try(fsel <- forward.sel(y, mem, adjR2thresh = R2adj, nperm = 999), TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- mem[, c(sign)]
      R2_W <- RsquareAdj(rda(y, MEM.FwdSel))$adj.r.squared
    } else { 
      R2_W <- NA
      pval <- 1      
    }
  } else { 
    R2_W <- NA
    pval <- 1      
  }
  list(pval = pval, R2_W = R2_W)
}

# Construction of a result matrix:                               
# ********************************
# For each B matrix, no A matrix (binary) and three A matrices tested.

results <- as.data.frame(matrix(nrow = 29, ncol = 2005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
                       paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), 
                                                                  sep = ""))
results[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed", "DBmax"), each = 4),
                 "DBMEM_PCNM")
results[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), times = 7), 
                 "1-(D/4t)^2")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"

# Sampling design:
design <- "clustered"   # "random" or "clustered"

style <- "B"                  # Either "B" or "W"

nperm <- 1000

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
thresh <- seq(lowlim, uplim, le = 10)   # 3 tested distances: thresh[1, 5, 9]
Y.listDB <- lapply(thresh[c(1, 5, 10)], dnearneigh, x = as.matrix(C), d1 = 0)

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

# Simulation procedure:
#######################

for (i in 1:nperm) {   
  
  set.seed(i)

    Y <- runif (nrow(C), min = 0, max = 20) ; ran <- "runif"  
    #   Y <- rnorm(nrow(C), mean = 0, sd = runif (1, 1, 3)) ; ran <- "rnorm" 
    #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
    #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"            

  # Significance test and MEM variable selection:
  # *********************************************
    for (q in 1:length(listW)) {
      test <- MEMfwd.test(Y, listW[[q]])
      results[q, i+5] <- test$pval
      results[q, i+1005] <- test$R2_W
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
