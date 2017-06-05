# ****************************************************************** #
### Connectivity parameters of the different connectivity matrices ###
# ****************************************************************** #

# Usefull packages and functions:
# *******************************

library(vegan)
library(adespatial)
library(spdep)

# Definition of the simulation parameters:
# ****************************************

# Sampling design:

design <- "random"    # Either "clustered" or "random"

nperm <- 1000

# Results matrices for the connectivity parameters of the different W matrices:
# *****************************************************************************
# Sera calculé quand on fait varier le design d'échantillonnage en maintenant
# fixe la variable réponse, dès lors que c'est la position des quadrats échantillonnés
# qui sera déterminant pour la connectivité (et non pas les valeurs dans les quadrats)
# Les paramètres ne dépendant que de la matrice de connectivité (la pondération,
# l'échelle ou l'intensité n'influencent pas).

connect_param <- as.data.frame(matrix(nrow = 8, ncol = 3007))
colnames(connect_param) <- c("Matrix B", "Median 0_prop", "sd 0_prop", 
                             "Median mean_nb_neigh", "sd mean_nb_neigh",
                             "Median sd_nb_neigh", "sd sd_nb_neigh",
                             paste("0_prop", c(1:1000), sep = ""),
                             paste("nb_neigh", c(1:1000), sep = ""))
connect_param[, 1] <- c("del", "gab", "rel", "mst", "DBmin", "DBmed","DBmax",
                        "column_count")
connect_param[8, ] <- c(1:ncol(connect_param))

# The MEM are built for a full grid (50 x 25 cells):
# **************************************************

xy <- expand.grid(x = seq(0.02, 1, 0.02), y = seq(0.04, 1, 0.04))
nb <- cell2nb(nrow = 25, ncol = 50, "queen")
nb2 <- nb2listw(nb)
MEM <- scores.listw(nb2, MEM.autocor = MEM_model)

# ***************************************************************** #
### I. The response remains unchanged and the sampling design varies:
# ***************************************************************** #

for (i in 1:nperm) {
  
  # Sampling scheme:
  # ****************
  
  if (design == "clustered") {
    C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
    set.seed(i)
    x1 <- runif(39, min = 0.04, max = 0.28)
    y1 <- runif(39, min = 0.34, max = 0.5)
    x2 <- runif(39, min = 0.38, max = 0.68)
    y2 <- runif(39, min = 0.3, max = 0.5)
    x3 <- runif(39, min = 0.6, max = 1)
    y3 <- runif(39, min = 0.02, max = 0.2)
    
    C[,1] <- rbind(x1, x2, x3)
    C[,2] <- rbind(y1, y2, y3)
    
    xy.d1 <- dist(C)
    
  } else {          # design = "random"
    
    C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
    set.seed(i)
    C[,1] <- runif(117, min = 0.02, max = 1)
    C[,2] <- runif(117, min = 0.02, max = 0.5) 
    
    xy.d1 <- dist(C)
  }
  
  # Attribute each sampled point to one of the grid cells:
  # ******************************************************
  
  grid.size <- 0.02
  tri <- c()
  for (k in 1:nrow(C)) {
    x <- floor((C[k, 1] * 50) / (grid.size * 50))
    y <- floor((C[k, 2] * 50) / (grid.size * 50))
    N <- y * 50 + x + 1
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
  C.list[[i]] <- C
  
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
  
  # Retrieval of the connectivity parameters:
  # *****************************************
  nb.list <- list(Y.del, Y.gab, Y.rel, Y.mst, Y.listDB[[1]], Y.listDB[[2]], 
                  Y.listDB[[3]])
  for (j in 1:length(nb.list)) {
    connect_param[j, 7+i] <- 1 - (sum(card(nb.list[[j]])) / length(nb.list[[j]])^2)
    connect_param[j, 1007+i] <- mean(card(nb.list[[j]])) 
    connect_param[j, 2007+i] <- sd(card(nb.list[[j]]))
  }
  
} # End of the simulation ('for') loop

# Summary of the results:
# ***********************

for (i in 1:nrow(connect_param)-1) {
  connect_param[i, 2] <- median(as.numeric(connect_param[i, c(8:(nperm + 7))]))
  connect_param[i, 3] <- sd(as.numeric(connect_param[i, c(8:(nperm + 7))]))
  connect_param[i, 4] <- median(as.numeric(connect_param[i, c(1008:(nperm + 1007))]))
  connect_param[i, 5] <- sd(as.numeric(connect_param[i, c(1008:(nperm + 1007))])) 
  connect_param[i, 6] <- median(as.numeric(connect_param[i, c(2008:(nperm + 2007))]))
  connect_param[i, 7] <- sd(as.numeric(connect_param[i, c(2008:(nperm + 2007))]))
}

# Output of the results:
# **********************
connect_file_name <- paste("Results_Connectivity_param", 
                           paste(design, ".txt", sep = ""), sep = "_")
write.table(connect_param, file = connect_file_name, sep = "\t")
