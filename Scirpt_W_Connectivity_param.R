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

xy <- expand.grid(x = seq(1, 50, 1), y = seq(1, 25, 1))
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
    x1 <- runif(39, min = 2, max = 14)
    y1 <- runif(39, min = 17, max = 25)
    x2 <- runif(39, min = 19, max = 34)
    y2 <- runif(39, min = 15, max = 25)
    x3 <- runif(39, min = 30, max = 50)
    y3 <- runif(39, min = 1, max = 10)
    
    C[, 1] <- rbind(x1, x2, x3)
    C[, 2] <- rbind(y1, y2, y3)
    
  } else {          # design = "random"
    
    C <- as.matrix(matrix(0, ncol = 2, nrow = 117))
    set.seed(i)
    C[, 1] <- runif(117, min = 1, max = 50)
    C[, 2] <- runif(117, min = 1, max = 25) 
    
  }
  
  # Attribute each sampled point to one of the grid cells:
  # ******************************************************
  
  # /25 = taille d'un côté du quadrat ; pour que la numérotation des quadrats se
  # fasse de gauche à droite en partant du coin inférieur gauche de la grille :
  # N <- y * (nb horizontal de quadrats) + x + 1 ; pour que la numérotation des 
  # quadrats se fasse de bas en haut en partant du coin inférieur gauche :
  # N <- x * (nb vertical de quadrats) + y + 1
  
  #  # Comme les coord. sont comprises entre 0.02 et 1 et 0.02 et 0.5, on les multiplie
  #  # par 50 (momentanément) pour attributer une cellule à chaque point.
  
  grid.size <- 1
  tri <- c()
  for (k in 1:nrow(C)) {
    x <- floor((C[k, 1]) / (grid.size))
    y <- floor((C[k, 2]) / (grid.size))
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
  xy.d1 <- dist(C)
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
