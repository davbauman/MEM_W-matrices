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

source("lmp.R")

# Definition of the simulation parameters:
# ****************************************

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"   # Either "positive" or "negative"

# Sampling design:
design <- "clustered"    # Either "clustered" or "random"

nperm <- 1000

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
resultsB_pop[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "R2_W", "pvalReal")
resultsB_pop[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/", "/")

resultsM_pop <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsM_pop) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_pop[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "R2_W", "pvalReal")
resultsM_pop[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/", "/")

resultsB_sub <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsB_sub) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsB_sub[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "R2_W", "pvalReal")
resultsB_sub[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/", "/")

resultsM_sub <- as.data.frame(matrix(nrow = 33, ncol = 4009))
colnames(resultsM_sub) <- c("Matrix B", "Matrix A", "Power", 
                          "Median deltaR2_pop", "sd deltaR2_pop", 
                          "Median deltaR2_sub", "sd deltaR2_sub",
                          "Median deltaR2_subpop", "sd deltaR2_subpop",
                           paste("pval", c(1:1000), sep = ""), 
                           paste("delR2_pop", c(1:1000), sep = ""),
                           paste("delR2_sub", c(1:1000), sep = ""),
                           paste("delR2_subpop", c(1:1000), sep = ""))
resultsM_sub[,1] <- c(rep(c("del", "gab", "rel", "mst", "DBmin", "DBmed",
                      "DBmax"), each = 4), "DB_PCNM", "R2_pop", "R2_sub",
                      "R2_W", "pvalReal")
resultsM_sub[,2] <- c(rep(c("None", "Linear", "Concave-down", "Concave-up"), 
                      times = 7), "1-(D/4t)^2", "/", "/", "/", "/")

