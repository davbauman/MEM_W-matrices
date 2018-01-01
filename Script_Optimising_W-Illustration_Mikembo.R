# Useful packages:
# ****************

library(vegan)
library(adespatial)
library(spdep)
library(usdm)

source("listw.candidates.R")
source("MEM.modsel.R")

##################
# Fix parameters #
##################

# Abondance totale minimum sous laquelle on ne considère pas les espèces : 
ab <- 20
# Nb of quadrat to be randomly sampled:
sample_size = 160

##############
# Data input #
##############

C <- read.table("spa_25kr.txt", h = T, sep = "\t", row.names = 1)

env <- read.table("env_25kr.txt", h = T, sep = "\t", row.names = 1)

spe <- read.table("spe_25kr.txt", h = T, sep = "\t", row.names = 1)
sum <- c() ; for(i in 1:ncol(spe)) sum <- c(sum, sum(spe[,i]))
spe <- spe[, which(sum >= ab)]

# Parcimonious environmental model selection:
# *******************************************

vif(env)
vifcor(env, th = 0.7)
env <- env[, -c(2,4,6,8,10,12,13,15:17,19:22,25,27,29:31)]   # Pour cor < 70 %
dim(env)

###################
# Result matrices #
###################

results <- as.data.frame(matrix(nrow = ncol(spe), ncol = 7))
row.names(results) <- colnames(spe)
colnames(results) <- c("pval dbMEM_PCNM", "R2adj dbMEM_PCNM", "pval random choice", 
                       "R2adj random choice", "pval Optim", "R2adj Optim", 
                       "Optim_name")

MEM_model = "positive"   # Sign of the eigenvectors we are interested in

# We sample the grid to obtain an irregular sampling design:

set.seed(1)
sample <- sample(size = sample_size, x = c(1:160), replace = F)
spe_sam <- spe[sample, ]
C_sam <- C[sample, ]

for (a in 1:ncol(spe)) {
  
  Y <- spe_sam[, a]
  
  ###################
  # I. MEM Analysis #
  ###################
  
dbMEM_listw <- listw.candidates(C_sam, style = 'B', PCNM = T, del = F, gab = F, rel = F, 
                                mst = F, DB = F)
dbMEM <- MEM.modsel(Y, dbMEM_listw[[1]])

class <- class(dbMEM)
if (class == "list") {
  results[a, 1] <- dbMEM$best$pval
  results[a, 2] <- dbMEM$best$R2adj
}

thresholds <- seq(give.thresh(dist(C_sam)), max(dist(C_sam)), le = 10)
all <- listw.candidates(C_sam, style = 'B', DBthresh = sample(1, x = thresholds))
random_listw <- all[[sample(size = 1, c(1:length(all)))]]
random <- MEM.modsel(Y, random_listw)

class <- class(random)
if (class == "list") {
  results[a, 3] <- random$best$pval
  results[a, 4] <- random$best$R2adj
}

candid <- listw.candidates(C_sam, style = 'B', del = F, rel = F, DB = F, bin = F, flin = F)
optim <- MEM.modsel(Y, candid)

class <- class(optim)
if (class == "list") {
  results[a, 5] <- optim$best$pval
  results[a, 6] <- optim$best$R2adj
  results[a, 7] <- optim$best$name
}

}

write.table(results, paste("results_", sample_size, ".txt", sep = ""), sep = "\t")

# Complete grid (160 quadrats):
# *****************************

results <- as.data.frame(matrix(nrow = ncol(spe), ncol = 2))
row.names(results) <- colnames(spe)
colnames(results) <- c("pval queen", "R2adj queen")

MEM_model = "positive"   # Sign of the eigenvectors we are interested in

#xy <- expand.grid(x = seq(1, 20, 1), y = seq(1, 8, 1))
nb <- cell2nb(nrow = 20, ncol = 8, "queen")
nb2 <- nb2listw(nb, style = "B")
MEM <- scores.listw(nb2, MEM.autocor = MEM_model)
s.value(xy, MEM[,1])

# Réorganisons les lignes de "spe" de façon à correspondre à celles de "C" et "MEM":
spa <- read.table("spa160_reference.txt", h = T, sep = "\t", row.names = 1)
rnames <- row.names(spa)
spe2 <- spe
row.names(spe2) <- c(1:nrow(spe2))
for (i in 1:length(rnames)) {
  n <- which(row.names(spe) == rnames[i])
  spe2[i, ] <- spe[n, ]
  row.names(spe2)[i] <- row.names(spe)[n]
}

for (a in 1:ncol(spe)) {
  
  Y <- spe2[, a]
  ref <- MEM.modsel(Y, nb2)
  
  class <- class(ref)
  if (class == "list") {
    results[a, 1] <- ref$best$pval
    results[a, 2] <- ref$best$R2adj
  }

}

write.table(results, "results_160.txt", sep = "\t")
