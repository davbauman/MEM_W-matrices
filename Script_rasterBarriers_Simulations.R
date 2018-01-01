## barr-rast est un raster crée par QGIS (fonction rasteriser après conversion d'un shapefile)
## rasteriser sur une grille de dimension grid.size x grid.size

library(adespatial)
library(spdep)
library(adegraphics)
library(raster)

grid.size <- 70

#setwd("/home/stephane/Documents/Collaborations/D_Bauman/Simul-barrieres/")
r2 = raster("barr-raster70")
r = raster()
crs(r) = '+proj=utm +zone=10' ## si on met rien, il considere une sphère (http://r-sig-geo.2731867.n2.nabble.com/gdistance-costDistance-with-barriers-td7587891.html)
bb = extent(0.5, grid.size + 0.5, 0.5, grid.size + 0.5) ## grille (centres des cellules)
extent(r) = bb
res(r) = 1
r[] = as.vector(t(as.matrix(r2)))
image(r, asp = 1)
r = r+1

## generate a full grid 
xy.complete = expand.grid(grid.size:1,1:grid.size)
xy.complete = xy.complete[,2:1]

load("simu-SD.RData")
## simu.pop contains the abundances of 30 species
## some plots to check the consistency between geographic coordinates
#s.value(xy.complete, as.vector(pop$p$habitat[1,,]), Sp = Sp1, pSp.col = c('grey', 'peru'), ppoint.cex = 0.1)

## subsampling the grid

buffer = 5
## exclude samples located in barriers
## closed to the borders to avoid edge effects
idx.selected <- which(as.vector(as.matrix(r))<10 & xy.complete[,1] < (grid.size - buffer) & xy.complete[,1]  > buffer & xy.complete[,2] < (grid.size - buffer) & xy.complete[,2]  > buffer)
xy.full <- xy.complete[idx.selected,]

## sample nsites
set.seed(12)
nsites = 100
idx.sample <- sample(idx.selected, nsites)
xy <- xy.complete[idx.sample,]
rownames(xy) = 1:nsites

par(mfrow = c(1,3))
image(r, asp = 1)
points(xy.complete, cex = 0.2)
image(r, asp = 1)
points(xy.full, cex = 0.2)
image(r, asp = 1)
points(xy , cex = 0.2)

## subsample the abundance data
species <- simu.pop[idx.sample,]

### compute MEM based on barriers
library(gdistance)
tr <- transition(r, transitionFunction=function(x) 1/mean(x), directions=8)
cout <- costDistance(tr,as.matrix(xy)) ## cout est une distance
par(mfrow = c(1, 1))
image(as.matrix(cout))
plot(raster(tr))

thr <- give.thresh(matdist = cout)
thr
mat01 <- ifelse(as.matrix(cout) > thr, 0, as.matrix(cout))
lw1 <- mat2listw(mat01)
me <- mem(lw1, MEM.autocor = "positive")
dim(me)
plot(lw1, xy, add = T)

Sp1 <- as(r, 'SpatialGridDataFrame')
s.label(xy, Sp=Sp1, nb = lw1, plabel.cex = 0.5, ppoint.cex = 2, pgrid.draw = F, pSp.col = c('grey', 'peru'))
s.value(xy, me[,1:9],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)
# Abundances of the first nine species in the sampled sites:
s.value(xy.complete, simu.pop[,1:9], ppoint.cex = 0.4 , Sp = Sp1, pSp.col = c('grey', 'peru'))

# Idée 1 : Créer matrice B avec Gabriel graph et pondération de sorte à
# ce que les distances supérieures à 'thresh' soit = 0. C'est une façon
# de secondairement éliminer les liens traversant les barrières.
#f1 <- function (D, dmax)     { 1 - (D/dmax) }
#dmax <- max(dist(xy))
#fun_barr <- function (x, th) {ifelse(x >= th, 0, f1(D = x, dmax))}
#Y.gab <- graph2nb(gabrielneigh(as.matrix(xy), nnmult = 5), sym = TRUE)
#listw.barr <- nb2listw(Y.gab, style = "B", glist = lapply(nbdists(Y.gab, as.matrix(xy)), 
#                                                          fun_bar, th = thr))
#mebarr <- scores.listw(listw.barr, MEM.autocor = "positive")
#RsquareAdj(rda(Y, mebarr))$adj.r.squared

## Idée 2 : Eliminer manuellement les liens n'ayant pas de sens avec edit():
Y.gab <- graph2nb(gabrielneigh(as.matrix(xy), nnmult = 5), sym = TRUE)
plot(raster(tr))
plot(Y.gab, xy, add = T)
# new.nb <- edit.nb(Y.gab, coords = xy) ### Ne marche pas dans RStudio
load("new.nb.RData")
f1 <- function (D, dmax)     { 1 - (D/dmax) }
dmax <- max(dist(xy))
w.newnb.binary <- nb2listw(new.nb, style = "B")
w.newnb.flin   <- nb2listw(new.nb, style = "B", glist = lapply(nbdists(new.nb, as.matrix(xy)),
                                                               f1, dmax = dmax))
### summary
### xy = coordinates for nsites
### species = abundances of 30 species for nsites

# Optimisation du choix de la matrice W:
# **************************************
source("listw.candidates.R")
source("MEM.modsel.R")

library(vegan)
Y <- decostand(species, "hellinger")

# Create a few contrasting W matrices:
candidates <- listw.candidates(xy, del = F, gab = T, rel = F, mst = F, DB = F, PCNM = T, 
                               bin = T, flin = T, fconcdown = F, fconcup = F)
length(candidates)
candidates[[length(candidates)+1]] <- w.newnb.binary
candidates[[length(candidates)+1]] <- w.newnb.flin
names(candidates)[c(length(candidates)-1, length(candidates))] <- c("barr_bin", "barr_lin")
names(candidates)

modsel <- MEM.modsel(Y, candidates, autocor = "positive")
modsel$all$R2adj

   # Visualiser les variables MEM d'une matrice W donnée (élément i de modsel$all$MEM.all):
i <- 4
s.value(xy, modsel$all$MEM.all[[i]][, 1:4], Sp = Sp1, symbol = 'circle', pgrid.draw = F, 
        pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)
   # Visualiser les axes RDA (Y, MEM) d'une des matrices W:
r2 <- RsquareAdj(rda(Y, modsel$all$MEM.all[[4]]))$adj.r.squared
fwd <- forward.sel(Y, modsel$all$MEM.all[[4]], adjR2thresh = r2)
sub.me <- modsel$all$MEM.all[[4]][, sort(fwd$order)]
rda <- rda(Y, sub.me)
sc <- scores(rda, choices = c(1:6), display = c("lc", scaling = 1))
s.value(xy, sc[, 1:6],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)


# Univariate analysis (species by species):
# *****************************************
modsel.list <- vector("list", ncol(species))
for (i in 1:ncol(species)) {
  modsel.list[[i]] <- MEM.modsel(species[, i], candidates, autocor = "positive")
}
signif <- c()
for (i in 1:length(modsel.list)) signif <- c(signif, 
                                             ifelse(class(modsel.list[[i]]) == "NULL", 0, 1))
which(signif == 1)
modsel.list[[29]]$all$R2adj   
# species n° 12, 13, and 18 --> Higher R2adj with 'barriers' than with classical listw
i <- 12
rda <- rda(species[, i], modsel.list[[i]]$best$MEM.select)
sc <- scores(rda, choices = 1, display = c("lc", scaling = 1))
library(adegraphics)
plot(raster(tr))
s.value(xy, sc, add.plot = TRUE)

s.value(xy, me[,1:9],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)
# Abundances of the first nine species in the sampled sites:
s.value(xy.complete, simu.pop[, c(12, 13, 19)], symbol = 'circle', ppoint.cex = 0.4 , 
        Sp = Sp1, pSp.col = c('grey', 'peru'))

### Construction d'une variable réponse sur base des MEM :
# ********************************************************
# ********************************************************
nspe <- 10
Y <- matrix(ncol = nspe, nrow = nrow(me))
for (i in 1:nspe) {
  set.seed(i)
  intensity <- sample(seq(0.5, 0.9, 0.05), 1)
  Y[, i] <- (intensity * scale(me[, i])) + ((1 - intensity) * scale(rnorm(nrow(me), 0, 1)))
}

s.value(xy, me[, 1:10],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)
s.value(xy, Y,  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)

# Create a few contrasting W matrices:
candidates <- listw.candidates(xy, del = F, gab = T, rel = F, mst = F, DB = F, PCNM = T, 
                               bin = T, flin = T, fconcdown = F, fconcup = F)
length(candidates)
candidates[[length(candidates)+1]] <- w.newnb.binary
candidates[[length(candidates)+1]] <- w.newnb.flin
names(candidates)[c(length(candidates)-1, length(candidates))] <- c("barr_bin", "barr_lin")
names(candidates)

modsel <- MEM.modsel(Y, candidates, autocor = "positive")
modsel$all$R2adj

rda <- rda(Y, modsel$best$MEM.select)
sc <- scores(rda, choices = c(1:6), display = c("lc", scaling = 1))
s.value(xy, sc[, 1:6],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)

i <- 4
s.value(xy, modsel$all$MEM.all[[i]][, 1:4], Sp = Sp1, symbol = 'circle', pgrid.draw = F, 
        pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)

r2 <- RsquareAdj(rda(Y, modsel$all$MEM.all[[4]]))$adj.r.squared
fwd <- forward.sel(Y, modsel$all$MEM.all[[4]], adjR2thresh = r2)
sub.me <- modsel$all$MEM.all[[4]][, sort(fwd$order)]
rda <- rda(Y, sub.me)
sc <- scores(rda, choices = c(1:6), display = c("lc", scaling = 1))
s.value(xy, sc[, 1:6],  Sp = Sp1, symbol = 'circle', pgrid.draw = F, pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)

#
rda.list <- lapply(modsel$all$MEM.all, function (x) rda(Y, x))
sc.list <- lapply(rda.list, function (x) scores(x, choices = 1, display = c("lc", scaling = 1)))
for (i in 1:length(sc.list)) {
  s.value(xy, sc.list[[i]], Sp = Sp1, symbol = 'circle', pgrid.draw = F, 
          pSp.col = c('grey', 'peru'), plegend.drawKey = FALSE)
}
