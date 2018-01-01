set.seed(5)

design <- "clustered"   # "random"

if (design == "clustered") {
  zones <- matrix(c(1:9, rep(c(0, 30, 60), times = 3), rep(c(0, 30, 60), each = 3)), ncol = 3)
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

cclus <- C
cran  <- C

par(mfrow = c(1, 2), mar = c(2, 3, 0, 0))
plot(cclus, xlab = "", ylab = "")
plot(cran, xlab = "", ylab = "")

##################################################################
# Fig. 2:
# Illustration of the connection schemes:
# ***************************************

nsam <- 25
set.seed(1)
grid <- expand.grid(seq(1, 20), seq(1, 20))
subgrid <- grid[sample(c(1:nrow(grid)), nsam), ]
plot(subgrid)

del <- tri2nb(subgrid)
gab <- graph2nb(gabrielneigh(as.matrix(subgrid), nnmult = 5), sym = TRUE)
rel <- graph2nb(relativeneigh(as.matrix(subgrid), nnmult = 5), sym = TRUE)
mst <- mst.nb(dist(subgrid))
db  <- dnearneigh(give.thresh(dist(subgrid)), x = as.matrix(subgrid), d1 = 0)
par(mfrow = c(2, 3), mar = c(0, 0, 0, 2))
plot(db, subgrid)
plot(del, subgrid)
plot(gab, subgrid)
plot(rel, subgrid)
plot(mst, subgrid)

###########################################################
# Fig. annexe présentant les fonctions de pondérations :
# ******************************************************

f1 <- function (D, dmax)     { 1 - (D/dmax) }       # Linear function
f2 <- function (D, dmax, y)  { 1 - (D/dmax)^y }     # Concave-down function
f3 <- function (D, y)        { 1 / D^y }            # Concave-up function

dist <- c(1:100)
dmax <- max(dist)
ycdown <- c(1:10)
ycup <- seq(0.1, 1, 0.1)

f2d <- vector("list", length(ycdown))
for (i in 1:length(ycdown)) {
  f2d[[i]] <- sapply(dist, function(x) f2(x, 100, ycdown[i]))
}

f3d <- vector("list", length(ycup))
for (i in 1:length(ycup)) {
  f3d[[i]] <- sapply(dist, function(x) f3(x, ycup[i]))
}

par(mar = c(2, 2, 1, 2))
plot(dist, f2d[[1]], type = "l", col = "red", lwd = 3)
for (i in 2:length(ycdown)) {
  if (i != 5) {
    lines(dist, f2d[[i]], col = "blue")
  } else lines(dist, f2d[[i]], col = "red", lwd = 3)
}
plot(dist, f2d[[1]], type = "n")
for (i in ycup) {
  k <- i * 10
  if (i != 0.5) {
    lines(dist, f3d[[k]], col = "blue")
  } else lines(dist, f3d[[k]], col = "red", lwd = 3)
}
