listf <- vector("list", length(y_exp))

for (i in 1:length(y_exp)) {
  listf[[i]] <- nb2listw(Y.DBMEM, style = style, 
                         glist = lapply(nbdists(Y.DBMEM, as.matrix(C)),
                                        f3, y = i, dmax = dmax))
}

# Représentation graphique des fonctions de pondération en fonction de la distance pour
# différentes valeurs d'exposants:
# *************************************************************************************
#c <- expand.grid(x = seq(1, 10000, 100), y = seq(1, 5000, 100))
#xy.d1 <- dist(c)   ; xy.d1big <- xy.d1
#lowlim <- give.thresh(xy.d1)
#dmaxbig <- max(xy.d1big)
#Y.DBMEMbig <- dnearneigh(x = as.matrix(c), d2 = lowlim, d1 = 0)
#Y.DBMEM.f3.list <- lapply(y_exp, function(x) test.W.R2(nb= Y.DBMEMbig, xy = c,
                                                       f = f3, y= x, dmax = dmaxbig,
                                                       MEM.autocor = MEM_model))

f2 <- function (D, dmax, y) {1 - (D/dmax)^y}      # Concave-down function
f3 <- function (D, dmax, y) {1 / (D/dmax)^y}      # Concave-up function

D <- c(1:round(dmax))
listf <- vector("list", length(y_exp))

f <- f2
####
# Si voulons visualiser ancienne f3:
      f <- function (D, y) {1 / (D)^y}

      for (j in 1:length(y_exp)) {
        vec <- c()
        for (i in 1:length(D)) {
          vec <- c(vec, f(D[i], y_exp[j]))
        }
        listf[[j]] <- vec
      }
####

for (j in 1:length(y_exp)) {
  vec <- c()
  for (i in 1:length(D)) {
    vec <- c(vec, f(D[i], dmax, y_exp[j]))
  }
  listf[[j]] <- vec
}

par(mar=c(4, 4, 1, 1))
plot(D, listf[[10]], type = "n")
for (i in 1:length(listf)) lines(D, listf[[i]])
for (i in 1:length(listf)) lines(D, listf[[i]], col = "blue")

# Min et Max pour chaque valeur de 'y':
for (i in 1:length(listf)) {
  print(min(listf[[i]]))
}
for (i in 1:length(listf)) {
  print(max(listf[[i]]))
}
