MEM.modsel <- function(x, candidates, autocor = c("positive", "negative", "all"), 
                       alpha_thresh = 0.05)
{
  
  library(vegan)       # Eliminer quand sera dans le package
  library(adespatial)  # Eliminer quand sera dans le package
  
  x <- as.data.frame(x)
  if (any(is.na(x)) | any(is.na(coord))) stop("NA entries in x or coord")
  if (nrow(x) != nrow(coord)) stop("different number of rows")
  
  # **********************************************************************************
  # The MEM.test function tests the significance of a W matrix while taking into
  # consideration the total number of W matrices tested in MEM.modsel. This total
  # number of tests is used to apply a correction to the p-value in order not to
  # inflate the type I error rate. If the tested W matrix is significant, a model
  # selection is performed using Blanchet et al.'s forward selection with two stopping 
  # criteria.
  MEM.test <- function (a = x, b, c = autocor, alpha = alpha_thresh)
  {
    pval <- anova.cca(rda(a, b), permutations = 10000)$Pr[1]
    if (c == "all") pval <- 1-(1-pval)^2   # Sidak correction for autocor = "all" 
    if (pval <= alpha) {  
      R2adj <- RsquareAdj(rda(a, b))$adj.r.squared
      class <- class(try(fsel <- forward.sel(a, b, adjR2thresh = R2adj, nperm = 999),
                         TRUE))
      if (class != "try-error") { 
        sign <- sort(fsel$order)
        MEM.select <- b[, c(sign)] 
        list(MEM.select = MEM.select, NbVar = length(sign), pval = pval, 
             R2adj = RsquareAdj(rda(a, MEM.select))$adj.r.squared, 
             AdjR2Cum = fsel$AdjR2Cum)
      } 
    } else return(NA)
  }                                                  # End of the MEM.test() function
  # **********************************************************************************
  
  # **********************************************************************************
  # Function aiming at choosing the weighting function parameter that maximises the 
  # R2adj. The function returns 1) the complete W matrix corresponding to the best 
  # value and 2) the index of the selected parameter value (e.g., 'y' varies from
  # 5 to 9 and y = 6 is selected, then the function returns y = 2, that is, the index).
  chooseparam <- function (x = x, lw) 
  {
    listw <- vector("list", length(lw))
    listR2 <- vector("numeric", length(lw))
    for (i in 1:length(listw)) {
      listw[[i]] <- scores.listw(lw[[i]], MEM.autocor = cor[h])
      listR2[i] <- RsquareAdj(rda(x, listw[[i]]))$adj.r.squared
    }
    best <- which.max(listR2)
    list(W = listw[[best]], y_index = best)
  }                                                # End of the chooseparam() function
  # **********************************************************************************
  
  autocor <- match.arg(autocor) 
  
  # Since the loop is entered only once if autocor = "positive" or "negative" 
  # but twice if autocor = "all" (once for the positive and once for the negative 
  # MEM models), we begin by setting up the loop.
  
  if (autocor != "all") { 
    k <- 1   # number of times the for loop will run
    if (autocor == "positive") { 
      cor <- c("positive", "negative") 
    }  else cor <- c("negative", "positive") 
  } else {
    k <- 2
    cor <- c("positive", "negative")
  }
  
  # We define the total nb of tests:
  nbtest <- length(candidates)
  
  for (h in 1:k) {
    
    # For model comparison and selection
    results <- as.data.frame(matrix(nrow = nbtest, ncol = 4))
    colnames(results) <- c("pval_sidak", "R2adj", "NbMEM", "y_index")     
    # List of the MEM.select matrices
    listMEM <- vector("list", nbtest)
    # and corresponding R2
    listR2 <- vector("list", nbtest)
    # List of 'MEM.test' results
    listtest <- vector("list", nbtest)
    # List of global W matrices
    listW <- vector("list", nbtest)
    # Vector of best parameter values (when severall compared):
    param <- rep("NA", nbtest)
    
    for (q in 1:nbtest) {
      # If length(candidates[[q]]) > 3, we have a list in the candidates list, 
      # meaning that severall W matrices were build on the basis of a single set of
      # connectivity and weighting matrices. This occurs when different parameters of
      # a weighting functions are tested. If this case, the length of candidates[[q]]
      # corresponds to the number of parameters tested. One parameter value has to be 
      # chosen based on the global R2adj: the corresponding W matrix is then tested 
      # using MEM.test(). 
      # If length == 3, we only have a listw object in candidates[[q]] and the length
      # of 3 corresponds to "style", "neighbours", and "weights". We only have one 
      # weighted list that can be directly tested.
      if (length(candidates[[q]]) == 3) {
        W <- scores.listw(candidates[[q]], MEM.autocor = cor[h])
        listW[[q]] <- W
        listtest[[q]] <- MEM.test(x, W)
      } else {
        bestparam <- chooseparam(x = x, lw = candidates[[q]])
        listW[[q]] <- bestparam$W
        param[q] <- bestparam$y_index
        listtest[[q]] <- MEM.test(x, W)
      }
    }
    # Save the results in order to compare them and choose the best model:
    for (i in 1:nbtest) {
      if (is.list(listtest[[i]]) == TRUE) {
        results[i, 1] <- listtest[[i]]$pval
        results[i, 2] <- listtest[[i]]$R2adj
        results[i, 3] <- listtest[[i]]$NbVar
        if (param[i] != "NA") results[i, 4] <- param[i]
        listMEM[[i]] <- listtest[[i]]$MEM.select
        listR2[[i]] <- listtest[[i]]$AdjR2Cum
      }
    }
    # Selection of the best W matrix (and best model within it):
    if (length(which(results[, 1] <= alpha_thresh)) > 0) {
      best <- which.max(results[, 2])
      lenlist <- c()   # Will help with the result output  
      if (autocor != "all") {
        lenlist <- c(lenlist, cor[h])
        L <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                  MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                  param_index = results[best, 4], pval = results[best, 1], 
                  R2adj = results[best, 2], NbVar = results[best, 3], 
                  bestw_index = best)
      } else {
        lenlist <- c(lenlist, cor[h])
        if (h == 1) {
          L1 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                     param_index = results[best, 4], pval = results[best, 1], 
                     R2adj = results[best, 2], NbVar = results[best, 3], 
                     bestw_index = best)
        } else {
          L2 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                     MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                     param_index = results[best, 4], pval = results[best, 1], 
                     R2adj = results[best, 2], NbVar = results[best, 3], 
                     bestw_index = best)
        }
      }
    }
  }    # End of the for loop
  
  # Output of the MEM.modsel function:
  if (length(lenlist) == 2) {
    cat("\n", "\n", "*****************************************************", "\n",
        "*****************************************************", "\n",
        "A best positive (corrected p-value = ", round(L1$pval, 5), ", R2adj of the", 
        "\n", "selected MEM variables = ", round(L1$R2adj, 3), 
        ") and best negative (corrected p-value = ", round(L2$pval, 5), ",", "\n", 
        "R2adj of the selected MEM variables = ", round(L2$R2adj, 3), ")", 
        "MEM models were selected.", "\n", 
        "The corresponding spatial weighting W matrices are ", 
        L1$name, " (parameter_index = ", L1$param_index, ")", "\n", " and ", 
        L2$name, " (parameter_index = ", L2$param_index, ")", ", respectively.", "\n",
        "*****************************************************", "\n",
        "*****************************************************", "\n",
        "The output of the function is a list of two lists (MEM.pos and MEM.neg).",
        sep = "")
    list(MEM.pos = L1, MEM.neg = L2)
  } else 
    if (length(lenlist) == 1) { 
      if (autocor == "all") if (lenlist == "positive") L <- L1 else L <- L2
      cat("\n", "\n", "*****************************************************", 
          "\n", "*****************************************************", "\n",
          "A best ", lenlist, " MEM model was selected (corrected p-value = ", 
          round(L$pval, 5), ", R2adj of the selected", "\n",  "MEM variables = ", 
          round(L$R2adj, 3), ").", "\n", 
          "The corresponding spatial weighting W matrix is ", 
          L$name, "\n", "(parameter_index = ", L$param_index, ").", "\n",
          "*****************************************************", "\n",
          "*****************************************************", "\n", sep = "")
      list(MEM.all = L$MEM.all, MEM.select = L$MEM.select, 
           MEM.AdjR2Cum = L$MEM.AdjR2Cum, name = L$name, param_index = L$param_index,
           pval = L$pval, R2adj = L$R2adj, NbVar = L$NbVar, bestw_index = L$bestw_index)
    } else cat("\n", "\n", "*****************************************************", 
               "\n", "*****************************************************", "\n",
               "No significant spatial structure was detected in the data.", "\n",
               "\n", sep = "")   
}

# Written by:
# David Bauman, in Bauman et al. (2017)
