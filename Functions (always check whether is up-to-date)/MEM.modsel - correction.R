# Function MEM.modsel:
# ********************
# The function requires the packages vegan, adespatial and spdep to be installed. 
# The MEM.modsel function will load them.
# INPUT:
# ******
# Arguments:
# - x is the vector or dataframe of response variable(s).
# - coord is the matrix or dataframe of X and Y spatial coordinates. It can also be a 
# vector in the case of a transect.
# - autocor corresponds to the type of spatial correlation modelled by the MEM 
# eigenfunctions to be considered: autocor is either "positive" (default argument), 
# "negative" or "all" --> positively, negatively or both types of MEM eigenfunctions
# are considered. In the latter case, the positive and negative models are built, 
# tested and selected separately (a global adjusted R2 cannot be obtained and therefore 
# cannot but used as second stopping criterion of the forward selection for saturated 
# models (n-1 variables here, if positive and negative MEM eigenfunctions considered 
# together).


# - style: A COMPLETER

# - PCNM: whether a distance-based MEM based on the PCNM criterion (Dray et al. 2006) 
# should be computed and tested.
# - del, gab, rel, mst: Connectivity matrices to be tested (Delaunay triangulation, 
# Gabriel's graph, relative neighbourhood graph, and minimum spanning tree, 
# respectively). The default arguments test all connectivity matrices.
# - weightfun: Whether weighting matrices should be added to the connectivity matrices 
# (del, gab, rel, mst). If weightfun = TRUE, three weighting functions can be tested: 
# a decreasing linear (flin = TRUE), a concave-down (fconcdown = TRUE), and a 
# concave-up function (fconcup = TRUE). 
# - ymax: The two latter functions are tested by default for a 'y' exponent parameter 
# taking the value 2 to 10 and 1 to 10, respectively. The maximum value of the 'y' 
# exponent can be set to a different value through the 'ymax' argument.
# If weightfun = FALSE, only the binary MEM based on the selected connectivity matrices
# are tested.
# - alpha_thresh: Significance threshold value.
# OUTPUT:
# *******
# If only a positive or negative spatial model is tested, or if both are tested but 
# only one was significant, the function returns a list containing all characteristics 
# of the best spatial model selected: the connectivity and weighting matrices used for 
# constructing the spatial weighting matrix (W), the 'y' exponent parameter (if a 
# concave-down or -up weighting function was used), the p-value of the global model 
# (corrected by a Sidak correction depending on the number of W matrices tested), 
# the adjusted R2 of the final model (after performing the forward selection with 
# double stopping criterion), the R2 of each MEM variable and the eigenvectors (i.e., 
# the selected MEM variables). The latter can directly be used as spatial explanatory 
# variables in statistical models.
# If both positive and negative spatial models are tested and if both are significant, 
# the function returns a list of two lists (MEM.pos and MEM.neg), each one containing 
# the same information as described above. 

MEM.modsel <- function(x, candidates, autocor = c("positive", "negative", "all"), 
                       alpha_thresh = 0.05)
  
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
   MEM.test <- function (a = x, b, c = autocor, d = nbtest, alpha = alpha_thresh)
   {
      pval <- anova.cca(rda(a, b), permutations = 10000)$Pr[1]
      pval <- 1-(1-pval)^d                   # Sidak correction (nb of MEM model tested) 
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

   # A multitest p-value correction is needed for controling the type-I error rate. 
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
                    R2adj = results[best, 2], NbVar = results[best, 3])
        } else {
          lenlist <- c(lenlist, cor[h])
          if (h == 1) {
            L1 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                       MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                       param_index = results[best, 4], pval = results[best, 1], 
                       R2adj = results[best, 2], NbVar = results[best, 3])
          } else {
            L2 <- list(MEM.all = listW[[best]], MEM.select = listMEM[[best]], 
                       MEM.AdjR2Cum = listR2[[best]], name = names(candidates)[best], 
                       param_index = results[best, 4], pval = results[best, 1], 
                       R2adj = results[best, 2], NbVar = results[best, 3])
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
            pval = L$pval, R2adj = L$R2adj, NbVar = L$NbVar)
     } else cat("\n", "\n", "*****************************************************", 
                "\n", "*****************************************************", "\n",
                "No significant spatial structure was detected in the data.", "\n",
                "\n", sep = "")   
}

# Written by:
# David Bauman, in Bauman et al. (2017)
