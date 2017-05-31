###################################
#### Figure Type I error rates ####
###################################

#### Load ggplot2 ####
library(ggplot2)

#### Input of the data ####
# Pour les points :
chic <- read.table("typeIerror.txt", h = T, sep = "\t")
str(chic)
chic$Connectivity <- factor(chic$Connectivity, levels = c("del", "gab", "rel", "mst", 
                                                          "DB"))
chic$Weighting <- factor(chic$Weighting, levels = c("bin", "Linear", "Concave-down", 
                                                    "Concave-up"))
chic$Design <- factor(chic$Design, levels = c("Regular", "Random"))
# Pour les barplots :
chic2 <- read.table("typeIerror_means.txt", h = T, sep = "\t")
chic2$Connectivity <- factor(chic2$Connectivity, levels = c("del", "gab", "rel", "mst", 
                                                          "DB"))
chic2$Design <- factor(chic2$Design, levels = c("Regular", "Random"))

# Barplots et points superposés :
# *******************************

b <- ggplot(chic2, aes(Connectivity, mean))
b <- b + geom_bar(stat = "identity", position = "stack", color = "black", 
                  fill = "gray80") + 
   facet_wrap(~Design, ncol =  2) 
b <- b + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                   element_line(colour = "white"))
b <- b + geom_point(data = chic, aes(Connectivity, TypeIerror, 
                                     color = factor(Weighting),
                                     shape = factor(Weighting)), size = 1.5)  
b <- b + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up"),
                            values = c(15:18))
b <- b + labs(x = "Connectivity matrix", y = "Type I error rate")
b <- b + theme(axis.title = element_text(size = 10.5))
b <- b + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
b <- b + expand_limits(y = 0)
b <- b + scale_color_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up", 
                                       "PCNM"), 
                            values = c("blue1", "black", "firebrick3", "forestgreen")) 
b <- b + geom_hline(yintercept = 0.05, linetype = 2) + 
    scale_y_continuous(breaks = seq(0.05, 0.35, by = 0.05))
b

###################################
#### Figure Accuracy (deltaR²) ####
###################################

# Pour les moyennes (barplots):
data <- read.table("power_accuracy.txt", h = T, sep = "\t")
str(data)
data$Connectivity <- factor(data$Connectivity, levels = c("del", "gab", "rel", "mst", 
                                                          "DB"))
data$Design <- factor(data$Design, levels = c("Regular", "Random"))
data$Scale <- factor(data$Scale, levels = c("Broad", "Medium", "Fine"))
# Pour les matrices A:
data_w <- read.table("power_accuracy_detail.txt", h = T, sep = "\t")
str(data_w)
data_w$Connectivity <- factor(data_w$Connectivity, levels = c("del", "gab", "rel", "mst", 
                                                              "DB"))
data_w$Design <- factor(data_w$Design, levels = c("Regular", "Random"))
data_w$Scale <- factor(data_w$Scale, levels = c("Broad", "Medium", "Fine"))
data_w$Weighting <- factor(data_w$Weighting, levels = c("Binary", "Linear",
                                                        "Concave-down", "Concave-up"))

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare
# les barres :
(g <- ggplot(data, aes(x = Connectivity, y = deltaR2)) 
  + facet_grid(Design~Scale))
(g <- g + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80"))
(g <- g + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white")))
(g <- g + labs(x = "Connectivity matrix", y = "Accuracy (ΔR²)"))
(g <- g + theme(axis.title = element_text(size = 10.5)))
(g <- g + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid")))
(g <- g + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2)))
(g <- g + geom_point(data = data_w, aes(Connectivity, deltaR2, 
                                  color = factor(Weighting),
                                  shape = factor(Weighting)), size = 1.5)) 
(g <- g + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up"),
                            values = c(15:18)))
(g <- g + scale_color_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", 
                                       "Concave-up"), 
                            values = c("blue1", "black", "firebrick3", "forestgreen")))

#### Power ####
###############

(w <- ggplot(data, aes(x = Connectivity, y = power)) 
 + facet_grid(Design~Scale))
(w <- w + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80"))
(w <- w + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white")))
(w <- w + labs(x = "Connectivity matrix", y = "Power"))
(w <- w + theme(axis.title = element_text(size = 10.5)))
(w <- w + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid")))
(w <- w + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1)))
(w <- w + geom_point(data = data_w, aes(Connectivity, power, 
                                        color = factor(Weighting),
                                        shape = factor(Weighting)), size = 1.5)) 
(w <- w + scale_shape_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", "Concave-up"),
                             values = c(15:18)))
(w <- w + scale_color_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", 
                                        "Concave-up"), 
                             values = c("blue1", "black", "firebrick3", "forestgreen")))

# Visualisation simultannée de la puissance et de la précision :
# **************************************************************
library(gridExtra)
grid.arrange(b, w, g, nrow = 3)

#############################
#############################
#### INFORMATION DIVERSE ####
#############################

#### A Default ggplot ####
(g <- ggplot(chic, aes(Connectivity, TypeIerror, color = factor(Weighting))))
(g <- g + geom_point())

#### Working with Axes ####
(g <- g + labs(x = "Connectivity matrix", y = "Type I error rate") + 
   theme(axis.title = element_text(size = 10.5)))
# To include 0 on the y-axis:
(g <- g + expand_limits(y=0))
# Pour avoir des axes x et y apparents (barre noire) :
(g <- g + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"), 
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid")))

#### Working with Legends ####

# Si on ne veut pas de légende :
#g + theme(legend.position = "none")
# Pas de titre à la légende :
# g + theme(legend.title = element_blank())
# Titre et redéfinir "labels" de la légende
(g <- g + scale_color_discrete("Weighting function:", labels=c("Binary", "Linear", 
                                                               "Concave-down", 
                                                               "Concave-up", "PCNM")))
#### Working with the background ####

# tout est dans theme(). theme(panel.background = element_rect(fill = "grey60")) permet
# de choisir la couleur du fond. panel.grid.major = element_line(colour = ) et 
# panel.grid.minor décident de la couleur et de la taille (si on ajoute ", size = " 
# dans element_line) des lignes de quadrillages principales et secondaires du graphique.

(g <- g + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white")))

#### Working with themes ####

library(ggthemes)
g + theme_tufte()   # Fond blanc, pas d'axes, ni lignes horizontales ou verticales

# - Utiliser theme_set() pour définir un thème et ne plus devoir le préciser.
# - Taper le nom d'un thème pour voir son code --> créer un nouveau thème en modifiant
# certains aspects et en enregistrant le résultat dans un nouveau thème (en créant une
# fonction), qu'on définit ensuite comme thème par défaut avec theme_set().
# - theme_Existant <- theme_update(...) --> theme_update permet de changer certains
# aspects d'un thème existant.
# - If we do not want to use our theme anymore: theme_set(theme_gray())

#### Working with Colors ####
## categorial variables

# Pour définir manuellement : g + scale_color_manual(values = c("dodgerblue4", ...))
# Pour utiliser palette existante :
g + scale_color_brewer(palette = "Set1")
# ou
library(ggthemes)
g + scale_color_tableau()  # entre les () : palette = "nom palette"
g + scale_color_tableau(palette = "colorblind10")
## continiuous variables
# utiliser + scale_color_continuous("nom de la variable continue")
# ou 
# + scale_color_gradient(low = "darkkhaki", high = "darkgreen", "Ozone:")
# ou encore
# scale_color_gradient2(midpoint = mid, low = "blue4", mid = "white", high = "red4", 
# "Ozone:")

#### Ajout de lignes de références au plot ####
# Utilisation de geom_abline(), geom_hline(), geom_vline()
(g <- g + geom_hline(yintercept = 0.05, linetype = 2) + 
   scale_y_continuous(breaks = seq(0.05, 0.35, by = 0.05)))
# Et changement des valeurs d'axes avec scale_x_continuous ou scale_y_continuous
