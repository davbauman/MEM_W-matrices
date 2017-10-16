###################################
#### Figure Type I error rates ####
###################################

#### Load ggplot2 ####
library(ggplot2)

#### Input of the data ####
# Pour les points :
chic <- read.table("fig_typeIerror_all_(details).txt", h = T, sep = "\t")
str(chic)
chic$Connectivity <- factor(chic$Connectivity, levels = c("del", "gab", "rel", "mst", "db"))
#chic$Connectivity <- factor(chic$Connectivity, levels = c("graph-based", "distance-based"))
chic$Weighting <- factor(chic$Weighting, levels = c("Binary", "Linear", "Concave-down", 
                                                    "Concave-up", "PCNM"))
chic$Design <- factor(chic$Design, levels = c("Clustered", "Random"))

# Pour les barplots :
chic2 <- read.table("fig_typeIerror_means_(details).txt", h = T, sep = "\t")
chic2$Connectivity <- factor(chic2$Connectivity, levels = c("del", "gab", "rel", "mst", "db"))
#chic2$Connectivity <- factor(chic2$Connectivity, levels = c("graph-based", "distance-based"))
chic2$Design <- factor(chic2$Design, levels = c("Clustered", "Random"))

# Barplots et points superposés :
# *******************************

b <- ggplot(chic2, aes(Connectivity, Mean))
b <- b + geom_bar(stat = "identity", position = "stack", color = "black", 
                  fill = "gray80") + facet_wrap(~Design, ncol =  2) 
b <- b + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                   element_line(colour = "white"))
b <- b + geom_point(data = chic, aes(Connectivity, typeIerror, 
                                     color = factor(Weighting),
                                     shape = factor(Weighting)), size = 1.5)  
b <- b + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up", 
                                       "PCNM"), values = c(15:19))
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
                            values = c("blue1", "black", "firebrick3", "forestgreen", 
                                       "darkorange2")) 
b <- b + geom_hline(yintercept = 0.05, linetype = 2) + 
    scale_y_continuous(breaks = seq(0.01, 0.06, by = 0.01))
b

###################################
#### Figure Accuracy (deltaR²) ####
###################################

# Pour les moyennes (barplots):
data <- read.table("fig_power_POP_means_styleB.txt", h = T, sep = "\t")
str(data)
data$Connectivity <- factor(data$Connectivity, levels = c("del", "gab", "rel", "mst", "db"))
data$Strength <- factor(data$Strength, levels = c("Strong", "Weak"))
data$Scale <- factor(data$Scale, levels = c("Broad", "Fine"))
data$Design <- factor(data$Design, levels = c("Clustered", "Random"))
# Pour les matrices A:
data_w <- read.table("fig_power_POP_all_styleB.txt", h = T, sep = "\t")
str(data_w)
data_w$Connectivity <- factor(data_w$Connectivity, levels = c("del", "gab", "rel", "mst", "db"))
data_w$Design <- factor(data_w$Design, levels = c("Clustered", "Random"))
data_w$Scale <- factor(data_w$Scale, levels = c("Broad", "Fine"))
data_w$Weighting <- factor(data_w$Weighting, levels = c("Binary", "Linear",
                                                        "Concave-down", "Concave-up", "PCNM"))

# On construit le graphique séparément pour les deux niveaux du facteur Strength:
# *******************************************************************************
strength <- "Strong"    # Strong or Weak
datasub <- subset(data, Strength == strength)
data_wsub <- subset(data_w, Strength == strength)

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare
# les barres :
g <- ggplot(datasub, aes(x = Connectivity, y = dR2sub)) + facet_grid(Design~Scale)
g <- g + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80")
g <- g + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
g <- g + labs(x = "Connectivity matrix", y = "Accuracy (ΔR²sub)")
g <- g + theme(axis.title = element_text(size = 10.5))
g <- g + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
g <- g + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2))
g <- g + geom_point(data = data_wsub, aes(Connectivity, dR2sub, 
                                  color = factor(Weighting),
                                  shape = factor(Weighting)), size = 1.5)
g <- g + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up", 
                                       "PCNM"), values = c(15:19))
(g <- g + scale_color_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", 
                                       "Concave-up", "PCNM"), 
                            values = c("blue1", "black", "firebrick3", "forestgreen", 
                                       "darkorange2")))

#### Power ####
###############

w <- ggplot(datasub, aes(x = Connectivity, y = Power)) + facet_grid(Design~Scale)
w <- w + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80")
w <- w + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white"))
w <- w + labs(x = "Connectivity matrix", y = "Power")
w <- w + theme(axis.title = element_text(size = 10.5))
w <- w + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"))
w <- w + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1))
w <- w + geom_point(data = data_wsub, aes(Connectivity, Power, 
                                        color = factor(Weighting),
                                        shape = factor(Weighting)), size = 1.5)
w <- w + scale_shape_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", "Concave-up",
                                        "PCNM"),
                             values = c(15:19))
(w <- w + scale_color_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", 
                                        "Concave-up", "PCNM"), 
                             values = c("blue1", "black", "firebrick3", "forestgreen", 
                                        "darkorange2")))

# Visualisation simultannée de la puissance et de la précision :
# **************************************************************
library(gridExtra)
grid.arrange(b, w, g, nrow = 3)

####################################
### Figure TypeIerror MEM.modsel ###
####################################

typeIopt <- read.table("fig_typeIerror_Opt.txt", h = T, sep = "\t")

k <- ggplot(typeIopt, aes(Correction, typeIerror))
k <- k + geom_bar(stat = "identity", position = "stack", color = "black", 
                  fill = "gray80") + facet_wrap(~Design, ncol =  2) 
k <- k + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
k <- k + labs(x = "global p-value correction", y = "Type I error rate")
k <- k + theme(axis.title = element_text(size = 10.5))
k <- k + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
k <- k + expand_limits(y = 0)
k <- k + geom_hline(yintercept = 0.05, linetype = 2) + 
  scale_y_continuous(breaks = seq(0.05, 0.35, by = 0.05))
k

#########################################
### Power and Accuracy - Optimisation ###
#########################################

data <- read.table("fig_power_Optim_POP_styleB.txt", h = T, sep = "\t")
str(data)
data$Strength <- factor(data$Strength, levels = c("Strong", "Weak"))
data$Scale <- factor(data$Scale, levels = c("Broad", "Fine"))
data$Design <- factor(data$Design, levels = c("Clustered", "Random"))
data$W_mat <- factor(data$W_mat, levels = c("Optimisation", "Random"))

# On construit le graphique séparément pour les deux niveaux du facteur Strength:
# *******************************************************************************
strength <- "Strong"    # Strong or Weak
datasub <- subset(data, Strength == strength)

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare les barres :
(q <- ggplot(datasub, aes(x = W_mat, y = dR2sub)) + facet_grid(Design~Scale))
(q <- q + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80"))
(q <- q + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white")))
(q <- q + labs(x = "Choice of the W matrix", y = "Accuracy (ΔR²sub)"))
(q <- q + theme(axis.title = element_text(size = 10.5)))
(q <- q + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid")))
(q <- q + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2)))
(q <- q + geom_errorbar(data = datasub, aes(ymin = dR2sub - sd, ymax = dR2sub + sd), 
                        width = .2, position = position_dodge(0.05)))

#### Power ####
###############

(p <- ggplot(datasub, aes(x = W_mat, y = Power)) 
 + facet_grid(Design~Scale))
(p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80"))
(p <- p + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white")))
(p <- p + labs(x = "Connectivity matrix", y = "Power"))
(p <- p + theme(axis.title = element_text(size = 10.5)))
(p <- p + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid")))
(p <- p + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1)))
(p <- p + geom_errorbar(data = datasub, aes(ymin = Power - sd, ymax = Power + sd), 
                        width = .2, position = position_dodge(0.05)))


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
