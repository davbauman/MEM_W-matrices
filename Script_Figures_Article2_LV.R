###################################
#### Figure Accuracy (deltaR²) ####
###################################

# Pour les moyennes (barplots):
data <- read.table("fig_power_SUB_means_styleB.txt", h = T, sep = "\t")
str(data)
data$Connectivity <- factor(data$Connectivity, levels = c("del", "gab", "rel", "mst", "db"))
data$Strength <- factor(data$Strength, levels = c("Strong", "Weak"))
data$Scale <- factor(data$Scale, levels = c("Broad", "Fine"))
data$Design <- factor(data$Design, levels = c("Clustered", "Random"))
# Pour les matrices A:
data_w <- read.table("fig_power_SUB_all_styleB.txt", h = T, sep = "\t")
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
g <- g + labs(x = "", y = "Accuracy (ΔR²sub)")
g <- g + theme(axis.title.x = element_blank())
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
                                        "darkorange2")) + theme(legend.position = "none"))

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
                                        "darkorange2")) + theme(legend.position = "none"))

# On construit le graphique séparément pour les deux niveaux du facteur Strength:
# *******************************************************************************
strength <- "Weak"    # Strong or Weak
datasub <- subset(data, Strength == strength)
data_wsub <- subset(data_w, Strength == strength)

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare
# les barres :
h <- ggplot(datasub, aes(x = Connectivity, y = dR2sub)) + facet_grid(Design~Scale)
h <- h + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                  fill = "gray80")
h <- h + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
#h <- h + labs(x = "", y = "")
h <- h + theme(axis.title = element_blank())
h <- h + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
h <- h + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2))
h <- h + geom_point(data = data_wsub, aes(Connectivity, dR2sub, 
                                          color = factor(Weighting),
                                          shape = factor(Weighting)), size = 1.5)
h <- h + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up", 
                                       "PCNM"), values = c(15:19))
(h <- h + scale_color_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", 
                                        "Concave-up", "PCNM"), 
                             values = c("blue1", "black", "firebrick3", "forestgreen", 
                                        "darkorange2")) + theme(legend.position = "none"))

#### Power ####
###############

x <- ggplot(datasub, aes(x = Connectivity, y = Power)) + facet_grid(Design~Scale)
x <- x + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                  fill = "gray80")
x <- x + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
x <- x + labs(x = "Connectivity matrix", y = "")
x <- x + theme(axis.title.y = element_blank())
x <- x + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
x <- x + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1))
x <- x + geom_point(data = data_wsub, aes(Connectivity, Power, 
                                          color = factor(Weighting),
                                          shape = factor(Weighting)), size = 1.5)
x <- x + scale_shape_manual("Weighting function:", 
                            labels = c("Binary", "Linear", "Concave-down", "Concave-up",
                                       "PCNM"),
                            values = c(15:19))
(x <- x + scale_color_manual("Weighting function:", 
                             labels = c("Binary", "Linear", "Concave-down", 
                                        "Concave-up", "PCNM"), 
                             values = c("blue1", "black", "firebrick3", "forestgreen", 
                                        "darkorange2")) + theme(legend.position = "none") )

# Visualisation simultannée de la puissance et de la précision :
# **************************************************************
library(gridExtra)
grid.arrange(g, h, w, x, nrow = 2)

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

data <- read.table("fig_power_Optim_SUB_styleB.txt", h = T, sep = "\t")
str(data)
data$Strength <- factor(data$Strength, levels = c("Strong", "Weak"))
data$Scale <- factor(data$Scale, levels = c("Broad", "Fine"))
data$Design <- factor(data$Design, levels = c("Clustered", "Random"))
data$W_mat <- factor(data$W_mat, levels = c("Optim", "Random"))

# On construit le graphique séparément pour les deux niveaux du facteur Strength:
# *******************************************************************************
strength <- "Strong"    # Strong or Weak
datasub <- subset(data, Strength == strength)

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare les barres :
q <- ggplot(datasub, aes(x = W_mat, y = dR2sub)) + facet_grid(Design~Scale)
q <- q + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80")
q <- q + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white"))
q <- q + labs(x = "Choice of the W matrix", y = "Accuracy (ΔR²sub)")
q <- q + theme(axis.title.x = element_blank())
q <- q + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"))
q <- q + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2))
(q <- q + geom_errorbar(data = datasub, aes(ymin = dR2sub - sd, ymax = dR2sub + sd), 
                        width = .2, position = position_dodge(0.05)))

#### Power ####
###############

p <- ggplot(datasub, aes(x = W_mat, y = Power)) + facet_grid(Design~Scale)
p <- p + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                   fill = "gray80")
p <- p + theme(panel.background = element_rect(fill = "white"), 
                panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                  element_line(colour = "white"))
p <- p + labs(x = "Choice of the W matrix", y = "Power")
p <- p + theme(axis.title = element_text(size = 10.5))
p <- p + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"),
                axis.line.y = element_line(color = "black", size = 0.5, 
                                           linetype = "solid"))
p <- p + geom_hline(yintercept = 0, linetype = 1) + 
    scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1))
(p <- p + geom_errorbar(data = datasub, aes(ymin = Power - sd, ymax = Power + sd), 
                        width = .2, position = position_dodge(0.05)))

# On construit le graphique séparément pour les deux niveaux du facteur Strength:
# *******************************************************************************
strength <- "Weak"    # Strong or Weak
datasub <- subset(data, Strength == strength)

# Dans aes, on précise ce qui sur x et y, et à partir de quel facteur on sépare les barres :
r <- ggplot(datasub, aes(x = W_mat, y = dR2sub)) + facet_grid(Design~Scale)
r <- r + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                  fill = "gray80")
r <- r + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
r <- r + labs(x = "", y = "")
r <- r + theme(axis.title = element_blank())
r <- r + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
r <- r + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(-0.5, 0.1, by = 0.1), 2))
(r <- r + geom_errorbar(data = datasub, aes(ymin = dR2sub - sd, ymax = dR2sub + sd), 
                        width = .2, position = position_dodge(0.05)))

#### Power ####
###############

o <- ggplot(datasub, aes(x = W_mat, y = Power)) + facet_grid(Design~Scale)
o <- o + geom_bar(stat = "identity", position = position_dodge(), color = "black",
                  fill = "gray80")
o <- o + theme(panel.background = element_rect(fill = "white"), 
               panel.grid.major = element_line(colour = "white"), panel.grid.minor = 
                 element_line(colour = "white"))
o <- o + labs(x = "Choice of the W matrix", y = "Power")
o <- o + theme(axis.title.y = element_blank())
o <- o + theme(axis.line.x = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"),
               axis.line.y = element_line(color = "black", size = 0.5, 
                                          linetype = "solid"))
o <- o + geom_hline(yintercept = 0, linetype = 1) + 
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.2), 1))
(o <- o + geom_errorbar(data = datasub, aes(ymin = Power - sd, ymax = Power + sd), 
                        width = .2, position = position_dodge(0.05)))

grid.arrange(q, r, p, o, nrow = 2)
