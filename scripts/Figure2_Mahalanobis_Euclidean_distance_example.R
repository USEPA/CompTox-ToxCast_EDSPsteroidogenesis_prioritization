#-----------------------------------------------------------------------------------#
#Figure 1. Illustrating the difference between Mahalanabos and Euclidean distance

# This figure is to show the difference in Mahalanobis and Euclidean distance metrics

rm(list = ls())

library(ggplot2)
library(ggthemes)
library(ellipse)
library(grid)
library(gridExtra)
library(mvtnorm)
library(ggpubr)
library(data.table)

#setwd("./Misc")

Cor <- matrix(c(1,.8,.8,1),nrow=2, ncol=2)
Cor0 <- matrix(c(1,0,0,1),nrow=2, ncol=2)
Scale <- c(0.08, 0.16)
Scale0 <- c(0.16,0.16)
Cov <- diag(Scale) %*% Cor %*% diag(Scale)
Cov0 <- diag(Scale) %*% Cor0 %*% diag(Scale)
iCov  <- solve(Cov)
MHD <- function(x1,x2,iCov) {
    (x2 - x1) %*% iCov %*% (x2 - x1)
}

ED <- function(x1,x2) {
    (x2 - x1) %*% (x2 - x1)
}

## x, y
z0 <- c(0,0)
z1 <- c(0.1, 0.2)
z2 <- c(0.2236068, 0.0) #z2 has same distance from z0 as z1

##generate sampled data
set.seed(125)
sample_data <- as.data.table(rmvnorm(mean = c(0, 0), sigma = Cov, n = 25))

pts1 <- as.data.frame(rbind(z0,z1,z2))
#pts1 <- sweep(pts1,2,c(-3,-3),"+")
colnames(pts1) <- c("x","y")
pts1T <- pts1
pts1T$labels=c("conc 1","conc 2","conc 3")
TL <- as.data.frame(ellipse(Cor, scale=Scale, centre=unlist(pts1[1,])))
p1 <- ggplot() +
  #geom_point(data=pts1T, aes(x=x,y=y, shape = factor(labels))) +
  geom_text(data=pts1T[2,], aes(x=x,y=y,label=labels), nudge_y = 0.015) +
  geom_text(data=pts1T[3,], aes(x=x,y=y,label=labels), nudge_x = 0.05) +
  geom_text(data=pts1T[1,], aes(x=x,y=y,label=labels), nudge_y = -0.01) +
  geom_path(data=TL, aes(x=x,y=y), linetype = "dashed") +
  geom_segment(aes(x=pts1$x[1], xend=pts1$x[2], y=pts1$y[1], yend=pts1$y[2]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_segment(aes(x=pts1$x[1], xend=pts1$x[3], y=pts1$y[1], yend=pts1$y[3]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_point(data = sample_data, aes(x = V1, y = V2), alpha = 0.2) +
  scale_x_continuous("Measured Hormone X (log uM)", limits = c(-0.4, 0.4)) +
  scale_y_continuous("Measured Hormone Y (log uM)", limits = c(-0.4, 0.4)) +
  theme_few() +
  #theme(aspect.ratio = 1)
  coord_fixed() +
  theme(plot.margin = margin(5,5,5,5, "mm"))


zz <- eigen(iCov)
tr0 <- t(zz$vectors %*% diag(sqrt(zz$values)))

pts2 <- as.data.frame(t(tr0 %*% t(data.matrix(pts1))))
sample_data2 <- as.data.frame(t(tr0 %*% t(data.matrix(sample_data))))
colnames(pts2) <- c("x","y")
pts2T <- pts2
pts2T$labels=c("conc 1","conc 2","conc 3")
TL2 <- as.data.frame(t(tr0 %*% t(data.matrix(TL))))
colnames(TL2) <- c("x","y")
p2 <- ggplot() +
  #geom_point(data=pts2T, aes(x=x,y=y, shape = factor(labels))) +
  geom_text(data=pts2T[2,], aes(x=x,y=y,label=labels), nudge_y = -0.1) +
  geom_text(data=pts2T[3,], aes(x=x,y=y,label=labels), nudge_x = -0.45) +
  geom_text(data=pts2T[1,], aes(x=x,y=y,label=labels), nudge_y = 0.2) +
  geom_path(data=TL2, aes(x=x,y=y), linetype = "dashed") +
  geom_segment(aes(x=pts2$x[1], xend=pts2$x[2], y=pts2$y[1], yend=pts2$y[2]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_segment(aes(x=pts2$x[1], xend=pts2$x[3], y=pts2$y[1], yend=pts2$y[3]),
               arrow=arrow(length=unit(0.125,"inches")), size = 0.6) +
  geom_point(data = sample_data2, aes(x = V1, y = V2), alpha = 0.2) +
  scale_x_continuous("Rotated and Scaled Axis 1", limits = c(-5.25, 3)) +
  scale_y_continuous("Rotated and Scaled Axis 2", limits = c(-5.25, 3)) +
  theme_few() +
  #theme(aspect.ratio = 1)
  coord_fixed() +
  theme(plot.margin = margin(5,5,5,5, "mm"))


pdf("./figures/Figure2_mahalanobis_example.pdf", width = 12, height =8)
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, align = "h", hjust = -1, vjust = 7, widths = c(1, 1), heights = c(1, 1))
dev.off()

png("./figures/Figure2_mahalanobis_example.png", res = 300, units = "in", width = 12, height =8)
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, align = "h", hjust = -1, vjust = 7, widths = c(1, 1), heights = c(1, 1))
dev.off()
