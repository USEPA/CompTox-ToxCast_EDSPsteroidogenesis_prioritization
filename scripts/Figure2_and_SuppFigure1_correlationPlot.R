#-----------------------------------------------------------------------------------#

#This script generates figure 2 (correlation plots) and supplemental figure 1

library(data.table)
library(boot)
library(stringr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)
library(MASS)
library(mvtnorm)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(corrplot)

rm(list = ls())

#setwd("/share/home/dhaggard/Desktop/Mahalanobis_Follow_Up") #change accordingly

#-----Load in original mahalanobis distances and the new simulated log_uMs
load("./RData/AllResps_outliersRemoved2018-10-02.RData")
load("./RData/simulated_Mahalanobis_distance_output_2018-10-03.RData")
load("./RData/subsample_Model_output2018-10-03.RData")
#load("./simulated_Mahalanobis_distance_output_originalCovariance_2018-05-10.RData")
#load("./Mahalanobis_distance_output_withSimulatedCovariances_2018-05-10.RData")

#-----Calculate Pooled covariance matrices for upper/lower subsamples of original data
Covs_original <- vector("list", length = 8)
names(Covs_original) <- names(Models_all_subsample)

Covs_original_upper <- vector("list", length = 8)
names(Covs_original_upper) <- names(Models_upper_subsample)

Covs_original_lower <- vector("list", length = 8)
names(Covs_original_lower) <- names(Models_lower_subsample)

#estimate block-level covariance
for(block in names(Models_all_subsample)){
  Covs_original[[block]] <- estVar(Models_all_subsample[[block]]$fit)
  Covs_original_upper[[block]] <- estVar(Models_upper_subsample[[block]]$fit)
  Covs_original_lower[[block]] <- estVar(Models_lower_subsample[[block]]$fit)
}

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(colnames(Covs_original[[1]]), rownames(Covs_original[[1]]), names(Models_all_subsample)))
for (blk in names(Models_all_subsample)) CovT[rownames(Covs_original[[blk]]),colnames(Covs_original[[blk]]),blk] <- Covs_original[[blk]]
CovTP_original <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(colnames(Covs_original[[1]]), rownames(Covs_original[[1]]), names(Models_all_subsample)))
for (blk in names(Models_upper_subsample)) CovT[rownames(Covs_original_upper[[blk]]),colnames(Covs_original_upper[[blk]]),blk] <- Covs_original_upper[[blk]]
CovTP_original_upper <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(colnames(Covs_original[[1]]), rownames(Covs_original[[1]]), names(Models_all_subsample)))
for (blk in names(Models_lower_subsample)) CovT[rownames(Covs_original_lower[[blk]]),colnames(Covs_original_lower[[blk]]),blk] <- Covs_original_lower[[blk]]
CovTP_original_lower <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

colnames(CovTP_original) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                              "E1", "E2")
rownames(CovTP_original) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                              "E1", "E2")

colnames(CovTP_original_upper) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                                    "E1", "E2")
rownames(CovTP_original_upper) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                                    "E1", "E2")

colnames(CovTP_original_lower) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                                    "E1", "E2")
rownames(CovTP_original_lower) <- c("CORT", "CORTIC", "11-DCORT", "ANDR", "DOC", "TESTO", "OH-PROG", "OH-PREG", "PROG",
                                    "E1", "E2")

#-----Generate correlation plots using ggplot2
#correlation matrices
cor_original <- round(cov2cor(CovTP_original), digits = 2)
cor_original_lower <- round(cov2cor(CovTP_original_lower), digits = 2)
cor_original_upper <- round(cov2cor(CovTP_original_upper), digits = 2)

cor_diag <- cor_original #to label diagonal

#cluster based on original and reorder
cor_original_order <- hclust(dist((1-cor_original)/2), method = "ward.D")

cor_original <- cor_original[cor_original_order$order, cor_original_order$order]
cor_original_lower <- cor_original_lower[cor_original_order$order, cor_original_order$order]
cor_original_upper <- cor_original_upper[cor_original_order$order, cor_original_order$order]

#remove upper and diagonal
cor_original[!lower.tri(cor_original)] <- NA
cor_original_lower[!lower.tri(cor_original_lower)] <- NA
cor_original_upper[!lower.tri(cor_original_upper)] <- NA

cor_diag[lower.tri(cor_diag)] <- NA
cor_diag[upper.tri(cor_diag)] <- NA

#melt matrices
cor_original_long <- data.table(melt(cor_original))
cor_original_lower_long <- data.table(melt(cor_original_lower))
cor_original_upper_long <- data.table(melt(cor_original_upper))

cor_diag_long <- data.table(melt(cor_diag))
cor_diag_long[value == 1, diag := as.character(Var1)]

#get diagonals of covariance matrices to show differences in variance
cov_original <- round(CovTP_original, 3)
cov_original_lower <- round(CovTP_original_lower, 3)
cov_original_upper <- round(CovTP_original_upper, 3)

cov_original[lower.tri(cov_original)] <- NA
cov_original[upper.tri(cov_original)] <- NA
cov_original_lower[lower.tri(cov_original_lower)] <- NA
cov_original_lower[upper.tri(cov_original_lower)] <- NA
cov_original_upper[lower.tri(cov_original_upper)] <- NA
cov_original_upper[upper.tri(cov_original_upper)] <- NA

#melt matrices
cov_original_long <- data.table(melt(cov_original))
cov_original_lower_long <- data.table(melt(cov_original_lower))
cov_original_upper_long <- data.table(melt(cov_original_upper))

#generate plot using clustered correlations of original and add annotations for original lower and upper

p <- ggplot(data = cor_original_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(data = cor_diag_long, aes(x = Var1, y = Var2, label = diag), size = 5, fontface = "bold", hjust = "right", nudge_x = 0.3, nudge_y = 0.55) +
  scale_fill_gradient2(low = "#ca0020", high = "#0571b0", mid = "white", midpoint = 0, limit = c(-1, 1),
                       name = "Correlation", na.value = "white") +
  theme_minimal() +
  geom_text(aes(x = Var1, y = Var2, label = value), nudge_y = 0, fontface = "bold") +
  geom_text(data = cor_original_lower_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0.25) + 
  geom_text(data = cor_original_upper_long, aes(x = Var1, y = Var2, label = value), nudge_y = -0.25, fontface = "italic") +
  geom_text(data = cov_original_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0, fontface = "bold") +
  geom_text(data = cov_original_lower_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0.25) +
  geom_text(data = cov_original_upper_long, aes(x = Var1, y = Var2, label = value), nudge_y = -0.25, fontface = "italic") +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.85),
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_colorbar(barwidth = 18, barheight = 1.5, title.position = "top", title.hjust = 0.5))

pdf("./supplement/SupplementalFigure1_merged_correlation_plot_originalCovs.pdf", width = 8, height = 8)
p
dev.off()

png("./supplement/SupplementalFigure1_merged_correlation_plot_originalCovs.png", width = 8, height = 8, units = "in", res = 300)
p
dev.off()

#-----Analysis using simulated covariance matrices
#-----Separate the two nested tables and merge lists into two data tables 
simulated_CovTP <- vector("list", length = 2000)
simulated_CovTP_upper <- vector("list", length = 2000)
simulated_CovTP_lower <- vector("list", length = 2000)

for(x in 1:2000){

  simulated_CovTP[[x]] <- simulated_mMd_output[[x]]$CovTP
  simulated_CovTP_upper[[x]] <- simulated_mMd_output[[x]]$CovTP_upper
  simulated_CovTP_lower[[x]] <- simulated_mMd_output[[x]]$CovTP_lower

}

#-----Generate corrplots or the different covariance types to compare them, use median or covariance elements

CovTP_original <- array(NA, dim=c(11, 11, 2000), dimnames=list(colnames(simulated_CovTP[[1]]), rownames(simulated_CovTP[[1]]), 1:2000))
for(i in 1:2000) CovTP_original[rownames(simulated_CovTP[[i]]),colnames(simulated_CovTP[[i]]),i] <- simulated_CovTP[[i]]
CovTP_original <- apply(CovTP_original, c(1, 2), mean)

CovTP_lower <- array(NA, dim=c(11, 11, 2000), dimnames=list(colnames(simulated_CovTP_lower[[1]]), rownames(simulated_CovTP_lower[[1]]), 1:2000))
for(i in 1:2000) CovTP_lower[rownames(simulated_CovTP_lower[[i]]),colnames(simulated_CovTP_lower[[i]]),i] <- simulated_CovTP_lower[[i]]
CovTP_lower <- apply(CovTP_lower, c(1, 2), mean)

CovTP_upper <- array(NA, dim=c(11, 11, 2000), dimnames=list(colnames(simulated_CovTP_upper[[1]]), rownames(simulated_CovTP_upper[[1]]), 1:2000))
for(i in 1:2000) CovTP_upper[rownames(simulated_CovTP_upper[[i]]),colnames(simulated_CovTP_upper[[i]]),i] <- simulated_CovTP_upper[[i]]
CovTP_upper <- apply(CovTP_upper, c(1, 2), mean)

colnames(CovTP_original) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                              "PROG", "TESTO")
rownames(CovTP_original) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                              "PROG", "TESTO")

colnames(CovTP_lower) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                           "PROG", "TESTO")
rownames(CovTP_lower) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                           "PROG", "TESTO")

colnames(CovTP_upper) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                           "PROG", "TESTO")
rownames(CovTP_upper) <- c("11-DCORT", "ANDR", "CORTIC", "CORT", "DOC", "E2", "E1", "OH-PREG", "OH-PROG",
                           "PROG", "TESTO")
#generate corrplots
# png("./figures/corr_cov_original.png", width = 16, height = 12, units = "in", res = 300)
# corrplot(corr = cov2cor(CovTP_original), method = "color", tl.col = "black", tl.pos = "d", diag = FALSE,
#          cl.align.text = "l", cl.cex = 1, tl.cex = 1, number.cex = 0.9, addgrid.col = "grey", cl.pos = "b",
#          order = "hclust", addCoef.col = "black")
# dev.off()
# 
# png("./figures/corr_cov_lower.png", width = 16, height = 12, units = "in", res = 300)
# corrplot(corr = cov2cor(CovTP_lower), method = "color", tl.col = "black", tl.pos = "d", diag = FALSE,
#          cl.align.text = "l", cl.cex = 1, tl.cex = 1, number.cex = 0.9, addgrid.col = "grey", cl.pos = "b",
#          order = "hclust", addCoef.col = "black")
# dev.off()
# 
# png("./figures/corr_cov_upper.png", width = 16, height = 12, units = "in", res = 300)
# corrplot(corr = cov2cor(CovTP_upper), method = "color", tl.col = "black", tl.pos = "d", diag = FALSE,
#          cl.align.text = "l", cl.cex = 1, tl.cex = 1, number.cex = 0.9, addgrid.col = "grey", cl.pos = "b",
#          order = "hclust", addCoef.col = "black")
# dev.off()

#-----Generate correlation plots using ggplot2
#correlation matrices
cor_original <- round(cov2cor(CovTP_original), digits = 2)
cor_lower <- round(cov2cor(CovTP_lower), digits = 2)
cor_upper <- round(cov2cor(CovTP_upper), digits = 2)

cor_diag <- cor_original #to label diagonal

#cluster based on original and reorder
cor_original_order <- hclust(dist((1-cor_original)/2), method = "ward.D")

cor_original <- cor_original[cor_original_order$order, cor_original_order$order]
cor_lower <- cor_lower[cor_original_order$order, cor_original_order$order]
cor_upper <- cor_upper[cor_original_order$order, cor_original_order$order]

#remove upper and diagonal
cor_original[!lower.tri(cor_original)] <- NA
cor_lower[!lower.tri(cor_lower)] <- NA
cor_upper[!lower.tri(cor_upper)] <- NA

cor_diag[lower.tri(cor_diag)] <- NA
cor_diag[upper.tri(cor_diag)] <- NA

#melt matrices
cor_original_long <- data.table(melt(cor_original))
cor_lower_long <- data.table(melt(cor_lower))
cor_upper_long <- data.table(melt(cor_upper))

cor_diag_long <- data.table(melt(cor_diag))
cor_diag_long[value == 1, diag := as.character(Var1)]

#get diagonals of covariance matrices to show differences in variance
cov_original <- round(CovTP_original, 3)
cov_lower <- round(CovTP_lower, 3)
cov_upper <- round(CovTP_upper, 3)

cov_original[lower.tri(cov_original)] <- NA
cov_original[upper.tri(cov_original)] <- NA
cov_lower[lower.tri(cov_lower)] <- NA
cov_lower[upper.tri(cov_lower)] <- NA
cov_upper[lower.tri(cov_upper)] <- NA
cov_upper[upper.tri(cov_upper)] <- NA

#melt matrices
cov_original_long <- data.table(melt(cov_original))
cov_lower_long <- data.table(melt(cov_lower))
cov_upper_long <- data.table(melt(cov_upper))

#generate plot using clustered correlations of original and add annotations for original lower and upper

p <- ggplot(data = cor_original_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(data = cor_diag_long, aes(x = Var1, y = Var2, label = diag), size = 5, fontface = "bold", hjust = "right", nudge_x = 0.3, nudge_y = 0.55) +
  scale_fill_gradient2(low = "#ca0020", high = "#0571b0", mid = "white", midpoint = 0, limit = c(-1, 1),
                       name = "Correlation", na.value = "white") +
  theme_minimal() +
  geom_text(aes(x = Var1, y = Var2, label = value), nudge_y = 0, fontface = "bold") +
  geom_text(data = cor_lower_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0.25) + 
  geom_text(data = cor_upper_long, aes(x = Var1, y = Var2, label = value), nudge_y = -0.25, fontface = "italic") +
  geom_text(data = cov_original_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0, fontface = "bold") +
  geom_text(data = cov_lower_long, aes(x = Var1, y = Var2, label = value), nudge_y = 0.25) +
  geom_text(data = cov_upper_long, aes(x = Var1, y = Var2, label = value), nudge_y = -0.25, fontface = "italic") +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.85),
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_colorbar(barwidth = 18, barheight = 1.5, title.position = "top", title.hjust = 0.5))


pdf("./figures/Figure2_merged_correlation_plot.pdf", width = 8, height = 8)
p
dev.off()

png("./figures/Figure2_merged_correlation_plot.png", width = 8, height = 8, units = "in", res = 300)
p
dev.off()

