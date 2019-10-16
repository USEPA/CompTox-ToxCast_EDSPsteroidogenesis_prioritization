#-----------------------------------------------------------------------------------#

#This script generates figure 4 and supp data 1#


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
load("./RData/subsample_mMd_output2018-10-03.RData")
#load("./simulated_Mahalanobis_distance_output_originalCovariance_2018-05-10.RData")
#load("./Mahalanobis_distance_output_withSimulatedCovariances_2018-05-10.RData")

#-----Separate the two nested tables and merge lists into two data tables 
simulated_mMds <- vector("list", length = 2000)


for(x in 1:2000){
  simulated_mMds[[x]] <- simulated_mMd_output[[x]]$all_dists
  
  simulated_mMds[[x]][, sim := x]
}

simulated_mMds <- do.call(rbind, simulated_mMds)

#-----Format data to plot distribution of mMds with varying covariance type and generate Supplemental Data 1
#remove original lower and upper
simulated_mMds <- simulated_mMds[, c(-11, -12)]

simulated_mMds_long <- melt(simulated_mMds, id.vars = colnames(simulated_mMds)[c(1:2, 4:8, 11)], value.name = "mMds",
                            variable.name = "simulation_set")

simulated_mMds_long[, simulation_set := factor(simulation_set, levels = c("D11P_lower", "D11P", "D11P_upper"))]

simulated_mMds_long[simulation_set == "D11P_lower", subsample := "lower"]
simulated_mMds_long[simulation_set == "D11P", subsample := "original"]
simulated_mMds_long[simulation_set == "D11P_upper", subsample := "upper"]
simulated_mMds_long[, subsample := factor(subsample, levels = c("lower", "original", "upper"))]

Dists_merged <- rbind(Dists_lower_subsample, Dists[, subsample := "original"], Dists_upper_subsample)
Dists_merged[, subsample := factor(subsample, levels = c("lower", "original", "upper"))]

#supplemendat data 
pdf("./supplement/SupplementalFigure2_mMdDists.pdf", width = 12, height = 8)
for(x in unique(simulated_mMds_long[, date_chnm_plate])){
  p <- ggplot(data = simulated_mMds_long[date_chnm_plate == x,], aes(x = factor(Conc), y = mMds, color = subsample)) +
    geom_boxplot(position = position_dodge()) +
    geom_point(data = Dists_merged[CA == x], aes(x = factor(Conc), y = D11P, color = subsample), position = position_dodge(0.75), size = 3, shape = 17) +
    theme_few() +
    ylab("mMd") +
    xlab("Concentration (uM)") +
    ggtitle(x)
  
  plot(p)
}
dev.off()

#grab plots for chemicals in the upper, lower, and a moderate chemical as examples for Figure 3
prochloraz <- "20130321_Prochloraz_Plate4"
eds <- "20170411_Ethylene dimethanesulfonate_Plate3"
butylparaben <- "20150610_Butylparaben_Plate12"

prochloraz_plot <- ggplot(data = simulated_mMds_long[date_chnm_plate == prochloraz,], aes(x = factor(Conc), y = mMds, color = subsample)) +
  geom_boxplot(position = position_dodge()) +
  geom_point(data = Dists_merged[CA == prochloraz], aes(x = factor(Conc), y = D11P, color = subsample), position = position_dodge(0.75), size = 3, shape = 17) +
  theme_few() +
  ylab("mMd") +
  xlab("Concentration (uM)") +
  theme(legend.position = "none") +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  ylim(c(0, 20))

butylparaben_plot <- ggplot(data = simulated_mMds_long[date_chnm_plate == butylparaben,], aes(x = factor(Conc), y = mMds, color = subsample)) +
  geom_boxplot(position = position_dodge()) +
  geom_point(data = Dists_merged[CA == butylparaben], aes(x = factor(Conc), y = D11P, color = subsample), position = position_dodge(0.75), size = 3, shape = 17) +
  theme_few() +
  ylab("mMd") +
  xlab("Concentration (uM)") +
  theme(legend.position = "none") +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  ylim(c(0, 20))

eds_plot <- ggplot(data = simulated_mMds_long[date_chnm_plate == eds,], aes(x = factor(Conc), y = mMds, color = subsample)) +
  geom_boxplot(position = position_dodge()) +
  geom_point(data = Dists_merged[CA == eds], aes(x = factor(Conc), y = D11P, color = subsample), position = position_dodge(0.75), size = 3, shape = 17) +
  theme_few() +
  ylab("mMd") +
  xlab("Concentration (uM)") +
  theme(legend.position = "none") +
  theme(axis.title = element_text(face = "bold", size = 14)) +
  ylim(c(0, 20))

#-----95% confidence intervals for mMds based on simulation
simulated_mMds_CI <- copy(simulated_mMds)

simulated_mMds_CI[, ID := paste(date_chnm_plate, Conc, conc_index, sep = "_")]

conf_intervals <- vector("list", length = 4499)
#names(conf_intervals) <- unique(simulated_mMds_CI[, ID])

interval_calc <- function(z){
  
  x <- unique(simulated_mMds_CI[z, ID])
  
  mD11P <- median(simulated_mMds_CI[ID == x, D11P])
  mD11P_upper <- median(simulated_mMds_CI[ID == x, D11P_upper])
  mD11P_lower <- median(simulated_mMds_CI[ID == x, D11P_lower])
  #mD11P_original_upper <- mean(simulated_mMds_CI[ID == x, D11P_original_upper])
  #mD11P_original_lower <- mean(simulated_mMds_CI[ID == x, D11P_original_lower])
  
  m_CI <- quantile(simulated_mMds_CI[ID == x, D11P], c(0.025, 0.975))
  m_upper_CI <- quantile(simulated_mMds_CI[ID == x, D11P_upper], c(0.025, 0.975))
  m_lower_CI <- quantile(simulated_mMds_CI[ID == x, D11P_lower], c(0.025, 0.975))
  #m_original_upper_CI <- quantile(simulated_mMds_CI[ID == x, D11P_original_upper], c(0.025, 0.975))
  #m_original_lower_CI <- quantile(simulated_mMds_CI[ID == x, D11P_original_lower], c(0.025, 0.975))
  
  intervals <- data.table(date_chnm_plate = paste(str_split_fixed(x, "_", 5)[,1], str_split_fixed(x, "_", 5)[,2], str_split_fixed(x, "_", 5)[,3], sep = "_"),
                          conc = str_split_fixed(x, "_", 5)[,4], conc_index = str_split_fixed(x, "_", 5)[,5],
                          mD11P = mD11P, mD11P_upCI = m_CI[2], mD11P_lowCI = m_CI[1],
                          mD11P_upper = mD11P_upper, mD11P_upper_upCI = m_upper_CI[2], mD11P_upper_lowCI = m_upper_CI[1],
                          mD11P_lower = mD11P_lower, mD11P_lower_upCI = m_lower_CI[2], mD11P_lower_lowCI = m_lower_CI[1]
                          #mD11P_original_upper = mD11P_original_upper, mD11P_original_upper_upCI = m_original_upper_CI[2], mD11P_original_upper_lowCI = m_original_upper_CI[1],
                          #mD11P_original_lower = mD11P_original_lower, mD11P_original_lower_upCI = m_original_lower_CI[2], mD11P_original_lower_lowCI = m_original_lower_CI[1]
  )
  return(intervals)
}

conf_intervals <- mclapply(1:4499, interval_calc, mc.cores = 30)

conf_intervals_merged <- do.call(rbind, conf_intervals)
conf_intervals_merged <- conf_intervals_merged[order(mD11P)]

colnames(Dists_all_subsample)[c(3, 4)] <- c("original_D11P", "date_chnm_plate")
Dists_all_subsample[, ID := paste(date_chnm_plate, Conc, sep = "_")]
conf_intervals_merged[, ID := paste(date_chnm_plate, conc, sep = "_")]

conf_intervals_merged <- merge(conf_intervals_merged, Dists_all_subsample[, c(3, 10)], by = "ID")
conf_intervals_merged[, ID := NULL]

conf_intervals_merged_long <- melt.data.table(conf_intervals_merged, id.vars = c("date_chnm_plate", "conc", "conc_index", "original_D11P"),
                                              value.name = c("median_D11P", "upper_CI", "lower_CI"), variable.name = "simulation_type",
                                              measure.vars = list(c("mD11P", "mD11P_lower", "mD11P_upper"),#, "mD11P_original_lower", "mD11P_original_upper"),
                                                                  c("mD11P_upCI", "mD11P_lower_upCI", "mD11P_upper_upCI"),# "mD11P_original_lower_upCI", "mD11P_original_upper_upCI"),
                                                                  c("mD11P_lowCI", "mD11P_lower_lowCI", "mD11P_upper_lowCI")))#, "mD11P_original_lower_lowCI", "mD11P_original_upper_lowCI")))

conf_intervals_merged_long[simulation_type == 1, simulation_type := "global"]
conf_intervals_merged_long[simulation_type == 2, simulation_type := "lower"]
conf_intervals_merged_long[simulation_type == 3, simulation_type := "upper"]
conf_intervals_merged_long[, simulation_type := factor(simulation_type, levels = c("lower", "global", "upper"))]
#conf_intervals_merged_long[simulation_type == 4, simulation_type := "original_lower"]
#conf_intervals_merged_long[simulation_type == 5, simulation_type := "original_upper"]

#plot for scatterplot of simulated mMd medians vs original mMds
scatterplot <- ggplot(data = conf_intervals_merged_long, aes(x = original_D11P, y = median_D11P, color = simulation_type)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), size = 1) +
  geom_point(size = 4, alpha = 0.75)+
  geom_abline(slope = 1, linetype = "dashed", alpha = 0.65) +
  theme_few() +
  ylim(c(0, 55)) +
  xlim(c(0, 55)) +
  xlab("Original mMd") +
  ylab("Median Simulated mMd") +
  scale_color_discrete(name = "Covariance Type") +
  theme(axis.title = element_text(face = "bold", size = 16),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))

#-----Generate final plot

p <- ggarrange(scatterplot, ggarrange(prochloraz_plot,
                                      butylparaben_plot,
                                      eds_plot, ncol = 3,
                                      labels = c("B", "C", "D")),
               nrow = 2, labels = "A")

pdf("./figures/Figure4_scatterplot_mMdBoxplots.pdf", width = 12, height = 10)
p
dev.off()

png("./figures/Figure4_scatterplot_mMdBoxplots.png", width = 12, height = 10, units = "in", res = 300)
p
dev.off()


