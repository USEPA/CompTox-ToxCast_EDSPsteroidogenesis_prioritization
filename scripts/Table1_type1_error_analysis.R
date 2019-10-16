##############################################
#This script generates the Type I error table#
##############################################

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
load("./RData/AllResps_outliersRemoved2017-08-09.RData")
load("./RData/simulated_Mahalanobis_distance_output_2018-08-14.RData")
load("./RData/subsample_Model_output2018-06-01.RData")
#load("./simulated_Mahalanobis_distance_output_originalCovariance_2018-05-10.RData")
#load("./Mahalanobis_distance_output_withSimulatedCovariances_2018-05-10.RData")

#-----Separate the two nested tables and merge lists into two data tables 
simulated_mMds <- vector("list", length = 2000)


for(x in 1:2000){
  simulated_mMds[[x]] <- simulated_mMd_output[[x]]$all_dists
  
  simulated_mMds[[x]][, sim := x]
}

simulated_mMds <- do.call(rbind, simulated_mMds)

#-----Type I error analysis
#Critical values for a chemical with 11 hormones and N = 2 is 1.642736 (0.01) and 1.503230 (0.05)
negative_maxmMd <- simulated_mMds[chnm == "negative",]
negative_maxmMd <- negative_maxmMd[, .(maxmMd = max(D11P), maxmMd_upper = max(D11P_upper), maxmMd_lower = max(D11P_lower)),
                                   by = c("date_chnm_plate", "sim", "chnm", "NH", "N")]

negative_maxmMd[maxmMd >= 1.642736, .N] #16 of 2000 (0.008)
negative_maxmMd[maxmMd >= 1.503230, .N] #38 of 2000 (0.019)
negative_maxmMd[maxmMd_lower >= 1.642736, .N] #13 of 2000 (0.0065)
negative_maxmMd[maxmMd_lower >= 1.503230, .N] #30 of 2000 (0.015)
negative_maxmMd[maxmMd_upper >= 1.642736, .N] #3 of 2000 (0.0055)
negative_maxmMd[maxmMd_upper >= 1.503230, .N] #10 of 2000 (0.0115)
