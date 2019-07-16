#-----------------------------------------------------------------------------------#

#This script generates figure 4 and also output the power analysis as a txt file#


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

#-----Separate the two nested tables and merge lists into two data tables 
simulated_mMds <- vector("list", length = 2000)


for(x in 1:2000){
  simulated_mMds[[x]] <- simulated_mMd_output[[x]]$all_dists
  
  simulated_mMds[[x]][, sim := x]
}

simulated_mMds <- do.call(rbind, simulated_mMds)


#-----Look at sensitivity of responses over the critical value for the simulated response data

sim_dat <- simulated_mMds[, date := str_split_fixed(date_chnm_plate, "_", 2)[,1]]
sim_dat <- sim_dat[date == "simulated"]

sim_dat_combined <- sim_dat[, .(mD11P = median(D11P), mD11P_lower = median(D11P_lower),
                                mD11P_upper = median(D11P_upper)), by = c("date_chnm_plate", "chnm", "nChem",
                                                                          "Conc", "conc_index")]
for(x in unique(sim_dat[, chnm])){
  for(y in unique(sim_dat[chnm == x, Conc])){
    sim_dat_combined[chnm == x & Conc == y, global := sim_dat[chnm == x & Conc == y & D11P >= 1.642736, .N]/2000]
    sim_dat_combined[chnm == x & Conc == y, lower := sim_dat[chnm == x & Conc == y & D11P_lower >= 1.642736, .N]/2000]
    sim_dat_combined[chnm == x & Conc == y, upper := sim_dat[chnm == x & Conc == y & D11P_upper >= 1.642736, .N]/2000]
  }
}

fwrite(x = sim_dat_combined, file = "./misc/dummyData_powerTable.txt", sep = "\t")

#plot power across response type
sim_dat_combined_long <- melt.data.table(sim_dat_combined, id.vars = colnames(sim_dat_combined)[1:8], variable.name = "type", value.name = "power")
#sim_dat_combined_long[, power := 1-power]
sim_dat_combined_long[, type := factor(type, levels = c("lower", "global", "upper"))]
sim_dat_combined_long[Conc == 1, Conc := 1.1]
sim_dat_combined_long[Conc == 2, Conc := 1.5]
sim_dat_combined_long[Conc == 3, Conc := 2]
sim_dat_combined_long[Conc == 4, Conc := 2.5]
sim_dat_combined_long[power == 1, power := 0.99]
sim_dat_combined_long[, chnm := factor(chnm, levels = c("allUp", "someUp", "gluccoUp", "mineralUp", "andrUp", "testoUp",
                                                        "androgensUp", "e1Up", "e2Up", "estrogensUp"))]

sim_dat_combined_long[, Conc_type := paste(type, Conc, sep = "_")]

p <- ggplot(data = sim_dat_combined_long[chnm != "negative" & Conc != 2.5 & Conc != 2,], aes(x = factor(chnm), y = power, color = factor(Conc_type, levels = c("lower_1.1", "global_1.1", "upper_1.1",
                                                                                                                                                               "lower_1.5", "global_1.5", "upper_1.5"))
                                                                                             , shape = factor(Conc_type, levels = c("lower_1.1", "global_1.1", "upper_1.1",
                                                                                                                                    "lower_1.5", "global_1.5", "upper_1.5")))) +
  geom_point(size = 4, alpha = 0.75)+
  theme_few() +
  geom_hline(yintercept=0.8, linetype = "dashed") +
  annotate("text", x = 9.15, y = 1, label = "*2 and 2.5 fold changes were >0.999", size = 3.75) +
  #ylim(c(0, 0.999)) +
  #ylab("Positive Rate") +
  xlab("") +
  scale_shape_manual(name = "Covariance Type (Fold Change)", labels = c("lower (1.1)", "global (1.1)", "upper (1.1)",
                                                                          "lower (1.5)", "global (1.5)", "upper (1.5)"),
                     values = c(19, 19, 19, 17, 17, 17)) +
  scale_color_manual(name = "Covariance Type (Fold Change)", labels = c("lower (1.1)", "global (1.1)", "upper (1.1)",
                                                                          "lower (1.5)", "global (1.5)", "upper (1.5)"),
                     values = c("#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")) +
  scale_y_continuous(name = "Power", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 0.99)) +
  scale_x_discrete(labels = c("All", "Some", "Glucocorticoids", "Mineralocorticoids", "Androstenedione", "Testosterone", 
                              "Androgens", "Estrone", "Estradiol", "Estrogens")) +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1, vjust = 1), axis.title.y = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12)) +
  guides(shape = guide_legend(ncol = 2), color = guide_legend(ncol = 2))

pdf("./figures/Figure4_powerAnalysis.pdf", width = 14, height = 8)
p
dev.off()

png("./figures/Figure4_powerAnalysis.png", width = 14, height = 8, units = "in", res = 300)
p
dev.off()

