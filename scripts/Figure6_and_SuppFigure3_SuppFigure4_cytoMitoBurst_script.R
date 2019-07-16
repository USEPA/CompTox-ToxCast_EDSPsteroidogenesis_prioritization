#-----------------------------------------------------------------------------------#

#Selectivity scoring of the Mahalanobis distance through
#incorporation of cyto and mito burst.
#Genrate figure 6 and supplemental figures 3 and 4


#devtools::install_git("http://dhaggard@bitbucket.zn.epa.gov/scm/tox/tcpl.git", branch = 'scitovation-updates')

library(tcpl)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(plyr)
library(stringr)
library(ggpubr)

rm(list = ls())

#setwd("L:/Lab/NCCT_ToxCast/Derik Haggard/Mahalanobis Follow Up") #change accordingly
tcplConf(user = "dhaggard", pass = "new_pass_dhaggard", db = "invitrodb", host = "mysql-res1.epa.gov", drvr = "MySQL")

#-----Load in Mahalanobis distances and mMd fits
load("./RData/AllResps_outliersRemoved2018-10-02.RData")
#load("./RData/mMdFits2018-10-03.RData")
load("./RData/mMD_AUC_Fits_2018-10-08.RData") #load in AUC Fits
#load("./RData/old/mMdFits.RData")

Fits <- copy(data.table(Fits_AUC))
colnames(Fits)[1] <- "date_chnm_plate"
Fits[, colnames(Fits)[2:12] := .(as.numeric(T), as.numeric(cc), as.numeric(d), as.numeric(BMD), as.numeric(MaxmMd),
                                 as.numeric(Scrit), as.numeric(cor), as.numeric(cor_pvalue), as.numeric(convergence),
                                 as.numeric(AUC_trap), as.numeric(AUC_stats))]
Fits[BMD >= 150, BMD := 1e03]
Fits[BMD <= 1e-04, BMD := 1e-04]

dats <- merge(Fits, Mahalanobis_dat[, c("date_chnm_plate", "casn")], by = "date_chnm_plate")
dats[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]

for(x in unique(dat[chnm != "DMSO", date_chnm_plate])){
  dats[date_chnm_plate == x, spid := unique(dat[date_chnm_plate == x, spid])]
}

dats <- dats[MaxmMd >= Scrit]

#-----Load all chemical IDs in HT-H295R in tcpl
chems <- tcplLoadChem()
chids <- chems[casn %in% unique(dat[, casn]),]

#-----Generate the cyto burst data for HT-H295R mc chemicals
cyto_H295R <- tcplCytoPt(chid = unique(chems[, chid]))
cyto_H295R <- cyto_H295R[chid %in% unique(chids[, chid])]

#-----Get all of the aeids for the mito-relevant assays and run tcplCytoPt on the mito aeids
#assays <- tcplLoadAeid()
#aeids <- assays[str_detect(aenm, "MitoM") | str_detect(aenm, "MMP_ratio") | str_detect(aenm, "MitoFxn") | str_detect(aenm, "hPBR") | str_detect(aenm, "rPBR"),]

#mito_H295R <- tcplCytoPt(chid = unique(chems[, chid]), aeid = unique(aeids[, aeid])) #try combinig cyto and mito burst aeids
#names(mito_H295R)[11:12] <- c("mito_pt", "mito_pt_um")
#mito_H295R <- mito_H295R[chid %in% unique(chids[, chid])]

#-----Add the HT-H295R MTT ACC values (conc at the 70% cutoff)
tcplLoadAeid(fld='asid', val=8)
mtt.mc5 <-tcplLoadData(lvl=5,fld='aeid', val=1663, type='mc')

for(x in unique(chems[spid %in% mtt.mc5[, spid], spid])){
  mtt.mc5[spid == x, chnm := unique(chems[spid == x, chnm])]
  mtt.mc5[spid == x, casn := unique(chems[spid == x, casn])]
}

mtt.mc5[, chnm := str_replace(chnm, "GSID_", "GSID")]
mtt.mc5 <- mtt.mc5[spid != "DMSO"]
mtt.mc5[spid == "TP0000884B12" | spid == "TP0000884A04" | spid == "TP0000883A04", chnm := "Triadimenol2"]
mtt.mc5 <- mtt.mc5[chnm != "Digoxigenin" & chnm != "Colchicine",]

mtt.mc5 <- mtt.mc5[spid %in% unique(dat[, spid])]

#-----Add cyto and mito burst values and top MTT concs to dats
distance_dat <- merge(dats, unique(mtt.mc5[, c("spid", "modl_ga", "modl_acc")]), by = "spid")
distance_dat <- merge(distance_dat, cyto_H295R[, c("casn", "cyto_pt_um", "lower_bnd_um")], by = "casn") 
#distance_dat <- merge(distance_dat, mito_H295R[, c("casn", "mito_pt_um")], by = "casn")

distance_dat[lower_bnd_um <= 1e-04, lower_bnd_um := 1e-04]
#distance_dat[mito_pt_um <= 1e-04, mito_pt_um := 1e-04]
distance_dat[is.na(modl_ga), modl_ga := 3]
distance_dat[is.na(modl_acc), modl_acc := 3]

distance_dat[, log_cyto_pt_um := log10(cyto_pt_um)]
distance_dat[, log_lower_bnd_um := log10(lower_bnd_um)]
#distance_dat[, log_mito_pt_um := log10(mito_pt_um)]
distance_dat[, log_BMD := log10(BMD)]

#-----Calculate the fold change compared to DMSO control to use for filtering
# dat[, date_plate := paste(date, plate, sep = "_")]
# dmso_dat <- dat[chnm == "DMSO", .(mean_um = mean(uM)), by = c("steroid", "date_chnm_plate", "date_plate")]
# 
# #for(x in unique(dat[, date_chnm_plate])){
#   for(y in unique(dat[, steroid])){
#     for(z in unique(dat[, date_plate])){
#       dat[steroid == y & date_plate == z, fold_change := uM/dmso_dat[steroid == y & date_plate == z, mean_um]]
#     }
#   }
# 
# dat[, log_fold_change := log(fold_change)]
# dat_max <- dat[chnm != "DMSO", .(max_log_FC = max(abs(log_fold_change))), by = c("date_chnm_plate")]

#-----Filter based on a maximum absolute fold change of 1.5 or larger
# distance_dat <- distance_dat[date_chnm_plate %in% dat_max[max_log_FC >= log(1.5), date_chnm_plate],]

#-----Add in bordelrine hit filter using the derived 0.33 standard deviation
#-----(i.e. multiplying the Scrit by 2*SD) to get the minimum maxmMd needed to be
#-----outisde the 95% confidence interval
distance_dat[, upper_bound := Scrit * exp(0.66)]
distance_dat[, borderline := ifelse(MaxmMd >= upper_bound, "No", "Yes")]

#-----Add selectivity metric
distance_dat[, selectivity := .(pmin(lower_bnd_um, 10^modl_acc) - BMD)]
distance_dat[, log_selectivity := .(pmin(log_lower_bnd_um, modl_acc) - log_BMD)]
#distance_dat[, selectivity := .(pmin(mito_pt_um, lower_bnd_um, 10^modl_acc) - BMD)]
#distance_dat[, log_selectivity := .(pmin(log_mito_pt_um, log_lower_bnd_um, modl_acc) - log_BMD)]

distance_dat[, spid := NULL]
#distance_dat[, selectivity := .(pmin(mito_pt_um, cyto_pt_um) - BMD)]
#distance_dat[, log_selectivity := .(pmin(log_mito_pt_um, log_cyto_pt_um) - log_BMD)]

#distance_dat[, viability_selectivity := .(pmin(mito_pt_um, cyto_pt_um) - max.c.by.chem)]
#distance_dat[, log_viability_selectivity := .(pmin(log_mito_pt_um, log_cyto_pt_um) - max.logc.by.chem)]

#-----Reformat data for plotting
distance_dat_long <- melt.data.table(distance_dat, id.vars = colnames(distance_dat)[c(1:6, 8:14, 22:25)],
                                     variable.name = "model_type", value.name = "model_val")

distance_dat_long[model_type != "BMD" & model_type != "log_BMD", MaxmMd := 3]
distance_dat_long[model_type != "BMD" & model_type != "log_BMD", borderline := "No"]

distance_dat_long_log <- distance_dat_long[str_detect(model_type, "log") | str_detect(model_type, "modl_acc"),]
distance_dat_long_log[model_val >= 3, model_val := 3]

distance_dat_long_no_log <- distance_dat_long[!str_detect(model_type, "log"),]

# #-----Generate plot
# p <- ggplot(data = distance_dat_long_log[model_type != "log_cyto_pt_um" & log_selectivity >= 0.5,], aes(x = model_val, y = reorder(date_chnm_plate, log_selectivity), group = model_type)) +
#   geom_point(aes(color = factor(model_type, levels = c("log_BMD", "modl_acc", "log_lower_bnd_um"))), alpha = 0.5, size = 4) +
#   geom_hline(yintercept = 201.5, linetype = "dashed", alpha = 0.5) + #422 samples have selectivity >= 0.5 log units
#   theme_few() +
#   #scale_x_log10() +
#   theme(axis.text.y = element_text(size = 5)) +
#   ylab("Samples") +
#   xlab(expression(log[10]~Concentration~(uM))) +
#   scale_size(name = "maxmMd Value", breaks = c(1, 2, 4, 8, 16, 32)) +
#   scale_color_manual(name = "Model Type", labels = c("BMD", expression("MTT"[acc]), expression("cyto"[burst])), values = c("#F8766D", "#7CAE00", "#00BFC4")) +
#   theme(axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"),
#         legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
#         legend.text.align = 0, axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#   #scale_shape(name = "Borderline") +
#   annotation_logticks(base = 10, sides = "b") +
#   scale_x_continuous(breaks = seq(-4, 4, 1))
# 
# pdf("./supplement/SupplementalFigure_selectivityPlot.pdf", width = 12, height = 20)
# p
# dev.off()
# 
# png("./supplement/SupplementalFigure_selectivityPlot.png", width = 12, height = 20, units = "in", res = 300)
# p
# dev.off()

#-----Look at chemical overlap in selectivity and AUC
Fits <- Fits[order(AUC_stats),]

AUC_selected <- distance_dat[log_selectivity >= 0.5,]

AUC_selected <- AUC_selected[order(AUC_stats),]

distance_dat_long_log <- distance_dat_long_log[order(log_selectivity),]

#selected_chems <- unique(distance_dat_long_log[, date_chnm_plate])[c(1:25, 599:623)] #top and bottom 25 chemicals
selected_chems <- unique(AUC_selected[, date_chnm_plate])[c(1:25, 404:428)] #top and bottom 25 chemicals
#selected_chems <- unique(AUC_selected[, date_chnm_plate])[c(1:20, 202:221, 403:422)] #top, middle, and bottom 20 chemicals

p <- ggplot(data = distance_dat_long_log[date_chnm_plate %in% selected_chems & model_type != "log_cyto_pt_um",], aes(x = model_val, y = reorder(date_chnm_plate, AUC_stats), group = model_type)) +
  geom_point(aes(color = factor(model_type, levels = c("log_BMD", "modl_acc", "log_lower_bnd_um"))), alpha = 0.75, size = 4) +
  theme_few() +
  #scale_x_log10() +
  theme(axis.text.y = element_text(size = 5)) +
  geom_hline(yintercept = 25.5, linetype = "dashed", alpha = 0.5) +
  #geom_hline(yintercept = 20.5, linetype = "dashed", alpha = 0.5) + #for top, mid, bottom 20
  #geom_hline(yintercept = 40.5, linetype = "dashed", alpha = 0.5) + #for top, mid, bottom 20
  ylab("") +
  xlab(expression(log[10]~Concentration~(uM))) +
  scale_size(name = "maxmMd Value", breaks = c(1, 2, 4, 8, 16, 32)) +
  scale_color_manual(name = "Model Type", labels = c("BMD", expression("MTT"[acc]), expression("cyto"[burst])), values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  theme(axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10), legend.position = "bottom",
        legend.text.align = 0, axis.text.y = element_text(size = 10)) +
  #scale_shape(name = "Borderline") +
  annotation_logticks(base = 10, sides = "b") +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  scale_y_discrete(labels = str_split_fixed(selected_chems, "_", 3)[,2])


#-----Use geom_tile to add annotation showing selectivity, maxmMd, and AUC for the selected chemicals

p1 <- ggplot(data = distance_dat[date_chnm_plate %in% selected_chems]) +
  geom_text(aes(label = sprintf("%.2f", round(log_selectivity, 2)), x = "A", y = reorder(date_chnm_plate, AUC_stats), hjust = "middle")) +
  geom_text(aes(label = sprintf("%.2f", round(AUC_stats, 2)), x = "B", y = reorder(date_chnm_plate, AUC_stats), hjust = "middle")) +
  geom_text(aes(label = sprintf("%.2f", round(MaxmMd, 2)), x = "C", y = reorder(date_chnm_plate, AUC_stats), hjust = "middle")) +
  theme_few() +
  theme(rect = element_blank(), line = element_blank(),axis.title = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, hjust = 0.5), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,2,0,0), "mm")) +
  scale_x_discrete(labels = c("Selectivity", "AUC", "maxmMd"), position = "top")

#-----Merge p and p1

p3 <- ggarrange(p, p1, common.legend = TRUE, widths = c(1, 0.2), align = "h", legend = "right")

pdf("./figures/Figure6_selectivityPlot.pdf", width = 15, height = 12)
p3
dev.off()

png("./figures/Figure6_selectivityPlot.png", width = 15, height = 12, units = "in", res = 300)
p3
dev.off()
# 
# pdf("./figures/test.pdf", width = 14, height = 12)
# p3
# dev.off()
# 
# png("./figures/test.png", width = 14, height = 12, units = "in", res = 300)
# p3
# dev.off()

#-----Generate plot of all selective chemicals, but rank by AUC_stats
p <- ggplot(data = distance_dat_long_log[model_type != "log_cyto_pt_um" & log_selectivity >= 0.5,], aes(x = model_val, y = reorder(date_chnm_plate, AUC_stats), group = model_type)) +
  geom_point(aes(color = factor(model_type, levels = c("log_BMD", "modl_acc", "log_lower_bnd_um"))), alpha = 0.5, size = 4) +
  #geom_hline(yintercept = 201.5, linetype = "dashed", alpha = 0.5) + #422 samples have selectivity >= 0.5 log units
  theme_few() +
  #scale_x_log10() +
  theme(axis.text.y = element_text(size = 5)) +
  ylab("Samples") +
  xlab(expression(log[10]~Concentration~(uM))) +
  scale_size(name = "maxmMd Value", breaks = c(1, 2, 4, 8, 16, 32)) +
  scale_color_manual(name = "Model Type", labels = c(expression("BMD", "MTT"[acc]), expression("cyto"[burst])), values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  theme(axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
        legend.text.align = 0, axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #scale_shape(name = "Borderline") +
  annotation_logticks(base = 10, sides = "b") +
  scale_x_continuous(breaks = seq(-4, 4, 1))

pdf("./supplement/SupplementalFigure3_selectivityPlot_AUC.pdf", width = 12, height = 20)
p
dev.off()

png("./supplement/SupplementalFigure3_selectivityPlot_AUC.png", width = 12, height = 20, units = "in", res = 300)
p
dev.off()

#-----Save output

save(distance_dat, distance_dat_long_log, file = paste("./RData/mMd_selectivityScoring_output_", Sys.Date(), ".RData", sep = ""))

#-----Generate supplemental table with AUC and other mMd metrics
supp_table <- distance_dat[, c(1, 2, 7, 21, 16, 19:20, 25, 13)]

fwrite(supp_table, "./supplement/SuppTable1_rankPrioritizationTable.txt", sep = "\t")

#-----Generate supplemental figure of maxmMd and AUC

summary(lm(AUC_stats ~ MaxmMd, data = distance_dat))

p1 <- ggplot(data = distance_dat, aes(x = MaxmMd, y = AUC_stats)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  theme_few() +
  xlab("maxmMd") +
  ylab("AUC")

p2 <- ggplot(data = distance_dat, aes(x = MaxmMd, y = log_selectivity)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  theme_few() +
  xlab("maxmMd") +
  ylab(expression('log'[10]*" Selectivity"))

pdf("./supplement/SupplementalFigure3_maxmMd_AUC_correlationPlot.pdf", width = 14, height = 10)
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2)
dev.off()

png("./supplement/SupplementalFigure3_maxmMd_AUC_correlationPlot.png", width = 14, height = 10, res = 300, units = "in")
ggarrange(p1, p2, labels = c("A", "B"), ncol = 2)
dev.off()

#only plot selectivity by sample
p <- ggplot(data = distance_dat, aes(x = log_selectivity, y = reorder(date_chnm_plate, log_selectivity))) +
  geom_point(aes(size = log(MaxmMd), shape = factor(borderline))) +
  theme_few() +
  #scale_x_log10() +
  theme(axis.text.y = element_text(size = 5))
#annotation_logticks(sides = "b")

p

