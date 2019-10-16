#-----------------------------------------------------------------------------------#

#Pattern analysis of HT-H295R data near the calculated BMD#
#Generates figure 6#


# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
library(tcpl)
library(data.table)
library(stringr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(fmsb)
library(reshape2)
library(NMF)
library(ggdendro)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

rm(list = ls())

#setwd("L:/Lab/NCCT_ToxCast/Derik Haggard/Mahalanobis Follow Up") #change accordingly

#-----Load in Mahalanobis distances and mMd fits
load("./RData/AllResps_outliersRemoved2018-10-02.RData")

#-----Load in Richard's reference chemical list
ref_list <- fread(file = "./misc/S12 Candidate reference chemicals after mode aggregation.txt", header = TRUE)

ref_list <- ref_list[target == "CYP19A1", ] #only CYP19A1 actives
ref_list[, name := tolower(name)]

#-----Load in HT-H295R ANOVA hitcalls to select chems that affected androgens or estrogens
anova_list <- fread("./supplement/SuppData2_Global_H295R_ANOVA_OECD_directionality_filtered_output_2018-08-09.txt", header = TRUE)
anova_list_wide <- as.data.table(t(anova_list[,-1]), row.names(colnames(anova_list)))
colnames(anova_list_wide) <- c("date_chnm_plate", anova_list[, steroid])

anova_list_wide[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
anova_list_wide[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]
anova_list_wide_filtered <- anova_list_wide[Estrone == -1 & Estradiol == -1, ]

#-----Cast dat_mean
dat_mean <- merge(dat_mean, unique(dat[, c(18, 3)]), by = "date_chnm_plate")
dat_mean[, conc_index := factor(conc_index, levels = c("CONC_DMSO", "CONC6", "CONC5", "CONC4", "CONC3", "CONC2", "CONC1"))]
steroids <- c("Androstenedione", "Testosterone", "Estrone", "Estradiol")
#dat_mean <- dat_mean[steroid %in% steroids] #only include steroids of interest
dat_mean[, date_plate := paste(date, str_split_fixed(date_chnm_plate, "_", 3)[,3], sep = "_")]
dat_mean[, chnm := tolower(chnm)]
#dat_mean <- dat_mean[chnm %in% unique(ref_list[, name]) | chnm == "dmso" | casn %in% unique(ref_list[, casrn])]
dat_mean <- dat_mean[date_chnm_plate %in% anova_list_wide_filtered[, date_chnm_plate] | chnm == "dmso",]

#remove steroid hormones reference chemicals
#dat_mean <- dat_mean[chnm != "5alpha-dihydrotestosterone" & chnm != "estrone" & chnm != "4-androstene-3,17-dione",]


#-----Calculate the log fold change between test chemical and dmso control
for(x in unique(dat_mean[chnm != "dmso", date_chnm_plate])){
  for(y in unique(dat_mean[ date_chnm_plate == x, steroid])){
    for(z in unique(dat_mean[date_chnm_plate == x & steroid == y, date_plate])){
      dat_mean[steroid == y & date_chnm_plate == x & date_plate == z,
                   log_fold_change := mean_log_uM - dat_mean[chnm == "dmso" & steroid == y & date_plate == z, mean_log_uM]]
    }
  }
}

dat_mean[, log2_fold_change := log2(exp(log_fold_change))]
#dat_mean[, steroid := factor(steroid, levels = c("Androstenedione", "Testosterone", "Estrone", "Estradiol"))]

dat_mean_wide <- dcast.data.table(data = dat_mean[chnm != "dmso", c(-6:-9, -11)], formula = date_chnm_plate + chnm + casn ~ steroid + conc_index,
                                      value.var = "log2_fold_change")

#-----Filter data based on -1.5 FC cutoff for E1 and E2
for(x in unique(dat_mean_wide[, date_chnm_plate])){
  dat_mean_wide[date_chnm_plate == x, FC_filter := ifelse(min(dat_mean_wide[date_chnm_plate == x, 34:45], na.rm = TRUE) <= -0.5849625, 1, 0)]
}

dat_mean_wide_filtered <- dat_mean_wide[FC_filter == 1,]
dat_mean_wide_filtered[, chnm := tolower(str_split_fixed(date_chnm_plate, "_", 3)[,2])]

dat_mean_wide_filtered <- dat_mean_wide_filtered[chnm != "dinp branched"]
#-----use hclust to define order of the heatmap rows
# dat_mean_wide <- dcast.data.table(data = dat_mean[chnm != "dmso", c(-6:-9, -11)], formula = date_chnm_plate + chnm + casn ~ steroid + conc_index,
#                                   value.var = "log_fold_change")
# 
# model <- hclust(dist(dat_mean_wide[, c(4:27)]), method = "ward.D")
# order <- model$order

#-----Make the same plot in ComplexHeatmap and add maxmMd row annotations
heat_maxmMd <- Mahalanobis_dat[date_chnm_plate %in% dat_mean_wide[, date_chnm_plate],]

heatmap_dat1 <- as.matrix(copy(dat_mean_wide_filtered[, c(10:15, 64:69, 40:45, 34:39)]))
heatmap_dat2 <- as.matrix(copy(dat_mean_wide_filtered[, c(58:63, 46:57, 28:33, 16:21, 4:9, 22:27)]))
rownames(heatmap_dat1) <- as.character(unlist(dat_mean_wide_filtered[, 2]))
rownames(heatmap_dat2) <- as.character(unlist(dat_mean_wide_filtered[, 2]))

#reference aromatase chemicals annotation
ref_chems <- dat_mean_wide_filtered[chnm %in% unique(ref_list[, name]) | chnm == "dmso" | casn %in% unique(ref_list[, casrn]), chnm]
lab_order <- which(rownames(heatmap_dat2) %in% ref_chems)
#labs <- str_split_fixed(rownames(heatmap_dat2)[lab_order], "_", 3)[,2]
labs <- rownames(heatmap_dat2)[lab_order]
labs <- c("Clotrimazole", "Chrysin", "Fadrozole hydrochloride", "Letrozole", "Metyrapone")

#define grouping of dendrogram for the split argument in Heatmap
dend <- hclust(dist(heatmap_dat1, ), method="ward.D")
dend_cut <- cutree(dend, h = 18)
dend_cut[dend_cut == 5] <- 3
dend_cut[dend_cut == 4] <- "Cluster 2"
dend_cut[dend_cut == 3] <- "Cluster 1"
dend_cut[dend_cut == 2] <- "Cluster 4"
dend_cut[dend_cut == 1] <- "Cluster 3"

ha_row <- rowAnnotation(df = data.frame(log10_maxmMd = log10(heat_maxmMd[, maxD11P])),
                  col = list(log10_maxmMd = colorRamp2(c(0, 1.5), c("white", "#b2182b"))), width = unit(5, "mm"),
                  annotation_legend_param = list(title_gp = gpar(fontsize = 14), title = bquote(paste("log"[10], " maxmMd"))))

h1_label <- HeatmapAnnotation(cn = function(index){
  grid.text("ANDR", x = 0.125, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("TESTO", x = 0.375, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("E1", x = 0.625, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("E2", x = 0.875, 1, just = "top", gp = gpar(fontsize = 13))
  })

h2_label <- HeatmapAnnotation(cn = function(index){
  grid.text("PROG", x = 1/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("OHPREG", x = 3/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("OHPROG", x = 5/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("DOC", x = 7/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("CORTIC", x = 9/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("11DCORT", x = 11/14, 1, just = "top", gp = gpar(fontsize = 13))
  grid.text("CORT", x = 13/14, 1, just = "top", gp = gpar(fontsize = 13))
})

heatmap1 <- Heatmap(name = "log2 FC", matrix = heatmap_dat1, cluster_columns = FALSE, col = colorRamp2(breaks = c(-2, 0, 2),
                                                                                   colors = c("#2166ac", "#f7f7f7", "#b2182b")),
                   show_row_names = FALSE, row_names_side = "left", show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
                   clustering_method_rows = "ward.D", bottom_annotation = h1_label, cluster_rows = TRUE,
                   column_title = "Aromatase-relevant", heatmap_legend_param = list(title = bquote(paste("log"[2], " fold-change")), title_gp = gpar(fontsize = 14,
                                                                                                    labels_gp = gpar(fontsize = 12))),
                   row_dend_width = unit(2, "cm"), row_dend_reorder = FALSE, split = dend_cut)

heatmap2 <- Heatmap(name = "heat2", matrix = heatmap_dat2, cluster_columns = FALSE, col = colorRamp2(breaks = c(-2, 0, 2),
                                                                                                    colors = c("#2166ac", "#f7f7f7", "#b2182b")),
                    show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE,
                    bottom_annotation =  h2_label, show_heatmap_legend = FALSE, clustering_method_rows = "ward.D",
                    column_title = "Rest of Pathway")
                    
row_annotation <- rowAnnotation(link = row_anno_link(at = lab_order, labels = labs, labels_gp = gpar(fontfontsize = 13)))

padding <- unit.c(unit(8, "mm"), unit(c(2, 4), "mm"), grobWidth(textGrob("long_text_string_for_padding_in_plot")) - unit(1, "cm"))

# ht_list1 <- ha_row + heatmap1 + row_annotation#+ heatmap2 #+ row_annotation
# 
# pdf("./figures/test_heatmap1_081618.pdf", width = 12, height = 8)
# draw(ht_list1, row_dend_side = "left", padding = padding)
# 
# decorate_heatmap_body("log2 FC", {
#   grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
# })
# 
# dev.off()
# 
# png("./figures/test_heatmap1_081618.png", width = 12, height = 8, units = "in", res = 300)
# draw(ht_list1, row_dend_side = "right", padding = padding)
# 
# decorate_heatmap_body("log2 FC", {
#   grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
# })
# 
# dev.off()

ht_list2 <- heatmap1 + heatmap2 + ha_row + row_annotation

pdf("./figures/Figure5_aromataseHeatmap.pdf", width = 22, height = 12)
draw(ht_list2, row_dend_side = "left", heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)

decorate_heatmap_body("log2 FC", slice = 1, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
#   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
#   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 2, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 3, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 4, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})

decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
#   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
#   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
#   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 1)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 2)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 3)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 4)
dev.off()

png("./figures/Figure6_aromataseHeatmap.png", width = 22, height = 12, res = 300, units = "in")
draw(ht_list2, row_dend_side = "left", heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)

decorate_heatmap_body("log2 FC", slice = 1, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 2, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 3, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})
decorate_heatmap_body("log2 FC", slice = 4, {
  grid.lines(c(0.25, 0.25), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.75, 0.75), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
})

decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 1)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 2)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 3)
decorate_heatmap_body("heat2", {
  grid.lines(c(0.1428, 0.1428), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.2857, 0.2857), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.4285, 0.4285), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.5714, 0.5714), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.7142, 0.7142), c(0, 1), gp = gpar(lwd = 2))
  grid.lines(c(0.8571, 0.8571), c(0, 1), gp = gpar(lwd = 2))
  #   grid.lines(c(0, 1), c(12/72, 12/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(20/72, 20/72), gp = gpar(lwd = 2, lty = 2))
  #   grid.lines(c(0, 1), c(47/72, 47/72), gp = gpar(lwd = 2, lty = 2))
}, slice = 4)
dev.off()


save(dat_mean_wide_filtered, file = "./misc/h295R_aromataseInhibitors.RData")
