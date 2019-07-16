#-----------------------------------------------------------------------------------#
# Tox21 Aromatase inhibition
# information from paper

# Katie Paul Friedman, paul-friedman.katie@epa.gov

# 27 Nov 2018
# Edited    
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Objectives

# Aromatase inhibition from Tox21 assay
# 10K screen is confounded by ERa-antagonists and cytotoxicity
# followup in Chen et al 2015 is analyzed here
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# Loading libraries
#-----------------------------------------------------------------------------------#

rm(list=ls())
library(data.table)
library(openxlsx)
library(plyr)
library(reshape2)
library(tcpl)
tcplConf(user='_dataminer', pass='pass', db='invitrodb', host='mysql-res1.epa.gov', drvr='MySQL')
library(VennDiagram) # capped at 5 areas
library(venn) # will do up to 7 areas

library(cowplot)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(gtable)
#library(ggExtra)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(GGally)
library(stringr)

#setwd("./additional_aromatase/") 
getwd()


#-----------------------------------------------------------------------------------#
# Loading files
#-----------------------------------------------------------------------------------#

tox21.aromatase <- fread('./additional_aromatase/tox21_putative_aromatase_inhibitors.csv') #50 putative aromatase inhibitors from screen
#Supp6 <- fread('./additional_aromatase/Supp6_Global_H295R_ANOVA_OECD_filtered_wide_output_2017-08-09.txt') # Haggard et al 2018
#Supp2 <- fread("./additional_aromatase/Supp2_Global_H295R_ANOVA_OECD_directionality_filtered_output_2018-08-09.txt")

tcplLoadAsid() #5 is NVS
nvs.aeids <- tcplLoadAeid(val=5, fld='asid') #aeid = 319 for NVS aromatase

nvs.arom <- tcplPrepOtpt(tcplLoadData(lvl=5,type='mc',fld='aeid', val=319))

#-----------------------------------------------------------------------------------#
#-----regenerate list of chemicals with significant decreases in E1 and/or E2 with 1.5 FC threshold
#-----------------------------------------------------------------------------------#

#-----Load in Mahalanobis distances and mMd fits
load("./RData/AllResps_outliersRemoved2018-10-02.RData")

#-----Load in HT-H295R ANOVA hitcalls to select chems that affected androgens or estrogens
anova_list <- fread("./supplement/Supp2_Global_H295R_ANOVA_OECD_directionality_filtered_output_2018-08-09.txt", header = TRUE)
anova_list_wide <- as.data.table(t(anova_list[,-1]), row.names(colnames(anova_list)))
colnames(anova_list_wide) <- c("date_chnm_plate", anova_list[, steroid])

anova_list_wide[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
anova_list_wide[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]
anova_list_wide_filtered <- copy(anova_list_wide)
#anova_list_wide_filtered <- anova_list_wide[Estrone == -1 | Estradiol == -1, ]

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

dat_mean_wide_filtered <- copy(dat_mean_wide)
#dat_mean_wide_filtered <- dat_mean_wide[FC_filter == 1,]
dat_mean_wide_filtered[, chnm := tolower(str_split_fixed(date_chnm_plate, "_", 3)[,2])]

dat_mean_wide_filtered <- dat_mean_wide_filtered[chnm != "dinp branched"]

anova_list_wide_filtered <- merge(anova_list_wide_filtered, dat_mean_wide_filtered[, c(1,3, 70)], by = "date_chnm_plate")

#-----------------------------------------------------------------------------------#
# comparison  
#-----------------------------------------------------------------------------------#
colnames(tox21.aromatase)
spids <- tox21.aromatase$spid # not sure casrn has been captured correctly, just proceed to get casns
chemdb <- tcplLoadChem(field='spid', val=spids)
casns <- chemdb$casn
chemdb.total <- tcplLoadChem(field='casn', val=casns)

hth295r.casn <- unique(anova_list_wide_filtered[,casn]) #104
nvs.casn <- unique(nvs.arom$casn) #248


length(setdiff(casns, hth295r.casn)) #29 of the 50 are not in the HT-H295R multi-conc set
50-29 #21, but this seems wrong because it looks like only 17 unique chems below

length(setdiff(casns, nvs.casn)) #32 of the 50 are not in the nvs multi-conc set
50-32 #18

colnames(anova_list_wide_filtered)

put.arom.h295r <- anova_list_wide_filtered[casn %in% casns]
length(unique(put.arom.h295r$casn)) #17 unique chems, but repeated chemicals
hth295r.e2 <- anova_list_wide_filtered[casn %in% casns & Estradiol==-1 & FC_filter == 1,] 
length(unique(hth295r.e2$casn)) #13 of 17
hth295r.est <- anova_list_wide_filtered[casn %in% casns & Estrone==-1 & FC_filter == 1] #20 of the 21 are shared estrone positives
length(unique(hth295r.est$casn)) #13 of 17
# we miss just propiconazole if we use a combination of e2 and estrone 

put.arom.nvs <- nvs.arom[casn %in% casns]
length(unique(put.arom.nvs$casn)) #14, all positive (triadimenol is repeated 3 times, and 2/3 samples are positive)

setdiff(nvs.casn, hth295r.casn)
setdiff(hth295r.casn, nvs.casn)

setdiff(put.arom.h295r$casn, put.arom.nvs$casn)
setdiff(put.arom.nvs$casn, put.arom.h295r$casn)
casn.intersect <- intersect(put.arom.h295r$casn, put.arom.nvs$casn) #12 casn

hth295r.intersect <- anova_list_wide_filtered[casn %in% casn.intersect]
nvs.intersect <- nvs.arom[casn %in% casn.intersect]

# venn diagram of performance in hth295r (positives from Tox21)
# venn diagram of NVS and hth295r intersect

# venn diagram with VennDiagram (limit at 5 areas for this package)
colnames(put.arom.h295r)
hth295r.venn <- unique(put.arom.h295r[,c('casn', 'date_chnm_plate', 'Estrone', 'Estradiol', 'FC_filter')]) 
hth295r.venn[, Estradiol := as.numeric(Estradiol)]
hth295r.venn[, Estrone := as.numeric(Estrone)]
hth295r.venn <- hth295r.venn[!(str_detect(date_chnm_plate, 'Triadimenol2')),] # reduce sample replication
hth295r.venn <- hth295r.venn[date_chnm_plate != "20140326_Triadimenol_Plate15",] # reduce sample replication, prefer positive result

area1 <- hth295r.venn[Estrone==-1 & FC_filter == 1, casn]

area2 <- hth295r.venn[Estradiol==-1 & FC_filter == 1, casn]

area3 <- hth295r.venn[,casn]

png("./misc/putative_aromataseInhibitors_venn.png", width = 10, height = 10, res = 300, units = "in")
venn(x = list(Estrone = area1, All = area3, Estradiol = area2),
                   lty = "dashed",
                   zcolor='white, green, yellow',
                   opacity = 0.2, cexil=1.5)
dev.off()

# now venn for overlap between hth295r and nvs aromatase

ht.h295r.venn <- hth295r.venn[casn %in% casn.intersect]
area4 <- hth295r.venn[Estrone==-1 & FC_filter == 1 | Estradiol==-1 & FC_filter == 1, casn]
area6 <- hth295r.venn[,casn]

nvs.intersect <- unique(nvs.intersect[,c('casn','chnm','hitc')])
nvs.venn <- nvs.intersect[!(chnm=='Triadimenol' & hitc==0)] # reduce sample replication, prefer positive

area7 <- nvs.venn[hitc==1,casn]
#area8 <- nvs.venn[hitc==0,casn] #there are none that are a hitc == 0

png("./misc/NVS_H295R_overlap_aroInhibitors_venn.png", width = 10, height = 10, res = 300, units = "in")
venn(x = list("HT-H295R Positives" = area4,
              "All HT-H295R"= area6, 
              "NVS Positives" = area7),
     lty = "dashed",
     zcolor='white, green, yellow',
     opacity = 0.2, cexil=1.5)
dev.off()


#-----Generate data table of the chemicals
aromatase_dat <- copy(hth295r.venn)

aromatase_dat <- unique(merge(aromatase_dat, Mahalanobis_dat[, c("date_chnm_plate", "maxD11P")], by = "date_chnm_plate"))

aromatase_dat[, NVS_hit := ifelse(casn %in% nvs.venn[, casn], "Yes", "No")]
aromatase_dat[, H295R_hit := ifelse(Estrone==-1 & FC_filter == 1 | Estradiol==-1 & FC_filter == 1, "Yes", "No")]
aromatase_dat[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]
aromatase_dat[, FC_filter := ifelse(FC_filter == 1, "Yes", "No")]

aromatase_dat <- aromatase_dat[order(chnm),]
aromatase_dat <- aromatase_dat[!(duplicated(chnm)),] #remove duplicates

aromatase_dat <- aromatase_dat[, c(9, 2, 7, 8, 3:5, 6)]

fwrite(aromatase_dat, file = "./misc/Table2_nvs_vs_hth295R.txt", sep = "\t")

# I think there is an error in how venn is handling duplicates





