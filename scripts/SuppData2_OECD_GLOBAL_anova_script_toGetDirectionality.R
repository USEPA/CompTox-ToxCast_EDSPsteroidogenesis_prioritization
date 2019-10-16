#-----------------------------------------------------------------------------------#

# This script with perform the ANOVA analysis fo the H295R data similar to Hecker et al., 2011.#
#but with added directionality that is used to define putative aromatase inhibitors#

#updated 09/21/2018

rm(list=ls())

library(data.table)
library(reshape2)
library(stringr)
library(tcpl)
library(multcomp)
library(ggplot2)
library(ggthemes)
library(gvlma)
library(gridExtra)
library(grid)

#setwd("./ANOVA/Global ANOVA")

#-------------------------------------------------------------------------------------------------#
#-----Load in H295R master data table
#-------------------------------------------------------------------------------------------------#

load("./RData/H295R_master_table_2017-08-08.RData")
#load("./Global_ANOVA_06202017.RData")

#-------------------------------------------------------------------------------------------------#
#-----Clean up data
#-------------------------------------------------------------------------------------------------#

dat0_combined[spid == "TP0000884B12" | spid == "TP0000884A04" | spid == "TP0000883A04", chnm := "Triadimenol2"]

dat0_combined[, date := str_split_fixed(apid, "\\.", n = 2)[,1]]
dat0_combined[, plate := str_split_fixed(apid, "\\.", n = 2)[,2]]
dat0_combined[, plate := str_replace(plate, "\\.", "")]
dat0_combined[, chnm := str_replace(chnm, "PharmaGSID_", "PharmaGSID")]

reform <- c("02Apr14" = "20140402",
            "02Apr2014" = "20140402",
            "03Sep2014" = "20140903",
            "09Apr2014" = "20140409",
            "10Jun2015" = "20150610",
            "20130320" = "20130320",
            "20130321" = "20130321",
            "26Mar2014" = "20140326",
            "04112017" = "20170411")

dat0_combined[, date := reform[date]]

dat0_combined[, date_chnm_plate := paste(date, chnm, plate, sep = "_")]
dat0_combined[, date_chnm_plate := str_replace(date_chnm_plate, "_NA_", "_DMSO_")]

#-------------------------------------------------------------------------------------------------#
#-----fix rounding errors for chemcial concentrations
#-------------------------------------------------------------------------------------------------#

dat0_combined[conc >= 99.9, conc := 100]

#fix odd Triclosan conc
dat0_combined[chnm == "Triclosan" & conc == 0.004115226, conc := 0.004] 

dat0_combined[, conc := signif(conc, digits = 2)]

#-------------------------------------------------------------------------------------------------#
#-----Drop DHEA and Pregnenolone
#-------------------------------------------------------------------------------------------------#

dat0_combined <- dat0_combined[steroid != "DHEA" & steroid != "Pregnenolone"]
 
#-------------------------------------------------------------------------#
#-----OPTIONAL: Filter the samples with residual outliers-----------------#
#-------------------------------------------------------------------------#
# dat0_combined[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
# 
# dat0_combined <- dat0_combined[!c(conc_date_chnm_plate %in% Residuals[maxSD == 1, conc_date_chnm_plate]),] #filter based on Residuals with the maxSD filter == 1
# dat0_combined[, conc_date_chnm_plate := NULL]

#-------------------------------------------------------------------------------------------------#
#-----factor steroid analytes
#-------------------------------------------------------------------------------------------------#

dat0_combined[,steroid := factor(steroid, levels = c("OH-Pregnenolone", "Progesterone", "OH-Progesterone",
                                                     "DOC", "Corticosterone", "11-deoxycortisol", "Cortisol",
                                                     "Androstenedione", "Testosterone", "Estrone", "Estradiol"))]

#-------------------------------------------------------------------------------------------------#
#-----FOR loop structure:
#1. for each casNO
##2. for each hormone
###3. for each apid (include the plate DMSO for each loop)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
#-----Output ANOVA data into a table
#-------------------------------------------------------------------------------------------------#

out <- vector("list", length = length(unique(dat0_combined[spid != "DMSO", date_chnm_plate])))
names(out) <- unique(dat0_combined[spid != "DMSO", date_chnm_plate])
out2 <- vector("list", length = length(unique(dat0_combined[spid != "DMSO", date_chnm_plate])))
names(out2) <- unique(dat0_combined[spid != "DMSO", date_chnm_plate])
for(x in unique(dat0_combined[order(chnm, steroid) & spid != "DMSO",date_chnm_plate])){
  steroid_list <- vector("list", length = 1+length(unique(dat0_combined[spid != "DMSO", steroid])))
  names(steroid_list) <- c("date_chnm_plate", as.character(levels(unique(dat0_combined[chnm != "DMSO" ,steroid]))))
  steroid_list2 <- vector("list", length = 1+length(unique(dat0_combined[spid != "DMSO", steroid])))
  names(steroid_list2) <- c("date_chnm_plate", as.character(levels(unique(dat0_combined[chnm != "DMSO" ,steroid]))))
  for(y in levels(unique(dat0_combined[,steroid]))){
    for(z in unique(dat0_combined[date_chnm_plate == x, apid])){
      fit_data <- dat0_combined[date_chnm_plate == x | spid == "DMSO",]
      check_data <- fit_data
      fit_data <- dat0_combined[date_chnm_plate == x & steroid == y | spid == "DMSO" & steroid == y,]
      fit_data <- fit_data[apid == z,]
      
#       if(length(unique(fit_data[, conc])) < 7){
#         empty_row <- fit_data[1,]
#         empty_row[,c("conc", "uM") := .(100000, -5)]
#         fit_data <- rbind(fit_data, empty_row)
#       }
      
      fit_data[, Fconc := factor(conc)]
      
      chemname <- unique(fit_data[wllt == "t", chnm])
      
      ft <- lm(uM ~ Fconc, data = fit_data)
      sum_ft <- summary(ft)
      
      if(!is.na(sum_ft$fstatistic[1])) {
        set.seed(123456789)
        anova_ft <- glht(ft, linfct = mcp(Fconc = "Dunnett"))
        sum_anova <- summary(anova_ft)
        if(length(unique(fit_data[, Fconc])) < length(unique(check_data[, conc]))){
          df <- c(sum_anova$test$pvalues, NA)
          names(df) <- c(sort(unique(fit_data[, conc]))[-1], unique(check_data[!(conc %in% fit_data[,conc]),conc]))
          df <- ifelse(df <= 0.05, 1, 0)
          de <- c(sum_anova$test$coefficients, NA) #slope estimates for directionality
          names(de) <- c(sort(unique(fit_data[, conc]))[-1], unique(check_data[!(conc %in% fit_data[,conc]),conc]))
          de <- ifelse(de > 0, 1, -1)
          de[df == 0] <- 0
        } else{
          df <- c(sum_anova$test$pvalues)
          names(df) <- c(sort(unique(fit_data[, conc]))[-1])
          df <- ifelse(df <= 0.05, 1, 0)
          de <- c(sum_anova$test$coefficients)
          names(de) <- c(sort(unique(fit_data[, conc]))[-1])
          de <- ifelse(de > 0, 1, -1)
          de[df == 0] <- 0 
        }
        steroid_list[[y]] <- df
        steroid_list2[[y]] <- de
      } else{
        df <- c(NA, NA, NA, NA, NA, NA)
        steroid_list[[y]] <- df
        de <- c(NA, NA, NA, NA, NA, NA)
        steroid_list2[[y]] <- de
      }
      steroid_list[[1]] <- x
      steroid_list2[[1]] <- x
      #steroid_list[[2]] <- data.frame(conc = sort(unique(fit_data[,Fconc]))[-1])
    }
  }
  steroid_list <- data.frame(as.character(steroid_list)[-1], row.names = names(steroid_list)[-1])
  colnames(steroid_list) <- x
  out[[x]] <- steroid_list
  steroid_list2 <- data.frame(as.character(steroid_list2)[-1], row.names = names(steroid_list2)[-1])
  colnames(steroid_list2) <- x
  out2[[x]] <- steroid_list2
  #break
}

ANOVA_out <- data.table(do.call("cbind", out))
row.names(ANOVA_out) <- levels(unique(dat0_combined[,steroid]))
ANOVA_out2 <- data.table(do.call("cbind", out2))
row.names(ANOVA_out2) <- levels(unique(dat0_combined[,steroid]))

t_ANOVA_out <- data.table(t(ANOVA_out), keep.rownames= TRUE)
colnames(t_ANOVA_out) <- c("date_chnm_plate", rownames(ANOVA_out))
t_ANOVA_out[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
t_ANOVA_out[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]
t_ANOVA_out2 <- data.table(t(ANOVA_out2), keep.rownames= TRUE)
colnames(t_ANOVA_out2) <- c("date_chnm_plate", rownames(ANOVA_out2))
t_ANOVA_out2[, date_chnm_plate := str_replace(date_chnm_plate, "PharmaGSID_", "PharmaGSID")]
t_ANOVA_out2[, chnm := str_split_fixed(date_chnm_plate, "_", 3)[,2]]

for(x in unique(t_ANOVA_out[, chnm])){
  t_ANOVA_out[chnm == x, casn := unique(dat0_combined[chnm == x, casn])]
  t_ANOVA_out2[chnm == x, casn := unique(dat0_combined[chnm == x, casn])]
}

#fwrite(t_ANOVA_out, file = paste0("Global_H295R_ANOVA_OECD_strings_filtered_output_", Sys.Date(), ".txt"), sep = "\t")
fwrite(t_ANOVA_out2, file = paste0("./supplement/Supp2_Global_H295R_ANOVA_OECD_directionality_filtered_output_", Sys.Date(), ".txt"), sep = "\t")
