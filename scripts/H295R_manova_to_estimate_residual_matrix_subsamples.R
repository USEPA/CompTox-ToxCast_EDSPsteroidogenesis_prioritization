#-----------------------------------------------------------------------------------#

#This code will generate the residual matrices from the H295R
#MC data for use in the Mahalanobis distance calculations
#It will subset the data as well to allow for comparison
#of the different covariance matrix partitions


#-----------------NOTES-----------------#
##The loaded H295R master data file has been corrected for LLOQ/sqrt(2) for ND/NQ values and all NR values were removed##
##Similarly, Forskolin and Prochloraz controls were removed##
##rvals have also already been converted to uM##
#---------------------------------------#

rm(list = ls())

library(data.table)
library(multcomp)
library(stringr)
library(reshape2)

#setwd("L:/Lab/NCCT_ToxCast/Derik Haggard/Mahalanobis Follow Up")

#-------------------------------------------------------------------------------------------------#
#-----Functions used in the script
#-------------------------------------------------------------------------------------------------#

Ngood <- function(X) {
  inx <- intersect(colnames(X), hormones_sub)
  x2 <- X[apply(X,1,function(x) all(!is.na(x[inx]))),]
  clF <- unique(x2$classFull)
  clF <- clF[-grep("^DMSO",clF)]
  length(clF)
}

delcols <- function(x) {
  CCols <- intersect(colnames(x), hormones_sub)
  lastNgood <- Ngood(x)
  lasti <- NA
  lastj <- NA
  for (i in CCols)
    for (j in CCols) {
      tmp <- Ngood(x[,-match(unique(c(i,j)),colnames(x)), with = FALSE])
      if (tmp > lastNgood) {
        lastNgood <- tmp
        lasti <- i
        lastj <- j
      }
    }
  list(ngood=lastNgood, dropcols=c(lasti, lastj))
}

#-------------------------------------------------------------------------------------------------#
#-----Load H295R master data file
#-------------------------------------------------------------------------------------------------#

load("./RData/AllResps_outliersRemoved2018-10-02.RData") #load first
load("./RData/H295R_master_table_2017-08-08.RData")

dat <- copy(dat0_combined)


#-----Give DMSO control a chnm


dat[spid == "DMSO", chnm := "DMSO"]


#-----fix Triadimenol to have two chnms reflecting different conc ranges


dat[spid == "TP0000884B12" | spid == "TP0000884A04" | spid == "TP0000883A04", chnm := "Triadimenol2"]


#-----Extract date component of each apid and reformat to a consistent format


dat[, date := str_split_fixed(apid, "\\.", n = 2)[,1]]
dat[, plate := str_split_fixed(apid, "\\.", n = 2)[,2]]
dat[, plate := str_replace(plate, "\\.", "")]

reform <- c("02Apr14" = "20140402",
            "02Apr2014" = "20140402",
            "03Sep2014" = "20140903",
            "09Apr2014" = "20140409",
            "10Jun2015" = "20150610",
            "20130320" = "20130320",
            "20130321" = "20130321",
            "26Mar2014" = "20140326",
            "04112017" = "20170411")

dat[, date := reform[date]]


#-----make chnm_date variable since the date will be the blocking variable for the MANOVA


dat[, chnm := str_replace(chnm, "PharmaGSID_", "PharmaGSID")]

dat[, chnm_date := paste0(chnm, "_", date)]

dat[, date_chnm_plate := paste(date, chnm, plate, sep = "_")]

dat[, chnm_plate := paste0(chnm, "_", plate)]


#-----log uM


dat[, log_uM := log(uM)]


#-----fix rounding errors for chemcial concentrations


dat[conc >= 99.9, conc := 100]

#fix odd Triclosan conc
dat[chnm == "Triclosan" & conc == 0.004115226, conc := 0.004] 

dat[, conc := signif(conc, digits = 2)]
dat[str_detect(date_chnm_plate, "DMSO"), conc_index := "CONC_DMSO"]

#fix cas
dat[, casn := str_replace(casn, "NOCAS_", "NOCAS")]
#setkey(dat, "date_chnm_plate")

#dat[, fixed_conc := signif(conc, digits = 2)]
#conc_count <- dat[, .(sum_conc = length(unique(conc)), sum_fixed_conc = length(unique(fixed_conc))), by = chnm]

#-------------------------------------------------------------------------------------------------#
#-----Select which hormones we want to use
#-------------------------------------------------------------------------------------------------#

hormones_all <- unique(dat[,steroid])

#remove Pregnenolone and DHEA
hormones_sub <- hormones_all[c(-8, -11)]

#First, select upper and lower third of chemicals to sample

Mahalanobis_dat[, date := str_split_fixed(date_chnm_plate, "_", n = 3)[,1]]

block_upper_third <- vector("list", length = length(unique(Mahalanobis_dat[, date])))
names(block_upper_third) <- unique(Mahalanobis_dat[, date])
block_lower_third <- vector("list", length = length(unique(Mahalanobis_dat[, date])))
names(block_lower_third) <- unique(Mahalanobis_dat[, date])

for(block in unique(Mahalanobis_dat[, date])){
  d <- Mahalanobis_dat[date == block]
  d <- d[order(maxD11P),]
  l <- length(d[, maxD11P])/2
  
  lower <- d[1:l,]
  upper <- d[c(max(length(d[, maxD11P]))-l+1):max(length(d[, maxD11P])),]
  
  block_lower_third[[block]] <- lower
  block_upper_third[[block]] <- upper
}

lower_third <- do.call(rbind, block_lower_third)
upper_third <- do.call(rbind, block_upper_third)

lower_third[, c("date", "plate") := .(str_split_fixed(date_chnm_plate, "_", n = 3)[,1], str_split_fixed(date_chnm_plate, "_", n = 3)[,3]),]
upper_third[, c("date", "plate") := .(str_split_fixed(date_chnm_plate, "_", n = 3)[,1], str_split_fixed(date_chnm_plate, "_", n = 3)[,3]),]

#subset into upper and lower and include plate-level DMSO controls

dat_upper <- copy(dat[date_chnm_plate %in% upper_third[, date_chnm_plate] | date_chnm_plate %in% paste(upper_third[, date], "DMSO", upper_third[, plate], sep = "_"),])
dat_lower <- copy(dat[date_chnm_plate %in% lower_third[, date_chnm_plate] | date_chnm_plate %in% paste(lower_third[, date], "DMSO", lower_third[, plate], sep = "_"),])

# bootstrap <- sort(sample(unique(dat_upper[, date_chnm_plate]), replace = TRUE)) #bootstrap method
# 
# test <- dat_upper[bootstrap, allow.cartesian = TRUE] #allows replicated rows in dat
# 
# for(x in unique(test[, date_chnm_plate])){
#   test[date_chnm_plate == x, ]
# }

#MANOVA on all data
dat_all <- copy(dat)


Blocks <- unique(dat_all$date)
Models_all_subsample <- vector("list", length=length(Blocks))
names(Models_all_subsample) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_all[, date])) {
  dt <- copy(dat_all[date == block, c("conc", "conc_index", "steroid", "chnm", "casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
  dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
  dt[, casn := NULL]
  chnm_apids <- unique(dt[,date_chnm_plate])
  chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
  chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
  #     for (i in seq_along(chnm_apids)) {
  #       chnmapid <- chnm_apids[i]
  #       dt[,paste0("Reduced_",i)] <-
  #         factor(ifelse(dt$date_chnm_plate != chnmapid & dt$wllt == "t", paste0(make.names(dt$date_chnm_plate), "_", dt$conc),
  #                       paste0("DMSO_n")))
  #     }
  #dt[, ID := paste0(coli, "_", rowi, "_", chnm)]
  dt[, ID := paste0(coli, "_", rowi, "_", chnm, "_", plate)]
  #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
  #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
  dtw <- dcast.data.table(dt[,c(-7, -8)], ID + ... ~ steroid, value.var = "log_uM")
  colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
  dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
  
  #Bootstrap dtw
#   dtw_list <- vector("list", length = length(unique(dtw[, date_chnm_plate])))
#   
#   for(x in unique(dtw[, date_chnm_plate])){
#     conc_sample <- vector("list", length = length(unique(dtw[date_chnm_plate == x, conc])))
#     names(conc_sample) <- unique(dtw[date_chnm_plate == x, conc_index])
#     for(y in unique(dtw[date_chnm_plate == x, conc_index])){
#       dtw_sample <- (dtw[date_chnm_plate == x & conc_index == y])
#       conc_sample[[y]] <- dtw_sample[sample(nrow(dtw_sample), replace = TRUE),]
#       #conc_sample[[y]][date_chnm_plate == x & conc_index == y, rep := 1:nrow(conc_sample[[y]][date_chnm_plate == x & conc_index == y,])]
#     }
#     dtw_list[[x]] <- do.call(rbind, conc_sample)
#    }
#   
#   dtw <- do.call(rbind, dtw_list)
  #dtw[, ID := paste(rep, date_chnm_plate, sep = "_"),]
  #dtw[, rep := NULL]
  
  ## For the Manova, drop columns that have an excessive number of NA's.
  ## Criterion is to maximize the number of chemxconcxapid levels.
  out <- delcols(dtw)
  if (any(!is.na(out$dropcols))) {
    dtw2 <- dtw[,-match(out$dropcols[!is.na(out$dropcols)], colnames(dtw)), with = FALSE]
    LHS <- paste0("cbind(",paste(paste0("`", intersect(hormones_sub, colnames(dtw2)),"`"), collapse=", "),")")
  } else {
    dtw2 <- dtw
    LHS <- paste0("cbind(",paste(paste0("`", hormones_sub, "`"), collapse=", "),")")
  }
  ## Initial test to find bad Residuals
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw2), silent=TRUE)
  
  #filter dtw to remove singleton samples and outlier paired samples with standard deviation of the Residuals > 1
  inx <- intersect(colnames(dtw2), hormones_sub)
  x2 <- dtw2[apply(dtw2,1,function(x) all(!is.na(x[inx]))),] #NAs are ommitted and allows for proper cbind of Residuals
  residual_list <- cbind(data.table(zz[[1]]$residuals), ID = x2[,ID], date_chnm_plate = x2[,date_chnm_plate], conc = x2[,conc])
  residual_list[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  
  #count up the number of replicates for each conc_date_chnm_plate
  for(x in unique(residual_list[, conc_date_chnm_plate])){
    residual_list[conc_date_chnm_plate == x, reps := residual_list[conc_date_chnm_plate == x, .N]]
  }
  
  residual_list_filtered <- residual_list[reps > 1,] #remove singleton samples
  
  #calculate standard deviation of the Residuals
  residual_diff_data <- residual_list_filtered[, lapply(.SD, function(x) sd(x)),
                                               by = .(conc_date_chnm_plate, date_chnm_plate),
                                               .SDcols = colnames(residual_list_filtered)[-c((length(residual_list_filtered)-4):length(residual_list_filtered))]]
  
  
  residual_diff_data[, maxSD := ifelse(apply(.SD, 1, function(x) max(x)) > 1, 1, 0), .SDcols = colnames(residual_diff_data)[3:length(residual_diff_data)]]
  
  Residuals[[block]] <- residual_diff_data
  
  #make new filtered dtw data table
  dtw3 <- copy(dtw2)
  dtw3[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  dtw3 <- dtw3[conc_date_chnm_plate %in% residual_diff_data[maxSD == 0, conc_date_chnm_plate],] #filter out bad data
  dtw3[, conc_date_chnm_plate := NULL]
  
  ## Real Tests
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw3), silent=TRUE)
  
  Models_all_subsample[[block]] <- list(fit = zz[[1]],
                              Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                              IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#MANOVA on upper data
dat_upper <- copy(dat_upper)


Blocks <- unique(dat_upper$date)
Models_upper_subsample <- vector("list", length=length(Blocks))
names(Models_upper_subsample) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_upper[, date])) {
  dt <- copy(dat_upper[date == block, c("conc", "conc_index", "steroid", "chnm", "casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
  dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
  dt[, casn := NULL]
  chnm_apids <- unique(dt[,date_chnm_plate])
  chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
  chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
  #     for (i in seq_along(chnm_apids)) {
  #       chnmapid <- chnm_apids[i]
  #       dt[,paste0("Reduced_",i)] <-
  #         factor(ifelse(dt$date_chnm_plate != chnmapid & dt$wllt == "t", paste0(make.names(dt$date_chnm_plate), "_", dt$conc),
  #                       paste0("DMSO_n")))
  #     }
  #dt[, ID := paste0(coli, "_", rowi, "_", chnm)]
  dt[, ID := paste0(coli, "_", rowi, "_", chnm, "_", plate)]
  #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
  #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
  dtw <- dcast.data.table(dt[,c(-7, -8)], ID + ... ~ steroid, value.var = "log_uM")
  colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
  dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
  
  ## For the Manova, drop columns that have an excessive number of NA's.
  ## Criterion is to maximize the number of chemxconcxapid levels.
  out <- delcols(dtw)
  if (any(!is.na(out$dropcols))) {
    dtw2 <- dtw[,-match(out$dropcols[!is.na(out$dropcols)], colnames(dtw)), with = FALSE]
    LHS <- paste0("cbind(",paste(paste0("`", intersect(hormones_sub, colnames(dtw2)),"`"), collapse=", "),")")
  } else {
    dtw2 <- dtw
    LHS <- paste0("cbind(",paste(paste0("`", hormones_sub, "`"), collapse=", "),")")
  }
  ## Initial test to find bad Residuals
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw2), silent=TRUE)
  
  #filter dtw to remove singleton samples and outlier paired samples with standard deviation of the Residuals > 1
  inx <- intersect(colnames(dtw2), hormones_sub)
  x2 <- dtw2[apply(dtw2,1,function(x) all(!is.na(x[inx]))),] #NAs are ommitted and allows for proper cbind of Residuals
  residual_list <- cbind(data.table(zz[[1]]$residuals), ID = x2[,ID], date_chnm_plate = x2[,date_chnm_plate], conc = x2[,conc])
  residual_list[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  
  #count up the number of replicates for each conc_date_chnm_plate
  for(x in unique(residual_list[, conc_date_chnm_plate])){
    residual_list[conc_date_chnm_plate == x, reps := residual_list[conc_date_chnm_plate == x, .N]]
  }
  
  residual_list_filtered <- residual_list[reps > 1,] #remove singleton samples
  
  #calculate standard deviation of the Residuals
  residual_diff_data <- residual_list_filtered[, lapply(.SD, function(x) sd(x)),
                                               by = .(conc_date_chnm_plate, date_chnm_plate),
                                               .SDcols = colnames(residual_list_filtered)[-c((length(residual_list_filtered)-4):length(residual_list_filtered))]]
  
  
  residual_diff_data[, maxSD := ifelse(apply(.SD, 1, function(x) max(x)) > 1, 1, 0), .SDcols = colnames(residual_diff_data)[3:length(residual_diff_data)]]
  
  Residuals[[block]] <- residual_diff_data
  
  #make new filtered dtw data table
  dtw3 <- copy(dtw2)
  dtw3[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  dtw3 <- dtw3[conc_date_chnm_plate %in% residual_diff_data[maxSD == 0, conc_date_chnm_plate],] #filter out bad data
  dtw3[, conc_date_chnm_plate := NULL]
  
  ## Real Tests
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw3), silent=TRUE)
  
  Models_upper_subsample[[block]] <- list(fit = zz[[1]],
                                Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#MANOVA on all data
dat_lower <- copy(dat_lower)


Blocks <- unique(dat_lower$date)
Models_lower_subsample <- vector("list", length=length(Blocks))
names(Models_lower_subsample) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_lower[, date])) {
  dt <- copy(dat_lower[date == block, c("conc", "conc_index", "steroid", "chnm","casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
  dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
  dt[, casn := NULL]
  chnm_apids <- unique(dt[,date_chnm_plate])
  chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
  chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
  #     for (i in seq_along(chnm_apids)) {
  #       chnmapid <- chnm_apids[i]
  #       dt[,paste0("Reduced_",i)] <-
  #         factor(ifelse(dt$date_chnm_plate != chnmapid & dt$wllt == "t", paste0(make.names(dt$date_chnm_plate), "_", dt$conc),
  #                       paste0("DMSO_n")))
  #     }
  #dt[, ID := paste0(coli, "_", rowi, "_", chnm)]
  dt[, ID := paste0(coli, "_", rowi, "_", chnm, "_", plate)]
  #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
  #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
  dtw <- dcast.data.table(dt[,c(-7, -8)], ID + ... ~ steroid, value.var = "log_uM")
  colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
  dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
  
  ## For the Manova, drop columns that have an excessive number of NA's.
  ## Criterion is to maximize the number of chemxconcxapid levels.
  out <- delcols(dtw)
  if (any(!is.na(out$dropcols))) {
    dtw2 <- dtw[,-match(out$dropcols[!is.na(out$dropcols)], colnames(dtw)), with = FALSE]
    LHS <- paste0("cbind(",paste(paste0("`", intersect(hormones_sub, colnames(dtw2)),"`"), collapse=", "),")")
  } else {
    dtw2 <- dtw
    LHS <- paste0("cbind(",paste(paste0("`", hormones_sub, "`"), collapse=", "),")")
  }
  ## Initial test to find bad Residuals
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw2), silent=TRUE)
  
  #filter dtw to remove singleton samples and outlier paired samples with standard deviation of the Residuals > 1
  inx <- intersect(colnames(dtw2), hormones_sub)
  x2 <- dtw2[apply(dtw2,1,function(x) all(!is.na(x[inx]))),] #NAs are ommitted and allows for proper cbind of Residuals
  residual_list <- cbind(data.table(zz[[1]]$residuals), ID = x2[,ID], date_chnm_plate = x2[,date_chnm_plate], conc = x2[,conc])
  residual_list[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  
  #count up the number of replicates for each conc_date_chnm_plate
  for(x in unique(residual_list[, conc_date_chnm_plate])){
    residual_list[conc_date_chnm_plate == x, reps := residual_list[conc_date_chnm_plate == x, .N]]
  }
  
  residual_list_filtered <- residual_list[reps > 1,] #remove singleton samples
  
  #calculate standard deviation of the Residuals
  residual_diff_data <- residual_list_filtered[, lapply(.SD, function(x) sd(x)),
                                               by = .(conc_date_chnm_plate, date_chnm_plate),
                                               .SDcols = colnames(residual_list_filtered)[-c((length(residual_list_filtered)-4):length(residual_list_filtered))]]
  
  
  residual_diff_data[, maxSD := ifelse(apply(.SD, 1, function(x) max(x)) > 1, 1, 0), .SDcols = colnames(residual_diff_data)[3:length(residual_diff_data)]]
  
  Residuals[[block]] <- residual_diff_data
  
  #make new filtered dtw data table
  dtw3 <- copy(dtw2)
  dtw3[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
  dtw3 <- dtw3[conc_date_chnm_plate %in% residual_diff_data[maxSD == 0, conc_date_chnm_plate],] #filter out bad data
  dtw3[, conc_date_chnm_plate := NULL]
  
  ## Real Tests
  zz <- list()
  Formula <- as.formula(paste0(LHS," ~ 0 + classFull"))
  zz[[1]] <- try(lm(Formula, data=dtw3), silent=TRUE)
  
  Models_lower_subsample[[block]] <- list(fit = zz[[1]],
                                Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#Save data
save(Models_all_subsample, Models_upper_subsample, Models_lower_subsample, file = paste0("./RData/subsample_Model_output", Sys.Date(), ".RData"))
