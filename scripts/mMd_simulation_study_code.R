#-----------------------------------------------------------------------------------#

#Code to resample data from multivariate normal distribution#
#This generates all of the simulated data used thorughout the paper#

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

rm(list = ls())

#setwd("/share/home/dhaggard/Desktop/Mahalanobis_Follow_Up") #adjust accordingly

#-------------------------------------------------------------------------#
#-----Functions needed for Manova
#-------------------------------------------------------------------------#
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

#-------------------------------------------------------------------------#
#-----Functions needed for mMd calculations
#-------------------------------------------------------------------------#
PT2maxmv  <-  function(c, K, p) {
  fn <- function(v, c, p, K) {
    exp((K - 1) * pchisq(c, df=p, ncp=v, log.p=TRUE)) *
      dchisq(v, p)
  }
  tmp <- integrate(fn, 0, Inf, c=c, p=p, K=K, rel.tol=1e-12)
  tmp$value
}

critval <- function(alpha, K, p) {
  if (K < 2 || p < 1) {
    list(root=NA)
  } else {
    ##    Find c in PT2maxmv such that PT2maxmv == 1 - alpha
    obfun <- function(c) {
      PT2maxmv(c, K=K, p=p) - (1 - alpha)
    }
    ## When c is small, obfun < 0, increases with c.
    ## Start with c=1
    cc <- 1
    Z <- obfun(cc)
    Niter <- 0
    while ((Z <- obfun(cc)) < 0 && ((Niter <- Niter + 1) < 100)) cc <- cc * 2
    if (Z > 0) {
      out <- uniroot(obfun, interval=c(cc/2,cc))
      out
    } else {
      list(root=NA)
    }
  }
}

#-------------------------------------------------------------------------#
#-----Define LLOQ values to aid simulation
#-------------------------------------------------------------------------#
MW <- numeric(13)
names(MW) <- c("11-deoxycortisol", "Androstenedione", "Corticosterone", 
               "Cortisol", "DHEA", "DOC", "Estradiol", "Estrone", 
               "OH-Pregnenolone", "OH-Progesterone", "Pregnenolone",
               "Progesterone", "Testosterone")
MW["Pregnenolone"] <- 316.4776
MW["Progesterone"] <- 314.46
MW["DOC"] <- 330.4611
MW["Corticosterone"] <- 346.461
MW["OH-Pregnenolone"] <- 332.477
MW["OH-Progesterone"] <- 330.4611
MW["11-deoxycortisol"] <- 346.4605
MW["Cortisol"] <- 362.460
MW["DHEA"] <- 288.424
MW["Androstenedione"] <- 286.409
MW["Estrone"] <- 270.366
MW["Testosterone"] <- 288.424
MW["Estradiol"] <- 272.382

LOQ <- data.frame(acsn = sort(unique(names(MW))),
                  LLOQ = numeric(13),
                  ULOQ = numeric(13))

LOQ[1,2:3] <- c(5,1000)
LOQ[2,2:3] <- c(1,200)
LOQ[3,2:3] <- c(0.5,100)
LOQ[4,2:3] <- c(0.5,100)
LOQ[5,2:3] <- c(3,600)
LOQ[6,2:3] <- c(0.5,100)
LOQ[7,2:3] <- c(0.03,6)
LOQ[8,2:3] <- c(0.03,6)
LOQ[9,2:3] <- c(5,1000)
LOQ[10,2:3] <- c(0.2,40)
LOQ[11,2:3] <- c(2,400)
LOQ[12,2:3] <- c(0.2,40)
LOQ[13,2:3] <- c(0.1,20)

#-----convert to uM
LOQuM <- LOQ[2]
rownames(LOQuM) <- as.character(LOQ$acsn)
LOQuM <- data.table(steroid = rownames(LOQuM), sweep(LOQuM, 1, 1/MW, "*")) #data with no detection will be LOQuM/sqrt(2)

LLOQuM <- as.numeric(t(LOQuM[, 2]))
names(LLOQuM) <- as.character(LOQ$acsn)

LLOQuM <- log(LLOQuM/sqrt(2))

#-------------------------------------------------------------------------#
#-----Step 1: Partition original data based on maxmMd and calc mean
#-----response by chemical
#-------------------------------------------------------------------------#
#load data
load("./RData/AllResps_outliersRemoved2018-10-02.RData")
load("./RData/updated_manova_output2018-02-06.RData")
load("./RData/CovTP0_2018-01-31.RData")
load("./RData/subsample_Model_output2018-10-03.RData")

#-----fix dat
dat <- copy(dat)

dat[spid == "DMSO", chnm := "DMSO"]
dat[spid == "DMSO", conc_index := "CONC_DMSO"]
dat[, conc_index := factor(conc_index, levels = c("CONC_DMSO", "CONC6", "CONC5", "CONC4", "CONC3", "CONC2", "CONC1"))]

dat[, conc := signif(conc, digits = 2)]

dat[, casn := str_replace(casn, "NOCAS_", "NOCAS")]

dat[, N := .N, by = .(date_chnm_plate, steroid, conc, conc_index)]

#filter out samples with bad residuals
dat[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]
dat <- dat[!c(conc_date_chnm_plate %in% Residuals[maxSD == 1, conc_date_chnm_plate]),] #filter based on Residuals with the maxSD filter == 1
dat[, conc_date_chnm_plate := NULL]

#-----average responses across replicates
dat_mean <- dat[, .(mean_log_uM = mean(log_uM)), by = .(date_chnm_plate, steroid, chnm, casn, conc, conc_index, wllt, date, N)]
dat_mean <- dat_mean[steroid != "DHEA" & steroid != "Pregnenolone",]

#remove estrone and estradiol from block 20140326
# dat_mean <- dat_mean[, filter := ifelse(date == "20140326" & steroid == "Estrone" | date == "20140326" & steroid == "Estradiol", 0, 1)]
# dat_mean <- dat_mean[filter == 1,]
# dat_mean[, filter := NULL]

for(x in unique(dat_mean[,date_chnm_plate])){
  for(y in unique(dat_mean[date_chnm_plate == x, conc_index])){
    if(any(dat_mean[date_chnm_plate == x & conc_index == y, N] == 1)){
      dat_mean[date_chnm_plate == x & conc_index == y, N := 1]
    }  else if(all(dat_mean[date_chnm_plate == x & conc_index == y, N] == 3)){
      dat_mean[date_chnm_plate == x & conc_index == y, N := 3]
    } else if(all(dat_mean[date_chnm_plate == x & conc_index == y, N] == 4)){
      dat_mean[date_chnm_plate == x & conc_index == y, N := 4]
    }
  }
}

dat_mean_wide <- dcast.data.table(dat_mean[,], formula = ... ~ steroid, value.var = "mean_log_uM")
dat_mean_wide[chnm == "DMSO", casn := "DMSO"]

dat_mean_wide[, ID := paste(date_chnm_plate, casn, conc, conc_index, sep = "_")]

#-----partition data into upper/mid/lower/dmso subsets based on maxmMd values from original analysis
Mahalanobis_dat[, date := str_split_fixed(date_chnm_plate, "_", n = 3)[,1]]

block_upper <- vector("list", length = length(unique(Mahalanobis_dat[, date])))
names(block_upper) <- unique(Mahalanobis_dat[, date])

block_mid <- vector("list", length = length(unique(Mahalanobis_dat[, date])))
names(block_mid) <- unique(Mahalanobis_dat[, date])

block_lower <- vector("list", length = length(unique(Mahalanobis_dat[, date])))
names(block_lower) <- unique(Mahalanobis_dat[, date])

for(block in unique(Mahalanobis_dat[, date])){
  d <- Mahalanobis_dat[date == block]
  d <- d[order(maxD11P),]
  l <- round(length(d[, maxD11P])/3)
  
  lower <- d[1:l,]
  mid <- d[(l+1):(2*l),]
  upper <- d[c((2*l)+1):max(length(d[, maxD11P])),]
  
  block_lower[[block]] <- lower
  block_mid[[block]] <- mid
  block_upper[[block]] <- upper
}

lower <- do.call(rbind, block_lower)
mid <- do.call(rbind, block_mid)
upper <- do.call(rbind, block_upper)

lower[, c("date", "plate") := .(str_split_fixed(date_chnm_plate, "_", n = 3)[,1], str_split_fixed(date_chnm_plate, "_", n = 3)[,3]),]
mid[, c("date", "plate") := .(str_split_fixed(date_chnm_plate, "_", n = 3)[,1], str_split_fixed(date_chnm_plate, "_", n = 3)[,3]),]
upper[, c("date", "plate") := .(str_split_fixed(date_chnm_plate, "_", n = 3)[,1], str_split_fixed(date_chnm_plate, "_", n = 3)[,3]),]

#-----subset data: high, low, and DMSO
dat_upper <- copy(dat[date_chnm_plate %in% upper[, date_chnm_plate],])
dat_mid <- copy(dat[date_chnm_plate %in% mid[, date_chnm_plate],])
dat_lower <- copy(dat[date_chnm_plate %in% lower[, date_chnm_plate],])
dat_dmso <- copy(dat[chnm == "DMSO",])

#-------------------------------------------------------------------------#
#-----Step 2: Rum manova on partitions to generate pooled covariance matrices
#-------------------------------------------------------------------------#

#-----select hormones we care about
hormones_all <- unique(dat[,steroid])
hormones_sub <- sort((hormones_all[c(-8, -11)]))

#-----Manova on upper subset
Blocks <- unique(dat_upper$date)
Models_upper <- vector("list", length=length(Blocks))
names(Models_upper) <- Blocks

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
  
  Models_upper[[block]] <- list(fit = zz[[1]],
                                Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#-----Manova on middle subset
Blocks <- unique(dat_mid$date)
Models_mid <- vector("list", length=length(Blocks))
names(Models_mid) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_mid[, date])) {
  dt <- copy(dat_mid[date == block, c("conc", "conc_index", "steroid", "chnm", "casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
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
  
  Models_mid[[block]] <- list(fit = zz[[1]],
                              Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                              IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#-----Manova on lower subset
Blocks <- unique(dat_lower$date)
Models_lower <- vector("list", length=length(Blocks))
names(Models_lower) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_lower[, date])) {
  dt <- copy(dat_lower[date == block, c("conc", "conc_index", "steroid", "chnm", "casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
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
  
  Models_lower[[block]] <- list(fit = zz[[1]],
                                Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#-----Manova on dmso subset
Blocks <- unique(dat_dmso$date)
Models_dmso <- vector("list", length=length(Blocks))
names(Models_dmso) <- Blocks

Residuals <- vector("list", length=length(Blocks))
names(Residuals) <- Blocks

for (block in unique(dat_dmso[, date])) {
  dt <- copy(dat_dmso[date == block, c("conc", "conc_index", "steroid", "chnm", "casn", "log_uM", "date", "coli", "rowi", "wllt","date_chnm_plate", "plate"), with = FALSE])
  
  dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
  dt[, casn := NULL]
  chnm_apids <- unique(dt[,date_chnm_plate])
  #chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
  #chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
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
  Formula <- as.formula(paste0(LHS," ~ 0 + date_chnm_plate"))
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
  Formula <- as.formula(paste0(LHS," ~ 0 + date_chnm_plate"))
  zz[[1]] <- try(lm(Formula, data=dtw3), silent=TRUE)
  
  Models_dmso[[block]] <- list(fit = zz[[1]],
                               Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                               IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
} 

#-----Calculate Pooled covariance matrices for subsamples
Covs_upper <- vector("list", length = 8)
names(Covs_upper) <- names(Models_upper)

Covs_mid <- vector("list", length = 8)
names(Covs_mid) <- names(Models_mid)

Covs_lower <- vector("list", length = 8)
names(Covs_lower) <- names(Models_lower)

Covs_dmso <- vector("list", length = 8)
names(Covs_dmso) <- names(Models_dmso)


#estimate block-level covariance
for(block in names(Models)){
  Covs_upper[[block]] <- estVar(Models_upper[[block]]$fit)
  Covs_mid[[block]] <- estVar(Models_mid[[block]]$fit)
  Covs_lower[[block]] <- estVar(Models_lower[[block]]$fit)
  Covs_dmso[[block]] <- estVar(Models_dmso[[block]]$fit)
}

#estimate pooled covariance matrices
CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_upper)))
for (blk in names(Models_upper)) CovT[rownames(Covs_upper[[blk]]),colnames(Covs_upper[[blk]]),blk] <- Covs_upper[[blk]]
CovTP_upper <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_mid)))
for (blk in names(Models_mid)) CovT[rownames(Covs_mid[[blk]]),colnames(Covs_mid[[blk]]),blk] <- Covs_mid[[blk]]
CovTP_mid <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_lower)))
for (blk in names(Models_lower)) CovT[rownames(Covs_lower[[blk]]),colnames(Covs_lower[[blk]]),blk] <- Covs_lower[[blk]]
CovTP_lower <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_dmso)))
for (blk in names(Models_dmso)) CovT[rownames(Covs_dmso[[blk]]),colnames(Covs_dmso[[blk]]),blk] <- Covs_dmso[[blk]]
CovTP_dmso <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

#-----Calculate Pooled covariance matrices for upper/lower subsamples of original data
Covs_original_upper <- vector("list", length = 8)
names(Covs_original_upper) <- names(Models_upper_subsample)

Covs_original_lower <- vector("list", length = 8)
names(Covs_original_lower) <- names(Models_lower_subsample)

#estimate block-level covariance
for(block in names(Models)){
  Covs_original_upper[[block]] <- estVar(Models_upper_subsample[[block]]$fit)
  Covs_original_lower[[block]] <- estVar(Models_lower_subsample[[block]]$fit)
}

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_upper_subsample)))
for (blk in names(Models_upper_subsample)) CovT[rownames(Covs_original_upper[[blk]]),colnames(Covs_original_upper[[blk]]),blk] <- Covs_original_upper[[blk]]
CovTP_original_upper <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_lower_subsample)))
for (blk in names(Models_lower_subsample)) CovT[rownames(Covs_original_lower[[blk]]),colnames(Covs_original_lower[[blk]]),blk] <- Covs_original_lower[[blk]]
CovTP_original_lower <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

#################################################################################################################################
#################################################################################################################################
#-----Begin loop!
#################################################################################################################################
#################################################################################################################################

#-----Define important R objects

simulated_mMd_output <- vector("list", length = 2000)

simulated_distance_function <- function(z){
  
  #-------------------------------------------------------------------------#
  #-----Step 3: MVN sampling
  #-------------------------------------------------------------------------#
  
  #-----Upper
  dat_mean_upper <- dat_mean_wide[date_chnm_plate %in% unique(dat_upper[, date_chnm_plate])]
  
  mvn_dat_upper <- vector("list", length = 8)
  names(mvn_dat_upper) <- unique(dat[, date])
  
  for(block in names(Models)){
    mvn_dat_upper_block <- vector("list", length = nrow(dat_mean_upper[date == block,]))
    names(mvn_dat_upper_block) <- dat_mean_upper[date == block, ID]
    
    for(sample in names(mvn_dat_upper_block)){
      means <- as.numeric(dat_mean_upper[ID == sample, c(9:19), with = FALSE])
      names(means) <- colnames(dat_mean_upper[ID == sample, c(9:19), with = FALSE])
      #naCols <- names(means[is.na(means) == TRUE])
      #means[is.na(means) == TRUE] <- 1000
      
      mvn <- data.table(ID = dat_mean_upper[ID == sample, ID],
                        rmvnorm(n = dat_mean_upper[ID == sample, N], mean = means, sigma = Covs_upper[[block]]),
                        rep = 1:dat_mean_upper[ID == sample, N], N = dat_mean_upper[ID == sample, N])
      #mvn[, c(naCols) := NA]
      mvn_dat_upper_block[[sample]] <- mvn
      
    }
    mvn_dat_upper[[block]] <- do.call(rbind, mvn_dat_upper_block)
    
    #-----clean up IDs and add date_chnm_plate, conc, and other info
    mvn_dat_upper[[block]][, date := str_split_fixed(ID, "_", n = 6)[,1], ]
    mvn_dat_upper[[block]][, chnm := str_split_fixed(ID, "_", n = 6)[,2], ]
    mvn_dat_upper[[block]][, plate := str_split_fixed(ID, "_", n = 6)[,3], ]
    mvn_dat_upper[[block]][, casn := str_split_fixed(ID, "_", n = 6)[,4], ]
    mvn_dat_upper[[block]][, conc := str_split_fixed(ID, "_", n = 6)[,5], ]
    mvn_dat_upper[[block]][, conc_index := str_split_fixed(ID, "_", n = 6)[,6], ]
    mvn_dat_upper[[block]][, date_chnm_plate := paste(str_split_fixed(ID, "_", n = 6)[,1], str_split_fixed(ID, "_", n = 6)[,2],
                                                      str_split_fixed(ID, "_", n = 6)[,3], sep = "_"), ]
    #mvn_dat_upper[[block]][ID == "DMSO_n",date_chnm_plate := paste(block, "DMSO", "Plate", sep = "_")]
    mvn_dat_upper[[block]][, ID := NULL]
    mvn_dat_upper[[block]][, N := .N, by = .(date_chnm_plate, conc, casn, date)]
  }
  
  mvn_dat_upper <- do.call(rbind, c(mvn_dat_upper, list(fill = TRUE)))
  
  #-----Mid
  dat_mean_mid <- dat_mean_wide[date_chnm_plate %in% unique(dat_mid[, date_chnm_plate])]
  
  mvn_dat_mid <- vector("list", length = 8)
  names(mvn_dat_mid) <- unique(dat[, date])
  
  for(block in names(Models)){
    mvn_dat_mid_block <- vector("list", length = nrow(dat_mean_mid[date == block,]))
    names(mvn_dat_mid_block) <- dat_mean_mid[date == block, ID]
    
    for(sample in names(mvn_dat_mid_block)){
      means <- as.numeric(dat_mean_mid[ID == sample, c(9:19), with = FALSE])
      names(means) <- colnames(dat_mean_mid[ID == sample, c(9:19), with = FALSE])
      #naCols <- names(means[is.na(means) == TRUE])
      #means[is.na(means) == TRUE] <- 1000
      
      mvn <- data.table(ID = dat_mean_mid[ID == sample, ID],
                        rmvnorm(n = dat_mean_mid[ID == sample, N], mean = means, sigma = Covs_mid[[block]]),
                        rep = 1:dat_mean_mid[ID == sample, N], N = dat_mean_mid[ID == sample, N])
      #mvn[, c(naCols) := NA]
      mvn_dat_mid_block[[sample]] <- mvn
      
    }
    mvn_dat_mid[[block]] <- do.call(rbind, mvn_dat_mid_block)
    
    #-----clean up IDs and add date_chnm_plate, conc, and other info
    mvn_dat_mid[[block]][, date := str_split_fixed(ID, "_", n = 6)[,1], ]
    mvn_dat_mid[[block]][, chnm := str_split_fixed(ID, "_", n = 6)[,2], ]
    mvn_dat_mid[[block]][, plate := str_split_fixed(ID, "_", n = 6)[,3], ]
    mvn_dat_mid[[block]][, casn := str_split_fixed(ID, "_", n = 6)[,4], ]
    mvn_dat_mid[[block]][, conc := str_split_fixed(ID, "_", n = 6)[,5], ]
    mvn_dat_mid[[block]][, conc_index := str_split_fixed(ID, "_", n = 6)[,6], ]
    mvn_dat_mid[[block]][, date_chnm_plate := paste(str_split_fixed(ID, "_", n = 6)[,1], str_split_fixed(ID, "_", n = 6)[,2],
                                                    str_split_fixed(ID, "_", n = 6)[,3], sep = "_"), ]
    #mvn_dat_mid[[block]][ID == "DMSO_n",date_chnm_plate := paste(block, "DMSO", "Plate", sep = "_")]
    mvn_dat_mid[[block]][, ID := NULL]
    mvn_dat_mid[[block]][, N := .N, by = .(date_chnm_plate, conc, casn, date)]
  }
  
  mvn_dat_mid <- do.call(rbind, c(mvn_dat_mid, list(fill = TRUE)))
  
  
  #-----Lower
  dat_mean_lower <- dat_mean_wide[date_chnm_plate %in% unique(dat_lower[, date_chnm_plate])]
  
  mvn_dat_lower <- vector("list", length = 8)
  names(mvn_dat_lower) <- unique(dat[, date])
  
  for(block in names(Models)){
    mvn_dat_lower_block <- vector("list", length = nrow(dat_mean_lower[date == block,]))
    names(mvn_dat_lower_block) <- dat_mean_lower[date == block, ID]
    
    for(sample in names(mvn_dat_lower_block)){
      means <- as.numeric(dat_mean_lower[ID == sample, c(9:19), with = FALSE])
      names(means) <- colnames(dat_mean_lower[ID == sample, c(9:19), with = FALSE])
      #naCols <- names(means[is.na(means) == TRUE])
      #means[is.na(means) == TRUE] <- 1000
      
      mvn <- data.table(ID = dat_mean_lower[ID == sample, ID],
                        rmvnorm(n = dat_mean_lower[ID == sample, N], mean = means, sigma = Covs_lower[[block]]),
                        rep = 1:dat_mean_lower[ID == sample, N])
      #mvn[, c(naCols) := NA]
      mvn_dat_lower_block[[sample]] <- mvn
      
    }
    mvn_dat_lower[[block]] <- do.call(rbind, mvn_dat_lower_block)
    
    #-----clean up IDs and add date_chnm_plate, conc, and other info
    mvn_dat_lower[[block]][, date := str_split_fixed(ID, "_", n = 6)[,1], ]
    mvn_dat_lower[[block]][, chnm := str_split_fixed(ID, "_", n = 6)[,2], ]
    mvn_dat_lower[[block]][, plate := str_split_fixed(ID, "_", n = 6)[,3], ]
    mvn_dat_lower[[block]][, casn := str_split_fixed(ID, "_", n = 6)[,4], ]
    mvn_dat_lower[[block]][, conc := str_split_fixed(ID, "_", n = 6)[,5], ]
    mvn_dat_lower[[block]][, conc_index := str_split_fixed(ID, "_", n = 6)[,6], ]
    mvn_dat_lower[[block]][, date_chnm_plate := paste(str_split_fixed(ID, "_", n = 6)[,1], str_split_fixed(ID, "_", n = 6)[,2],
                                                      str_split_fixed(ID, "_", n = 6)[,3], sep = "_"), ]
    mvn_dat_lower[[block]][, ID := NULL]
    mvn_dat_lower[[block]][, N := .N, by = .(date_chnm_plate, conc, casn, date)]
  }
  
  mvn_dat_lower <- do.call(rbind, c(mvn_dat_lower, list(fill = TRUE)))
  
  #-----DMSO
  dat_mean_dmso <- dat_mean_wide[date_chnm_plate %in% unique(dat_dmso[, date_chnm_plate])]
  
  mvn_dat_dmso <- vector("list", length = 8)
  names(mvn_dat_dmso) <- unique(dat[, date])
  
  for(block in names(Models)){
    mvn_dat_dmso_block <- vector("list", length = nrow(dat_mean_dmso[date == block,]))
    names(mvn_dat_dmso_block) <- dat_mean_dmso[date == block, ID]
    
    for(sample in names(mvn_dat_dmso_block)){
      means <- as.numeric(dat_mean_dmso[ID == sample, c(9:19), with = FALSE])
      names(means) <- colnames(dat_mean_dmso[ID == sample, c(9:19), with = FALSE])
      #naCols <- names(means[is.na(means) == TRUE])
      #means[is.na(means) == TRUE] <- 1000
      
      mvn <- data.table(ID = dat_mean_dmso[ID == sample, ID],
                        rmvnorm(n = dat_mean_dmso[ID == sample, N], mean = means, sigma = Covs_dmso[[block]]),
                        rep = 1:dat_mean_dmso[ID == sample, N])
      #mvn[, c(naCols) := NA]
      mvn_dat_dmso_block[[sample]] <- mvn
      
    }
    mvn_dat_dmso[[block]] <- do.call(rbind, mvn_dat_dmso_block)
    
    #-----clean up IDs and add date_chnm_plate, conc, and other info
    mvn_dat_dmso[[block]][, date := str_split_fixed(ID, "_", n = 6)[,1], ]
    mvn_dat_dmso[[block]][, chnm := str_split_fixed(ID, "_", n = 6)[,2], ]
    mvn_dat_dmso[[block]][, plate := str_split_fixed(ID, "_", n = 6)[,3], ]
    mvn_dat_dmso[[block]][, casn := "DMSO", ]
    mvn_dat_dmso[[block]][, conc := str_split_fixed(ID, "_", n = 6)[,5], ]
    mvn_dat_dmso[[block]][, conc_index := str_split_fixed(ID, "_", n = 6)[,6], ]
    mvn_dat_dmso[[block]][, date_chnm_plate := paste(str_split_fixed(ID, "_", n = 6)[,1], str_split_fixed(ID, "_", n = 6)[,2],
                                                     str_split_fixed(ID, "_", n = 6)[,3], sep = "_"), ]
    #mvn_dat_dmso[[block]][ID == "DMSO_n",date_chnm_plate := paste(block, "DMSO", "Plate", sep = "_")]
    mvn_dat_dmso[[block]][, ID := NULL]
    mvn_dat_dmso[[block]][, N := .N, by = .(date_chnm_plate, conc, casn, date)]
  }
  
  mvn_dat_dmso <- do.call(rbind, c(mvn_dat_dmso, list(fill = TRUE)))
  
  #-----Generate fake negative data based on one DMSO value, six concentrations
  dat_mean_negative <- dat_mean_wide[date_chnm_plate == "20170411_DMSO_Plate1"]
  dat_mean_negative[, c("date_chnm_plate", "date", "chnm", "ID") := .("simulated_DMSO_Plate1",
                                                                      "simulated", "DMSO",
                                                                      "simulated_DMSO_Plate1_DMSO_0_CONC_DMSO")]
  
  dat_mean_negative <- rbind(dat_mean_negative, dat_mean_negative, dat_mean_negative,
                             dat_mean_negative, dat_mean_negative, dat_mean_negative,
                             dat_mean_negative)
  
  dat_mean_negative[2:7, c("date_chnm_plate", "date", "chnm") := .("simulated_negative_Plate1",
                                                                   "simulated", "negatve")]
  dat_mean_negative[2:7, c("conc", "casn", "conc_index") := .(1:6, "negative", c("CONC1", "CONC2", "CONC3", "CONC4", "CONC5", "CONC6"))]
  dat_mean_negative[chnm != "DMSO", wllt := "t"]
  dat_mean_negative[2:7, ID := paste(date_chnm_plate, casn, conc, conc_index, sep = "_")]
  
  #add in fake response data 
  fake_sim <- fread("./misc/simulated_data.txt")
  fake_sim[chnm == "DMSO", wllt := "n"]
  
  fake_sim <- fake_sim[, c(1:2, 19, 3:4, 16:18, 7, 8, 6, 5, 9, 15, 14, 12, 11, 13, 10)]
  fake_sim[, ID := paste(date_chnm_plate, casn, conc, conc_index, sep = "_")]
  
  #merge negative and fake response data
  dat_mean_negative <- rbind(dat_mean_negative, fake_sim)

  
  #make simulated negative dat
  mvn_dat_negative <- list()
  
  for(x in 1:nrow(dat_mean_negative)){
    means <- as.numeric(dat_mean_negative[x, c(9:19), with = FALSE])
    names(means) <- colnames(dat_mean_negative[x, c(9:19), with = FALSE])
    
    if (str_detect(dat_mean_negative[x, date_chnm_plate],  "Plate1")){
      mvn <- data.table(ID = dat_mean_negative[x, ID],
                        rmvnorm(n = dat_mean_negative[x, N], mean = means, sigma = Covs_dmso[["20170411"]]),
                        rep = 1:dat_mean_negative[x, N])
    } else if (dat_mean_negative[x, chnm] == "DSMO"){
      mvn <- data.table(ID = dat_mean_negative[x, ID],
                        rmvnorm(n = dat_mean_negative[x, N], mean = means, sigma = CovTP_dmso),
                        rep = 1:dat_mean_negative[x, N])
    } else{
      mvn <- data.table(ID = dat_mean_negative[x, ID],
                        rmvnorm(n = dat_mean_negative[x, N], mean = means, sigma = estVar(Models[["20170411"]]$fit)),
                        rep = 1:dat_mean_negative[x, N])
    }
    #mvn[, c(naCols) := NA]
    mvn_dat_negative[[x]] <- mvn
  }
  
  mvn_dat_negative <- do.call(rbind, mvn_dat_negative)
  
  #-----clean up IDs and add date_chnm_plate, conc, and other info
  mvn_dat_negative[, date := str_split_fixed(ID, "_", n = 6)[,1], ]
  mvn_dat_negative[, chnm := str_split_fixed(ID, "_", n = 6)[,2], ]
  mvn_dat_negative[, casn := str_split_fixed(ID, "_", n = 6)[,4], ]
  mvn_dat_negative[, plate := str_split_fixed(ID, "_", n = 6)[,3], ]
  mvn_dat_negative[, conc := str_split_fixed(ID, "_", n = 6)[,5], ]
  mvn_dat_negative[, conc_index := str_split_fixed(ID, "_", n = 6)[,6], ]
  mvn_dat_negative[, date_chnm_plate := paste(str_split_fixed(ID, "_", n = 6)[,1], str_split_fixed(ID, "_", n = 6)[,2],
                                                   str_split_fixed(ID, "_", n = 6)[,3], sep = "_"), ]
  #mvn_dat_negative[ID == "DMSO_n",date_chnm_plate := paste(block, "DMSO", "Plate", sep = "_")]
  mvn_dat_negative[, ID := NULL]
  mvn_dat_negative[, N := .N, by = .(date_chnm_plate, conc, casn, date)]
  
  #-------------------------------------------------------------------------#
  #-----Step 4: Run manova on global simulated data for covariance matrix
  #-------------------------------------------------------------------------#
  
  #-----Merge simulated data 
  mvn_sim_dat <- rbind(mvn_dat_upper, mvn_dat_mid, mvn_dat_lower, mvn_dat_dmso, mvn_dat_negative)
  mvn_sim_dat[, wllt := ifelse(chnm == "DMSO", "n", "t")]
  
  mvn_sim_dat[, ID := paste(date_chnm_plate, casn, conc, conc_index, sep = "_")]
  mvn_sim_dat <- mvn_sim_dat[, c(hormones_sub, colnames(mvn_sim_dat)[12:22]), with = FALSE]
  
  mvn_sim_dat_long <- melt(mvn_sim_dat, id.vars = colnames(mvn_sim_dat)[c(19, 13:18, 12, 20, 21, 22)], value.name = "log_uM",
                           variable.name = "steroid") 
  
  mvn_sim_dat_long <- na.omit(mvn_sim_dat_long) #remove NA values
  
  dat_mean[chnm == "DMSO", casn := "DMSO"]
  dat_mean[, ID := paste(date_chnm_plate, casn, conc, conc_index, sep = "_")]
  
  #-----Fix non-detect values since they should have no variance
    for(x in unique(mvn_sim_dat_long[date != "simulated", ID])){
      for(y in unique(mvn_sim_dat_long[ID == x, steroid])){
        if(dat_mean[ID == x & steroid == y, mean_log_uM] == LLOQuM[y]){
          mvn_sim_dat_long[ID == x & steroid == y,  log_uM := LLOQuM[y]]
  
        }
      }
    }
    
    for(y in unique(mvn_sim_dat_long[, steroid])){ #if simulated data is below LLOQ set to LLOQ
      mvn_sim_dat_long[steroid == y, log_uM := ifelse(log_uM < LLOQuM[y], LLOQuM[y], log_uM)]
    }
  
  #-----Manova on global simulated data
  Blocks <- unique(mvn_sim_dat_long[date != "simulated", date])
  Models_sim_all <- vector("list", length=length(Blocks))
  names(Models_sim_all) <- Blocks
  
  Residuals <- vector("list", length=length(Blocks))
  names(Residuals) <- Blocks
  
  for (block in unique(mvn_sim_dat_long[date != "simulated", date])) {
    dt <- copy(mvn_sim_dat_long[date == block, c("conc", "conc_index", "steroid", "chnm", "rep", "casn", "log_uM", "date", "wllt","date_chnm_plate", "plate"), with = FALSE])
    
    dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
    dt[, casn := NULL]
    chnm_apids <- unique(dt[,date_chnm_plate])
    chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
    chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
    
    dt[, ID := paste0(rep, "_", chnm, "_", plate)]
    #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
    #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
    dtw <- dcast.data.table(dt[,c(-5)], ID + ... ~ steroid, value.var = "log_uM")
    colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
    #dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
    
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
    
    Models_sim_all[[block]] <- list(fit = zz[[1]],
                                    Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                    IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
  } 
  
  #-----Estimate pooled covariance matrix from Models_sim_all
  Covs_sim_all <- vector("list", length = 8)
  names(Covs_sim_all) <- names(Models_sim_all)
  
  for(block in names(Models_sim_all)){
    Covs_sim_all[[block]] <- estVar(Models_sim_all[[block]]$fit)
  }
  
  CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_upper)))
  for (blk in names(Models_sim_all)) CovT[rownames(Covs_sim_all[[blk]]),colnames(Covs_sim_all[[blk]]),blk] <- Covs_sim_all[[blk]]
  CovTP_sim_all <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances
  #CovTP_sim_all <- CovTP0 #test with original covariance
  
  #-------------------------------------------------------------------------#
  #-----Step 5: Calculate mean Mahalanobis and max mean Mahalanobis distance 
  #-----on global simulated data
  #-------------------------------------------------------------------------#
  
  #-----Generate average response of simluated data and define some variables
  mvn_sim_dat_mean <- mvn_sim_dat_long[, .(mean_log_uM = mean(log_uM)), by = .(date_chnm_plate, steroid, chnm, conc, conc_index, wllt, date, N)]
  
  CA <- unique(mvn_sim_dat_mean$date_chnm_plate)
  CA <- CA[!str_detect(CA, "DMSO")]
  
  #-----Mean Mahalanobis distance calculation
  out <- vector("list", length=length(CA))
  names(out) <- CA
  SSnms <- paste("N",hormones_sub,sep="_")
  for (ca in CA) {
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca,]
    #hormones_sub <- as.character(unique(tdta[, steroid]))
    con <- str_split_fixed(ca, "_", 3)[,3]
    con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
    blk <- unique(tdta[,date])
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
    #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
    tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
    tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
    tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
    tdta[, conc := as.numeric(conc)]
    indx <- order(tdta$conc)
    tdta <- tdta[indx,]
    colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
    ## Hormones are in columns 7:ncol(tdta)
    ## Drop hormones which are NA in the bottom row
    refrow <- 1
    maxconc <- nrow(tdta)
    y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
    y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
    y1 <- sweep(y1,2,y0,"-")
    y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
    NH <- length(which(hormones_sub %in% colnames(y1)))
    CI11P <- solve(CovTP_sim_all[which(rownames(CovTP_sim_all) %in% colnames(y1)),
                                 which(colnames(CovTP_sim_all) %in% colnames(y1))])
    SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
    maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
    out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                            conc_index = tdta$conc_index[-1],
                            #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                            D11P = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                            #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                            #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                            CA = rep(ca, nrow(y1)),
                            row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                            NH = NH,
                            N = maxsqrtN[-1]^2)
  }
  
  Dists_mvn_sim <- data.table(do.call("rbind", out))
  Dists_mvn_sim[, chnm := str_split_fixed(CA, "_", 3)[,2]]
  Dists_mvn_sim$nChem <- make.names(Dists_mvn_sim$chnm)
  
  
  Mahalanobis_sim_dat <- Dists_mvn_sim[, .(maxD11P = max(D11P)), by = .(CA, chnm, nChem)]
  for(x in unique(Mahalanobis_sim_dat[, CA])){
    Mahalanobis_sim_dat[CA == x, c("Conc", "conc_index", "NH", "N") := .(Dists_mvn_sim[CA == x & D11P == Mahalanobis_sim_dat[CA == x, maxD11P], Conc],
                                                                         Dists_mvn_sim[CA == x & D11P == Mahalanobis_sim_dat[CA == x, maxD11P], conc_index],
                                                                         Dists_mvn_sim[CA == x & D11P == Mahalanobis_sim_dat[CA == x, maxD11P], NH],
                                                                         Dists_mvn_sim[CA == x & D11P == Mahalanobis_sim_dat[CA == x, maxD11P], N])]
  }
  
  colnames(Mahalanobis_sim_dat)[1] <- "date_chnm_plate"
  
  #-------------------------------------------------------------------------#
  #-----Step 6: Re-partition the simulated data based on new maxmMd values
  #-------------------------------------------------------------------------#
  
  #-----partition data into upper/lower subsets based on maxmMd values from original analysis
  Mahalanobis_sim_dat[, date := str_split_fixed(date_chnm_plate, "_", n = 3)[,1]]
  
  block_upper_sim <- vector("list", length = length(unique(Mahalanobis_sim_dat[date != "simulated", date])))
  names(block_upper_sim) <- unique(Mahalanobis_sim_dat[date != "simulated", date])
  block_lower_sim <- vector("list", length = length(unique(Mahalanobis_sim_dat[date != "simulated", date])))
  names(block_lower_sim) <- unique(Mahalanobis_sim_dat[date != "simulated", date])
  
  for(block in unique(Mahalanobis_sim_dat[date != "simulated", date])){
    d <- Mahalanobis_sim_dat[date == block]
    d <- d[order(maxD11P),]
    l <- round(length(d[, maxD11P])/2)
    
    lower <- d[1:l,]
    upper <- d[c(l+1):max(length(d[, maxD11P])),]
    
    block_lower_sim[[block]] <- lower
    block_upper_sim[[block]] <- upper
  }
  
  lower_sim <- do.call(rbind, block_lower_sim)
  upper_sim <- do.call(rbind, block_upper_sim)
  
  lower_sim[, c("plate") := str_split_fixed(date_chnm_plate, "_", n = 3)[,3],]
  upper_sim[, c("plate") := str_split_fixed(date_chnm_plate, "_", n = 3)[,3],]
  
  #-----subset data: high and low
  dat_upper_sim <- copy(mvn_sim_dat_long[date_chnm_plate %in% upper_sim[, date_chnm_plate] | date == "simulated" | date_chnm_plate %in% paste(upper_sim[, date], "DMSO", upper_sim[, plate], sep = "_"),])
  dat_lower_sim <- copy(mvn_sim_dat_long[date_chnm_plate %in% lower_sim[, date_chnm_plate] | date == "simulated" | date_chnm_plate %in% paste(lower_sim[, date], "DMSO", lower_sim[, plate], sep = "_"),])
  
  #-------------------------------------------------------------------------#
  #-----Step 7: Rum manova on upper and lower simulated partitions to generate pooled covariance matrices
  #-------------------------------------------------------------------------#
  
  #-----Upper simulated manova
  Blocks <- unique(dat_upper_sim[date != "simulated", date])
  Models_sim_upper <- vector("list", length=length(Blocks))
  names(Models_sim_upper) <- Blocks
  
  Residuals <- vector("list", length=length(Blocks))
  names(Residuals) <- Blocks
  
  for (block in unique(dat_upper_sim[date != "simulated", date])) {
    dt <- copy(dat_upper_sim[date == block, c("conc", "conc_index", "steroid", "chnm", "rep", "casn", "log_uM", "date", "wllt","date_chnm_plate", "plate"), with = FALSE])
    
    dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
    dt[, casn := NULL]
    chnm_apids <- unique(dt[,date_chnm_plate])
    chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
    chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
    
    dt[, ID := paste0(rep, "_", chnm, "_", plate)]
    #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
    #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
    dtw <- dcast.data.table(dt[,c(-5)], ID + ... ~ steroid, value.var = "log_uM")
    colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
    #dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
    
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
    
    Models_sim_upper[[block]] <- list(fit = zz[[1]],
                                      Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                      IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
  } 
  
  #-----Lower simulated manova
  Blocks <- unique(dat_lower_sim[date != "simulated", date])
  Models_sim_lower <- vector("list", length=length(Blocks))
  names(Models_sim_lower) <- Blocks
  
  Residuals <- vector("list", length=length(Blocks))
  names(Residuals) <- Blocks
  
  for (block in unique(dat_lower_sim[date != "simulated", date])) {
    dt <- copy(dat_lower_sim[date == block, c("conc", "conc_index", "steroid", "chnm", "rep", "casn", "log_uM", "date", "wllt","date_chnm_plate", "plate"), with = FALSE])
    
    dt$classFull <- factor(ifelse(dt$wllt == "t", paste0(make.names(dt[,date_chnm_plate]),"_",dt[, casn],"_",dt[,conc], "_", dt[,conc_index]), paste0("DMSO_",dt[,wllt])))
    dt[, casn := NULL]
    chnm_apids <- unique(dt[,date_chnm_plate])
    chnm_apids <- c(chnm_apids, "DMSO") #so it doesn't become empty vector
    chnm_apids <- chnm_apids[-grep("DMSO", chnm_apids)]
    
    dt[, ID := paste0(rep, "_", chnm, "_", plate)]
    #dtw <- reshape(dt, v.names="log_uM", timevar="steroid", idvar="ID", direction="wide", drop=c("coli","rowi"))
    #dtw <- dcast.data.table(dt[,c(-6, -7, -10)], ID + ... ~ steroid, value.var = "log_uM", fun.aggregate = mean)
    dtw <- dcast.data.table(dt[,c(-5)], ID + ... ~ steroid, value.var = "log_uM")
    colnames(dtw) <- gsub("^log_uM\\.","",colnames(dtw))
    #dtw <- dtw[,-match(c("DHEA","Pregnenolone"), colnames(dtw)), with = FALSE]
    
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
    
    Models_sim_lower[[block]] <- list(fit = zz[[1]],
                                      Hormones = if (any(!is.na(out$dropcols))) setdiff(names(hormones_all), c(out$dropcols, "Pregnenolone","DHEA")) else setdiff(names(hormones_all), c("Pregnenolone","DHEA")),
                                      IDs = data.table(dtw3[, c("date_chnm_plate", "conc_index", "conc")], ID = rownames(dtw3)))
  } 
  
  #-----Estimate pooled covariance matrix from upper and lower simulated manova output
  Covs_sim_upper <- vector("list", length = 8)
  names(Covs_sim_upper) <- names(Models_sim_upper)
  
  Covs_sim_lower <- vector("list", length = 8)
  names(Covs_sim_lower) <- names(Models_sim_lower)
  
  for(block in names(Models_sim_upper)){
    Covs_sim_upper[[block]] <- estVar(Models_sim_upper[[block]]$fit)
    Covs_sim_lower[[block]] <- estVar(Models_sim_lower[[block]]$fit)
  }
  
  CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_sim_upper)))
  for (blk in names(Models_sim_upper)) CovT[rownames(Covs_sim_upper[[blk]]),colnames(Covs_sim_upper[[blk]]),blk] <- Covs_sim_upper[[blk]]
  CovTP_sim_upper <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances
  
  CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models_sim_lower)))
  for (blk in names(Models_sim_lower)) CovT[rownames(Covs_sim_lower[[blk]]),colnames(Covs_sim_lower[[blk]]),blk] <- Covs_sim_lower[[blk]]
  CovTP_sim_lower <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances
  
  #-------------------------------------------------------------------------#
  #-----Step 8: Recalculate mean Mahalanobis distance using the upper and lower
  #-----simulated covariance matrices and also determine the maxmMd for each
  #-------------------------------------------------------------------------#
  
  #-----Mean Mahalanobis distance calculation for upper covariance (use mvn_sim_dat_mean to compare with global simulated Mahalanobis distances)
  out <- vector("list", length=length(CA))
  names(out) <- CA
  SSnms <- paste("N",hormones_sub,sep="_")
  for (ca in CA) {
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca,]
    #hormones_sub <- as.character(unique(tdta[, steroid]))
    con <- str_split_fixed(ca, "_", 3)[,3]
    con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
    blk <- unique(tdta[,date])
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
    #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
    tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
    tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
    tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
    tdta[, conc := as.numeric(conc)]
    indx <- order(tdta$conc)
    tdta <- tdta[indx,]
    colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
    ## Hormones are in columns 7:ncol(tdta)
    ## Drop hormones which are NA in the bottom row
    refrow <- 1
    maxconc <- nrow(tdta)
    y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
    y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
    y1 <- sweep(y1,2,y0,"-")
    y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
    NH <- length(which(hormones_sub %in% colnames(y1)))
    CI11P <- solve(CovTP_sim_upper[which(rownames(CovTP_sim_upper) %in% colnames(y1)),
                                   which(colnames(CovTP_sim_upper) %in% colnames(y1))])
    SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
    maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
    out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                            conc_index = tdta$conc_index[-1],
                            #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                            D11P_upper = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                            #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                            #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                            CA = rep(ca, nrow(y1)),
                            row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                            NH = NH,
                            N = maxsqrtN[-1]^2)
  }
  
  Dists_mvn_sim_upper <- data.table(do.call("rbind", out))
  Dists_mvn_sim_upper[, chnm := str_split_fixed(CA, "_", 3)[,2]]
  Dists_mvn_sim_upper$nChem <- make.names(Dists_mvn_sim_upper$chnm)
  
  maxmMd_upper <- Dists_mvn_sim_upper[, .(maxD11P_upper = max(D11P_upper)), by = .(CA, chnm, nChem)]
  
  #-----Mean Mahalanobis distance calculation for lower covariance (use mvn_sim_dat_mean to compare with global simulated Mahalanobis distances)
  out <- vector("list", length=length(CA))
  names(out) <- CA
  SSnms <- paste("N",hormones_sub,sep="_")
  for (ca in CA) {
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca,]
    #hormones_sub <- as.character(unique(tdta[, steroid]))
    con <- str_split_fixed(ca, "_", 3)[,3]
    con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
    blk <- unique(tdta[,date])
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
    #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
    tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
    tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
    tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
    tdta[, conc := as.numeric(conc)]
    indx <- order(tdta$conc)
    tdta <- tdta[indx,]
    colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
    ## Hormones are in columns 7:ncol(tdta)
    ## Drop hormones which are NA in the bottom row
    refrow <- 1
    maxconc <- nrow(tdta)
    y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
    y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
    y1 <- sweep(y1,2,y0,"-")
    y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
    NH <- length(which(hormones_sub %in% colnames(y1)))
    CI11P <- solve(CovTP_sim_lower[which(rownames(CovTP_sim_lower) %in% colnames(y1)),
                                   which(colnames(CovTP_sim_lower) %in% colnames(y1))])
    SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
    maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
    out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                            conc_index = tdta$conc_index[-1],
                            #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                            D11P_lower = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                            #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                            #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                            CA = rep(ca, nrow(y1)),
                            row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                            NH = NH,
                            N = maxsqrtN[-1]^2)
  }
  
  Dists_mvn_sim_lower <- data.table(do.call("rbind", out))
  Dists_mvn_sim_lower[, chnm := str_split_fixed(CA, "_", 3)[,2]]
  Dists_mvn_sim_lower$nChem <- make.names(Dists_mvn_sim_lower$chnm)
  
  maxmMd_lower <- Dists_mvn_sim_lower[, .(maxD11P_lower = max(D11P_lower)), by = .(CA, chnm, nChem)]

  #-----Mean Mahalanobis distance calculation for upper covariance based on original data (use mvn_sim_dat_mean to compare with global simulated Mahalanobis distances)
  out <- vector("list", length=length(CA))
  names(out) <- CA
  SSnms <- paste("N",hormones_sub,sep="_")
  for (ca in CA) {
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca,]
    #hormones_sub <- as.character(unique(tdta[, steroid]))
    con <- str_split_fixed(ca, "_", 3)[,3]
    con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
    blk <- unique(tdta[,date])
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
    #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
    tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
    tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
    tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
    tdta[, conc := as.numeric(conc)]
    indx <- order(tdta$conc)
    tdta <- tdta[indx,]
    colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
    ## Hormones are in columns 7:ncol(tdta)
    ## Drop hormones which are NA in the bottom row
    refrow <- 1
    maxconc <- nrow(tdta)
    y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
    y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
    y1 <- sweep(y1,2,y0,"-")
    y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
    NH <- length(which(hormones_sub %in% colnames(y1)))
    CI11P <- solve(CovTP_original_upper[which(rownames(CovTP_original_upper) %in% colnames(y1)),
                                   which(colnames(CovTP_original_upper) %in% colnames(y1))])
    SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
    maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
    out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                            conc_index = tdta$conc_index[-1],
                            #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                            D11P_original_upper = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                            #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                            #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                            CA = rep(ca, nrow(y1)),
                            row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                            NH = NH,
                            N = maxsqrtN[-1]^2)
  }
  
  Dists_upper <- data.table(do.call("rbind", out))
  Dists_upper[, chnm := str_split_fixed(CA, "_", 3)[,2]]
  Dists_upper$nChem <- make.names(Dists_upper$chnm)
  
  maxmMd_original_upper <- Dists_upper[, .(maxD11P_original_upper = max(D11P_original_upper)), by = .(CA, chnm, nChem)]

  #-----Mean Mahalanobis distance calculation for lower covariance based on original data (use mvn_sim_dat_mean to compare with global simulated Mahalanobis distances)
  out <- vector("list", length=length(CA))
  names(out) <- CA
  SSnms <- paste("N",hormones_sub,sep="_")
  for (ca in CA) {
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca,]
    #hormones_sub <- as.character(unique(tdta[, steroid]))
    con <- str_split_fixed(ca, "_", 3)[,3]
    con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
    blk <- unique(tdta[,date])
    tdta <- mvn_sim_dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
    #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
    tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
    tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
    tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
    tdta[, conc := as.numeric(conc)]
    indx <- order(tdta$conc)
    tdta <- tdta[indx,]
    colnames(tdta) <- gsub("^mean_log_uM_","",colnames(tdta))
    ## Hormones are in columns 7:ncol(tdta)
    ## Drop hormones which are NA in the bottom row
    refrow <- 1
    maxconc <- nrow(tdta)
    y0 <- data.matrix(tdta[1,c(hormones_sub), with = FALSE])
    y1 <- data.matrix(tdta[-1,c(hormones_sub), with = FALSE])
    y1 <- sweep(y1,2,y0,"-")
    y1 <- y1[, apply(y1, 2, function(x) !any(is.na(x)))]
    NH <- length(which(hormones_sub %in% colnames(y1)))
    CI11P <- solve(CovTP_original_lower[which(rownames(CovTP_original_lower) %in% colnames(y1)),
                               which(colnames(CovTP_original_lower) %in% colnames(y1))])
    SSmeasured <- colnames(tdta)[colnames(tdta) %in% SSnms]
    maxsqrtN <- sqrt(apply(tdta[, SSmeasured, with = FALSE],1,max, na.rm=TRUE))
    out[[ca]] <- data.frame(Conc = tdta$conc[-1],
                            conc_index = tdta$conc_index[-1],
                            #D11 = apply(y1,1, function(x) sqrt((x %*% CI11[[blk]] %*% x)/NH)),
                            D11P_original_lower = apply(y1,1, function(x) sqrt((x %*%  CI11P %*% x)/NH)),
                            #sD11 = apply(y1,1, function(x) sqrt(sum(x^2 * VI11[[blk]])/NH)),
                            #sD11P = apply(y1,1, function(x) sqrt(sum(x^2 *  VI11P[which(names(VI11P) %in% colnames(y1))])/NH)),
                            CA = rep(ca, nrow(y1)),
                            row.names=paste(rep(ca, nrow(y1)),tdta$conc[-1],sep="_X_"),
                            NH = NH,
                            N = maxsqrtN[-1]^2)
  }
  
  Dists_lower <- data.table(do.call("rbind", out))
  Dists_lower[, chnm := str_split_fixed(CA, "_", 3)[,2]]
  Dists_lower$nChem <- make.names(Dists_lower$chnm)
  
  maxmMd_original_lower <- Dists_lower[, .(maxD11P_original_lower = max(D11P_original_lower)), by = .(CA, chnm, nChem)]
  
  #-------------------------------------------------------------------------#
  #-----Step 9: Merge mean Mahalanobis distances and maxmMd data tables for
  #-----the simulated and store in global loop variable
  #-------------------------------------------------------------------------#
  
  colnames(Dists_mvn_sim)[4] <- "date_chnm_plate"
  colnames(Dists_mvn_sim_upper)[4] <- "date_chnm_plate"
  colnames(Dists_mvn_sim_lower)[4] <- "date_chnm_plate"
  
  colnames(Dists_upper)[4] <- "date_chnm_plate"
  colnames(Dists_lower)[4] <- "date_chnm_plate"
  
  Dists_mvn_merged <- copy(Dists_mvn_sim)
  
  for(x in unique(Dists_mvn_merged[, date_chnm_plate])){
    for(y in unique(Dists_mvn_merged[date_chnm_plate == x, Conc])){
      Dists_mvn_merged[date_chnm_plate == x & Conc == y, c("D11P_upper", "D11P_lower", "D11P_original_upper", "D11P_original_lower") := .(Dists_mvn_sim_upper[date_chnm_plate == x & Conc == y, D11P_upper],
                                                                                                                                          Dists_mvn_sim_lower[date_chnm_plate == x & Conc == y, D11P_lower],
                                                                                                                                          Dists_upper[date_chnm_plate == x & Conc == y, D11P_original_upper],
                                                                                                                                          Dists_lower[date_chnm_plate == x & Conc == y, D11P_original_lower])] 
    }
  }
  #Dists_mvn_merged[, sim := i] #keep track of loop number
  
  Mahalanobis_sim_dat_merged <- copy(Mahalanobis_sim_dat)
  
  for(x in unique(Mahalanobis_sim_dat_merged[, date_chnm_plate])){
    Mahalanobis_sim_dat_merged[date_chnm_plate == x, c("maxD11P_upper", "maxD11P_lower", "maxD11P_original_upper", "maxD11P_original_lower") := .(maxmMd_upper[CA == x, maxD11P_upper],
                                                                                                                                                  maxmMd_lower[CA == x, maxD11P_lower],
                                                                                                                                                  maxmMd_original_upper[CA == x, maxD11P_original_upper],
                                                                                                                                                  maxmMd_original_lower[CA == x, maxD11P_original_lower])]
  }
  #Mahalanobis_sim_dat_merged[, sim := i] #keep track of loop number
  
  #-----Combine Dists_mvn_merged and Mahalanobis_sim_dat_merged into a list to return to the function call
  returned_list <- list(all_dists = Dists_mvn_merged, maxmMds = Mahalanobis_sim_dat_merged, sim_dat = mvn_sim_dat, CovTP = CovTP_sim_all, CovTP_upper = CovTP_sim_upper, CovTP_lower = CovTP_sim_lower)
  
  return(returned_list)
  #-----Store simulated Mahalanobis results
  # simulated_mMds[[iteration]] <- Dists_mvn_merged
  # simulated_maxmMds[[iteration]] <- Mahalanobis_sim_dat_merged
  # 
  # break #kill loop to test
}


simulated_mMd_output <- mclapply(X = 1:2000, FUN = simulated_distance_function, mc.cores = 34) #now use mclapply()!!!!!!!!!!!

save(simulated_mMd_output, file = paste0("./RData/simulated_Mahalanobis_distance_output_", Sys.Date(), ".RData"))
#save(simulated_mMd_output, file = paste0("./simulated_Mahalanobis_distance_output_originalCovariance_", Sys.Date(), ".RData"))
