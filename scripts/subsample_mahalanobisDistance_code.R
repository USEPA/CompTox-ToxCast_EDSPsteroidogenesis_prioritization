#-----------------------------------------------------------------------------------#

#This code will generate mMds using the three
#different covariance matrix types based on the
#original published HT-H295R data set. These 
#values are used to help generate Supp Figure 1
#and the filled triangels in Figure 3 and Supp Figure 2


library(data.table)
library(boot)
library(stringr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)
library(MASS)
library(mvtnorm)

rm(list = ls())

#setwd("/share/home/dhaggard/Desktop/Mahalanobis Follow Up")

#-------------------------------------------------------------------------#
#-----Load in MANOVA output and original mMd RData file to get mMd distribution
#-------------------------------------------------------------------------#
load("./RData/subsample_Model_output2018-10-03.RData")
load("./RData/AllResps_outliersRemoved2018-10-02.RData")
load("./RData/updated_manova_output2018-02-06.RData")

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

#-------------------------------------------------------------------------------------------------#
#-----Clean up dat
#-------------------------------------------------------------------------------------------------#

dat[spid == "DMSO", chnm := "DMSO"]

dat[spid == "DMSO", date_chnm_plate := str_replace(date_chnm_plate, "NA", "DMSO")]

dat <- dat[steroid != "DHEA" & steroid != "Pregnenolone"]

#-------------------------------------------------------------------------#
#-----OPTIONAL: Filter the samples with residual outliers-----------------#
#-------------------------------------------------------------------------#

dat[, conc_date_chnm_plate := paste(conc, date_chnm_plate, sep = "_")]

dat <- dat[!c(conc_date_chnm_plate %in% Residuals[maxSD == 1, conc_date_chnm_plate]),] #filter based on Residuals with the maxSD filter == 1
dat[, conc_date_chnm_plate := NULL]

#-------------------------------------------------------------------------------------------------#
#-----Subset dat_mean to only include ~100 test chemicals and the simulated data and 
#-----verage the uM values between replicates across test chemicals
#-------------------------------------------------------------------------------------------------#
#OECD reference chemcials
# OECD <- c("1912-24-9", "17804-35-2", "94-26-8", "66575-29-9", "2212-67-1", "26027-38-3", "67747-09-5", "125-84-8", "4672-49-5",
#           "112809-51-5", "13647-35-3", "80-05-7", "51-28-5", "17230-88-5", "117-81-7", "60-51-5", "60168-88-9", "98319-26-7",
#           "13311-84-7", "446-72-0", "65277-42-1", "84371-65-3", "51-03-6", "1610-18-0", "52-01-7", "1330-78-5")
# 
# others <- c("Fadrozole hydrochloride", "Benzyl butyl phthalate", "EPN", "Epoxiconazole", "Flumetsulam", "Metyrapone", 
#             "Diethyl phthalate", "Isophorone diisocyanate", "Nelivaptan", "1,4-Dihydroxy-2-naphthoic acid",
#             "2-Methoxy-5-nitroaniline", "Allethrin", "Captan", "Cerivastatin sodium", "Clotrimazole",
#             "2-Methoxy-4-nitroaniline", "3,3'-Dimethylbenzidine")
# 
# dat_sub <- copy(dat[casn %in% OECD | chnm %in% others | chnm == "DMSO",])
dat_sub <- copy(dat)
dat_sub[, N := .N, by = .(date_chnm_plate, steroid, chnm, conc, conc_index)]

dat_mean <- dat_sub[, .(mean_log_uM = mean(log_uM)), by = .(date_chnm_plate, steroid, chnm, conc, conc_index, wllt, date, N, casn)]

#-------------------------------------------------------------------------------------------------#
#-----Define some variables for mMd calculation
#-------------------------------------------------------------------------------------------------#
CA <- unique(dat_mean$date_chnm_plate)
CA <- CA[!str_detect(CA, "DMSO")]

hormones_all <- colnames(Pvalues)[5:17]
hormones_sub <- hormones_all[c(-8, -11)]

#-------------------------------------------------------------------------------------------------#
#-----all data
#-------------------------------------------------------------------------------------------------#

Covs <- vector("list", length = 8)
names(Covs) <- names(Models)

for(block in names(Models_all_subsample)){
  nms <- colnames(Models_all_subsample[[block]]$fit$residuals)
  hormones <- length(nms)
  Covs[[block]] <- estVar(Models_all_subsample[[block]]$fit)
  rownames(Covs[[block]]) <- nms
  colnames(Covs[[block]]) <- nms
}

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models)))
for (blk in names(Models_all_subsample)) CovT[rownames(Covs[[blk]]),colnames(Covs[[blk]]),blk] <- Covs[[blk]]
CovTP <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

nms <- hormones_sub

out <- vector("list", length=length(CA))
names(out) <- CA
SSnms <- paste("N",nms,sep="_")
for (ca in CA) {
  tdta <- dat_mean[date_chnm_plate == ca,]
  con <- str_split_fixed(ca, "_", 3)[,3]
  con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
  blk <- unique(tdta[,date])
  tdta <- dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
  #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
  tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
  tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
  tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
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
  CI11P <- solve(CovTP[which(rownames(CovTP) %in% colnames(y1)),
                       which(colnames(CovTP) %in% colnames(y1))])
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
Dists_all_subsample <- data.table(do.call("rbind", out))
Dists_all_subsample[, chnm := str_split_fixed(CA, "_", 3)[,2]]
Dists_all_subsample$nChem <- make.names(Dists_all_subsample$chnm)
rownames(Dists_all_subsample) <- as.character(1:nrow(Dists_all_subsample))

#-------------------------------------------------------------------------------------------------#
#-----upper data
#-------------------------------------------------------------------------------------------------#

Covs <- vector("list", length = 8)
names(Covs) <- names(Models)

for(block in names(Models_upper_subsample)){
  nms <- colnames(Models_upper_subsample[[block]]$fit$residuals)
  hormones <- length(nms)
  Covs[[block]] <- estVar(Models_upper_subsample[[block]]$fit)
  rownames(Covs[[block]]) <- nms
  colnames(Covs[[block]]) <- nms
}

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models)))
for (blk in names(Models_upper_subsample)) CovT[rownames(Covs[[blk]]),colnames(Covs[[blk]]),blk] <- Covs[[blk]]
CovTP <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

nms <- hormones_sub

out <- vector("list", length=length(CA))
names(out) <- CA
SSnms <- paste("N",nms,sep="_")
for (ca in CA) {
  tdta <- dat_mean[date_chnm_plate == ca,]
  con <- str_split_fixed(ca, "_", 3)[,3]
  con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
  blk <- unique(tdta[,date])
  tdta <- dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
  #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
  tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
  tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
  tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
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
  CI11P <- solve(CovTP[which(rownames(CovTP) %in% colnames(y1)),
                       which(colnames(CovTP) %in% colnames(y1))])
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
Dists_upper_subsample <- data.table(do.call("rbind", out))
Dists_upper_subsample[, chnm := str_split_fixed(CA, "_", 3)[,2]]
Dists_upper_subsample$nChem <- make.names(Dists_upper_subsample$chnm)
rownames(Dists_upper_subsample) <- as.character(1:nrow(Dists_upper_subsample))

#-------------------------------------------------------------------------------------------------#
#-----lower data
#-------------------------------------------------------------------------------------------------#

Covs <- vector("list", length = 8)
names(Covs) <- names(Models)

for(block in names(Models_lower_subsample)){
  nms <- colnames(Models_lower_subsample[[block]]$fit$residuals)
  hormones <- length(nms)
  Covs[[block]] <- estVar(Models_lower_subsample[[block]]$fit)
  rownames(Covs[[block]]) <- nms
  colnames(Covs[[block]]) <- nms
}

CovT <- array(NA, dim=c(11, 11, 8), dimnames=list(hormones_sub, hormones_sub, names(Models)))
for (blk in names(Models_lower_subsample)) CovT[rownames(Covs[[blk]]),colnames(Covs[[blk]]),blk] <- Covs[[blk]]
CovTP <- apply(CovT, c(1, 2), mean, na.rm=TRUE) #make pooled matrix based on block subset of covariances

nms <- hormones_sub

out <- vector("list", length=length(CA))
names(out) <- CA
SSnms <- paste("N",nms,sep="_")
for (ca in CA) {
  tdta <- dat_mean[date_chnm_plate == ca,]
  con <- str_split_fixed(ca, "_", 3)[,3]
  con <- paste0(unique(tdta[,date]), "_", "DMSO_", con)
  blk <- unique(tdta[,date])
  tdta <- dat_mean[date_chnm_plate == ca | date == blk & date_chnm_plate == con]
  #tdta <- tdta[tdta$acsn %in% colnames(Covs[[blk]]),]
  tdta[chnm == "DMSO", conc_index := "CONC_DMSO"]
  tdta$id <- paste0(tdta$date_chnm_plate,"_X_",tdta$conc)
  tdta <- dcast.data.table(tdta, ... ~ steroid, value.var = c("mean_log_uM", "N"))
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
  CI11P <- solve(CovTP[which(rownames(CovTP) %in% colnames(y1)),
                       which(colnames(CovTP) %in% colnames(y1))])
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
Dists_lower_subsample <- data.table(do.call("rbind", out))
Dists_lower_subsample[, chnm := str_split_fixed(CA, "_", 3)[,2]]
Dists_lower_subsample$nChem <- make.names(Dists_lower_subsample$chnm)
rownames(Dists_lower_subsample) <- as.character(1:nrow(Dists_lower_subsample))

#-------------------------------------------------------------------------------------------------#
#-----original distances
#-------------------------------------------------------------------------------------------------#

Dists_original <- Dists[CA %in% Dists_all_subsample[, CA],]

#-------------------------------------------------------------------------------------------------#
#-----merge and plot
#-------------------------------------------------------------------------------------------------#
Dists_all_subsample[, subsample := "all"]
Dists_upper_subsample[, subsample := "upper"]
Dists_lower_subsample[, subsample := "lower"]
Dists_original[, subsample := "original"]

save(Dists_all_subsample, Dists_upper_subsample, Dists_lower_subsample, file = paste0("./RData/subsample_mMd_output", Sys.Date(), ".RData"))


#-----Plot everything
# Dists_merged <- rbind(Dists_all_subsample[!str_detect(CA, "simulated"),], Dists_upper_subsample[!str_detect(CA, "simulated"),],
#                       Dists_lower_subsample[!str_detect(CA, "simulated"),], Dists_original[!str_detect(CA, "simulated"),])
# 
# pdf("new_merged_subsample_mMd.pdf", width = 12, height = 8)
# for(x in unique(Dists_merged[, CA])){
#   p <- ggplot(data = Dists_merged[CA == x,], aes(x = factor(Conc), y = D11P, color = subsample)) +
#     geom_point(position = position_dodge(.25), size = 3) +
#     #geom_jitter(position = position_jitterdodge()) +
#     #geom_point(aes(x = factor(Conc), y = original_D11P, color = "original")) +
#     theme_few() +
#     ggtitle(x) +
#     xlab("Concentration (uM)") +
#     ylab("mMd") +
#     scale_color_brewer(palette = "Dark2", name = "subsample")
#   plot(p)
# }
# dev.off()
