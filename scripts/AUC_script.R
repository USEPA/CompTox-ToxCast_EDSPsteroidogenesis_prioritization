#-----------------------------------------------------------------------------------#

#Calculate the Area Under the Curve (AUC) of the BMD fits to determine efficacy of mMd responses#


library(tcpl)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(plyr)
library(stringr)
library(pracma)
library(stats)

rm(list = ls())

#setwd("L:/Lab/NCCT_ToxCast/Derik Haggard/Mahalanobis Follow Up") #change accordingly

#-----Load in Mahalanobis distances and mMd fits
load("./RData/AllResps_outliersRemoved2018-10-02.RData")
load("./RData/mMdFits2018-10-03.RData")


#-----Define logistic function
logistic4par  <- function(lx, T, cc, d) {
  (1 + (cc - 1) / (1 + exp(-d * log(lx) + T)))
}

logistic4par_integrate  <- function(lx) {
  (1 + (cc - 1) / (1 + exp(-d * (lx) + T)))
}

#-----Cacluate AUC for hill fits
Fits <- data.table(Fits)
Fits[, colnames(Fits)[2:10] := .(as.numeric(T), as.numeric(cc), as.numeric(d), as.numeric(BMD), as.numeric(MaxmMd),
                                 as.numeric(Scrit), as.numeric(cor), as.numeric(cor_pvalue), as.numeric(convergence))]

#Calculate AUC and add to Fits
for(x in unique(Mahalanobis_dat[, date_chnm_plate])){
  minC <- min(dat[date_chnm_plate == x, conc])
  maxC <- max(dat[date_chnm_plate == x, conc])
  
  #conc vector
  concs <- seq(minC, maxC, length.out = 1000)
  
  #function params
  T <- Fits[Name == x, T]
  cc <- Fits[Name == x, cc]
  d <- Fits[Name == x, d]
  
  output <- c(1:1000)
  for(i in 1:1000){
    output[i] <- logistic4par(lx = concs[i], T = T, cc = cc, d = d)
  }
  
  Fits[Name == x, c("AUC_trap", "AUC_stats") := .(trapz(log(concs), output), integrate(logistic4par_integrate, log(minC), log(maxC), subdivisions = 1000)[[1]])]
}

#-----Plot Curve fits along with original mMd values
# pdf("./figures/mMd_hill_fits.pdf", width = 12, height = 8)
# for(x in unique(Mahalanobis_dat[, date_chnm_plate])){
#   minC <- min(dat[date_chnm_plate == x, conc])
#   maxC <- max(dat[date_chnm_plate == x, conc])
#   
#   #conc vector
#   concs <- seq(minC, maxC, length.out = 1000)
#   
#   #function params
#   T <- Fits[Name == x, T]
#   cc <- Fits[Name == x, cc]
#   d <- Fits[Name == x, d]
#   
#   output <- c(1:1000)
#   for(i in 1:1000){
#     output[i] <- logistic4par(lx = concs[i], T = T, cc = cc, d = d)
#   }
#   
#   dat_int <- data.table(concs = concs, mMd = output)
#   
#   p <- ggplot(data = dat_int, aes(x = concs, y = mMd)) +
#     geom_line() +
#     geom_point(data = Dists[CA == x,], aes(x = Conc, y = D11P)) +
#     geom_hline(yintercept = Mahalanobis_dat[date_chnm_plate == x, Scrit01], linetype = "dashed") +
#     theme_few() +
#     scale_x_log10() +
#     annotation_logticks(base = 10, sides = "b") +
#     ggtitle(label = x)
#   
#   plot(p)
# }
# dev.off()

View(Fits[order(AUC_stats)])

#-----Plot histogram of AUC curves
# p <- ggplot(data = Fits, aes(x = reorder(Name, AUC_stats), y = AUC_stats)) +
#   geom_bar(stat = "identity", width = 0.5, color = "black") +
#   theme_few() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   scale_x_discrete(expand = c(0,0)) +
#   ylab("AUC") +
#   xlab("Sample")
# 
# pdf("./figures/mmD_AUC_histogram.pdf", width = 12, height = 8)
# p
# dev.off()

#save output with AUC estiamtes
Fits_AUC <- Fits

save(Fits_AUC, file = paste("./RData/mMD_AUC_Fits_", Sys.Date() ,".RData", sep = ""))
