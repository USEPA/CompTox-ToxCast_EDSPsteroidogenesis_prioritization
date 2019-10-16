#-----------------------------------------------------------------------------------#

# Code to run 4-parameter Hill models on mean Mahalanobis Distance values to calculate the BMD (conce where
# mMd goes above the critical limit).


rm(list = ls())

library(parallel)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(GGally)
library(stringr)
library(data.table)
library(ggthemes)
library(dplyr)
library(cowplot)


#setwd("./Mahalanobis Distance") #change accordingly

#-------------------------------------------------------------------------------------------------#
#-----Functions used
#-------------------------------------------------------------------------------------------------#

logistic4par  <- function(lx, T, cc, d) {
    (1 + (cc - 1) / (1 + exp(-d * lx + T)))
}

BMD <- function(Z, T, cc, d) {
    lD <- unlist((log((cc - 1) / (Z - 1) - 1) - T) / (- d))
    if (Z > cc | d < 0.01) Inf else if (is.finite(lD)) exp(lD) else lD
}

dofit <- function(Dsnm) {
    of <- function(parms) {
        T <- parms[1]
        cc <- exp(parms[2])
        d <- exp(parms[3])
        sum((DS[[Dsnm]]$y - logistic4par(DS[[Dsnm]]$lx, T, cc, d))^2)
    }
    Spcor <- cor(DS[[Dsnm]]$y, DS[[Dsnm]]$lx, method="spearman")
    Spcortest <- cor.test(DS[[Dsnm]]$y, DS[[Dsnm]]$lx, method="spearman")$p.value
    out <- optim(c(0, 1, 0), of, control=list(maxit=2000))
    Parms <- c(T = out$par[1], cc = exp(out$par[2]), d = exp(out$par[3]))
    list(Parms=Parms, BMD=unname(BMD(DS[[Dsnm]]$Scrit01, Parms[1], Parms[2], Parms[3])),
         MaxmMd = DS[[Dsnm]]$MaxmMd, Scrit = DS[[Dsnm]]$Scrit01,
         cor=Spcor, cor_pvalue=Spcortest, 
         convergence = out$convergence, Name = Dsnm)
}

doplots <- function(zz) {
    pdta1 <- data.frame(Concentration = exp(DS[[zz[["Name"]]]]$lx),
                        mMd = DS[[zz[["Name"]]]]$y)
    xmin <- if(is.finite(zz$BMD) & zz$BMD > 0) min(c(DS[[zz[["Name"]]]]$lx, log(zz$BMD))) else min(DS[[zz[["Name"]]]]$lx)
    pdta2 <- data.frame(x = exp(z <- seq(xmin,
                                         max(DS[[zz[["Name"]]]]$lx), length=300)),
                        y = logistic4par(z, zz$Parms["T"], zz$Parms["cc"], zz$Parms["d"]))
    p <- ggplot() +
        geom_point(data=pdta1, aes(x = Concentration, y=mMd)) +
        geom_line(data=pdta2, aes(x=x, y=y)) +
        annotate(geom="segment", x=min(pdta2$x), xend=zz$BMD, y=zz$Scrit, yend=zz$Scrit,
                 lty=2, color=brewer.pal(5, "YlGnBu")[3]) +
        annotate(geom="segment", x=zz$BMD, xend=zz$BMD, y=zz$Scrit, yend=0,
                 color=brewer.pal(5, "YlGnBu")[3]) +
        scale_x_log10() +
        ggtitle(zz[["Name"]])
    p
}



#-------------------------------------------------------------------------------------------------#
#-----Load output from mahalanobis_distance_calculation_and_Supp9.R 
#-------------------------------------------------------------------------------------------------#

load("./RData/AllResps_outliersRemoved2018-10-02.RData")

## split out data 

ChemNames <- unique(Mahalanobis_dat$date_chnm_plate)

DS <- lapply(ChemNames, function(nm) {
    ds <- Dists[Dists$CA == nm,]
    list(lx = log(ds$Conc),
         y = ds$D11P,
         Scrit01 = Mahalanobis_dat[Mahalanobis_dat$date_chnm_plate == nm,"Scrit01"],
         MaxmMd = Mahalanobis_dat[Mahalanobis_dat$date_chnm_plate == nm,"maxD11P"]
         )
    })
names(DS) <- ChemNames

## Delete the objects we don't need (to reduce the amount of stuff copied over to the child
## processes
rm(CritLim, dat, dat_mean, Dists, Mahalanobis_dat, Residuals, Resps)
gc()

## -----------------------------------------------------------------
## Do the fits and construct plots

out <- mclapply(ChemNames, function(nm) {
    zz <- try(dofit(nm))
    p <- if (!is(zz, "try-error")) {
        doplots(zz)
    } else NA
        
    list(Name = nm,
         Fit = zz,
         Plot = p)
}, mc.cores=35)

Fits <- do.call(rbind,
                lapply(out,
                       function(x) {
                           cbind(data.frame(Name=x$Name,
                                            as.data.frame(matrix(c(x$Fit$Parms,BMD=x$Fit$BMD,
                                                                   MaxmMd=x$Fit$MaxmMd,
                                                                   Scrit=x$Fit$Scrit,
                                                                   cor=x$Fit$cor,
                                                                   cor_pvalue=x$Fit$cor_pvalue,
                                                                   convergence=x$Fit$convergence),
                                                                 nrow=1,
                                                                 dimnames=list(NULL, c(names(x$Fit$Parms),
                                                                                       "BMD","MaxmMd","Scrit","cor","cor_pvalue",
                                                                                       "convergence"))))))}))

save(Fits, file=paste0("./RData/mMdFits", Sys.Date(), ".RData"))
