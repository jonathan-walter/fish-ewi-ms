## Apply EWI statistics to region-by-sampling-method time series

library(lubridate)
library(tidyr)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat_raw <- read.csv("../data/fish_abund_cleaned.csv")
spp <- unique(dat_raw$Taxa)
dat <- dat_raw[dat_raw$Region != "",]
dat$RegionMethod <- as.factor(paste(dat$Region, dat$Method, sep="-"))
dat$Year <- as.factor(year(dat$Date))
RegionMethod <- unique(dat$RegionMethod)
length(unique(dat$RegionMethod))
years <- as.numeric(levels(dat$Year))

### Organize into a list of matrices --------------------------------------------------------------

sppmats <- list()

for(ii in 1:length(spp)){
  tmp <- dat[dat$Taxa == spp[ii],]
  tmp <- aggregate(CPUE ~ RegionMethod + Year, data=tmp, FUN=mean)
  tmp <- pivot_wider(tmp, id_cols="RegionMethod", names_from="Year", values_from="CPUE", 
                     id_expand=TRUE, names_expand=TRUE)
  tmpmat <- as.matrix(tmp[,-1])
  rownames(tmpmat) <- sort(RegionMethod)
  sppmats[[paste0(spp[ii])]] <- tmpmat
}
rm(tmp, tmpmat)

saveRDS(sppmats, "../data/spp_abund_for_ewi.RDS")

### Calculate EWIs on time series -----------------------------------------------------------------

# Parameters for data filtering, moving windows
wwidth <- 5 #width of moving window (years)
dmin.all <- wwidth/length(years) # fraction of data needed for a time series to be analyzed -- currently arbitrarily low
dmin.w <- wwidth-1 #years needed in a moving window to analyze it.
dmin.space <- 3 #region-method combs needed to compute spatial CV and spatial synchrony

cv <- function(x){ #helper function for computing CV
  return(sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
}

ar1 <- function(x1,x2){ #helper function for computing lag-1 autocorrelation
  return(cor(x1, x2, use="pairwise.complete.obs"))
}

bkshft <- function(x){ #helper function used to shift indexing for lag-1 autocorrelation
  x_bkshft <- x-1
  x_bkshft[x_bkshft <=0] <- NA
  return(x_bkshft)
}

sppewi <- list()

for(ii in 1:length(sppmats)){ #start loop over species matrices
  mat <- sppmats[[ii]] #select species
  
  out <- list(cv.time=matrix(NA, nrow(mat), ncol(mat)) #initialize output
              ,cv.space=matrix(NA, 1, ncol(mat))
              ,ar1=matrix(NA, nrow(mat), ncol(mat))
              ,spatsync=matrix(NA, 1, ncol=ncol(mat)))
  
  rownames(out$cv.time) <- sort(RegionMethod)
  rownames(out$ar1) <- sort(RegionMethod)
  colnames(out$cv.time) <- years
  colnames(out$cv.space) <- years
  colnames(out$ar1) <- years
  colnames(out$spatsync) <- years
  
  for(jj in 1:nrow(mat)){ #start loop over region-method combs (rows)
    if(mean(!is.na(mat[jj,])) < dmin.all){
      next
    }
    else{
      # calculate windowed temporal CV
      for(kk in wwidth:ncol(mat)){
        if(sum(!is.na(mat[jj,(kk-wwidth+1):kk])) >= dmin.w){
          out$cv.time[jj,kk] <- cv(mat[jj,(kk-wwidth+1):kk])
        }
      }
      # calculate windowed lag-1 correlation
      for(kk in wwidth:ncol(mat)){
        if(sum(!is.na(mat[jj,(kk-wwidth+1):kk])) >= dmin.w){
          out$ar1[jj,kk] <- ar1(mat[jj,(kk-wwidth+1):kk], mat[jj, bkshft((kk-wwidth+1):kk)])
        }
      }
    }
  }
  #calculate spatial CV
  for(kk in 1:ncol(mat)){
    if(sum(!is.na(mat[,kk])) >= dmin.space){
      out$cv.space[1,kk] <- cv(mat[,kk])
    }
  }
  #calculate spatial synchrony
  for(kk in wwidth:ncol(mat)){
    nacheck <- apply(mat[,(kk-wwidth+1):kk], MARGIN=1, FUN=function(x){sum(!is.na(x))})
    matsub <- mat[nacheck >= dmin.w,(kk-wwidth+1):kk]
    if(!is.null(nrow(matsub))){
      if(nrow(matsub) >= dmin.space){
        out$spatsync[1,kk] <- mean(cor(matsub, method="spearman", use="pairwise.complete.obs")[lower.tri(cor(matsub, method="spearman", use="pairwise.complete.obs"))])
      }
    }
  }
  
  sppewi[[paste0(spp[ii])]] <- out
}

saveRDS(sppewi, "../data/ewi_abund_wwidth5.rds")


### Recompute EWIs using a different window width -------------------------------------------------

# Parameters for data filtering, moving windows
wwidth <- 7 #width of moving window (years)
dmin.all <- wwidth/length(years) # fraction of data needed for a time series to be analyzed -- currently arbitrarily low
dmin.w <- wwidth-2 #years needed in a moving window to analyze it.
dmin.space <- 5 #region-method combs needed to compute spatial CV

cv <- function(x){ #helper function for computing CV
  return(sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
}

ar1 <- function(x1,x2){ #helper function for computing lag-1 autocorrelation
  return(cor(x1, x2, use="pairwise.complete.obs"))
}

bkshft <- function(x){ #helper function used to shift indexing for lag-1 autocorrelation
  x_bkshft <- x-1
  x_bkshft[x_bkshft <=0] <- NA
  return(x_bkshft)
}

sppewi <- list()

for(ii in 1:length(sppmats)){ #start loop over species matrices
  mat <- sppmats[[ii]] #select species
  
  out <- list(cv.time=matrix(NA, nrow(mat), ncol(mat)) #initialize output
              ,cv.space=matrix(NA, 1, ncol(mat))
              ,ar1=matrix(NA, nrow(mat), ncol(mat))
              ,spatsync=matrix(NA, 1, ncol=ncol(mat)))
  
  rownames(out$cv.time) <- sort(RegionMethod)
  rownames(out$ar1) <- sort(RegionMethod)
  colnames(out$cv.time) <- years
  colnames(out$cv.space) <- years
  colnames(out$ar1) <- years
  colnames(out$spatsync) <- years
  
  for(jj in 1:nrow(mat)){ #start loop over region-method combs (rows)
    if(mean(!is.na(mat[jj,])) < dmin.all){
      next
    }
    else{
      # calculate windowed temporal CV
      for(kk in wwidth:ncol(mat)){
        if(sum(!is.na(mat[jj,(kk-wwidth+1):kk])) >= dmin.w){
          out$cv.time[jj,kk] <- cv(mat[jj,(kk-wwidth+1):kk])
        }
      }
      # calculate windowed lag-1 correlation
      for(kk in wwidth:ncol(mat)){
        if(sum(!is.na(mat[jj,(kk-wwidth+1):kk])) >= dmin.w){
          out$ar1[jj,kk] <- ar1(mat[jj,(kk-wwidth+1):kk], mat[jj, bkshft((kk-wwidth+1):kk)])
        }
      }
    }
  }
  #calculate spatial CV
  for(kk in 1:ncol(mat)){
    if(sum(!is.na(mat[,kk])) >= dmin.space){
      out$cv.space[1,kk] <- cv(mat[,kk])
    }
  }
  #calculate spatial synchrony
  for(kk in wwidth:ncol(mat)){
    nacheck <- apply(mat[,(kk-wwidth+1):kk], MARGIN=1, FUN=function(x){sum(!is.na(x))})
    matsub <- mat[nacheck >= dmin.w,(kk-wwidth+1):kk]
    if(!is.null(nrow(matsub))){
      if(nrow(matsub) >= dmin.space){
        out$spatsync[1,kk] <- mean(cor(matsub, method="spearman", use="pairwise.complete.obs")[lower.tri(cor(matsub, method="spearman", use="pairwise.complete.obs"))])
      }
    }
  }
  
  sppewi[[paste0(spp[ii])]] <- out
}


saveRDS(sppewi, "../data/ewi_abund_wwidth7.rds")
