## Calculate trends on ewi metrics and develop risk and certainty metrics

library(lubridate)
library(tidyr)
library(trend)
library(viridis)
library(lmerTest)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sppmats <- readRDS("../data/spp_abund_for_ewi.rds")
sppewi5 <- readRDS("../data/ewi_abund_wwidth5.rds")
#sppewi7 <- readRDS("../data/ewi_abund_wwidth7.rds") # replace to use different window width
taxa <- names(sppewi5)
sppInfo <- read.csv("../data/species_details.csv")
sunits <- rownames(sppewi5$`Acanthogobius flavimanus`$cv.time)



## Order sampling units by region (bay to upriver) and method ------------------
sunits.sort <- data.frame(sunits=sunits, region=NA, method=NA)
for(ii in 1:length(sunits)){
  tmp <- unlist(strsplit(sunits[ii], split="-"))
  sunits.sort$region[ii] <- tmp[1]
  sunits.sort$method[ii] <- tmp[2]
  rm(tmp)
}

sunits.sort$region <- factor(sunits.sort$region, 
                             levels=c("South Bay","Central Bay", "San Pablo", "Napa River", "Suisun Bay and Marsh",
                                      "Confluence", "South", "North", "SJ Upriver", "Sac Upriver"))
sunits.sort <- sunits.sort[order(sunits.sort$region, sunits.sort$method),]

sunits.sort$prettyname <- c("South Bay-MWT","Central Bay-BS","Central Bay-MWT",
                            "San Pablo Bay-BS","San Pablo Bay-MWT","Napa River-MWT",          
                            "Suisun Bay-MWT","Confluence-BS","Confluence-MWT",
                            "South Delta-BS","South Delta-MWT","North Delta-BS",                  
                            "North Delta-MWT","San Joaquin R.-BS","Sacramento R.-BS")

## Order species according to salinity zone, water column position, and taxonomic family
sppInfo$Salinity <- factor(sppInfo$Salinity, 
                              levels=c("Marine","Marine/Brackish","Marine/Brackish/Freshwater",
                                       "Marine/Brackish/Freshwater (anadromous)","Brackish","Brackish/Freshwater",
                                       "Brackish/Freshwater (anadromous)","Freshwater"))
sppInfo$Water.Column <- factor(sppInfo$Water.Column,
                                  levels=c("Demersal","Benthopelagic","Pelagic"))
sppInfo <- sppInfo[order(sppInfo$Salinity, sppInfo$Water.Column, sppInfo$Family),]
sppInfo <- sppInfo[sppInfo$Scientific.name %in% taxa,]


## Remove outlier EWIs ----------------------------------------------------------------------------

zmax <- 3.5 #+/-

for(ii in 1:length(sppewi5)){
  tmp <- sppewi5[[ii]]
  tmp$cv.space[scale(c(tmp$cv.space)) > zmax | scale(c(tmp$cv.space)) < -zmax] <- NA
  tmp$spatsync[scale(c(tmp$spatsync)) > zmax | scale(c(tmp$spatsync)) < -zmax] <- NA
  for(jj in 1:length(sunits)){
    tmp$cv.time[jj, scale(c(tmp$cv.time[jj,])) > zmax | scale(c(tmp$cv.time[jj,])) < -zmax] <- NA
    tmp$ar1[jj, scale(c(tmp$ar1[jj,])) > zmax | scale(c(tmp$ar1[jj,])) < -zmax] <- NA
  }
  sppewi5[[ii]] <- tmp
}

## Make some plots of EWIs and species abundances -------------------------------------------------

ymin <- 5 #number of years with ewi data
min.span <- 10 #minimum span of years with ewi data

years <- 1980:2023

# Commented out are some preliminary plots to visualize time series
# pdf("../figures/ewi_vs_time_allspp_rough.pdf", onefile=TRUE)
# 
# for(spp in taxa){
#   
#   sppmat <- sppmats[[spp]]
#   tmpewi <- sppewi5[[spp]]
#   
#   for(ii in 1:nrow(tmpewi$cv.time)){
#     if(all(is.na(tmpewi$cv.time[ii,]))){next}
#     if(sum(!is.na(tmpewi$cv.time[ii,])) >= ymin &
#        (max(which(!is.na(tmpewi$cv.time[ii,]))) - min(which(!is.na(tmpewi$cv.time[ii,])))) + 1 >= min.span 
#     ){
#       par(mfrow=c(2,1), mar=c(2,4.1,3,1))
#       plot(years, sppmat[ii,], type="b", xlab="", ylab="CPUE")
#       mtext(sppInfo$Common.name[sppInfo$Scientific.name==spp], line=1.5)
#       mtext(rownames(sppmat)[ii], line=0.5)
#       
#       plot(years, scale(tmpewi$cv.time[ii,]), type="b", xlab="", ylab="EWI (z-scored)", ylim=c(-3,3), pch=16)
#       lines(years, scale(tmpewi$ar1[ii,]), type="b", col="red", pch=16)
#       #lines(years, scale(c(tmpewi$cv.space)), type="b", col="blue")
#       lines(years, scale(c(tmpewi$spatsync)), type="b", col="blue", pch=16)
#       legend("top", col=c("black","red","blue"), legend=c("cv_time","ar1","spatsync"), ncol=4, lwd=1, bty="n")
#     }
#   }
# }
# 
# dev.off()


## Selected plots for manuscript ------------------------------------------------------------------


pdf("../figures/ewi_examples.pdf", width=6.5, height=6.5)

par(mfrow=c(3,1), mar=c(2,3.5,3.5,3.5), mgp=c(2.4,1,0), oma=c(1,0,0,0))

#Delta smelt, Confluence, MWT
plot(years, sppmats$`Hypomesus transpacificus`[rownames(sppmats[[1]])=="Confluence-Midwater trawl",], 
     type="b", lwd=1.5, xlab="", ylab="CPUE")
par(new=T)
plot(years, sppewi5$`Hypomesus transpacificus`$ar1[rownames(sppmats[[1]])=="Confluence-Midwater trawl",],
     yaxt="n", type="b", col="red", ylab="", xlab="")
axis(4, at=pretty(sppewi5$`Hypomesus transpacificus`$ar1[rownames(sppmats[[1]])=="Confluence-Midwater trawl",]),
     col="red", col.axis = "red")
text(par("usr")[1]+0.02*diff(par("usr")[1:2]), par("usr")[4]-0.05*diff(par("usr")[3:4]), "a)")

mtext("Delta smelt", line=1.6, cex=0.9)
mtext("Confluence - Midwater trawl", line=0.4, cex=0.9)
mtext("Lag-1 autocorrelation", side=4, line=2.25, cex=2/3)

#Striped bass, Confluence, MWT
plot(years, sppmats$`Morone saxatilis`[rownames(sppmats[[1]])=="Confluence-Midwater trawl",], 
     type="b", lwd=1.5, xlab="", ylab="CPUE")
par(new=T)
plot(years, sppewi5$`Morone saxatilis`$cv.time[rownames(sppmats[[1]])=="Confluence-Midwater trawl",],
     yaxt="n", type="b", col="red", ylab="", xlab="")
axis(4, at=pretty(sppewi5$`Monone saxatilis`$cv.time[rownames(sppmats[[1]])=="Confluence-Midwater trawl",]),
     col="red", col.axis = "red")
text(par("usr")[1]+0.02*diff(par("usr")[1:2]), par("usr")[4]-0.05*diff(par("usr")[3:4]), "b)")

mtext("Striped bass", line=1.6, cex=0.9)
mtext("Confluence - Midwater trawl", line=0.4, cex=0.9)
mtext("Temporal CV", side=4, line=2.25, cex=2/3)

#Mississippi Silverside, San Pablo Bay, MWT
plot(years, sppmats$`Menidia audens`[rownames(sppmats[[1]])=="Confluence-Beach seine",], 
     type="b", lwd=1.5, xlab="", ylab="CPUE")
par(new=T)
plot(years, sppewi5$`Menidia audens`$spatsync,
     yaxt="n", type="b", col="red", ylab="", xlab="")
axis(4, at=pretty(sppewi5$`Menidia audens`$spatsync),
     col="red", col.axis = "red")
text(par("usr")[1]+0.02*diff(par("usr")[1:2]), par("usr")[4]-0.05*diff(par("usr")[3:4]), "c)")

mtext("Mississippi silverside", line=1.6, cex=0.9)
mtext("Confluence - Beach seine", line=0.4, cex=0.9)
mtext("Spatial sync.", side=4, line=2.25, cex=2/3)

mtext("Year", side=1, outer=TRUE, cex=2/3)

dev.off()




## Trend in (windowed) temporal CV ----------------------------------------------------------------

trends_cv.time <- list()
ymin <- 5 #number of years with ewi data
min.span <- 10 #minimum span of years with ewi data

for(ii in 1:length(sppewi5)){
  
  tmpewi <- sppewi5[[ii]]$cv.time
  tmpout <- data.frame(Unit = rownames(tmpewi),
                       n = NA,
                       slope = NA,
                       slope_se = NA,
                       t_value = NA,
                       p_value = NA)
  for(jj in 1:nrow(tmpewi)){
    if(all(is.na(tmpewi[jj,]))){next}
    if(sum(!is.na(tmpewi[jj,])) >= ymin &
       (max(which(!is.na(tmpewi[jj,]))) - min(which(!is.na(tmpewi[jj,])))) + 1 >= min.span 
    ){
      slopetest <- summary(lm(tmpewi[jj,] ~ as.numeric(colnames(tmpewi))))
      tmpout[jj,2] <- sum(!is.na(tmpewi[jj,]))
      tmpout[jj,3:6] <- slopetest$coefficients[2,]
    }
  }
  trends_cv.time[[paste0(taxa[ii])]] <- tmpout
}



## Trend in (windowed) lag-1 autocorrelation ------------------------------------------------------

trends_ar1 <- list()
ymin <- 5 #number of years with ewi data
min.span <- 10 #minimum span of years with ewi data

for(ii in 1:length(sppewi5)){
  
  tmpewi <- sppewi5[[ii]]$ar1
  tmpout <- data.frame(Unit = rownames(tmpewi),
                       n = NA,
                       slope = NA,
                       slope_se = NA,
                       t_value = NA,
                       p_value = NA)
  for(jj in 1:nrow(tmpewi)){
    if(all(is.na(tmpewi[jj,]))){next}
    if(sum(!is.na(tmpewi[jj,])) >= ymin &
       (max(which(!is.na(tmpewi[jj,]))) - min(which(!is.na(tmpewi[jj,])))) + 1 >= min.span 
    ){
      slopetest <- summary(lm(tmpewi[jj,] ~ as.numeric(colnames(tmpewi))))
      tmpout[jj,2] <- sum(!is.na(tmpewi[jj,]))
      tmpout[jj,3:6] <- slopetest$coefficients[2,]
    }
  }
  trends_ar1[[paste0(taxa[ii])]] <- tmpout
}


## Trend in spatial synchrony ---------------------------------------------------------------------

trends_spatsync <- data.frame(Taxa=taxa,
                              n = NA,
                              slope = NA,
                              slope_se = NA,
                              t_value = NA,
                              p_value = NA)

ymin <- 5 #number of years with ewi data
min.span <- 10 #minimum span of years with ewi data

for(ii in 1:length(sppewi5)){
  
  tmpewi <- sppewi5[[ii]]$spatsync
  
  if(all(is.na(tmpewi))){next}
  if(sum(!is.na(tmpewi)) >= ymin &
     (max(which(!is.na(tmpewi))) - min(which(!is.na(tmpewi)))) + 1 >= min.span 
  ){
    slopetest <- summary(lm(c(tmpewi) ~ as.numeric(colnames(tmpewi))))
    trends_spatsync[ii,2] <- sum(!is.na(tmpewi[1,]))
    trends_spatsync[ii, 3:6] <- slopetest$coefficients[2,]
  }
}


## Make matrices of key values for plotting -----------------------------------------------

## trends in temporal cv (t values)
tval_trends_cv.time <- matrix(NA, nrow=length(taxa), ncol=length(sunits))
for(ii in 1:length(taxa)){
  tmp <- trends_cv.time[[ii]]
  tval_trends_cv.time[ii,] <- tmp$t_value
}

## trends in lag-1 correlation (t values)
tval_trends_ar1 <- matrix(NA, nrow=length(taxa), ncol=length(sunits))
for(ii in 1:length(taxa)){
  tmp <- trends_ar1[[ii]]
  tval_trends_ar1[ii,] <- tmp$t_value
}

#Exclude species that have NAs for all EWI trends
dropTaxa <- taxa[which(apply(tval_trends_ar1, MARGIN=1, function(x){all(is.na(x))}))]

tval_trends_cv.space <- trends_cv.space$t_value[!trends_cv.space$Taxa %in% dropTaxa]
tval_trends_ar1 <- tval_trends_ar1[!taxa %in% dropTaxa,]
tval_trends_cv.time <- tval_trends_cv.time[!taxa %in% dropTaxa,]
tval_trends_spatsync <- trends_spatsync$t_value[!trends_cv.space$Taxa %in% dropTaxa]

sppInfo.sub <- sppInfo[!sppInfo$Scientific.name %in% dropTaxa, ]
taxa.sub <- taxa[!taxa %in% dropTaxa]

permute.sunits <- match(sunits.sort$sunits, sunits)
permute.taxa <- match(sppInfo.sub$Scientific.name, taxa.sub)

## Do some basic plotting ------------------------------------------------------------------


## trends in spatial cv
pdf("../figures/ewi_spatsync_heatmap.pdf", width=6.5, height=3.25)
layout(matrix(1:2, nrow=2), heights=c(0.15,0.7))
par(mar=c(1.2,1,0.5,1))
image(matrix(1:50), col=turbo(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(-9,-4.5,0,4.5,9))
par(mar=c(10.5,1,1,1))
image(x=1:length(taxa.sub), z=matrix(tval_trends_spatsync[permute.taxa]), col=turbo(50), yaxt="n", xaxt="n",
      xlab="", zlim=c(-10.5, 10.5))
axis(1, at=1:length(taxa.sub), labels=sppInfo.sub$Common.name, las=2)
dev.off()


pdf("../figures/ewi_timecv_heatmap.pdf", width=6.5, height=6.5)
layout(matrix(1:2, nrow=2), height=c(0.1,0.9))
par(mar=c(1.2,8,1,1))
image(matrix(1:50), col=turbo(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(-17,-8.5,0,8.5,17))
par(mar=c(8,8,1,1))
image(x=1:length(taxa.sub), y=1:length(sunits), z=tval_trends_cv.time[permute.taxa,permute.sunits], col=turbo(50), yaxt="n", xaxt="n",
      xlab="", zlim=c(-10, 10), ylab="")
axis(1, at=1:length(taxa.sub), labels=sppInfo.sub$Common.name, las=2, cex.axis=0.78)
axis(2, at=1:length(sunits), labels=sunits.sort$prettyname, las=2, cex.axis=0.8)
dev.off()

pdf("../figures/ewi_ar1_heatmap.pdf", width=6.5, height=6.5)
layout(matrix(1:2, nrow=2), height=c(0.1,0.9))
par(mar=c(1.2,8,1,1))
image(matrix(1:50), col=turbo(50), yaxt="n", xaxt="n")
axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(-15,-7.5,0,7.5,15))
par(mar=c(8,8,1,1))
image(x=1:length(taxa.sub), y=1:length(sunits), z=tval_trends_ar1[permute.taxa,permute.sunits], col=turbo(50), yaxt="n", xaxt="n",
      xlab="", zlim=c(-10, 10), ylab="")
axis(1, at=1:length(taxa.sub), labels=sppInfo.sub$Common.name, las=2, cex.axis=0.9)
axis(2, at=1:length(sunits), labels=sunits.sort$prettyname, las=2, cex.axis=0.9)
dev.off()


scale01 <- function(x){
  x1 <- x-min(x, na.rm=TRUE)
  return(x1/max(x1, na.rm=TRUE))
}

risk_cv.space <- scale01(tval_trends_cv.space)
risk_cv.time <- scale01(apply(tval_trends_cv.time, MARGIN=1, FUN=mean, na.rm=TRUE))
risk_ar1 <- scale01(apply(tval_trends_ar1, MARGIN=1, FUN=mean, na.rm=TRUE))
risk_spatsync <- scale01(tval_trends_spatsync)


risk_comb <- scale01(colMeans(rbind(risk_spatsync, risk_cv.time, risk_ar1), na.rm=TRUE))
risk_matrix <- cbind(risk_spatsync, risk_cv.time, risk_ar1, risk_comb)

pdf("../figures/risk_heatmap.pdf", width=7, height=4)
layout(matrix(1:2, nrow=2), height=c(0.18,0.9))
par(mar=c(1.2,7.2,1,1))
image(matrix(1:50), col=rev(rocket(50)), yaxt="n", xaxt="n")
axis(1)
par(mar=c(9.2,7.2,1,1))
image(x=1:length(taxa.sub), y=1:4, z=risk_matrix[permute.taxa,], col=rev(rocket(50)), yaxt="n", xaxt="n",
      xlab="", zlim=c(0, 1), ylab="")
axis(1, at=1:length(taxa.sub), labels=sppInfo.sub$Common.name, las=2, cex.axis=0.9)
axis(2, at=1:4, labels=c("Spatial Sync.", "Temporal CV", "Lag-1 autocorr.", "Combined"), las=2, cex.axis=0.9)
abline(h=3.5, lwd=3, col="white")
dev.off()


ewi_matrix <- cbind(tval_trends_spatsync,
                    apply(tval_trends_cv.time, MARGIN=1, FUN=mean, na.rm=TRUE),
                    apply(tval_trends_ar1, MARGIN=1, FUN=mean, na.rm=TRUE))
colnames(ewi_matrix) <- c("trend_spatsync","trend_cv.time","trend_ar1")
ewi_sign_matrix <- ewi_matrix
ewi_sign_matrix[ewi_matrix<0] <- -1
ewi_sign_matrix[ewi_matrix>0] <- 1

ewi_agreement <- rep(0, length(taxa.sub))
for(ii in 1:nrow(ewi_sign_matrix)){
  ewi_agreement[ii] <- max(sum(ewi_sign_matrix[ii,]==1, na.rm=TRUE), 
                       sum(ewi_sign_matrix[ii,]==-1, na.rm=TRUE)) -1
}

ewi_agreement_scale <- ewi_agreement/(ncol(ewi_sign_matrix)-1)

cor(ewi_matrix, use="pairwise.complete.obs")


## Confidence metric ------------------------------------------------------------------------------

cv <- function(x){
  return(sd(x, na.rm=TRUE)/abs(mean(x, na.rm=TRUE)))
}

cv_trend_ar1 <- apply(tval_trends_ar1, MARGIN=1, FUN=sd, na.rm=TRUE)
cv_trend_ar1[is.na(cv_trend_ar1)] <- max(cv_trend_ar1, na.rm=TRUE) #assign NAs to the max
cv_trend_cv.time <- apply(tval_trends_cv.time, MARGIN=1, FUN=sd, na.rm=TRUE)
cv_trend_cv.time[is.na(cv_trend_cv.time)] <- max(cv_trend_cv.time, na.rm=TRUE)

n_sunits <- apply(tval_trends_ar1, MARGIN=1, FUN=function(x){sum(!is.na(x))})

cv_trend_ar1_scale <- scale01(-cv_trend_ar1)
cv_trend_cv.time_scale <- scale01(-cv_trend_cv.time)
n_sunits_scale <- scale01(n_sunits)

abund_dat <- read.csv("../data/fish_abund_cleaned.csv")
nobs <- aggregate(Count ~ Taxa, data=abund_dat, FUN="length")
nobs <- nobs[!nobs$Taxa %in% dropTaxa,]
nobs_scale <- scale01(log10(nobs$Count))

confidence_comb <- scale01(colMeans(rbind(ewi_agreement_scale, cv_trend_ar1_scale,
                                               cv_trend_cv.time_scale, n_sunits_scale, 
                                          nobs_scale)))

confidence_matrix <- cbind(nobs_scale, ewi_agreement_scale, cv_trend_cv.time_scale, 
                           cv_trend_ar1_scale, confidence_comb)

pdf("../figures/confidence_heatmap.pdf", width=7.2, height=4)
layout(matrix(1:2, nrow=2), height=c(0.18,0.9))
par(mar=c(1.2,8.5,1,1))
image(matrix(1:50), col=viridis(50), yaxt="n", xaxt="n")
axis(1)
par(mar=c(9.2,8.5,1,1))
image(x=1:length(taxa.sub), y=1:5, z=confidence_matrix[permute.taxa,], col=viridis(50), yaxt="n", xaxt="n",
      xlab="", zlim=c(0, 1), ylab="")
axis(1, at=1:length(taxa.sub), labels=sppInfo.sub$Common.name, las=2, cex.axis=0.9)
axis(2, at=1:5, labels=c("Sample number", "EWI trend agreement", "cv trend spat. var.",
                         "ar1 trend spat. var.", "Combined"), las=2, cex.axis=0.9)
abline(h=4.5, lwd=3, col="white")
dev.off()


#quadrant 1 - high risk, high confidence
sppInfo.sub$Common.name[risk_comb[permute.taxa]>=median(risk_comb) & confidence_comb[permute.taxa]>=median(confidence_comb)]
#quadrant 2 - high risk, low confidence
sppInfo.sub$Common.name[risk_comb[permute.taxa]>=median(risk_comb) & confidence_comb[permute.taxa]<median(confidence_comb)]
#quadrant 3 - low risk, low confidence
sppInfo.sub$Common.name[risk_comb[permute.taxa]<median(risk_comb) & confidence_comb[permute.taxa]<median(confidence_comb)]
#quadrant 4 - low risk, high confidence
sppInfo.sub$Common.name[risk_comb[permute.taxa]<median(risk_comb) & confidence_comb[permute.taxa]>=median(confidence_comb)]

cor.test(risk_comb, confidence_comb, use="pairwise.complete.obs")

pdf("../figures/risk_confidence_biplot.pdf", width=6.5, height=4)
par(mfrow=c(1,2), mar=c(3.1,3.1,0.6,0.6), mgp=c(2,0.7,0))
plot(confidence_comb, risk_comb, xlab="Relative confidence", ylab="Relative EWI score", pch=16, cex=0.75)
abline(h=median(risk_comb), col="grey")
abline(v=median(confidence_comb), col="grey")
text(x=1, y=median(risk_comb) + 0.025, "I", cex=0.9)
text(x=0, y=median(risk_comb) + 0.025, "II", cex=0.9)
text(x=0, y=median(risk_comb) - 0.025, "III", cex=0.9)
text(x=1, y=median(risk_comb) - 0.025, "IV", cex=0.9)
text(0.01,1,"a)")

par(mgp=c(1,0.7,0))
plot(NA, NA, ylim=c(0,1), xlim=c(0,1), xlab="Relative confidence", ylab="Relative EWI score",
     xaxt="n", yaxt="n")
abline(h=0.5, col="grey")
abline(v=0.5, col="grey")
text(x=par("usr")[1]+0.75*diff(par("usr")[1:2]), 
     y=par("usr")[3]+0.75*diff(par("usr")[3:4]), 
     "Shiner perch\nTopsmelt\nYellowfin goby\nMississippi silverside\nThreadfin shad\nAmerican shad\nRed shiner",
     cex=0.7)
text(x=par("usr")[1]+0.25*diff(par("usr")[1:2]),
     y=par("usr")[3]+0.75*diff(par("usr")[3:4]), 
     "White croaker\nSardine\nPlainfin midshipman\nBay pipefish\nShimofuri goby\nStarry flounder\nWhite sturgeon\nTule perch",
     cex=0.7)
text(x=par("usr")[1]+0.25*diff(par("usr")[1:2]), 
     y=par("usr")[3]+0.25*diff(par("usr")[3:4]), 
     "Pacific staghorn sculpin\nPacific pompano\nChinook salmon\nDelta smelt\nChannel catfish\nFathead minnow",
     cex=0.7)
text(x=par("usr")[1]+0.75*diff(par("usr")[1:2]), 
     y=par("usr")[3]+0.25*diff(par("usr")[3:4]), 
     "Jacksmelt\nNorthern anchovy\nStriped bass\nLongfin smelt\nSacramento splittail\nWhite catfish\nBluegill",
     cex=0.7)
text(0.01,1,"b)")
dev.off()


cor.test(confidence_comb, risk_comb)

## Calculate abundance trends and compare to risk -------------------------------------------------


abund_dat <- abund_dat[abund_dat$Region != "",]
abund_dat$RegionMethod <- as.factor(paste(abund_dat$Region, abund_dat$Method, sep="-"))
abund_dat <- abund_dat[abund_dat$Taxa %in% taxa.sub,]
abund_dat$Year <- year(abund_dat$Date)

betas <- rep(NA, length(taxa.sub))
tvals <- rep(NA, length(taxa.sub))

for(ii in 1:length(taxa.sub)){
  
  tmp <- abund_dat[abund_dat$Taxa==taxa.sub[ii],]
  mod <- lmer(log(CPUE) ~ Year + (1|RegionMethod), data=tmp)
  modsum <- summary(mod)
  hist(resid(mod), main=taxa.sub[ii])
  betas[ii] <- modsum$coefficients[2,1]
  tvals[ii] <- modsum$coefficients[2,4]
  
}

plot(risk_comb, tvals)
cor.test(tvals, risk_comb, method="pearson")


## Native/non-native status -----------------------------------------------------------------------

boxplot(risk_comb[permute.taxa] ~ sppInfo.sub$Native, ylab="Relative EWI score", xlab="", 
        names=c("Non-native","Native"))

summary(aov(risk_comb[permute.taxa] ~ sppInfo.sub$Native))

pdf("../figures/summarize_ewi_by_spp.pdf", width=6.5, height=3.25)
par(mfrow=c(1,2), mar=c(4.1,4.1,1,1))
plot(tvals, risk_comb, xlab="Population trend (t-value)", ylab="Relative EWI score", pch=16)
#mtext("Pearson r = 0.07, p = 0.71", line=0.1)
text(par("usr")[1] + 0.05*diff(par("usr")[1:2]),
     par("usr")[4] - 0.05*diff(par("usr")[3:4]), "a)")
boxplot(risk_comb[permute.taxa] ~ sppInfo.sub$Native, ylab="Relative EWI score", xlab="", 
        names=c("Non-native","Native"))
text(par("usr")[1] + 0.05*diff(par("usr")[1:2]),
     par("usr")[4] - 0.05*diff(par("usr")[3:4]), "b)")
dev.off()