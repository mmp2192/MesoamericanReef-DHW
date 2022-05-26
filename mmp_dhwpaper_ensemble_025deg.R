##############################################################################
# mmp_dhwpaper_ensemble_025deg
#
# SUMMARY: This script creates ensemble line plots from CMIP5 and CMIP6
# sea surface temperature data for the WWF Mesoamerican Reef project
# and corresponding paper with GISS Model-E. This version uses metrics
# calculated from 0.25-degree resolution SST data.
#
# -- Creates a line plot showing the seasonal cycle of DHW for each
# month of the year from the ensemble of GCMs for the 2050-2059 period
# with envelopes showing range of GCM outputs for each ssp245 and ssp585.
#
# -- Creates a line plot showing the seasonal cycle of DHW for each
# month of the year from MERRA-2 as an average across all baseline years
# (1980-2009) and for 2005 alone.
#
# -- Creates a line plot showing the maximum monthly delta value from
# each GCM for all scenarios and decades.
#
# -- Creates a line plot showing the average maximum annual DHW value
# for the baseline using MERRA-2 and for each future decade with envelopes
# showing range of GCM outputs for each ssp245 and ssp585, plus a version
# with extended MERRA-2 data out to 2019.
#
#
#
# 						author: Meridel Phillips
# 						created: August 3, 2021
#
##############################################################################

#### required libraries

library(pracma)
library(MASS)
library(raster)
library(gdata)
library(R.matlab)
library(R.utils)
library(ggplot2)
library(RNetCDF)
library(ncdf4)
library(PerformanceAnalytics)
library(chron)
library(RColorBrewer)
library(maps)
library(rworldmap)
library(mapdata)
library(ggmap)
library(maptools)
library(lubridate)
library(dplyr)
library(stringr)
library(reshape2)
library(rgdal)
library(fields)
library(zoo)
library(timeSeries)
library(broom)

rm(list = ls(all = T))

#### define root directory and set working directory

rootDir <- "/Users/mmphill2/"
setwd(paste(rootDir, "datasets/MarineHeatWaves/", sep = ""))

slices = c('2010s','2020s','2030s','2040s','2050s','2060s','2070s','2080s','2090s')
slicestart = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090)
sliceend = c(2019, 2029, 2039, 2049, 2059, 2069, 2079, 2089, 2099)

allgcm = c('ACCESS1-0','GISS-E2-R','NorESM1-M','ACCESS-CM2','CanESM5','GFDL-CM4','GISS_E213','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM')
allgcmsave = c('ACCESS1-0','GISS-E2-R','NorESM1-M','ACCESS-CM2','CanESM5','GFDL-CM4','GISS-E2-1-G','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM')

allensemble = c('r1i1p1','r6i1p1','r1i1p1','r1i1p1f1','r1i1p1f1','r1i1p1f1','none','r1i1p1f1','r1i1p1f1','r1i1p1f1')

allscen = c('rcp45','rcp85')

#####################################################################
#### average seasonal cycle of DHW in MERRA-2

allyears = c(1980:2019)
merraregionaldhw = array(data = NA, dim = c(length(allyears),105))
for (xyear in c(1:length(allyears))){
   thisyear = allyears[xyear]
   merraregionaldhwyymat = readMat(paste("MERRA2_RegionalDHW_HalfWeekly_", thisyear, ".mat",sep=""))
   data = merraregionaldhwyymat$merraregionaldhw
   merraregionaldhw[xyear,] = data
}

merraregionaldhw[is.na(merraregionaldhw)] = 0

# keep this to 1980-2009
merrayear = merraregionaldhw[1:30,]

plotdhw = apply(merrayear,2,mean)

meangcmyear = array(data = NA, dim = c(length(allgcm),length(allscen),30,105))

for (xgcm in c(1:length(allgcm))){
	  thisgcm = allgcm[xgcm]
    if (xgcm<4){
       thisscen = "rcp45"
    }
    if (xgcm>3){
       thisscen = "ssp245"
    }
    for (xyear in c(1:30)){
	     gcmregionaldhwyymat = readMat(paste(thisgcm, "_RegionalDHW_HalfWeekly_", thisscen, "_2050-2059_yy", xyear, ".mat", sep = ""))
       gcmregionaldhw = gcmregionaldhwyymat$regionaldhw
       meangcmyear[xgcm,1,xyear,] = gcmregionaldhw
    }

    if (xgcm<4){
       thisscen = "rcp85"
    }
    if (xgcm>3){
       thisscen = "ssp585"
    }
    for (xyear in c(1:30)){
	     gcmregionaldhwyymat = readMat(paste(thisgcm, "_RegionalDHW_HalfWeekly_", thisscen, "_2050-2059_yy", xyear, ".mat", sep = ""))
       gcmregionaldhw = gcmregionaldhwyymat$regionaldhw
       meangcmyear[xgcm,2,xyear,] = gcmregionaldhw
    }

}

meangcmyear[is.na(meangcmyear)] = 0
allgcmyear = meangcmyear

allgcm45 = allgcmyear[ ,1, , ]
allgcm85 = allgcmyear[ ,2, , ]
plotgcm45 = apply(allgcm45,c(1,3),mean)
plotgcm85 = apply(allgcm85,c(1,3),mean)
low45 = apply(plotgcm45,2,min)
mid45 = apply(plotgcm45,2,median)
high45 = apply(plotgcm45,2,max)
low85 = apply(plotgcm85,2,min)
mid85 = apply(plotgcm85,2,median)
high85 = apply(plotgcm85,2,max)

seasonaldf <- data.frame(years=c(1:105), hist=NA, rcp45=NA, rcp45low=NA, rcp45high=NA, rcp85=NA, rcp85low=NA, rcp85high=NA)
seasonaldf[,2] = plotdhw
seasonaldf[,3] = mid45
seasonaldf[,4] = low45
seasonaldf[,5] = high45
seasonaldf[,6] = mid85
seasonaldf[,7] = low85
seasonaldf[,8] = high85

# plot all years and scenarios with envelope
ggplot(seasonaldf, aes(x=seasonaldf$years, y = value, color = variable)) +
geom_line(aes(y = seasonaldf$hist, col = "MERRA-2 (1980-2009)"),size=2) +
geom_line(aes(y = seasonaldf$rcp45, col = "RCP4.5/SSP245 (2050-2059)")) +
geom_line(aes(y = seasonaldf$rcp85, col = "RCP8.5/SSP585 (2050-2059)")) +
scale_colour_manual(values = c("MERRA-2 (1980-2009)"="black", "RCP4.5/SSP245 (2050-2059)"="darkorange", "RCP8.5/SSP585 (2050-2059)"="red3")) +
geom_ribbon(aes(x=seasonaldf$years, ymin = seasonaldf$rcp45low, ymax = seasonaldf$rcp45high, fill = "RCP4.5/SSP245 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
geom_ribbon(aes(x=seasonaldf$years, ymin = seasonaldf$rcp85low, ymax = seasonaldf$rcp85high, fill = "RCP8.5/SSP585 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
scale_fill_manual(values = c("RCP4.5/SSP245 Min-Max"="orange", "RCP8.5/SSP585 Min-Max"="coral")) +
scale_x_continuous(breaks=seq(4.333,104,8.666), labels=c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"), limits = c(0,104), expand = c(0,0)) +
labs(y = "degree heating weeks (degrees C)", x = " ", title = "Average Seasonal Cycle of DHW", color = NULL) +
theme_bw() +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), legend.position = c(0.36,0.8), legend.background = element_rect(color = "black", fill = "white", linetype="solid"), axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18)) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 31)) +
geom_hline(yintercept=4, linetype="dashed",color="grey",size=1) +
  geom_hline(yintercept=8, linetype="dashed",color="grey",size=1) +
guides(fill=F)
ggsave("/Users/mmphill2/Desktop/SeasonalDHW_AllGCM_MAR.tiff", units = "in", width = 10, height = 7, dpi = 600)


#####
#### only mean seasonal DHW from MERRA-2

pdf(paste("/Users/mmphill2/Desktop/SeasonalDHW_MERRA2_MAR.pdf", sep=""), width = 12, height = 10, pointsize=12)
plot(1:105,plotdhw,'l', lty=1, lwd=3, col="black", xlim=c(1,105), xaxt="n", xlab = " ", ylab= "degree heating weeks (degrees C)", main = "MERRA-2 Mean Half-Weekly DHW 1980-2009",cex.axis=1.5,cex.lab=1.5)
axis(1,at=seq(4.333,104,8.666),labels=c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"),cex.axis=1.3,cex.lab=1.3)
dev.off()

#####
#### only 2005 half-weekly DHW from MERRA-2

plotdhw = merrayear[26,1:105]

plot(1:105,plotdhw,'l', lty=1, lwd=2, col="red3", xlim=c(1,105), xaxt="n", xlab = " ", ylab= "degree heating weeks (degrees C)", main = "MERRA-2 Half-Weekly DHW for 2005")
axis(1,at=seq(4.333,104,8.666),labels=c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"))

#####################################################################
#### maximum monthly delta value for each decade, model and scenario

maxdeltas = array(data=NA, dim=c(length(allgcm),length(slices),2))
deltaid = array(data=NA, dim=c(length(allgcm),length(slices),2))

for (xgcm in c(1:length(allgcm))){
	  thisgcm = allgcm[xgcm]
    if (xgcm<4){
       allscen = c("rcp45", "rcp85")
    }
    if (xgcm>3){
       allscen = c("ssp245", "ssp585")
    }

    regionaldeltas = array(data=NA, dim=c(length(slices),length(allscen),12))

    for (xscen in c(1:length(allscen))){
        thisscen = allscen[xscen]
        for (tt in c(1:length(slices))){
           thisslice = slices[tt]
           data = readMat(paste(thisgcm, "_", thisscen, "_", thisslice, "_RegionalDeltas.mat",sep=""));
           data = data$regionaldeltas
           regionaldeltas[tt,xscen,] = data
        }
    }

    for (xscen in c(1:length(allscen))){
        for (tt in c(1:length(slices))){

             md = max(regionaldeltas[tt,xscen,],na.rm=T)
             mdix = which(regionaldeltas[tt,xscen,] == md)

             maxdeltas[xgcm,tt,xscen] = md
             deltaid[xgcm,tt,xscen] = mdix
        }
    }
}


### plot maximum deltas

colorbar <- colorRampPalette(c("blue", "red3", "darkgreen", "dodgerblue", "purple", "darkorange2", "violetred1", "turquoise4", "plum4", "limegreen"))(n = 10)

pdf(paste("/Users/mmphill2/Desktop/MaximumDeltas_AllGCM.pdf", sep=""), width = 12, height = 10, pointsize=12)
for (xgcm in c(1:length(allgcm))){
    plot(1:length(slices),maxdeltas[xgcm, ,1],'l', lty=2, lwd=1.2, col=colorbar[xgcm], ylim=c(0.25,5), xlim=c(1.1,8.9),xaxt="n",xlab = " ",ylab= " ")
    par(new=T)
    plot(1:length(slices),maxdeltas[xgcm, ,2],'l', lty=1, lwd=1.2, col=colorbar[xgcm], ylim=c(0.25,5), xlim=c(1.1,8.9),xaxt="n",xlab=" ", ylab = " ")
    par(new=T)
}
plot(1:length(slices),c(1:length(slices))*NA,'l', lty=2, lwd=1.2, col="black",ylim=c(0.25,5), xlim=c(1.1,8.9),xaxt="n",xlab = " ",ylab = " ")
par(new=T)
plot(1:length(slices),c(1:length(slices))*NA,'l', lty=1, lwd=1.2, col="black",ylim=c(0.25,5), xlim=c(1.1,8.9),xaxt="n",xlab = " ",ylab = "maximum monthly regional delta value from baseline mean",main = "Maximum Monthly Delta Value For MAR Region", cex.lab=1.3)
axis(1,at=c(1.1,2,3,4,5,6,7,8,8.9),labels=c("2010s","2020s","2030s","2040s","2050s","2060s","2070s","2080s","2090s"),cex.axis=1.3,cex.lab=1.3)
legend(1.2, 5.1, legend = c('ACCESS1-0','GISS-E2-R','NorESM1-M','ACCESS-CM2','CanESM5','GFDL-CM4','GISS-E2.1-G','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM','RCP45/SSP245','RCP85/SSP585'), col = c("blue", "red3", "darkgreen", "dodgerblue", "purple", "darkorange2", "violetred1", "turquoise4", "plum4", "limegreen","black","black"),lty=c(1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,2,1.2),cex=1.1,ncol=2)
dev.off()

#####################################################################
#### line plots of regional annual DHW from MERRA-2 out into future
#### single point for each slice/decade for given model and scenario

merraregionaldhwyymat = readMat("MERRA2_RegionalDHW_AnnualMax.mat")
merraregionaldhwyy = merraregionaldhwyymat$merraregionaldhwyy
allgcmdhwyy = array(data=NA, dim = c(length(allgcm),length(slices),length(allscen),30))

for (xgcm in c(1:length(allgcm))){
    thisgcm = allgcm[xgcm]
    if (xgcm<4){
        allscen = c("rcp45", "rcp85")
    }
    if (xgcm>3){
        allscen = c("ssp245", "ssp585")
    }
    for (xscen in c(1:length(allscen))){
       thisscen = allscen[xscen]
       for (tt in c(1:length(slices))){
          thisslice = slices[tt]
          thisslicestart = slicestart[tt]
          thissliceend = sliceend[tt]

          for (xyear in c(1:30)){
              gcmregionaldhwyymat = readMat(paste(thisgcm, "_RegionalDHW_HalfWeekly_", thisscen, "_", thisslicestart, "-", thissliceend, "_yy", xyear, ".mat", sep = ""))
              gcmregionaldhw = gcmregionaldhwyymat$regionaldhw
              allgcmdhwyy[xgcm,tt,xscen,xyear] = max(gcmregionaldhw,na.rm=T)
           }
        }

    }
}

meangcmdhw = apply(allgcmdhwyy,c(1,2,3),mean)

modelmean = apply(meangcmdhw,c(2,3),mean)
modelmax = apply(meangcmdhw,c(2,3),max)
modelmin = apply(meangcmdhw,c(2,3),min)

merraregionaldhwyy[is.na(merraregionaldhwyy)] = 0
merrasmoothed = array(data=NA, dim=c(1,4))
merrasmoothed[1] = mean(merraregionaldhwyy[1:10],na.rm=T)
merrasmoothed[2] = mean(merraregionaldhwyy[11:20],na.rm=T)
merrasmoothed[3] = mean(merraregionaldhwyy[21:30],na.rm=T)
merrasmoothed[4] = mean(merraregionaldhwyy[31:40],na.rm=T)

df <- data.frame(years=seq(1985,2100,10), hist=NA, rcp45=NA, rcp45low=NA, rcp45high=NA, rcp85=NA, rcp85low=NA, rcp85high=NA)
df[1:3,2] = merrasmoothed[1,1:3]
df[1:3,4] = merrasmoothed[1,1:3]
df[1:3,5] = merrasmoothed[1,1:3]
df[1:3,7] = merrasmoothed[1,1:3]
df[1:3,8] = merrasmoothed[1,1:3]
df[4:12,3] = modelmean[,1]
#df$rcp45interp <- with(df, interp1(years, rcp45, years, "linear"))
df[4:12,4] = modelmin[,1]
#df$rcp45lowinterp <- with(df, interp1(years, rcp45low, years, "linear"))
df[4:12,5] = modelmax[,1]
#df$rcp45highinterp <- with(df, interp1(years, rcp45high, years, "linear"))
df[4:12,6] = modelmean[,2]
#df$rcp85interp <- with(df, interp1(years, rcp85, years, "linear"))
df[4:12,7] = modelmin[,2]
#df$rcp85lowinterp <- with(df, interp1(years, rcp85low, years, "linear"))
df[4:12,8] = modelmax[,2]
#df$rcp85highinterp <- with(df, interp1(years, rcp85high, years, "linear"))

# plot all years and scenarios with envelope
ggplot(df, aes(x=df$years, y = value, color = variable)) +
geom_line(aes(y = df$hist, col = "MERRA-2")) +
geom_line(aes(y = df$rcp45, col = "RCP4.5/SSP245")) +
geom_line(aes(y = df$rcp85, col = "RCP8.5/SSP585")) +
scale_colour_manual(values = c("MERRA-2"="black", "RCP4.5/SSP245"="darkorange", "RCP8.5/SSP585"="red3")) +
geom_ribbon(aes(x=df$years, ymin = df$rcp45low, ymax = df$rcp45high, fill = "RCP4.5/SSP245 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
geom_ribbon(aes(x=df$years, ymin = df$rcp85low, ymax = df$rcp85high, fill = "RCP8.5/SSP585 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
scale_fill_manual(values = c("RCP4.5/SSP245 Min-Max"="orange", "RCP8.5/SSP585 Min-Max"="coral")) +
scale_x_continuous(breaks=seq(1980,2090,20), limits = c(1980,2095), expand = c(0,0)) +
labs(y = "Annual DHW (degrees C)", x = " ", title = "Annual DHW for MERRA-2 Baseline + Future GCMs (MAR Region)", color = NULL) +
theme_bw() +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), legend.position = c(0.2,0.8), legend.background = element_rect(color = "black", fill = "white", linetype="solid"), axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18)) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
geom_hline(yintercept=4, linetype="dashed",color="grey",size=1) +
  geom_hline(yintercept=8, linetype="dashed",color="grey",size=1) +
guides(fill=F)
ggsave("/Users/mmphill2/Desktop/AnnualDHW_AllGCM_MAR.tiff", units = "in", width = 10, height = 7, dpi = 600)

## version with fourth MERRA point

df <- data.frame(years=seq(1985,2100,10), hist=NA, rcp45=NA, rcp45low=NA, rcp45high=NA, rcp85=NA, rcp85low=NA, rcp85high=NA)
df[1:4,2] = merrasmoothed[1,1:4]
df[1:3,4] = merrasmoothed[1,1:3]
df[1:3,5] = merrasmoothed[1,1:3]
df[1:3,7] = merrasmoothed[1,1:3]
df[1:3,8] = merrasmoothed[1,1:3]
df[4:12,3] = modelmean[,1]
#df$rcp45interp <- with(df, interp1(years, rcp45, years, "linear"))
df[4:12,4] = modelmin[,1]
#df$rcp45lowinterp <- with(df, interp1(years, rcp45low, years, "linear"))
df[4:12,5] = modelmax[,1]
#df$rcp45highinterp <- with(df, interp1(years, rcp45high, years, "linear"))
df[4:12,6] = modelmean[,2]
#df$rcp85interp <- with(df, interp1(years, rcp85, years, "linear"))
df[4:12,7] = modelmin[,2]
#df$rcp85lowinterp <- with(df, interp1(years, rcp85low, years, "linear"))
df[4:12,8] = modelmax[,2]
#df$rcp85highinterp <- with(df, interp1(years, rcp85high, years, "linear"))

# plot all years and scenarios with envelope
ggplot(df, aes(x=df$years, y = value, color = variable)) +
geom_line(aes(y = df$hist, col = "MERRA-2")) +
geom_line(aes(y = df$rcp45, col = "RCP4.5/SSP245")) +
geom_line(aes(y = df$rcp85, col = "RCP8.5/SSP585")) +
scale_colour_manual(values = c("MERRA-2"="black", "RCP4.5/SSP245"="darkorange", "RCP8.5/SSP585"="red3")) +
geom_ribbon(aes(x=df$years, ymin = df$rcp45low, ymax = df$rcp45high, fill = "RCP4.5/SSP245 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
geom_ribbon(aes(x=df$years, ymin = df$rcp85low, ymax = df$rcp85high, fill = "RCP8.5/SSP585 Min-Max"), alpha = 0.25, inherit.aes=FALSE) +
scale_fill_manual(values = c("RCP4.5/SSP245 Min-Max"="orange", "RCP8.5/SSP585 Min-Max"="coral")) +
scale_x_continuous(breaks=seq(1980,2090,20), limits = c(1980,2095), expand = c(0,0)) +
labs(y = "Annual DHW (degrees C)", x = " ", title = "Annual DHW for MERRA-2 Baseline + Future GCMs (MAR Region)", color = NULL) +
theme_bw() +
theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5), legend.position = c(0.2,0.8), legend.background = element_rect(color = "black", fill = "white", linetype="solid"), axis.text=element_text(size=18), axis.title=element_text(size=18), legend.text=element_text(size=18)) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
geom_hline(yintercept=4, linetype="dashed",color="grey",size=1) +
  geom_hline(yintercept=8, linetype="dashed",color="grey",size=1) +
guides(fill=F)
ggsave("/Users/mmphill2/Desktop/AnnualDHW_AllGCM_MAR_extendedHIST.tiff", units = "in", width = 9, height = 6, dpi = 300)
