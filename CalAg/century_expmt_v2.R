library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final/CalAg'
dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/long term studies'
rm_missing_plots <- TRUE #plot 8-3 will be removed since missing from 8-3
list.files(dataDir)
ltras_cn_1993 <- read.csv(file.path(dataDir, 'ltras_1993.csv'), stringsAsFactors = FALSE)
length(unique(ltras_cn_1993$Plot)) #53 plots
unique(ltras_cn_1993$Depth)
# ltras_cn_1993$Plot <- sub('-', '_', ltras_cn_1993$Plot)
tapply(ltras_cn_1993$PctC, ltras_cn_1993$Plot, function(x) length(x)) #all have five depths except 1-2
ltras_cn_2012 <- read.csv(file.path(dataDir, 'ltras_2012.csv'), stringsAsFactors = FALSE)
if(rm_missing_plots){
  ltras_cn_2012 <- ltras_cn_2012[-which(ltras_cn_2012$Plot=='8_3'), ]
  ltras_cn_2012 <- ltras_cn_2012[-which(ltras_cn_2012$Plot=='1_2' & ltras_cn_2012$Upper_Depth>=30), ]
  ltras_cn_2012$PctC[ltras_cn_2012$Plot=='1_8'& ltras_cn_2012$Depth=="000-015"] <- NA
}
length(unique(ltras_cn_2012$Plot)) #53 plots with edits
unique(ltras_cn_2012$Depth)
ltras_cn_2012[!(ltras_cn_2012$Plot %in% ltras_cn_1993$Plot), ] #8_3 missing from 1993
ltras_cn_1993[ltras_cn_1993$Plot=='1_2',]
ltras_cn_2012[ltras_cn_2012$Plot=='1_2',]
ltras_cn_1993[ltras_cn_1993$Plot=='1_8',]
ltras_cn_2012[ltras_cn_2012$Plot=='1_8',]
trtmt_list <- unique(ltras_cn_2012$Treatment)
trtmt_list

for(i in seq_along(trtmt_list)) {
  print(trtmt_list[i])
  print(unique(ltras_cn_1993$Plot[ltras_cn_1993$Treatment==trtmt_list[i]]))
}

#0-15 soil %C
soilC_0_15_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='000-015'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='000-015'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='000-015'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='000-015'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_1993)
soilC_0_15_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_0_15_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_0_15_1993$Treatment)

#15-30 soil C
soilC_15_30_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='015-030'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='015-030'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='015-030'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='015-030'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_1993)
soilC_15_30_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_15_30_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_15_30_1993$Treatment)

#30-60 soil C
soilC_30_60_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='030-060'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='030-060'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='030-060'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='030-060'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_1993)
soilC_30_60_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_30_60_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_30_60_1993$Treatment)

#60-100 soil C
soilC_60_100_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='060-100'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='060-100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='060-100'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='060-100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_1993)
soilC_60_100_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_60_100_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_60_100_1993$Treatment)

#100-150 soil C
soilC_100_150_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='100-150'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='100-150'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$PctC[ltras_cn_1993$Depth=='100-150'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='100-150'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_100_150_1993)
soilC_100_150_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_100_150_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_100_150_1993$Treatment)

#same for 2012 data
#0-15 soil %C
soilC_0_15_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='000-015'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='000-015'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='000-015'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='000-015'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_2012)
soilC_0_15_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_0_15_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_0_15_2012$Treatment)

#15-30 soil C
soilC_15_30_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='015-030'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='015-030'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='015-030'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='015-030'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_2012)
soilC_15_30_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_15_30_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_15_30_2012$Treatment)

#30-60 soil C
soilC_30_60_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='030-060'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='030-060'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='030-060'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='030-060'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_2012)
soilC_30_60_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_30_60_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_30_60_2012$Treatment)

#60-100 soil C
soilC_60_100_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='060-100'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='060-100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='060-100'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='060-100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_2012)
soilC_60_100_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_60_100_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_60_100_2012$Treatment)

#100-150 soil C
soilC_100_150_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='100-150'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='100-150'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$PctC[ltras_cn_2012$Depth=='100-150'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='100-150'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_100_150_2012)
soilC_100_150_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_100_150_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_100_150_2012$Treatment)

#consolidate 1993 data
soilC_1993 <- soilC_0_15_1993
soilC_1993 <- soilC_1993[,c(1,3,2)]
colnames(soilC_1993)[3] <- 'pctC_0_15'
soilC_1993$pctC_15_30 <- soilC_15_30_1993$soilC[match(soilC_1993$Plot, soilC_15_30_1993$Plot)]
soilC_1993$pctC_30_60 <- soilC_30_60_1993$soilC[match(soilC_1993$Plot, soilC_30_60_1993$Plot)]
soilC_1993$pctC_60_100 <- soilC_60_100_1993$soilC[match(soilC_1993$Plot, soilC_60_100_1993$Plot)]
soilC_1993$pctC_100_150 <- soilC_100_150_1993$soilC[match(soilC_1993$Plot, soilC_100_150_1993$Plot)]
lapply(soilC_1993[,3:7], summary)
lapply(soilC_1993[,3:7], hist)
lapply(soilC_1993[,3:7], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], mean))
soilC_1993_means <- do.call(cbind, lapply(soilC_1993[,3:7], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], mean)))
soilC_1993_sd <- do.call(cbind, lapply(soilC_1993[,3:7], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], sd)))


#consolidate 2012 data
soilC_2012 <- soilC_0_15_2012
soilC_2012 <- soilC_2012[,c(1,3,2)]
colnames(soilC_2012)[3] <- 'pctC_0_15'
soilC_2012$pctC_15_30 <- soilC_15_30_2012$soilC[match(soilC_2012$Plot, soilC_15_30_2012$Plot)]
soilC_2012$pctC_30_60 <- soilC_30_60_2012$soilC[match(soilC_2012$Plot, soilC_30_60_2012$Plot)]
soilC_2012$pctC_60_100 <- soilC_60_100_2012$soilC[match(soilC_2012$Plot, soilC_60_100_2012$Plot)]
soilC_2012$pctC_100_150 <- soilC_100_150_2012$soilC[match(soilC_2012$Plot, soilC_100_150_2012$Plot)]
lapply(soilC_2012[,3:7], summary)
lapply(soilC_2012[,3:7], hist)
lapply(soilC_2012[,3:7], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], mean))
soilC_2012_means <- do.call(cbind, lapply(soilC_2012[,3:7], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], mean)))
soilC_2012_sd <- do.call(cbind, lapply(soilC_2012[,3:7], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], sd)))
soilC_2012_means
soilC_2012_sd
View(soilC_2012)
depths <- c(7.5,22.5,45,80,125)

#maize vs. tomato comparisons
cex.setting <- 0.6
# tiff(file = file.path(FiguresDir, 'soilC_change_ltras.tif'), family = 'Times New Roman', width = 4.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
# par(mar=c(3.5,3.5,0.5,0.5))
# plot(soilC_1993_means['OMT',], -depths, type='b', xlim=c(0.33,1.4), col='tan4', yaxt='n', ylab='', xlab='', lty=2, pch=1, cex=cex.setting)
# axis(side = 2, at=-c(0,15,30,60,100), labels=c(0,15,30,60,100))
# mtext('Depth (cm)', 2, 2.25)
# mtext('Soil organic carbon (%)', 1, 2.25)
# lines(soilC_2012_means['OMT',], -depths, type='b', xlim=c(0.1,1.4), col='tan4', pch=16, cex=cex.setting)
# lines(soilC_1993_means['CMT',], -depths, type='b', col='deepskyblue', pch=2, lty=2, cex=cex.setting)
# lines(soilC_2012_means['CMT',], -depths, type='b', pch=17, col='deepskyblue', cex=cex.setting)
# lines(soilC_1993_means['LMT',], -depths, type='b', col='chartreuse4', pch=0, lty=2, cex=cex.setting)
# lines(soilC_2012_means['LMT',], -depths, type='b', pch=15, col='chartreuse4', cex=cex.setting)
# lines(soilC_1993_means['CWT',], -depths, type='b', col='gold', pch=6, lty=2, cex=cex.setting)
# lines(soilC_2012_means['CWT',], -depths, type='b', pch=25, col='gold', bg='gold', cex=cex.setting)
# text(x=0.35, y=-15, labels = 'a')
# legend(x='bottomright', legend=c("OMT ('93)", "OMT ('12)", "LMT ('93)", "LMT ('12)", "CMT ('93)", "CMT ('12)", "CWT ('93)", "CWT ('12)"), lty=rep(c(2,1), 4), col=c('tan4', 'tan4', 'chartreuse4', 'chartreuse4', 'deepskyblue', 'deepskyblue', 'gold', 'gold'), pch=c(1,16,0,15,2,17,6,25), pt.bg = c(rep(NA, 7), 'gold'), pt.cex = cex.setting)
# dev.off()

#version 2 of maize-tomato systems less wheat-tomato
tiff(file = file.path(FiguresDir, 'soilC_change_ltras_fig_a.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(3.5,3.5,0.1,0.1))
plot(soilC_1993_means['OMT',], -depths, type='b', xlim=c(0.33,1.4), col='tan4', yaxt='n', ylab='', xlab='', lty=2, pch=1, cex=cex.setting)
axis(side = 2, at=-c(0,15,30,60,100), labels=c(0,15,30,60,100))
mtext('Depth (cm)', 2, 2.25)
mtext('Soil organic carbon (%)', 1, 2.25)
lines(soilC_2012_means['OMT',], -depths, type='b', xlim=c(0.1,1.4), col='tan4', pch=16, cex=cex.setting)
lines(soilC_1993_means['CMT',], -depths, type='b', col='deepskyblue', pch=2, lty=2, cex=cex.setting)
lines(soilC_2012_means['CMT',], -depths, type='b', pch=17, col='deepskyblue', cex=cex.setting)
lines(soilC_1993_means['LMT',], -depths, type='b', col='chartreuse4', pch=0, lty=2, cex=cex.setting)
lines(soilC_2012_means['LMT',], -depths, type='b', pch=15, col='chartreuse4', cex=cex.setting)
points(y=-7.5,x=soilC_2012_means[row.names(soilC_2012_means)=='OMT', 1], cex=1.1)
points(y=-7.5,x=soilC_2012_means[row.names(soilC_2012_means)=='LMT', 1], cex=1.1)
points(y=-22.5,x=soilC_2012_means[row.names(soilC_2012_means)=='OMT', 2], cex=1.1)
text(x=0.4, y=-15, labels = 'a')
points(x=0.82,y=-80, pch=1, cex=1.1)
text(x=0.93,y=-75, labels="significant change\nafter 19 yrs.", adj=c(0,1), cex=0.85)
legend(x='bottomright', legend=c("ORG: yr 0", "ORG: yr 19", "CONV: yr 0", "CONV: yr 19", "CONV+WCC: yr 0", "CONV+WCC: yr 19"), lty=rep(c(2,1), 3), col=c('tan4', 'tan4', 'deepskyblue', 'deepskyblue', 'chartreuse4', 'chartreuse4'), pch=c(1,16,2,17,0,15), pt.cex = cex.setting, bty = 'n', cex=0.85)
dev.off()
#csv file convention c("OMT: yr 0", "OMT: yr 19", "LMT: yr 0", "LMT: yr 19", "CMT: yr 0", "CMT: yr 19")

#all wheat system comparisons
# tiff(file = file.path(FiguresDir, 'soilC_change_ltras_wheat.tif'), family = 'Times New Roman', width = 4.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
# par(mar=c(3.5,3.5,0.5,0.5))
# plot(soilC_1993_means['RWF',], -depths, type='b', xlim=c(0.35,1.1), col='springgreen4', yaxt='n', ylab='', xlab='', lty=2, pch=1, cex=cex.setting)
# axis(side = 2, at=-c(0,15,30,60,100), labels=c(0,15,30,60,100))
# mtext('Depth (cm)', 2, 2.25)
# mtext('Soil organic carbon (%)', 1, 2.25)
# lines(soilC_2012_means['RWF',], -depths, type='b', xlim=c(0.1,1.4), col='springgreen4', pch=16, cex=cex.setting)
# lines(soilC_1993_means['IWF',], -depths, type='b', col='royalblue3', pch=2, lty=2, cex=cex.setting)
# lines(soilC_2012_means['IWF',], -depths, type='b', pch=17, col='royalblue3', cex=cex.setting)
# lines(soilC_1993_means['IWC',], -depths, type='b', col='purple3', pch=0, lty=2, cex=cex.setting)
# lines(soilC_2012_means['IWC',], -depths, type='b', pch=15, col='purple3', cex=cex.setting)
# lines(soilC_1993_means['RWC',], -depths, type='b', col='red3', pch=6, lty=2, cex=cex.setting)
# lines(soilC_2012_means['RWC',], -depths, type='b', pch=25, col='red3', bg='red3', cex=cex.setting)
# lines(soilC_1993_means['RWL',], -depths, type='b', col='darkorange1', pch=5, lty=2, cex=cex.setting)
# lines(soilC_2012_means['RWL',], -depths, type='b', pch=23, col='darkorange1', bg='darkorange1', cex=cex.setting)
# text(x=0.4, y=-15, labels = 'b')
# legend(x='bottomright', legend=c("RWF ('93)", "RWF ('12)", "IWF ('93)", "IWF ('12)", "IWC ('93)", "IWC ('12)", "RWC ('93)", "RWC ('12)", "RWL ('93)", "RWL ('12)"), lty=rep(c(2,1), 5), col=c('springgreen4', 'springgreen4', 'royalblue3', 'royalblue3', 'purple3', 'purple3', 'red3', 'red3', 'darkorange1', 'darkorange1'), pch=c(1,16,0,15,2,17,6,25,5,23), pt.bg = c(rep(NA, 7), 'red3', NA, 'darkorange1'), pt.cex = cex.setting)
# dev.off()

#rainfed wheat system comparisons, fig b
tiff(file = file.path(FiguresDir, 'soilC_change_ltras_wheat_rainfed_fig_b.tif'), family = 'Times New Roman', width = 2.75, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(3.5,0.1,0.1,0.1))
plot(soilC_1993_means['RWF',], -depths, type='b', xlim=c(0.35,1.05), col='springgreen4', yaxt='n', ylab='', xlab='', lty=2, pch=1, cex=cex.setting)
#axis(side = 2, at=-c(0,15,30,60,100), labels=c(0,15,30,60,100))
#mtext('Depth (cm)', 2, 2.25)
mtext('Soil organic carbon (%)', 1, 2.25)
lines(soilC_2012_means['RWF',], -depths, type='b', xlim=c(0.1,1.4), col='springgreen4', pch=16, cex=cex.setting)
lines(soilC_1993_means['RWC',], -depths, type='b', col='red3', pch=6, lty=2, cex=cex.setting)
lines(soilC_2012_means['RWC',], -depths, type='b', pch=25, col='red3', bg='red3', cex=cex.setting)
lines(soilC_1993_means['RWL',], -depths, type='b', col='darkorange1', pch=5, lty=2, cex=cex.setting)
lines(soilC_2012_means['RWL',], -depths, type='b', pch=23, col='darkorange1', bg='darkorange1', cex=cex.setting)
points(y=-125,x=soilC_2012_means[row.names(soilC_2012_means)=='RWL', 5], cex=1.1)
points(y=-125,x=soilC_2012_means[row.names(soilC_2012_means)=='RWF', 5], cex=1.1)
text(x=0.4, y=-15, labels = 'b')
legend(x='bottomright', legend=c("RWF: yr 0", "RWF: yr 19", "RWF+N: yr 0", "RWF+N: yr 19",  "RWF+WCC: yr 0", "RWF+WCC: yr 19"), lty=rep(c(2,1), 5), col=c('red3', 'red3', 'springgreen4', 'springgreen4', 'darkorange1', 'darkorange1'), pch=c(6,25,1,16,5,23), pt.bg = c(NA, 'red3', rep(NA, 3), 'darkorange1'), pt.cex = cex.setting, bty = 'n', cex=0.85)
dev.off()
#csv convention for treatment abbreviations: c("RWF: yr 0", "RWF: yr 19", "RWC: yr 0", "RWC: yr 19", "RWL: yr 0", "RWL: yr 19")

tiff(file = file.path(FiguresDir, 'soilC_change_ltras_wheat_irrigated_fig_c.tif'), family = 'Times New Roman', width = 2.75, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(3.5,0.1,0.1,0.1))
plot(soilC_1993_means['IWF',], -depths, type='b', xlim=c(0.35,1.05), col='royalblue4', yaxt='n', ylab='', xlab='', lty=2, pch=2, cex=cex.setting)
#axis(side = 2, at=-c(0,15,30,60,100), labels=c(0,15,30,60,100))
#mtext('Depth (cm)', 2, 2.25)
mtext('Soil organic carbon (%)', 1, 2.25)
lines(soilC_2012_means['IWF',], -depths, type='b', pch=17, col='royalblue3', cex=cex.setting)
lines(soilC_1993_means['IWC',], -depths, type='b', col='purple3', pch=0, lty=2, cex=cex.setting)
lines(soilC_2012_means['IWC',], -depths, type='b', pch=15, col='purple3', cex=cex.setting)
points(y=-7.5,x=soilC_2012_means[row.names(soilC_2012_means)=='IWF', 1], cex=1.1)
points(y=-125,x=soilC_2012_means[row.names(soilC_2012_means)=='IWF', 5], cex=1.1)
text(x=0.4, y=-15, labels = 'c')
legend(x='bottomright', legend=c("IWF: yr 0", "IWF: yr 19", "IWF+N: yr 0", "IWF+N: yr 19"), lty=rep(c(2,1), 3), col=c('purple3', 'purple3', 'royalblue3', 'royalblue3'), pch=c(0,15,2,17), pt.cex = cex.setting, bty = 'n', cex=0.85)
dev.off()
#lines(soilC_1993_means['CWT',], -depths, type='b', col='gold', pch=1, lty=2, cex=cex.setting)
#lines(soilC_2012_means['CWT',], -depths, type='b', pch=21, col='gold', bg='gold', cex=cex.setting)
#, "IWT+N: yr 0", "IWT+N: yr 19"
#, 'gold', 'gold'
#,1,21
#, pt.bg = c(rep(NA,5), 'gold')


#t-tests to try to reproduce tautges et al. 2019
deltaC_t.test <- function(df_93, df_12, depth_class) {
  trtments <- unique(df_93$Treatment)
  results <- data.frame(Treatments=trtments, meanC_93=NA, sdevC_93=NA, meanC_12=NA, sdevC_12=NA, DeltaC=NA, low.CI.95=NA, high.CI.95=NA, p.value=NA, stringsAsFactors = FALSE)
  for (i in seq_along(trtments)) {
    result <- t.test(df_12$soilC[df_12$Treatment==trtments[i]], df_93$soilC[df_93$Treatment==trtments[i]], paired = TRUE)
    results$meanC_93[i] <- mean(df_93$soilC[df_93$Treatment==trtments[i]])
    results$sdevC_93[i] <- sd(df_93$soilC[df_93$Treatment==trtments[i]])
    results$meanC_12[i] <- mean(df_12$soilC[df_12$Treatment==trtments[i]])
    results$sdevC_12[i] <- sd(df_12$soilC[df_12$Treatment==trtments[i]])
    results$DeltaC[i] <- result$estimate
    results$low.CI.95[i] <- result$conf.int[1]
    results$high.CI.95[i] <- result$conf.int[2]
    results$p.value[i] <- result$p.value
  }
  write.csv(results, file.path(dataDir, 'CenturyResults/t.tests', paste0('t.test_', depth_class, '.csv')), row.names=FALSE)
}
deltaC_t.test(df_93 = soilC_0_15_1993, df_12 = soilC_0_15_2012, depth_class = '0_15')
deltaC_t.test(df_93 = soilC_15_30_1993, df_12 = soilC_15_30_2012, depth_class = '15_30')
deltaC_t.test(df_93 = soilC_30_60_1993, df_12 = soilC_30_60_2012, depth_class = '30_60')
deltaC_t.test(df_93 = soilC_60_100_1993, df_12 = soilC_60_100_2012, depth_class = '60_100')
deltaC_t.test(df_93 = soilC_100_150_1993, df_12 = soilC_100_150_2012, depth_class = '100_150')


#rainfed wheat
t.test(soilC_100_150_1993$soilC[soilC_100_150_1993$Treatment=='RWL'], soilC_100_150_2012$soilC[soilC_100_150_2012$Treatment=='RWL'], paired = TRUE) #0.002 (0.31 if not paired)
t.test(soilC_100_150_1993$soilC[soilC_100_150_1993$Treatment=='RWF'], soilC_100_150_2012$soilC[soilC_100_150_2012$Treatment=='RWF'], paired = TRUE) #0.05 (0.36 if not paired)
t.test(soilC_100_150_1993$soilC[soilC_100_150_1993$Treatment=='RWC'], soilC_100_150_2012$soilC[soilC_100_150_2012$Treatment=='RWC'], paired = TRUE) #0.55 (0.62 if not paired)

#create deltaC data.frame for inspection
dim(soilC_1993)
dim(soilC_2012)
all(soilC_1993$Plot==soilC_2012$Plot)
all(soilC_1993$Treatment==soilC_2012$Treatment)
deltaC_19yrs <- soilC_2012[,3:7] - soilC_1993[,3:7]
deltaC_19yrs$Plot <- row.names(deltaC_19yrs)
deltaC_19yrs$Treatment <- soilC_1993$Treatment
View(deltaC_19yrs)
lapply(deltaC_19yrs[,1:5], summary)

#delta C by depth class
all(soilC_0_15_1993$Plot==soilC_0_15_2012$Plot)
all(soilC_0_15_1993$Treatment==soilC_0_15_2012$Treatment)
deltasoilC_0_15 <- data.frame(Plot=soilC_0_15_1993$Plot, Treatment=soilC_0_15_1993$Treatment, deltaC=soilC_0_15_2012$soilC - soilC_0_15_1993$soilC, stringsAsFactors = FALSE)
tapply(deltasoilC_0_15$deltaC*10, deltasoilC_0_15$Treatment, mean) #RWC has five reps
summary(aov(deltaC ~ Treatment, data = deltasoilC_0_15))
TukeyHSD(aov(deltaC ~ Treatment, data = deltasoilC_0_15))

#RWC (aka RWF in Tautges et al.) shows large, non-sig decrease in upper 15 cm
cbind(soilC_0_15_1993[soilC_0_15_1993$Treatment=='RWC',c(1,2)], soilC_0_15_2012[soilC_0_15_1993$Treatment=='RWC',c(1,2)]) #plot 1-8 is suspect

#0-30 cm content plot

