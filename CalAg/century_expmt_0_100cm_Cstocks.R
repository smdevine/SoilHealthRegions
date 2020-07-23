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
  ltras_cn_1993 <-  ltras_cn_1993[-which(ltras_cn_1993$Plot=='1_8'& ltras_cn_1993$Depth=="000-015"), ]
  ltras_cn_2012 <-  ltras_cn_2012[-which(ltras_cn_2012$Plot=='1_8'& ltras_cn_2012$Depth=="000-015"), ]
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

#0-15 C stocks
soilC_0_15_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='000-015'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='000-015'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='000-015'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='000-015'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_1993)
soilC_0_15_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_0_15_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_0_15_1993$Treatment)


#15-30 soil C
soilC_15_30_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='015-030'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='015-030'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='015-030'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='015-030'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_1993)
soilC_15_30_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_15_30_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_15_30_1993$Treatment)

#30-60 soil C
soilC_30_60_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='030-060'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='030-060'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='030-060'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='030-060'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_1993)
soilC_30_60_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_30_60_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_30_60_1993$Treatment)

#60-100 soil C
soilC_60_100_1993 <- data.frame(Plot=names(tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='060-100'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='060-100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$C[ltras_cn_1993$Depth=='060-100'], ltras_cn_1993$Plot[ltras_cn_1993$Depth=='060-100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_1993)
soilC_60_100_1993$Treatment <- ltras_cn_1993$Treatment[match(soilC_60_100_1993$Plot, ltras_cn_1993$Plot)]
table(soilC_60_100_1993$Treatment)

#same for 2012 data
#0-15 soil %C
soilC_0_15_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='000-015'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='000-015'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='000-015'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='000-015'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_2012)
soilC_0_15_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_0_15_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_0_15_2012$Treatment)

#15-30 soil C
soilC_15_30_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='015-030'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='015-030'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='015-030'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='015-030'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_2012)
soilC_15_30_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_15_30_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_15_30_2012$Treatment)

#30-60 soil C
soilC_30_60_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='030-060'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='030-060'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='030-060'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='030-060'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_2012)
soilC_30_60_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_30_60_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_30_60_2012$Treatment)

#60-100 soil C
soilC_60_100_2012 <- data.frame(Plot=names(tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='060-100'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='060-100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$C[ltras_cn_2012$Depth=='060-100'], ltras_cn_2012$Plot[ltras_cn_2012$Depth=='060-100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_2012)
soilC_60_100_2012$Treatment <- ltras_cn_2012$Treatment[match(soilC_60_100_2012$Plot, ltras_cn_2012$Plot)]
table(soilC_60_100_2012$Treatment)


#consolidate 1993 data
soilC_1993 <- soilC_0_15_1993
soilC_1993 <- soilC_1993[,c(1,3,2)]
colnames(soilC_1993)[3] <- 'C_0_15'
soilC_1993$C_15_30 <- soilC_15_30_1993$soilC[match(soilC_1993$Plot, soilC_15_30_1993$Plot)]
soilC_1993$C_30_60 <- soilC_30_60_1993$soilC[match(soilC_1993$Plot, soilC_30_60_1993$Plot)]
soilC_1993$C_60_100 <- soilC_60_100_1993$soilC[match(soilC_1993$Plot, soilC_60_100_1993$Plot)]
# soilC_1993$C_100_150 <- soilC_100_150_1993$soilC[match(soilC_1993$Plot, soilC_100_150_1993$Plot)]
lapply(soilC_1993[,3:6], summary)
lapply(soilC_1993[,3:6], hist)
lapply(soilC_1993[,3:6], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], mean))
soilC_1993_means <- do.call(cbind, lapply(soilC_1993[,3:6], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], mean)))
soilC_1993_sd <- do.call(cbind, lapply(soilC_1993[,3:6], function(x) tapply(x[!is.na(x)], soilC_1993$Treatment[!is.na(x)], sd)))
View(soilC_1993)

#consolidate 2012 data
soilC_2012 <- soilC_0_15_2012
soilC_2012 <- soilC_2012[,c(1,3,2)]
colnames(soilC_2012)[3] <- 'C_0_15'
soilC_2012$C_15_30 <- soilC_15_30_2012$soilC[match(soilC_2012$Plot, soilC_15_30_2012$Plot)]
soilC_2012$C_30_60 <- soilC_30_60_2012$soilC[match(soilC_2012$Plot, soilC_30_60_2012$Plot)]
soilC_2012$C_60_100 <- soilC_60_100_2012$soilC[match(soilC_2012$Plot, soilC_60_100_2012$Plot)]
# soilC_2012$C_100_150 <- soilC_100_150_2012$soilC[match(soilC_2012$Plot, soilC_100_150_2012$Plot)]
lapply(soilC_2012[,3:6], summary)
lapply(soilC_2012[,3:6], hist)
lapply(soilC_2012[,3:6], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], mean))
soilC_2012_means <- do.call(cbind, lapply(soilC_2012[,3:6], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], mean)))
soilC_2012_sd <- do.call(cbind, lapply(soilC_2012[,3:6], function(x) tapply(x[!is.na(x)], soilC_2012$Treatment[!is.na(x)], sd)))
soilC_2012_means
soilC_2012_sd
View(soilC_2012)
depths <- c(7.5,22.5,45,80)

#0-30 cm stocks
dim(soilC_0_15_1993)
dim(soilC_15_30_1993)
dim(soilC_15_30_1993[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot,])
dim(soilC_0_15_2012)
dim(soilC_15_30_2012)
all(soilC_0_15_1993$Plot==soilC_15_30_1993$Plot[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot])
soilC_0_30_1993 <- soilC_15_30_1993[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot,]
soilC_0_30_1993$soilC <- soilC_0_15_1993$soilC + soilC_15_30_1993$soilC[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot]
tapply(soilC_0_30_1993$soilC, soilC_0_30_1993$Treatment, summary)
tapply(soilC_0_30_1993$soilC, soilC_0_30_1993$Treatment, function(x) length(x))

all(soilC_0_15_2012$Plot==soilC_15_30_2012$Plot[soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot])
soilC_0_30_2012 <- soilC_15_30_2012[soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot,]
soilC_0_30_2012$soilC <- soilC_0_15_2012$soilC + soilC_15_30_2012$soilC[soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot]
all(soilC_0_30_2012$Plot==soilC_0_30_2012$Plot)
tapply(soilC_0_30_2012$soilC, soilC_0_30_2012$Treatment, summary)
tapply(soilC_0_30_2012$soilC, soilC_0_30_2012$Treatment, function(x) length(x))

soilC_0_30_deltaC <- soilC_0_30_2012
soilC_0_30_deltaC$soilC <- soilC_0_30_2012$soilC - soilC_0_30_1993$soilC
colnames(soilC_0_30_deltaC)[2] <- 'deltaC_19yrs'
View(soilC_0_30_deltaC)

tapply(soilC_0_30_deltaC$deltaC_19yrs, soilC_0_30_deltaC$Treatment, summary)
tapply(soilC_0_30_deltaC$deltaC_19yrs, soilC_0_30_deltaC$Treatment, function(x) print(x[order(x)]))
tapply(soilC_0_30_deltaC$deltaC_19yrs, soilC_0_30_deltaC$Treatment, function(x) length(x))

soilC_0_30_1993means <- tapply(soilC_0_30_1993$soilC, soilC_0_30_1993$Treatment, mean)
soilC_0_30_2012means <- tapply(soilC_0_30_2012$soilC, soilC_0_30_2012$Treatment, mean)
soilC_0_30_1993SEs <- tapply(soilC_0_30_1993$soilC, soilC_0_30_1993$Treatment, function(x) sd(x)/sqrt(length(x)))
soilC_0_30_2012SEs <- tapply(soilC_0_30_2012$soilC, soilC_0_30_2012$Treatment, function(x) sd(x)/sqrt(length(x)))
soilC_0_30_1993means
soilC_0_30_2012means

#maize vs. tomato comparisons
soilC_0_30_MT_means <- cbind(year_0=soilC_0_30_1993means, year_19=soilC_0_30_2012means)
row.names(soilC_0_30_MT_means) <- names(soilC_0_30_1993means)
soilC_0_30_MT_means <- soilC_0_30_MT_means[c(1,5:6),]

soilC_0_30_MT_SEs <- cbind(year_0=soilC_0_30_1993SEs, year_19=soilC_0_30_2012SEs)
row.names(soilC_0_30_MT_SEs) <- names(soilC_0_30_1993SEs)
soilC_0_30_MT_SEs <- soilC_0_30_MT_SEs[c(1,5:6),]
tiff(file = file.path(FiguresDir, 'soilC_0_30change_Century.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,0.5,0.25))
century_bar_dims <- barplot(soilC_0_30_MT_means, beside=TRUE, legend.text=c('CONV', 'CONV-WCC', 'ORG'), ylab='', ylim=c(0,55), col=c('deepskyblue', 'chartreuse4', 'tan4'), names.arg = c('Year 0', 'Year 19'), args.legend=list(x='topleft', cex=0.9, bty='n', ncol=2, inset=c(0,0.05)), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
mtext(expression('0-30 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=century_bar_dims[1,], y0=soilC_0_30_MT_means[1,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_30_MT_means[1,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[1,], lwd = 1)
segments(x0=century_bar_dims[2,], y0=soilC_0_30_MT_means[2,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_30_MT_means[2,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[2,], lwd = 1)
segments(x0=century_bar_dims[3,], y0=soilC_0_30_MT_means[3,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_30_MT_means[3,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[3,], lwd = 1)
# segments(x0=century_bar_dims[4,], y0=soilC_0_30_MT_means[4,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_30_MT_means[4,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[4,], lwd = 1)
arrows(x0=century_bar_dims[1,], y0=soilC_0_30_MT_means[1,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_30_MT_means[1,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[1,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[2,], y0=soilC_0_30_MT_means[2,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_30_MT_means[2,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[3,], y0=soilC_0_30_MT_means[3,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_30_MT_means[3,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[3,], lwd = 1, angle = 90, code = 3, length = 0.05)
# arrows(x0=century_bar_dims[4,], y0=soilC_0_30_MT_means[4,] + qt(0.975,df=5) * soilC_0_30_MT_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_30_MT_means[4,] + qt(0.025,df=5) * soilC_0_30_MT_SEs[4,], lwd = 1, angle = 90, code = 3, length = 0.05)
dev.off()

#wheat comparisons
soilC_0_30_wheat_means <- cbind(year_0=soilC_0_30_1993means, year_19=soilC_0_30_2012means)
row.names(soilC_0_30_wheat_means) <- names(soilC_0_30_1993means)
soilC_0_30_wheat_means <- soilC_0_30_wheat_means[c(3:4,7:9),]

soilC_0_30_wheat_SEs <- cbind(year_0=soilC_0_30_1993SEs, year_19=soilC_0_30_2012SEs)
row.names(soilC_0_30_wheat_SEs) <- names(soilC_0_30_1993SEs)
soilC_0_30_wheat_SEs <- soilC_0_30_wheat_SEs[c(3:4,7:9),]

tiff(file = file.path(FiguresDir, 'soilC_0_30change_wheat_Century.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,0.5,0.25))
century_bar_dims <- barplot(soilC_0_30_wheat_means, beside=TRUE, legend.text=c('IWF', 'IWF+N', 'RWF', 'RWF+N', 'RWF+WCC'), ylab='', ylim=c(0,48), col=c('purple3', 'royalblue3', 'red3', 'springgreen4', 'darkorange1'), names.arg = c('Year 0', 'Year 19'), args.legend=list(x='topleft', cex=0.9, bty='n', ncol=3, inset=c(-0.25,-0.04)), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
mtext(expression('0-30 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=century_bar_dims[1,], y0=soilC_0_30_wheat_means[1,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_30_wheat_means[1,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[1,], lwd = 1)
segments(x0=century_bar_dims[2,], y0=soilC_0_30_wheat_means[2,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_30_wheat_means[2,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[2,], lwd = 1)
segments(x0=century_bar_dims[3,], y0=soilC_0_30_wheat_means[3,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_30_wheat_means[3,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[3,], lwd = 1)
segments(x0=century_bar_dims[4,], y0=soilC_0_30_wheat_means[4,] + qt(0.975,df=3) * soilC_0_30_wheat_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_30_wheat_means[4,] + qt(0.025,df=3) * soilC_0_30_wheat_SEs[4,], lwd = 1) #RWC treatment only had 4 available reps; hence df=3
segments(x0=century_bar_dims[5,], y0=soilC_0_30_wheat_means[5,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[5,], x1=century_bar_dims[5,], y1=soilC_0_30_wheat_means[5,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[5,], lwd = 1)
arrows(x0=century_bar_dims[1,], y0=soilC_0_30_wheat_means[1,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_30_wheat_means[1,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[1,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[2,], y0=soilC_0_30_wheat_means[2,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_30_wheat_means[2,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[3,], y0=soilC_0_30_wheat_means[3,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_30_wheat_means[3,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[3,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[4,], y0=soilC_0_30_wheat_means[4,] + qt(0.975,df=3) * soilC_0_30_wheat_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_30_wheat_means[4,] + qt(0.025,df=3) * soilC_0_30_wheat_SEs[4,], lwd = 1, angle = 90, code = 3, length = 0.05) #RWC treatment only had 4 available reps; hence df=3
arrows(x0=century_bar_dims[5,], y0=soilC_0_30_wheat_means[5,] + qt(0.975,df=5) * soilC_0_30_wheat_SEs[5,], x1=century_bar_dims[5,], y1=soilC_0_30_wheat_means[5,] + qt(0.025,df=5) * soilC_0_30_wheat_SEs[5,], lwd = 1, angle = 90, code = 3, length = 0.05)
dev.off()

#0-100 cm stocks
dim(soilC_0_15_1993)
dim(soilC_15_30_1993)
dim(soilC_30_60_1993)
dim(soilC_60_100_1993)

dim(soilC_0_15_2012)
dim(soilC_15_30_2012)
dim(soilC_30_60_2012)
dim(soilC_60_100_2012)

all(soilC_0_15_1993$Plot==soilC_15_30_1993$Plot[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot])
all(soilC_0_15_1993$Plot==soilC_30_60_1993$Plot)
all(soilC_0_15_1993$Plot==soilC_60_100_1993$Plot)
all(soilC_30_60_1993$Plot==soilC_60_100_1993$Plot)
sum(soilC_0_15_1993$Plot %in% soilC_30_60_1993$Plot)
sum(soilC_15_30_1993$Plot %in% soilC_30_60_1993$Plot)
sum(soilC_15_30_1993$Plot %in% soilC_30_60_1993$Plot & soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot)

soilC_0_100_1993 <- soilC_15_30_1993[soilC_15_30_1993$Plot %in% soilC_0_15_1993$Plot & soilC_15_30_1993$Plot %in% soilC_30_60_1993$Plot,]
dim(soilC_0_100_1993)
soilC_0_100_1993$soilC <- soilC_0_15_1993$soilC[match(soilC_0_100_1993$Plot, soilC_0_15_1993$Plot)] + soilC_15_30_1993$soilC[match(soilC_0_100_1993$Plot, soilC_15_30_1993$Plot)] + soilC_30_60_1993$soilC[match(soilC_0_100_1993$Plot, soilC_30_60_1993$Plot)] + soilC_60_100_1993$soilC[match(soilC_0_100_1993$Plot, soilC_60_100_1993$Plot)]
tapply(soilC_0_100_1993$soilC, soilC_0_100_1993$Treatment, summary)
tapply(soilC_0_100_1993$soilC, soilC_0_100_1993$Treatment, function(x) length(x))

all(soilC_0_15_2012$Plot==soilC_15_30_2012$Plot[soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot])
all(soilC_0_15_2012$Plot==soilC_30_60_2012$Plot)
all(soilC_0_15_2012$Plot==soilC_60_100_2012$Plot)
all(soilC_30_60_2012$Plot==soilC_60_100_2012$Plot)
sum(soilC_0_15_2012$Plot %in% soilC_30_60_2012$Plot)
sum(soilC_15_30_2012$Plot %in% soilC_30_60_2012$Plot)
sum(soilC_15_30_2012$Plot %in% soilC_30_60_2012$Plot & soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot)

soilC_0_100_2012 <- soilC_15_30_2012[soilC_15_30_2012$Plot %in% soilC_0_15_2012$Plot & soilC_15_30_2012$Plot %in% soilC_30_60_2012$Plot,]
dim(soilC_0_100_2012)
all(soilC_0_100_1993$Plot==soilC_0_100_2012$Plot)
soilC_0_100_2012$soilC <- soilC_0_15_2012$soilC[match(soilC_0_100_2012$Plot, soilC_0_15_2012$Plot)] + soilC_15_30_2012$soilC[match(soilC_0_100_2012$Plot, soilC_15_30_2012$Plot)] + soilC_30_60_2012$soilC[match(soilC_0_100_2012$Plot, soilC_30_60_2012$Plot)] + soilC_60_100_2012$soilC[match(soilC_0_100_2012$Plot, soilC_60_100_2012$Plot)]
tapply(soilC_0_100_2012$soilC, soilC_0_100_2012$Treatment, summary)
tapply(soilC_0_100_2012$soilC, soilC_0_100_2012$Treatment, function(x) length(x))

soilC_0_100_deltaC <- soilC_0_100_2012
soilC_0_100_deltaC$soilC <- soilC_0_100_2012$soilC - soilC_0_100_1993$soilC
colnames(soilC_0_100_deltaC)[2] <- 'deltaC_19yrs'
View(soilC_0_100_deltaC)

tapply(soilC_0_100_deltaC$deltaC_19yrs, soilC_0_100_deltaC$Treatment, summary)
tapply(soilC_0_100_deltaC$deltaC_19yrs, soilC_0_100_deltaC$Treatment, function(x) print(x[order(x)]))
tapply(soilC_0_100_deltaC$deltaC_19yrs, soilC_0_100_deltaC$Treatment, function(x) length(x))
soilC_0_100_deltaC[soilC_0_100_deltaC$Treatment=='LMT',]
View(soilC_60_100_1993)
View(soilC_60_100_2012) #plot 2-4 had unrealistic decline at 60-100 cm
ltras_cn_1993[ltras_cn_1993$Plot=='2_4' & ltras_cn_1993$Upper_Depth==60,]
ltras_cn_2012[ltras_cn_2012$Plot=='2_4' & ltras_cn_2012$Upper_Depth==60,]
ltras_cn_2012[ltras_cn_2012$Treatment=='LMT' & ltras_cn_2012$Upper_Depth==60,]
ltras_cn_2012$PctC[ltras_cn_2012$Upper_Depth==60][order(ltras_cn_2012$PctC[ltras_cn_2012$Upper_Depth==60])]
#get rid of plot 2-4
soilC_0_100_1993 <- soilC_0_100_1993[!soilC_0_100_1993$Plot=='2_4',]
soilC_0_100_2012 <- soilC_0_100_2012[!soilC_0_100_2012$Plot=='2_4',]
all(soilC_0_100_1993$Plot==soilC_0_100_2012$Plot)

soilC_0_100_1993means <- tapply(soilC_0_100_1993$soilC, soilC_0_100_1993$Treatment, mean)
soilC_0_100_2012means <- tapply(soilC_0_100_2012$soilC, soilC_0_100_2012$Treatment, mean)
soilC_0_100_1993SEs <- tapply(soilC_0_100_1993$soilC, soilC_0_100_1993$Treatment, function(x) sd(x)/sqrt(length(x)))
soilC_0_100_2012SEs <- tapply(soilC_0_100_2012$soilC, soilC_0_100_2012$Treatment, function(x) sd(x)/sqrt(length(x)))
soilC_0_100_1993means
soilC_0_100_2012means
soilC_0_100_1993SEs
soilC_0_100_2012SEs

#maize vs. tomato comparisons
soilC_0_100_MT_means <- cbind(year_0=soilC_0_100_1993means, year_19=soilC_0_100_2012means)
row.names(soilC_0_100_MT_means) <- names(soilC_0_100_1993means)
soilC_0_100_MT_means <- soilC_0_100_MT_means[c(1,5:6),]
soilC_0_100_MT_means

soilC_0_100_MT_SEs <- cbind(year_0=soilC_0_100_1993SEs, year_19=soilC_0_100_2012SEs)
row.names(soilC_0_100_MT_SEs) <- names(soilC_0_100_1993SEs)
soilC_0_100_MT_SEs <- soilC_0_100_MT_SEs[c(1,5:6),]

tiff(file = file.path(FiguresDir, 'soilC_0_100change_Century.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,0.5,0.25))
century_bar_dims <- barplot(soilC_0_100_MT_means, beside=TRUE, legend.text=c('CONV', 'CONV-WCC', 'ORG'), ylab='', ylim=c(0,155), col=c('deepskyblue', 'chartreuse4', 'tan4'), names.arg = c('Year 0', 'Year 19'), args.legend=list(x='topleft', cex=0.9, bty='n', ncol=3, inset=c(-0.15,-0.04)), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
mtext(expression('0-100 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=century_bar_dims[1,], y0=soilC_0_100_MT_means[1,] + qt(0.975,df=5) * soilC_0_100_MT_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_100_MT_means[1,] + qt(0.025,df=5) * soilC_0_100_MT_SEs[1,], lwd = 1)
segments(x0=century_bar_dims[2,], y0=soilC_0_100_MT_means[2,] + qt(0.975,df=4) * soilC_0_100_MT_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_100_MT_means[2,] + qt(0.025,df=4) * soilC_0_100_MT_SEs[2,], lwd = 1)
segments(x0=century_bar_dims[3,], y0=soilC_0_100_MT_means[3,] + qt(0.975,df=4) * soilC_0_100_MT_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_100_MT_means[3,] + qt(0.025,df=4) * soilC_0_100_MT_SEs[3,], lwd = 1)
# segments(x0=century_bar_dims[4,], y0=soilC_0_100_MT_means[4,] + qt(0.975,df=5) * soilC_0_100_MT_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_100_MT_means[4,] + qt(0.025,df=5) * soilC_0_100_MT_SEs[4,], lwd = 1)
arrows(x0=century_bar_dims[1,], y0=soilC_0_100_MT_means[1,] + qt(0.975,df=5) * soilC_0_100_MT_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_100_MT_means[1,] + qt(0.025,df=5) * soilC_0_100_MT_SEs[1,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[2,], y0=soilC_0_100_MT_means[2,] + qt(0.975,df=4) * soilC_0_100_MT_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_100_MT_means[2,] + qt(0.025,df=4) * soilC_0_100_MT_SEs[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[3,], y0=soilC_0_100_MT_means[3,] + qt(0.975,df=4) * soilC_0_100_MT_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_100_MT_means[3,] + qt(0.025,df=4) * soilC_0_100_MT_SEs[3,], lwd = 1, angle = 90, code = 3, length = 0.05)
# arrows(x0=century_bar_dims[4,], y0=soilC_0_100_MT_means[4,] + qt(0.975,df=5) * soilC_0_100_MT_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_100_MT_means[4,] + qt(0.025,df=5) * soilC_0_100_MT_SEs[4,], lwd = 1, angle = 90, code = 3, length = 0.05)
dev.off()

#wheat comparisons
soilC_0_100_wheat_means <- cbind(year_0=soilC_0_100_1993means, year_19=soilC_0_100_2012means)
row.names(soilC_0_100_wheat_means) <- names(soilC_0_100_1993means)
soilC_0_100_wheat_means <- soilC_0_100_wheat_means[c(3:4,7:9),]

soilC_0_100_wheat_SEs <- cbind(year_0=soilC_0_100_1993SEs, year_19=soilC_0_100_2012SEs)
row.names(soilC_0_100_wheat_SEs) <- names(soilC_0_100_1993SEs)
soilC_0_100_wheat_SEs <- soilC_0_100_wheat_SEs[c(3:4,7:9),]

tiff(file = file.path(FiguresDir, 'soilC_0_100change_wheat_Century.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,3,0.25))
century_bar_dims <- barplot(soilC_0_100_wheat_means, beside=TRUE, legend.text=c('IWF', 'IWF+N', 'RWF', 'RWF+N', 'RWF+WCC'), ylab='', ylim=c(0,140), col=c('purple3', 'royalblue3', 'red3', 'springgreen4', 'darkorange1'), names.arg = c('Year 0', 'Year 19'), args.legend=list(x='topleft', cex=0.9, bty='n', ncol=3, inset=c(-0.25,-0.2)), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
mtext(expression('0-100 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=century_bar_dims[1,], y0=soilC_0_100_wheat_means[1,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_100_wheat_means[1,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[1,], lwd = 1)
segments(x0=century_bar_dims[2,], y0=soilC_0_100_wheat_means[2,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_100_wheat_means[2,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[2,], lwd = 1)
segments(x0=century_bar_dims[3,], y0=soilC_0_100_wheat_means[3,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_100_wheat_means[3,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[3,], lwd = 1)
segments(x0=century_bar_dims[4,], y0=soilC_0_100_wheat_means[4,] + qt(0.975,df=3) * soilC_0_100_wheat_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_100_wheat_means[4,] + qt(0.025,df=3) * soilC_0_100_wheat_SEs[4,], lwd = 1) #RWC treatment only had 4 available reps; hence df=3
segments(x0=century_bar_dims[5,], y0=soilC_0_100_wheat_means[5,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[5,], x1=century_bar_dims[5,], y1=soilC_0_100_wheat_means[5,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[5,], lwd = 1)
arrows(x0=century_bar_dims[1,], y0=soilC_0_100_wheat_means[1,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[1,], x1=century_bar_dims[1,], y1=soilC_0_100_wheat_means[1,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[1,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[2,], y0=soilC_0_100_wheat_means[2,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[2,], x1=century_bar_dims[2,], y1=soilC_0_100_wheat_means[2,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[2,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[3,], y0=soilC_0_100_wheat_means[3,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[3,], x1=century_bar_dims[3,], y1=soilC_0_100_wheat_means[3,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[3,], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=century_bar_dims[4,], y0=soilC_0_100_wheat_means[4,] + qt(0.975,df=3) * soilC_0_100_wheat_SEs[4,], x1=century_bar_dims[4,], y1=soilC_0_100_wheat_means[4,] + qt(0.025,df=3) * soilC_0_100_wheat_SEs[4,], lwd = 1, angle = 90, code = 3, length = 0.05) #RWC treatment only had 4 available reps; hence df=3
arrows(x0=century_bar_dims[5,], y0=soilC_0_100_wheat_means[5,] + qt(0.975,df=5) * soilC_0_100_wheat_SEs[5,], x1=century_bar_dims[5,], y1=soilC_0_100_wheat_means[5,] + qt(0.025,df=5) * soilC_0_100_wheat_SEs[5,], lwd = 1, angle = 90, code = 3, length = 0.05)
dev.off()

#t-tests to try to reproduce tautges et al. 2019
deltaC_t.test <- function(df_93, df_12, depth_class) {
  if(depth_class=='60_100') {
    df_93 <- df_93[!df_93$Plot=='2_4',]
    df_12 <- df_12[!df_12$Plot=='2_4',]
  }
  trtments <- unique(df_93$Treatment)
  results <- data.frame(Treatments=trtments, meanC_93=NA, sdevC_93=NA, meanC_12=NA, sdevC_12=NA, DeltaC=NA, low.CI.95=NA, high.CI.95=NA, p.value=NA, stringsAsFactors = FALSE)
  for (i in seq_along(trtments)) {
    if(!all(df_12$Plot[df_12$Treatment==trtments[i]]==df_93$Plot[df_93$Treatment==trtments[i]])) {
      stop(print(paste('Plots do not align for', trtments[i])))
    }
    result <- t.test(df_12$soilC[df_12$Treatment==trtments[i]], df_93$soilC[df_93$Treatment==trtments[i]], paired = TRUE)
    results$meanC_93[i] <- mean(df_93$soilC[df_93$Treatment==trtments[i]], na.rm = TRUE)
    results$sdevC_93[i] <- sd(df_93$soilC[df_93$Treatment==trtments[i]], na.rm = TRUE)
    results$meanC_12[i] <- mean(df_12$soilC[df_12$Treatment==trtments[i]], na.rm = TRUE)
    results$sdevC_12[i] <- sd(df_12$soilC[df_12$Treatment==trtments[i]], na.rm = TRUE)
    results$DeltaC[i] <- result$estimate
    results$low.CI.95[i] <- result$conf.int[1]
    results$high.CI.95[i] <- result$conf.int[2]
    results$p.value[i] <- result$p.value
  }
  write.csv(results, file.path(dataDir, 'CenturyResults/t.tests/C_stocks', paste0('t.test_', depth_class, '_Cstocks.csv')), row.names=FALSE)
}
deltaC_t.test(df_93 = soilC_0_15_1993, df_12 = soilC_0_15_2012, depth_class = '0_15')
deltaC_t.test(df_93 = soilC_15_30_1993, df_12 = soilC_15_30_2012, depth_class = '15_30')
deltaC_t.test(df_93 = soilC_30_60_1993, df_12 = soilC_30_60_2012, depth_class = '30_60')
deltaC_t.test(df_93 = soilC_60_100_1993, df_12 = soilC_60_100_2012, depth_class = '60_100')
deltaC_t.test(df_93 = soilC_0_100_1993, df_12 = soilC_0_100_2012, depth_class = '0_100')
deltaC_t.test(df_93 = soilC_0_30_1993, df_12 = soilC_0_30_2012, depth_class = '0_30')