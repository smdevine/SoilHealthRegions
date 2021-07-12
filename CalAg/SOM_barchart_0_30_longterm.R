library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after blue <- rgb(0, 0, 1, alpha=0.5)updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/Figures'
dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/long term studies'

SOM_data <- read.csv(file.path(dataDir, 'long_term_0_30_SOM.csv'), stringsAsFactors = FALSE)
SOM_data$bar_color <- ifelse(SOM_data$Location=='LTRAS', 'tan4', ifelse(SOM_data$Location=='SOCS', 'lightgoldenrod', ifelse(SOM_data$Location=='WSREC', 'tan4', NA)))

# transparent <- rgb(0,0,0,0)
cex_labels <- 0.8
tiff(file = file.path(FiguresDir, 'SOM_0_30_change_LT.tif'), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(3,3.5,0.1,0.1))
bar_dims <- barplot(SOM_data$SOM_start, beside=TRUE, space=c(0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), col=SOM_data$bar_color, ylab='', ylim=c(0,2.25), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
legend(x=16, y=1.8, legend = c('region 1', 'region 2'),  col=c('lightgoldenrod', 'tan4'), lty=2, lwd=2, xjust=0, yjust=0, cex = 0.9, title = 'KSSL median SOM') #
barplot(SOM_data$SOM_end, beside=TRUE, space=c(0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), border=TRUE, angle = 45, density = 10, col='black', ylab='', ylim=c(0,2.25), axes=FALSE, cex.names=1, add = TRUE)
mtext('0-12" (0-30 cm) soil organic matter (%)', side = 2, line = 2.25)
mtext(SOM_data$line1, side=1, line=-0.15, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line2, side=1, line=0.5, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line3, side=1, line=1.15, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line4, side=1, line=1.8, at=bar_dims, cex=cex_labels)
abline(h=1.65, col='tan4', lty=2, lwd=2)
abline(h=1.07, col='lightgoldenrod', lty=2, lwd=2)
text(x=2.3, y=1.75, 'WSREC')
text(x=7.4, y=2.05, 'SOCS')
text(x=14.8, y=1.75, 'LTRAS')
dev.off()
bar_dims

cex_labels <- 0.8
cairo_pdf(file = file.path(FiguresDir, 'PDF versions', 'Figure_4_FINAL.pdf'), pointsize = 12, family = 'Times New Roman', width = 6, height = 4)
par(mar=c(3,3.5,0.1,0.1))
bar_dims <- barplot(SOM_data$SOM_start, beside=TRUE, space=c(0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), col=SOM_data$bar_color, ylab='', ylim=c(0,2.25), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
legend(x=15.5, y=1.75, legend = c('region 1', 'region 2'),  col=c('lightgoldenrod', 'tan4'), lty=2, lwd=2, xjust=0, yjust=0, cex = 0.9, title = 'KSSL median SOM', bty='n') #
barplot(SOM_data$SOM_end, beside=TRUE, space=c(0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1), border=TRUE, angle = 45, density = 10, col='black', ylab='', ylim=c(0,2.25), axes=FALSE, cex.names=1, add = TRUE)
mtext('0-12" (0-30 cm) soil organic matter (%)', side = 2, line = 2.25)
mtext(SOM_data$line1, side=1, line=-0.15, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line2, side=1, line=0.5, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line3, side=1, line=1.15, at=bar_dims, cex=cex_labels)
mtext(SOM_data$line4, side=1, line=1.8, at=bar_dims, cex=cex_labels)
abline(h=1.65, col='tan4', lty=2, lwd=2)
abline(h=1.07, col='lightgoldenrod', lty=2, lwd=2)
text(x=2.3, y=1.75, 'WSREC')
text(x=7.4, y=2.05, 'SOCS')
text(x=14.8, y=1.75, 'LTRAS')
dev.off()