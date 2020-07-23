library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final/CalAg'
dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/long term studies'
list.files(dataDir)
SJV_df <- read.csv(file.path(dataDir, "SJV_westside_tillage_CC.csv"), stringsAsFactors = FALSE)
head(SJV_df)

SJV_means <- SJV_df[,2:4]
SJV_means_by_yr <- cbind(SJV_means[1:4,3], SJV_means[5:8,3], SJV_means[9:12,3])
row.names(SJV_means_by_yr) <- SJV_means$Treatment[1:4]
colnames(SJV_means_by_yr) <- c('Year 0', 'Year 8', 'Year 13')
SJV_SEs <- SJV_df[,c(2:3,5)]
SJV_SEs_by_yr <- cbind(SJV_SEs[1:4,3], SJV_SEs[5:8,3])
row.names(SJV_SEs_by_yr) <- SJV_SEs$Treatment[1:4]
colnames(SJV_SEs_by_yr) <- c('Year 0', 'Year 8')

tiff(file = file.path(FiguresDir, 'soilC_change_SJV_tillage_CC.tif'), family = 'Times New Roman', width = 3.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,0.5,0.25))
sjv_bar_dims <- barplot(SJV_means_by_yr, beside=TRUE, legend.text=c('NTCC', 'STCC', 'NTNO', 'STNO'), col=c('forestgreen', 'chartreuse3', 'orange3', 'orange'), ylab='', ylim=c(0,40), args.legend=list(x='topleft', cex=0.9, bty='n', ncol=2, inset=c(0,0.05)), cex.names=1) #inset=c(0,-0.35) to place legend outside plot
mtext(expression('0-30 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=sjv_bar_dims[1,1:2], y0=SJV_means_by_yr[1,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[1,1:2], x1=sjv_bar_dims[1,1:2], y1=SJV_means_by_yr[1,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[1,1:2], lwd = 1)
segments(x0=sjv_bar_dims[2,1:2], y0=SJV_means_by_yr[2,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[2,1:2], x1=sjv_bar_dims[2,1:2], y1=SJV_means_by_yr[2,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[2,1:2], lwd = 1)
segments(x0=sjv_bar_dims[3,1:2], y0=SJV_means_by_yr[3,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[3,1:2], x1=sjv_bar_dims[3,1:2], y1=SJV_means_by_yr[3,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[3,1:2], lwd = 1)
segments(x0=sjv_bar_dims[4,1:2], y0=SJV_means_by_yr[4,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[4,1:2], x1=sjv_bar_dims[4,1:2], y1=SJV_means_by_yr[4,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[4,1:2], lwd = 1)
arrows(x0=sjv_bar_dims[1,1:2], y0=SJV_means_by_yr[1,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[1,1:2], x1=sjv_bar_dims[1,1:2], y1=SJV_means_by_yr[1,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[1,1:2], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=sjv_bar_dims[2,1:2], y0=SJV_means_by_yr[2,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[2,1:2], x1=sjv_bar_dims[2,1:2], y1=SJV_means_by_yr[2,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[2,1:2], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=sjv_bar_dims[3,1:2], y0=SJV_means_by_yr[3,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[3,1:2], x1=sjv_bar_dims[3,1:2], y1=SJV_means_by_yr[3,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[3,1:2], lwd = 1, angle = 90, code = 3, length = 0.05)
arrows(x0=sjv_bar_dims[4,1:2], y0=SJV_means_by_yr[4,1:2] + qt(0.975,df=7) * SJV_SEs_by_yr[4,1:2], x1=sjv_bar_dims[4,1:2], y1=SJV_means_by_yr[4,1:2] + qt(0.025,df=7) * SJV_SEs_by_yr[4,1:2], lwd = 1, angle = 90, code = 3, length = 0.05)
dev.off()

#Salinas plot
salinas_df <- read.csv(file.path(dataDir, 'Salinas_compost_CC.csv'), stringsAsFactors = FALSE)
salinas_means <- salinas_df[,2:4]
salinas_means_by_year <- cbind(salinas_means[1:5,3], salinas_means[6:10,3], salinas_means[11:15,3], salinas_means[16:20,3], salinas_means[21:25,3], salinas_means[41:45,3])
row.names(salinas_means_by_year) <- salinas_df$Treatment[1:5]
colnames(salinas_means_by_year) <- c('Year 0', 'Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 8')
salinas_means_by_year <- salinas_means_by_year[-c(4:5),]
salinas_means_by_year

salinas_SEs <- salinas_df[,c(2:3,5)]
salinas_SEs_by_year <- cbind(salinas_SEs[1:5,3], salinas_SEs[6:10,3], salinas_SEs[11:15,3], salinas_SEs[16:20,3], salinas_SEs[21:25,3], salinas_SEs[41:45,3])
row.names(salinas_SEs_by_year) <- salinas_df$Treatment[1:5]
colnames(salinas_SEs_by_year) <- c('Year 0', 'Year 1', 'Year 2', 'Year 3', 'Year 4', 'Year 8')
salinas_SEs_by_year <- salinas_SEs_by_year[-c(4:5),]

tiff(file = file.path(FiguresDir, 'soilC_change_salinas_compost_CC.tif'), family = 'Times New Roman', width = 4.5, height = 4, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(2.5,3.5,0.25,0.25))
salinas_bar_dims <- barplot(salinas_means_by_year, beside=TRUE, names.arg = NULL, legend.text=c('Organic system 1', 'Organic system 2', 'Organic system 3'), col = c('green2', 'green4', 'tan2'), ylab='', ylim=c(0,63), args.legend=list(x='topright', cex=0.9, bty='n', ncol=1, inset=c(0,0)), cex.names=1) #)c('4th yr. rye-legume CC', 'Compost + 4th yr. rye-legume CC', 'Compost + ann. rye-legume CC', 'Compost + ann. rye CC') title="Salinas intensive vegetable system experiment"
mtext(expression('0-30 cm soil organic carbon (Mg '~ha^-1*')'), side = 2, line = 2.25)
segments(x0=salinas_bar_dims[1,], y0=salinas_means_by_year[1,] + qt(0.975,df=3) * salinas_SEs_by_year[1,], x1=salinas_bar_dims[1,], y1=salinas_means_by_year[1,] + qt(0.025,df=3) * salinas_SEs_by_year[1,], lwd = 1)
segments(x0=salinas_bar_dims[2,], y0=salinas_means_by_year[2,] + qt(0.975,df=3) * salinas_SEs_by_year[2,], x1=salinas_bar_dims[2,], y1=salinas_means_by_year[2,] + qt(0.025,df=3) * salinas_SEs_by_year[2,], lwd = 1)
segments(x0=salinas_bar_dims[3,], y0=salinas_means_by_year[3,] + qt(0.975,df=3) * salinas_SEs_by_year[3,], x1=salinas_bar_dims[3,], y1=salinas_means_by_year[3,] + qt(0.025,df=3) * salinas_SEs_by_year[3,], lwd = 1)
# segments(x0=salinas_bar_dims[4,], y0=salinas_means_by_year[4,] + qt(0.975,df=3) * salinas_SEs_by_year[4,], x1=salinas_bar_dims[4,], y1=salinas_means_by_year[4,] + qt(0.025,df=3) * salinas_SEs_by_year[4,], lwd = 1)
arrows(x0=salinas_bar_dims[1,], y0=salinas_means_by_year[1,] + qt(0.975,df=3) * salinas_SEs_by_year[1,], x1=salinas_bar_dims[1,], y1=salinas_means_by_year[1,] + qt(0.025,df=3) * salinas_SEs_by_year[1,], lwd = 1, angle = 90, code = 3, length = 0.04)
arrows(x0=salinas_bar_dims[2,], y0=salinas_means_by_year[2,] + qt(0.975,df=3) * salinas_SEs_by_year[2,], x1=salinas_bar_dims[2,], y1=salinas_means_by_year[2,] + qt(0.025,df=3) * salinas_SEs_by_year[2,], lwd = 1, angle = 90, code = 3, length = 0.04)
arrows(x0=salinas_bar_dims[3,], y0=salinas_means_by_year[3,] + qt(0.975,df=3) * salinas_SEs_by_year[3,], x1=salinas_bar_dims[3,], y1=salinas_means_by_year[3,] + qt(0.025,df=3) * salinas_SEs_by_year[3,], lwd = 1, angle = 90, code = 3, length = 0.04)
# arrows(x0=salinas_bar_dims[4,], y0=salinas_means_by_year[4,] + qt(0.975,df=3) * salinas_SEs_by_year[4,], x1=salinas_bar_dims[4,], y1=salinas_means_by_year[4,] + qt(0.025,df=3) * salinas_SEs_by_year[4,], lwd = 1, angle = 90, code = 3, length = 0.04)
dev.off()
