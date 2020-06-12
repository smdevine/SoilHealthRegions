#this modified to work with aggregated dataset from ssurgo_calag_aggregate.R
#TO-DO
#(1) run cluster analysis on Salinas only [DONE]
#(2) summarize classes for 4 and 5 
#(3) log transform om and ksat [DONE]
#(4) identify outliers
laptop <- TRUE
library(vioplot)
library(raster)
library(corrplot)
library(cluster)
library(factoextra)
library(fpc)
library(fmsb)
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
}

cluster_results <- read.csv(file.path(dataDir, 'FINAL results', 'cluster_analysis_summary.csv'))
#v2 of Fig1a
tiff(file = file.path(FiguresDir, 'FINAL', 'validation plots', 'kmeans_comparison_6.8.20.tif'), family = 'Times New Roman', width = 6.5, height = 3.5, pointsize = 11, units = 'in', res=800, compression='lzw')
par(mar=c(4, 4, 0.5, 8))
plot(2:20, cluster_results$btwnSS_totSS[2:20], type='b', xlab='', ylab='', cex.axis=1, cex.lab=1, xaxt='n', ylim=c(20,82), cex=0.8)
lines(2:20, cluster_results$CH_index[2:20]/30, type='b', col='blue', cex=0.8)
lines(2:20, cluster_results$asw[2:20]*200, type='b', col='red', cex=0.8)
lines(2:20, cluster_results$gap_stat[2:20]*50, type='b', col='grey', cex=0.8)
mtext(text = 'Number of soil health regions (cluster size in conceptual model)', side=1, line=2.5)
mtext(text = 'SSURGO variability captured by clusters (%)', side=2, line=2.5, at=45)
axis(side=4, at=c(20,40,60,80), labels = FALSE, line=0.25, col='blue', col.ticks = 'blue')
mtext(text=c(20*30, 40*30, 60*30, 80*30), side = 4, at=c(20,40,60,80), col='blue', line=0.75)
mtext(text = 'Calinski Harabasz index', side=4, line=1.5, at=50, col = 'blue')
axis(side=4, at=c(30,50,70), labels = FALSE, line=2.75, col='red', col.ticks = 'red')
mtext(text=c(30/200, 50/200, 70/200), side = 4, at=c(30,50,70), col='red', line=3.25)
mtext(text = 'Average silhouette width', side=4, at=50, col='red', line=4) 
axis(side=4, at=c(40,50,60), labels = FALSE, line=5.25, col='grey', col.ticks = 'grey')
mtext(text=c(40/50, 50/50, 60/50), side = 4, at=c(40,50,60), col='grey', line=5.75)
mtext(text = 'Gap statistic', side=4, at=50, col='grey', line=6.5)
axis(side=1, at=seq(from=2, to=20, by=2))
#text(2:20, cluster_results$btwnSS_totSS[2:20], labels=as.character(2:20), pos=3, offset=0.5)
text(x=2, y=80, 'a', adj=c(0,0))
dev.off()

#v3
CH_index_rescale <- mean(cluster_results$tot.withinss[2:20]/cluster_results$CH_index[2:20])
asw_rescale <- mean(cluster_results$tot.withinss[2:20]/cluster_results$asw[2:20])
gap_rescale <- mean(cluster_results$tot.withinss[2:20]/cluster_results$gap_stat[2:20])
tiff(file = file.path(FiguresDir, 'FINAL', 'validation plots', 'kmeans_comparison_v3.tif'), family = 'Times New Roman', width = 6.5, height = 3.5, pointsize = 11, units = 'in', res=800, compression='lzw')
par(mar=c(4, 4, 0.5, 8))
plot(2:20, cluster_results$tot.withinss[2:20], type='b', xlab='', ylab='', cex.axis=1, cex.lab=1, xaxt='n', cex=0.8)
lines(2:20, cluster_results$CH_index[2:20]*CH_index_rescale, type='b', col='blue', cex=0.8)
lines(2:20, cluster_results$asw[2:20]*asw_rescale, type='b', col='red', cex=0.8)
lines(2:20, cluster_results$gap_stat[2:20]*gap_rescale, type='b', col='grey', cex=0.8)
mtext(text = 'Number of soil health regions (cluster size in conceptual model)', side=1, line=2.5)
mtext(text = 'Total within-cluster sum of squares', side=2, line=2.5, at=22500)
axis(side=4, at=c(800, 1200, 1600, 2000)*CH_index_rescale, labels = FALSE, line=0.25, col='blue', col.ticks = 'blue')
mtext(text=c(800, 1200, 1600, 2000), side = 4, at=c(800, 1200, 1600, 2000)*CH_index_rescale, col='blue', line=0.75)
mtext(text = 'Calinski Harabasz index', side=4, line=1.5, at=1400*CH_index_rescale, col = 'blue')
axis(side=4, at=c(0.3,0.25,0.2)*asw_rescale, labels = FALSE, line=2.75, col='red', col.ticks = 'red')
mtext(text=c(0.3,0.25,0.2), side = 4, at=c(0.3,0.25,0.2)*asw_rescale, col='red', line=3.25)
mtext(text = 'Average silhouette width', side=4, at=0.25*asw_rescale, col='red', line=4) 
axis(side=4, at=c(0.8,1,1.2)*gap_rescale, labels = FALSE, line=5.25, col='grey', col.ticks = 'grey')
mtext(text=c(0.8,1,1.2), side = 4, at=c(0.8,1,1.2)*gap_rescale, col='grey', line=5.75)
mtext(text = 'Gap statistic', side=4, at=gap_rescale, col='grey', line=6.5)
axis(side=1, at=seq(from=2, to=20, by=2))
#text(2:20, cluster_results$btwnSS_totSS[2:20], labels=as.character(2:20), pos=3, offset=0.5)
text(x=20, y=30000, 'a', adj=c(0,0))
dev.off()
