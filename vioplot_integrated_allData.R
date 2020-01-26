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
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}
mar_settings <- c(4, 4.5, 1, 1)
om_to_oc <- 1.72
#produced in ssurgo_calag_cluster_v2.R
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'v2 results', 'valley30cm_by_mukey_cluster.csv'), stringsAsFactors = FALSE)
order_lgnd_9 <- c(1,3,6,5,9,4,7,8,2)

kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_NArm.csv'), stringsAsFactors = FALSE)
match(kssl_points_30cm$cluster_9, order_lgnd_9)
kssl_points_30cm$xdim_vioplot <- match(kssl_points_30cm$cluster_9, order_lgnd_9)

kerri_points_30cm <- read.csv(file.path(kerriDir, 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
kerri_points_30cm$xdim_vioplot <- match(kerri_points_30cm$cluster_9, order_lgnd_9)

vioplot_mod_clus9_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, cdfa_pts, cdfa_varname, legend_plot, legendloc, legend_cex) {
  plot_order2 <- (1:9)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', 'validation plots', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_9==plot_order2[1]], times=round(df$area_ac[df$cluster_9==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[2]], times=round(df$area_ac[df$cluster_9==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[3]], times=round(df$area_ac[df$cluster_9==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[4]], times=round(df$area_ac[df$cluster_9==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[5]], times=round(df$area_ac[df$cluster_9==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[6]], times=round(df$area_ac[df$cluster_9==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[7]], times=round(df$area_ac[df$cluster_9==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[8]], times=round(df$area_ac[df$cluster_9==plot_order2[8]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[9]], times=round(df$area_ac[df$cluster_9==plot_order2[9]]/area_fact, 0)), col=c('lightgoldenrod', 'violetred', 'tan2', 'gray', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  points(x=kssl_df$xdim_vioplot[!is.na(kssl_df[[kssl_varname]])]-0.1, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.7, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_9, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.1, y=kssl_means$mean, pch=8, cex=0.9, col='orange')
  points(x=cdfa_pts$xdim_vioplot[!is.na(cdfa_pts[[cdfa_varname]])]+0.1, y=cdfa_pts[[cdfa_varname]][!is.na(cdfa_pts[[cdfa_varname]])], cex=0.7, pch=4, col='black')
  cdfa_means <- data.frame(mean=tapply(cdfa_pts[[cdfa_varname]], cdfa_pts$cluster_9, mean, na.rm=TRUE))
  cdfa_means$xdim_vioplot <- match(row.names(cdfa_means), plot_order)
  points(x=cdfa_means$xdim_vioplot+0.1, y=cdfa_means$mean, pch=8, cex=0.9, col='darkblue')
  if(legend_plot) {
    legend(x=legendloc, legend=c('SSURGO violin plots, area-weighted', 'KSSL data', 'KSSL mean', 'CDFA data', 'CDFA mean'), pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), cex=legend_cex)
  }
  dev.off()
}
vioplot_mod_clus9_validation(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,80), plot_order = order_lgnd_9, area_fact = 10, ylab='Clay (%)', fname='class9_clay_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.9, kssl_varname = 'clay_30cm', cdfa_varname = 'clay_30cm')
vioplot_mod_clus9_validation(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,80), plot_order = order_lgnd_9, area_fact = 10, ylab='Clay (%)', fname='class9_clay_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.9)
vioplot_mod_clus9_validation(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_9, area_fact = 10, ylab='soil pH', fname='class9_pH_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "pH_H2O_30cm", cdfa_varname = 'pH_H2O_30cm')
vioplot_mod_clus9_validation(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_9, area_fact = 10, ylab='Organic matter (%)', fname='class9_om_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "om_30cm", cdfa_varname = 'totC_30cm')




points(x=kssl_points_30cm$xdim_vioplot[!is.na(kssl_points_30cm$clay_30cm)]-0.1, y=kssl_points_30cm$clay_30cm[!is.na(kssl_points_30cm$clay_30cm)], cex=0.7, pch=1)
points(x=kssl_clay_means$xdim_vioplot-0.1, y=kssl_clay_means$mean_clay, pch=8, cex=0.9, c='orange')
points(x=kerri_points_30cm$xdim_vioplot[!is.na(kerri_points_30cm$clay_30cm)]+0.1, y=kerri_points_30cm$clay_30cm[!is.na(kerri_points_30cm$clay_30cm)], cex=0.7, pch=4)
kerri_clay_means <- data.frame(mean_clay=tapply(kerri_points_30cm$clay_30cm, kerri_points_30cm$cluster_9, mean, na.rm=TRUE))
kerri_clay_means$xdim_vioplot <- match(row.names(kerri_clay_means), order_lgnd_9)
points(x=kerri_clay_means$xdim_vioplot+0.1, y=kerri_clay_means$mean_clay, pch=8, cex=0.9, c='darkblue')
