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
crit_pH <- 7.8
#produced in ssurgo_calag_cluster_FINAL.R
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)

#updated 4/7/20 b/c of set.seed issue
clus_7_colors <- c('deepskyblue', 'gold', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1')
order_lgnd_7 <- c(4,5,2,3,7,1,6)

kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_FINAL.csv'), stringsAsFactors = FALSE) #oc and kgSOC_m2 updated 2/11/20
#replace OC with totC when pH sufficiently low
kssl_points_30cm$oc_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)] <- kssl_points_30cm$c_tot_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)]
sum(is.na(kssl_points_30cm$oc)) #NAs reduced from 88 to 56
#replace estimated om with oc times assumption (see 'aov_soil_properties_KSSL_CDFA.R' for details comparing KSSL OM est with this)
kssl_points_30cm$om_30cm <- kssl_points_30cm$oc_30cm * om_to_oc

# kssl_points_30cm$xdim_vioplot9 <- match(kssl_points_30cm$cluster_9, order_lgnd_9)
kssl_points_30cm$xdim_vioplot7 <- match(kssl_points_30cm$cluster_7, order_lgnd_7)
# kssl_points_30cm$xdim_vioplot6 <- match(kssl_points_30cm$cluster_6, order_lgnd_6)
# kssl_points_30cm$xdim_vioplot5 <- match(kssl_points_30cm$cluster_5, order_lgnd_5)
kssl_points_30cm$lep_30cm <- kssl_points_30cm$lep_30cm*100 #to match SSURGO scale

kerri_points_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
# kerri_points_30cm$xdim_vioplot9 <- match(kerri_points_30cm$cluster_9, order_lgnd_9)
kerri_points_30cm$xdim_vioplot7 <- match(kerri_points_30cm$cluster_7, order_lgnd_7)
# kerri_points_30cm$xdim_vioplot6 <- match(kerri_points_30cm$cluster_6, order_lgnd_6)
# kerri_points_30cm$xdim_vioplot5 <- match(kerri_points_30cm$cluster_5, order_lgnd_5)
colnames(kerri_points_30cm)

#7-region plots
table(kerri_points_30cm$cluster_7)
table(kssl_points_30cm$cluster_7)
table(kssl_points_30cm$cluster_7[!is.na(kssl_points_30cm$kgOrg.m2_30cm)])
table(kssl_points_30cm$cluster_7[!is.na(kssl_points_30cm$om_30cm)])
plot(x=rep(1, length(kssl_points_30cm$pH_H2O_30cm[kssl_points_30cm$xdim_vioplot7==4])), y=kssl_points_30cm$pH_H2O_30cm[kssl_points_30cm$xdim_vioplot7==4])

vioplot_mod_clus7_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, cdfa_pts, cdfa_varname, legend_plot, legendloc, legend_cex, sig_labels, fig_label) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'FINAL', 'validation plots', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.1, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.7, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.1, y=kssl_means$mean, pch=8, cex=0.9, col='orange')
  points(x=cdfa_pts$xdim_vioplot7[!is.na(cdfa_pts[[cdfa_varname]])]+0.1, y=cdfa_pts[[cdfa_varname]][!is.na(cdfa_pts[[cdfa_varname]])], cex=0.7, pch=4, col='black')
  cdfa_means <- data.frame(mean=tapply(cdfa_pts[[cdfa_varname]], cdfa_pts$cluster_7, mean, na.rm=TRUE))
  cdfa_means$xdim_vioplot <- match(row.names(cdfa_means), plot_order)
  points(x=cdfa_means$xdim_vioplot+0.1, y=cdfa_means$mean, pch=8, cex=0.9, col='darkblue')
  text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5)
  legend('topright', fig_label, bty='n', inset=0.005)
  if(legend_plot) {
    legend(x=legendloc, legend=c('SSURGO violin plots, area-weighted', 'KSSL data', 'KSSL mean', 'CDFA data', 'CDFA mean'), pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), cex=legend_cex)
  }
  dev.off()
}
vioplot_mod_clus7_validation(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(-2.5,80), plot_order = order_lgnd_7, area_fact = 10, ylab='Clay (%)', fname='class7_clay_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.9, kssl_varname = 'clay_30cm', cdfa_varname = 'clay_30cm', sig_labels = c('A', 'B', 'A', 'AB', 'A', 'C', 'D'), fig_label = 'a')

vioplot_mod_clus7_validation(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(4.9,10.1), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Soil pH'[H2O]), fname='class7_pH_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "pH_H2O_30cm", cdfa_varname = 'pH_H2O_30cm', sig_labels = c('CD', 'BC', 'A', 'AB', 'E', 'E', 'D'), fig_label = 'g')

vioplot_mod_clus7_validation(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(-0.2,12), plot_order = order_lgnd_7, area_fact = 10, ylab='Organic matter (%)', fname='class7_om_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "om_30cm", cdfa_varname = 'totC_30cm', sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'b')

vioplot_mod_clus7_validation(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.79,1.85), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class7_bd_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "bd_13b_30cm", cdfa_varname = 'bd_30cm', sig_labels = c('AB', 'A', 'B', 'AB', 'AB', 'AB', 'A'), fig_label = 'h')

vioplot_mod_clus7_validation(valley30cm_by_mukey, 'kgOrg.m2_30cm', ylim_vioplot = c(-0.15,12), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('0-30 cm SOC (kg m'^-2*')'), fname='class7_kgSOC_0_30cm_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "kgOrg.m2_30cm", cdfa_varname = 'kgOrg.m2_30cm', sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'c')

#KSSL only validation for 7-region model
vioplot_mod_clus7_KSSL_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, legend_plot, legendloc, legend_cex, sig_labels, fig_label) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'FINAL', 'validation plots', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.1, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.7, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.1, y=kssl_means$mean, pch=8, cex=0.9, col='orange')
  text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5)
  legend('topright', fig_label, bty='n', inset=0.005)
  if(legend_plot) {
    legend(x=legendloc, legend=c('SSURGO violin plots, area-weighted', 'KSSL data', 'KSSL mean'), pch=c(NA,1,8), col = c(NA, 'black', 'orange'), cex=legend_cex)
  }
  dev.off()
}

vioplot_mod_clus7_KSSL_validation(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(-1,60), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100 g'^-1*'soil)'), fname='class7_CEC_vioplots_KSSL_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = 'cec_7_30cm', sig_labels = c('A', 'B', 'A', 'AB', 'A', 'C', 'C'), fig_label = 'd')

vioplot_mod_clus7_KSSL_validation(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class7_AWC_vioplots_KSSL_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = 'awc_30cm', sig_labels = NA, fig_label='f') #no significant contrasts


vioplot_mod_clus7_KSSL_validation(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(-2,50), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Electrical conductivity (mmho cm'^-1*')'), fname='class7_EC_vioplots_KSSL_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = 'ec_30cm', sig_labels = c(rep('A', 3), 'AB', 'A', 'B', 'A'), fig_label = 'e')

vioplot_mod_clus7_KSSL_validation(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(-0.5,17), plot_order = order_lgnd_7, area_fact = 10, ylab='Linear extensibility (%)', fname='class7_lep_vioplots_KSSL_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = 'lep_30cm', sig_labels = c('A', 'BC', 'A', 'AB', 'AB', 'CD', 'D'), fig_label = 'f')

#no validation data plots
vioplot_mod_clus7 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, fig_label) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'FINAL', 'validation plots', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  legend('topright', fig_label, bty='n', inset=0.005)
  dev.off()
}
vioplot_mod_clus7(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_7, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class7_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1), fig_label = 'j')
valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus7(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class7_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1), fig_label = 'i')
vioplot_mod_clus7(valley30cm_by_mukey, 'storiemn', ylim_vioplot = c(0,100), plot_order = order_lgnd_7, area_fact = 10, ylab='Storie index', fname='class7_storie_vioplots.tif', mar=c(3.5, 4.25, 1, 1), fig_label = 'k')
