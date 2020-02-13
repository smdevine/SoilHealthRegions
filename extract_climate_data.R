laptop <- TRUE
library(raster)
library(rgeos)
library(prism)
library(vioplot)
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  prismDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/PRISM'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
}
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'for cluster analysis', 'valley30cm_by_mukey_final_v2.csv'), stringsAsFactors = FALSE)
sum(valley30cm_by_mukey$area_ac) #13034096
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'))
sum(valley_mu_shp_30cm$area_ac) #13873110 because some eliminated from analysis
sum(valley_mu_shp_30cm$area_ac[is.na(valley_mu_shp_30cm$cluster_7)]) #839014.1 acres
centroids_mu <- gCentroid(valley_mu_shp_30cm, byid = TRUE)
centroids.sp.df <- SpatialPointsDataFrame(centroids_mu, data.frame(valley_mu_shp_30cm))

#read in PRISM data
options(prism.path=prismDir)
get_prism_normals(type='ppt', resolution = '800m', annual=TRUE, keepZip = TRUE)
get_prism_normals(type='tmean', resolution = '800m', annual=TRUE, keepZip = TRUE)

list.dirs(prismDir)
list.dirs(prismDir)[grepl('ppt', list.dirs(prismDir))]
list.dirs(prismDir)[grepl('tmean', list.dirs(prismDir))]
prism_ppt <- raster(file.path(list.dirs(prismDir)[grepl('ppt', list.dirs(prismDir))], 'PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'))
prism_tmean <- raster(file.path(list.dirs(prismDir)[grepl('tmean', list.dirs(prismDir))], 'PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'))
centroids.sp.df_NAD83 <- spTransform(centroids.sp.df, crs(prism_ppt))
centroids.sp.df_NAD83 <-centroids.sp.df_NAD83[centroids.sp.df_NAD83$mukey %in% valley30cm_by_mukey$mukey,]
sum(centroids.sp.df_NAD83$area_ac) #13,034,096
centroids.sp.df_NAD83$annual.P <- extract( prism_ppt, centroids.sp.df_NAD83)
centroids.sp.df_NAD83$annual.T <- extract( prism_tmean, centroids.sp.df_NAD83)

valley30cm_climate_df <- as.data.frame(centroids.sp.df_NAD83)

vioplot_mod_clus7_climate <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, legend_plot, legendloc, legend_cex) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', 'climate', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
clus_7_colors <- c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3')
order_lgnd_7 <- c(4,6,1,7,3,2,5)
vioplot_mod_clus7_climate(valley30cm_climate_df, 'annual.P', ylim_vioplot = c(0,1400), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Annual precipitation (mm yr'^-1*')'), fname='class7_annual_P.tif', mar=c(3.5, 4.25, 1, 1))
