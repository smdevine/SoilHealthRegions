library(raster)
laptop <- TRUE
ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl/shapefiles'
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

list.files(ksslDir)
kssl_shp <- shapefile(file.path(ksslDir, 'kssl_CA.shp'))
crs(kssl_shp)
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'))
crs(valley_mu_shp_30cm)
kssl_shp_WGS84 <- spTransform(kssl_shp, crs(valley_mu_shp_30cm))
test <- extract(valley_mu_shp_30cm, kssl_shp_WGS84, df=TRUE)
