library(prism)
library(raster)
library(rgeos)
library(rgdal)
library(aqp)
library(vioplot)
library(extrafont)
library(extrafontdb)
loadfonts(device = 'win')
laptop <- FALSE

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
  prismDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/PRISM'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  prismDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/PRISM'
}
ca_ta <- showP4(showWKT("+init=epsg:3310"))
spatialCIMISdir <- 'D:/Dissertation/Allowable_Depletion/SpatialCIMIS'
clus_7_colors <- c('deepskyblue', 'gold', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1')
order_lgnd_7 <- c(4,5,2,3,7,1,6)
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')
list.files(ksslDir)
#add climate info
kssl_SHR_shp <- shapefile(file.path(ksslDir, 'shapefiles', 'kssl_SHR7.shp'))
kssl_SHR_shp #369 features in +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0

#prism data
prism_ppt <- raster(file.path(list.dirs(prismDir)[grepl('ppt', list.dirs(prismDir))], 'PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'))
prism_tmean <- raster(file.path(list.dirs(prismDir)[grepl('tmean', list.dirs(prismDir))], 'PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'))
crs(prism_ppt) == crs(kssl_SHR_shp)
kssl_SHR_NAD83 <- spTransform(kssl_SHR_shp, crs(prism_ppt))

kssl_SHR_NAD83$annual.P <- extract(prism_ppt, kssl_SHR_NAD83)
kssl_SHR_NAD83$annual.T <- extract(prism_tmean, kssl_SHR_NAD83)
summary(kssl_SHR_NAD83$annual.P)
summary(kssl_SHR_NAD83$annual.T)
kssl_climate <- as.data.frame(kssl_SHR_NAD83)
write.csv(kssl_climate, file.path(ksslDir, 'kssl_climate_SHRpts.csv'), row.names = FALSE)

#spatial CIMIS
spatialCIMIS <- raster(file.path(spatialCIMISdir, 'U2/2004', 'U220040101.tif'))
crs(spatialCIMIS)
kssl_SHR_CA_TA <- spTransform(kssl_SHR_NAD83, crs(spatialCIMIS))
CIMIScellnumbers <- as.integer(cellFromXY(object=spatialCIMIS, xy=kssl_SHR_CA_TA))
length(unique(CIMIScellnumbers)) #284 CIMIS cells need sampling
CIMIScellunique_df <- data.frame(CIMIS_cells=unique(CIMIScellnumbers))
write.csv(CIMIScellunique_df, file.path(mainDir, 'soil health', 'CIMIS', 'CIMIS_cells_unique_KSSL.csv'), row.names = FALSE)
kssl_SHR_CA_TA$CIMIScell <- as.integer(cellFromXY(object=spatialCIMIS, xy=kssl_SHR_CA_TA))
sum(is.na(kssl_SHR_CA_TA$CIMIScell))
ETo <- read.csv(file.path(mainDir, 'soil health', 'CIMIS', 'ETo_daily_2004_2018_QCpass.csv'), stringsAsFactors = FALSE)
ETo_annual_sum <- lapply(ETo[,6:ncol(ETo)], function(x) {
  tapply(x, ETo$year, sum)
})
ETo_annual_sum <- as.data.frame(ETo_annual_sum)
ETo_annual_mean <- apply(ETo_annual_sum, 2, mean)
kssl_SHR_CA_TA$CIMIScell <- paste0('cell_', kssl_SHR_CA_TA$CIMIScell)
sum(!kssl_SHR_CA_TA$CIMIScell %in% names(ETo_annual_mean)) #40 are missing from dataset
kssl_SHR_CA_TA$annual.ETo <- ETo_annual_mean[match(kssl_SHR_CA_TA$CIMIScell, names(ETo_annual_mean))]
summary(kssl_SHR_CA_TA$annual.ETo) #40 NA
plot(kssl_SHR_CA_TA[is.na(kssl_SHR_CA_TA$annual.ETo),])
