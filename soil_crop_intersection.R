laptop <- FALSE
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
  LandIQDir <- 'D:/Dissertation/Allowable_Depletion/LandIQ.crops'
}
list.files(file.path(dataDir, 'shapefiles with data'))
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'))
crs(valley_mu_shp_30cm)
sum(valley_mu_shp_30cm$area_ac) #13873110
sum(valley_mu_shp_30cm$area_ac[is.na(valley_mu_shp_30cm$cluster_9)]) #839014.1
list.files(LandIQDir)
valley_mu_shp_30cm <- valley_mu_shp_30cm[!is.na(valley_mu_shp_30cm$cluster_9),]
sum(valley_mu_shp_30cm$area_ac)
shapefile(valley_mu_shp_30cm, file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster_noNA.shp'))

crops <- shapefile(file.path(LandIQDir, 'i15_Crop_Mapping_2014_Final_LandIQonAtlas.shp'))
crs(crops)

#performed intersection between crops and soils in ArcGIS Desktop 10.5 on 1/17/20 using Analysis:Overlay:Intersection
valley_crops_mu <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster_crops.shp')) #539297 little polygons
valley_crops_mu$area_ac <- area(valley_crops_mu) / 10000 * 2.47105
sum(valley_crops_mu$area_ac) #8354961
length(unique(valley_crops_mu$Crop2014)) #47 crops
unique(valley_crops_mu$Crop2014)

crops_by_cluster <- tapply(valley_crops_mu$area_ac, list(valley_crops_mu$Crop2014, valley_crops_mu$cluster_9), sum)
dim(crops_by_cluster)
write.csv(crops_by_cluster, file.path(dataDir, 'crops_by_cluster.csv'))
