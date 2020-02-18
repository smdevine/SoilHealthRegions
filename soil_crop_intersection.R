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
unique(crops$Crop2014)
sum(crops$Acres) #14202468

#performed intersection between crops and soils in ArcGIS Desktop 10.5 on 1/17/20 using Analysis:Overlay:Intersection
valley_crops_mu <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster_crops.shp')) #539297 little polygons
valley_crops_mu$area_ac <- area(valley_crops_mu) / 10000 * 2.47105
sum(valley_crops_mu$area_ac) #8354961
length(unique(valley_crops_mu$Crop2014)) #47 crops
unique(valley_crops_mu$Crop2014)

#9 cluster summary
crops_by_cluster <- tapply(valley_crops_mu$area_ac, list(valley_crops_mu$Crop2014, valley_crops_mu$cluster_9), sum)
dim(crops_by_cluster)
write.csv(crops_by_cluster, file.path(dataDir, 'crops_by_cluster.csv'))

#7 cluster summary
crops_by_cluster7 <- as.data.frame(tapply(valley_crops_mu$area_ac, list(valley_crops_mu$Crop2014, valley_crops_mu$cluster_7), sum, na.rm=TRUE))
clus_7_names <- c('3. Coarse w/pans', '6. Fine saline-sodic', '5. Coarse saline-sodic', '1. Coarse w/no restrictions', '7. Fine shrink-swell', '2. Loamy w/no restrictions', '4. Loamy w/pans')
order_lgnd_7 <- c(4,6,1,7,3,2,5)
colnames(crops_by_cluster7) <- clus_7_names
crops_by_cluster7 <- crops_by_cluster7[,order_lgnd_7]


row.names(crops_by_cluster7)
apply(crops_by_cluster7, 1, function(x) round(sum(x, na.rm=TRUE), 0))

crops_simplifed <- c('Alfalfa and pasture', 'Almonds', 'Deciduous orchard', 'Deciduous orchard', 'Legumes', 'Small fruit and vegetables', 'Small fruit and vegetables', 'Deciduous orchard', 'Citrus', 'Small fruit and vegetables', 'Cereals and grain hay', 'Cotton', 'Specialty crops', 'Specialty crops', 'Grapes', 'Specialty crops', 'Idle', 'Specialty crops', 'Small fruit and vegetables', 'Managed wetland', 'Small fruit and vegetables', 'Deciduous orchard', 'Cereals and grain hay', 'Cereals and grain hay', 'Alfalfa and pasture', 'Specialty crops', 'Small fruit and vegetables', 'Alfalfa and pasture', 'Olives', 'Small fruit and vegetables', 'Deciduous orchard', 'Deciduous orchard', 'Small fruit and vegetables', 'Pistachios', 'Deciduous orchard', 'Deciduous orchard', 'Small fruit and vegetables', 'Rice', 'Oil crops', 'Small fruit and vegetables', 'Oil crops', 'Tomatoes', 'Urban', 'Walnuts', 'Cereals and grain hay', 'Rice', 'Young perennials')
crops_simplifed2 <- c('Alfalfa and pasture', 'Orchard', 'Orchard', 'Orchard', 'Annual field crops', 'Small fruit and vegetables', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Specialty crops', 'Specialty crops', 'Grapes', 'Specialty crops', 'Idle', 'Specialty crops', 'Small fruit and vegetables', 'Managed wetland', 'Small fruit and vegetables', 'Orchard', 'Annual field crops', 'Annual field crops', 'Alfalfa and pasture', 'Specialty crops', 'Small fruit and vegetables', 'Alfalfa and pasture', 'Orchard', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Orchard', 'Orchard', 'Orchard', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Small fruit and vegetables', 'Annual field crops', 'Annual field crops', 'Urban', 'Orchard', 'Annual field crops', 'Annual field crops', 'Young perennials')
length(crops_simplifed)
length(crops_simplifed2)
unique(crops_simplifed) #20
unique(crops_simplifed2) #10

crops_by_cluster7$crop_class <- crops_simplifed
crops_by_cluster7$crop_class2 <- crops_simplifed2
write.csv(crops_by_cluster7, file.path(dataDir, 'crops', 'crops_by_cluster_7SHR.csv'), row.names = TRUE)

valley_crops_mu$CropClass <- crops_by_cluster7$crop_class2[match(valley_crops_mu$Crop2014, row.names(crops_by_cluster7))]

cropClass_by_cluster7 <- as.data.frame(tapply(valley_crops_mu$area_ac, list(valley_crops_mu$CropClass, valley_crops_mu$cluster_7), sum, na.rm=TRUE))
colnames(cropClass_by_cluster7) <- clus_7_names
cropClass_by_cluster7 <- cropClass_by_cluster7[,order_lgnd_7]
write.csv(cropClass_by_cluster7, file.path(dataDir, 'crops', 'cropClass_by_cluster_7SHR.csv'), row.names = TRUE)
