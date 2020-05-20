#re-ran 5/4/20 with updated SHR shapefile
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
list.files(file.path(dataDir, 'FINAL results', 'shapefiles with data'))
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster.shp')) #updated file path on 5/4/20
crs(valley_mu_shp_30cm)
sum(valley_mu_shp_30cm$area_ac) #13873110
sum(valley_mu_shp_30cm$area_ac[is.na(valley_mu_shp_30cm$cluster_7)]) #839014.1
list.files(LandIQDir)
valley_mu_shp_30cm <- valley_mu_shp_30cm[!is.na(valley_mu_shp_30cm$cluster_7),]
sum(valley_mu_shp_30cm$area_ac) #13034096
sum(area(valley_mu_shp_30cm) / 10000 * 2.47105)
shapefile(valley_mu_shp_30cm, file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster_noNA.shp'))

crops <- shapefile(file.path(LandIQDir, 'i15_Crop_Mapping_2014_Final_LandIQonAtlas.shp'))
crs(crops)
unique(crops$Crop2014)
sum(crops$Acres) #14202468
sum(area(crops) / 10000 * 2.47105) #14202446


#performed intersection between crops and soils in ArcGIS Desktop 10.5 on 1/17/20 using Analysis:Overlay:Intersection
#performed intersection between crops and soils in ArcGIS Desktop 10.5 on 5/4/20 using Analysis:Overlay:Intersection
valley_crops_mu <- shapefile(file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster_crops.shp')) #539297 little polygons
valley_crops_mu$area_ac <- area(valley_crops_mu) / 10000 * 2.47105
sum(valley_crops_mu$area_ac) #8354961
length(unique(valley_crops_mu$Crop2014)) #47 crops
unique(valley_crops_mu$Crop2014)
sum(valley_crops_mu$area_ac[valley_crops_mu$Crop2014=='Almonds']) #1092400
sum(crops$Acres[crops$Crop2014=='Almonds']) #1127946
sum(valley_crops_mu$area_ac[valley_crops_mu$Crop2014=='Walnuts']) #349641
sum(crops$Acres[crops$Crop2014=='Walnuts']) #370363

#9 cluster summary
# crops_by_cluster <- tapply(valley_crops_mu$area_ac, list(valley_crops_mu$Crop2014, valley_crops_mu$cluster_9), sum)
# dim(crops_by_cluster)
# write.csv(crops_by_cluster, file.path(dataDir, 'crops_by_cluster.csv'))

#7 cluster summary
crops_by_cluster7 <- as.data.frame(tapply(valley_crops_mu$area_ac, list(valley_crops_mu$Crop2014, valley_crops_mu$cluster_7), sum, na.rm=TRUE))
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic') #updated 5/4/20
order_lgnd_7 <- c(4,5,2,3,7,1,6)
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
write.csv(crops_by_cluster7, file.path(dataDir, 'FINAL results', 'crops', 'crops_by_cluster_7SHR.csv'), row.names = TRUE)

valley_crops_mu$CropClass <- crops_by_cluster7$crop_class2[match(valley_crops_mu$Crop2014, row.names(crops_by_cluster7))]

cropClass_by_cluster7 <- as.data.frame(tapply(valley_crops_mu$area_ac, list(valley_crops_mu$CropClass, valley_crops_mu$cluster_7), sum, na.rm=TRUE))
colnames(cropClass_by_cluster7) <- clus_7_names
cropClass_by_cluster7 <- cropClass_by_cluster7[,order_lgnd_7]
write.csv(cropClass_by_cluster7, file.path(dataDir, 'FINAL results', 'crops', 'cropClass_by_cluster_7SHR.csv'), row.names = TRUE)


  
#SAGBI index
SAGBIdir <- 'C:/Users/smdevine/Desktop/SpatialData/SAGBI'
list.files(SAGBIdir)
SAGBMI_unmod <- shapefile(file.path(SAGBIdir, 'sagbi_unmod.shp'))
names(SAGBMI_unmod)
SAGBI_unmod_df <- as.data.frame(SAGBMI_unmod)
mukeys_unique <- unique(SAGBI_unmod_df$mukey)
SAGBI_unmod_by_mukey <- SAGBI_unmod_df[match(mukeys_unique, SAGBI_unmod_df$mukey),]

SAGBI_mod <- shapefile(file.path(SAGBIdir, 'sagbi_mod.shp'))
SAGBI_mod_df <- as.data.frame(SAGBI_mod)
mukeys_unique <- unique(SAGBI_mod_df$mukey)
SAGBI_mod_by_mukey <- SAGBI_mod_df[match(mukeys_unique, SAGBI_mod_df$mukey),]

#add SAGBI info to valley_crops_mu
valley_crops_mu$SAGBI_unmod <- SAGBI_unmod_by_mukey$sagbi[match(valley_crops_mu$mukey, SAGBI_unmod_by_mukey$mukey)]
valley_crops_mu$SAGBI_unmod_class <- SAGBI_unmod_by_mukey$rat_grp[match(valley_crops_mu$mukey, SAGBI_unmod_by_mukey$mukey)]
valley_crops_mu$SAGBI_mod <- SAGBI_mod_by_mukey$sagbi[match(valley_crops_mu$mukey, SAGBI_mod_by_mukey$mukey)]
valley_crops_mu$SAGBI_mod_class <- SAGBI_mod_by_mukey$rat_grp[match(valley_crops_mu$mukey, SAGBI_mod_by_mukey$mukey)]

#summaries by crops (almonds)
almonds_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Almonds',])
dim(almonds_df) #80789 rows
length(unique(almonds_df$mukey)) #1557 unique mukeys
almonds_by_SHR7 <- as.data.frame(t(tapply(almonds_df$area_ac, almonds_df$cluster_7, sum)))
colnames(almonds_by_SHR7) <- clus_7_names
almonds_by_SHR7 <- almonds_by_SHR7[ ,order_lgnd_7]
almonds_by_SHR7
almonds_by_sagbi_unmod <- as.data.frame(t(tapply(almonds_df$area_ac, almonds_df$SAGBI_unmod_class, sum)))
almonds_by_sagbi_unmod
almonds_by_sagbi_mod <- as.data.frame(t(tapply(almonds_df$area_ac, almonds_df$SAGBI_mod_class, sum)))
almonds_by_sagbi_mod
round(almonds_by_sagbi_mod/sum(almonds_by_sagbi_mod), 2) #74% in excellent-to-moderately good (81% through moderately poor)

grapes_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Grapes',])
dim(grapes_df) #89612 rows
length(unique(grapes_df$mukey)) #1799 unique mukeys
grapes_by_SHR7 <- as.data.frame(t(tapply(grapes_df$area_ac, grapes_df$cluster_7, sum)))
colnames(grapes_by_SHR7) <- clus_7_names
grapes_by_SHR7 <- grapes_by_SHR7[ ,order_lgnd_7]
grapes_by_SHR7
grapes_by_sagbi_unmod <- as.data.frame(t(tapply(grapes_df$area_ac, grapes_df$SAGBI_unmod_class, sum)))
grapes_by_sagbi_unmod
grapes_by_sagbi_mod <- as.data.frame(t(tapply(grapes_df$area_ac, grapes_df$SAGBI_mod_class, sum)))
grapes_by_sagbi_mod
round(grapes_by_sagbi_mod/sum(grapes_by_sagbi_mod),2) #73% in Excellent to Moderately Good (80% to Moderately Poor)

walnuts_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Walnuts',])
dim(walnuts_df) #29151 rows
length(unique(walnuts_df$mukey)) #1259 unique mukeys
walnuts_by_SHR7 <- as.data.frame(t(tapply(walnuts_df$area_ac, walnuts_df$cluster_7, sum)))
colnames(walnuts_by_SHR7) <- clus_7_names
walnuts_by_SHR7 <- walnuts_by_SHR7[ ,order_lgnd_7]
walnuts_by_SHR7
walnuts_by_sagbi_unmod <- as.data.frame(t(tapply(walnuts_df$area_ac, walnuts_df$SAGBI_unmod_class, sum)))
walnuts_by_sagbi_unmod
round(walnuts_by_sagbi_unmod/sum(walnuts_by_sagbi_unmod), 2)
walnuts_by_sagbi_mod <- as.data.frame(t(tapply(walnuts_df$area_ac, walnuts_df$SAGBI_mod_class, sum)))
walnuts_by_sagbi_mod
round(walnuts_by_sagbi_mod/sum(walnuts_by_sagbi_mod), 2) #64% in Excellent-to-Moderately Good (72% to Moderately Poor)

tomatoes_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Tomatoes',])
dim(tomatoes_df) #10525 rows
length(unique(tomatoes_df$mukey)) #674 unique mukeys
tomatoes_by_SHR7 <- as.data.frame(t(tapply(tomatoes_df$area_ac, tomatoes_df$cluster_7, sum)))
colnames(tomatoes_by_SHR7) <- clus_7_names
tomatoes_by_SHR7 <- tomatoes_by_SHR7[ ,order_lgnd_7]
tomatoes_by_SHR7
tomatoes_by_sagbi_unmod <- as.data.frame(t(tapply(tomatoes_df$area_ac, tomatoes_df$SAGBI_unmod_class, sum)))
tomatoes_by_sagbi_unmod
round(tomatoes_by_sagbi_unmod/sum(tomatoes_by_sagbi_unmod), 2) #only 37% in excellent ot moderately good (unmodified)
tomatoes_by_sagbi_mod <- as.data.frame(t(tapply(tomatoes_df$area_ac, tomatoes_df$SAGBI_mod_class, sum)))
tomatoes_by_sagbi_mod
round(tomatoes_by_sagbi_mod/sum(tomatoes_by_sagbi_mod), 2)

corn_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Corn, Sorghum and Sudan',])
dim(corn_df) #32471 rows
length(unique(corn_df$mukey)) #1239 unique mukeys
corn_by_SHR7 <- as.data.frame(t(tapply(corn_df$area_ac, corn_df$cluster_7, sum)))
colnames(corn_by_SHR7) <- clus_7_names
corn_by_SHR7 <- corn_by_SHR7[ ,order_lgnd_7]
corn_by_SHR7
corn_by_sagbi_unmod <- as.data.frame(t(tapply(corn_df$area_ac, corn_df$SAGBI_unmod_class, sum)))
corn_by_sagbi_unmod
round(corn_by_sagbi_unmod/sum(corn_by_sagbi_unmod), 2) #32% in Excellent to Moderately Good

wheat_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Wheat',])
dim(wheat_df) #29151 rows
length(unique(wheat_df$mukey)) #1259 unique mukeys
wheat_by_SHR7 <- as.data.frame(t(tapply(wheat_df$area_ac, wheat_df$cluster_7, sum)))
colnames(wheat_by_SHR7) <- clus_7_names
wheat_by_SHR7 <- wheat_by_SHR7[ ,order_lgnd_7]
wheat_by_SHR7
wheat_by_sagbi_unmod <- as.data.frame(t(tapply(wheat_df$area_ac, wheat_df$SAGBI_unmod_class, sum)))
wheat_by_sagbi_unmod
round(wheat_by_sagbi_unmod/sum(wheat_by_sagbi_unmod), 2) #31% in Excellent to Moderately Good
wheat_by_sagbi_mod <- as.data.frame(t(tapply(wheat_df$area_ac, wheat_df$SAGBI_mod_class, sum)))
wheat_by_sagbi_mod
round(wheat_by_sagbi_mod/sum(wheat_by_sagbi_mod), 2)

alfalfa_df <- as.data.frame(valley_crops_mu[valley_crops_mu$Crop2014=='Alfalfa and Alfalfa Mixtures',])
dim(alfalfa_df) #29151 rows
length(unique(alfalfa_df$mukey)) #1259 unique mukeys
alfalfa_by_SHR7 <- as.data.frame(t(tapply(alfalfa_df$area_ac, alfalfa_df$cluster_7, sum)))
colnames(alfalfa_by_SHR7) <- clus_7_names
alfalfa_by_SHR7 <- alfalfa_by_SHR7[ ,order_lgnd_7]
alfalfa_by_SHR7
alfalfa_by_sagbi_unmod <- as.data.frame(t(tapply(alfalfa_df$area_ac, alfalfa_df$SAGBI_unmod_class, sum)))
alfalfa_by_sagbi_unmod
round(alfalfa_by_sagbi_unmod/sum(alfalfa_by_sagbi_unmod), 2) #29% in Excellent to Moderately Good
alfalfa_by_sagbi_mod <- as.data.frame(t(tapply(alfalfa_df$area_ac, alfalfa_df$SAGBI_mod_class, sum)))
alfalfa_by_sagbi_mod
round(alfalfa_by_sagbi_mod/sum(alfalfa_by_sagbi_mod), 2)

sagbi_by_SHR7 <- tapply(valley_crops_mu$SAGBI_unmod, valley_crops_mu$cluster_7, mean, na.rm=TRUE)
names(sagbi_by_SHR7) <- clus_7_names
sagbi_by_SHR7 <- sagbi_by_SHR7[order_lgnd_7]

