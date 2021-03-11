library(raster)
dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final'
list.files(file.path(dataDir, 'FINAL results', 'shapefiles with data'))
# shr_order_shp <- shapefile(file.path('C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final/FINAL results/shapefiles with data', 'valley_30cm_cluster_SoilOrder.shp'))
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)
colnames(valley30cm_by_mukey)
dom_order_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', "soil survey facts", 'dom_order_by_mukey.csv'), stringsAsFactors = FALSE)
valley30cm_by_mukey$dom_order <- dom_order_by_mukey$dom_order[match(valley30cm_by_mukey$mukey, dom_order_by_mukey$mukey)]
colnames(valley30cm_by_mukey)
valley30cm_by_mukey <- valley30cm_by_mukey[,-c(10:19,26:34,47)]

# clus_3_names <- c('1. Loamy w/pans', '3. Shrink-swell', '2. Coarse salt-affected')
# clus_4_names <- c('1. Coarse w/no restrictions', '2. Loamy w/pans', '3. Loamy saline-sodic', '4. Fine shrink-swell')
# clus_5_names <- c('3. Loamy w/pans', '4. Loamy saline-sodic', '2. Loamy w/no restrictions', '5. Fine shrink-swell', '1. Coarse w/no restrictions')
clus_7_names <- c('6. Fine salt-affected', '3. Low OM with restrictive horizons', '4. High OM with restrictive horizons', '1. Coarse with no restrictions', '2. Loamy with no restrictions', '7. Shrink-swell', '5. Coarse-loamy salt-affected')
valley30cm_by_mukey$SHR7name <- clus_7_names[valley30cm_by_mukey$cluster_7]
tapply(valley30cm_by_mukey$clay_30cm, valley30cm_by_mukey$SHR7name, mean)
tapply(valley30cm_by_mukey$MnRs_dep, valley30cm_by_mukey$SHR7name, mean)
tapply(valley30cm_by_mukey$pH_30cm, valley30cm_by_mukey$SHR7name, mean)
# clus_9_names <- c('4. Coarse w/pans', '5. Loamy w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions',  '7. Coarse saline-sodic', '2. Coarse w/no restrictions', '8. Fine saline-sodic', '6. Loamy w/pans (high OM)', '9. Fine shrink-swell') #order corrected 4/7/20
# clus_8_names <- c('6. Coarse saline-sodic', '7. Fine saline-sodic', '5. Loamy w/pans', '2. Coarse w/no restrictions', '8. Fine shrink-swell', '4. Coarse w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions') #order corrected 4/7/20
# clus_6_names <- c('3. Coarse w/pans', '4. Loamy w/pans', '5. Loamy saline-sodic', '1. Coarse w/no restrictions', '6. Fine shrink-swell', '2. Loamy w/no restrictions') #order corrected 4/7/20
write.csv(valley30cm_by_mukey, file.path(dataDir, 'FINAL results', 'share', 'valley30cm_by_mukey_cluster_FINAL_to_share.csv'), row.names = FALSE)

#simple shapefile to share
valley_30cm_shp <- shapefile(file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster.shp'))
valley_30cm_shp$SHR7name <- clus_7_names[valley_30cm_shp$cluster_7]
area_by_shr <- tapply(valley_30cm_shp$area_ac, valley_30cm_shp$SHR7name, sum)
area_by_shr / 2.47105 #checks against sssaj table 2
names(valley_30cm_shp)
unique(valley_30cm_shp$arsymbl)
unique(valley_30cm_shp$spatlvr)
valley_simp_shp <- valley_30cm_shp[,c(1,5:6,8:9,19:20,56:63,65:66,70:ncol(valley_30cm_shp))]
names(valley_simp_shp)[9:21] <- c('clay', 'silt', 'sand', 'OM', 'CECpH7', 'BD', 'frags', 'EC', 'pH_H2O', 'LEP', 'ksat', 'SOC', 'AWC')
names(valley_simp_shp)
lapply(as.data.frame(valley_simp_shp)[,8:21], summary)
valley_simp_shp$dom_order <- dom_order_by_mukey$dom_order[match(valley_simp_shp$mukey, dom_order_by_mukey$mukey)]
area_by_order <- tapply(valley_simp_shp$area_ac, valley_simp_shp$dom_order, sum)
area_by_order / 2.47105 #checks against sssaj table 2
shapefile(valley_simp_shp, file.path(dataDir, 'FINAL results', 'share', 'CA_SHR7_30cm.shp'))
writeClipboard(names(valley_simp_shp))
