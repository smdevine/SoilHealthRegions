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
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
  dissertationDir <- 'D:/Dissertation/Allowable_Depletion/dissertation/v2.results'
}
mar_settings <- c(4, 4.5, 1, 1)
om_to_oc <- 1.72
crit_pH <- 7.8
#produced in ssurgo_calag_cluster_FINAL.R
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')
valley30cm_by_mukey$SHR7name <- clus_7_names[valley30cm_by_mukey$cluster_7]
#produced for chp 1 of dissertation
list.files(dissertationDir)
model_shp_0.5mAD30 <- shapefile(file.path(dissertationDir, 'shapefiles', 'results_0.5mAD30.shp'))
model_shp_0.5mAD30$mukey <- as.integer(model_shp_0.5mAD30$mukey)
model_shp_1.0mAD50 <- shapefile(file.path(dissertationDir, 'shapefiles', 'results_1.0mAD50.shp'))
model_shp_1.0mAD50$mukey <- as.integer(model_shp_1.0mAD50$mukey)
model_shp_2.0mAD50 <- shapefile(file.path(dissertationDir, 'shapefiles', 'results_2.0mAD50.shp'))
model_shp_2.0mAD50$mukey <- as.integer(model_shp_2.0mAD50$mukey)

length(unique(model_shp_0.5mAD30$mukey)) #4840
sum(unique(model_shp_0.5mAD30$mukey) %in% valley30cm_by_mukey$mukey)

model_shp_0.5mAD30$SHR7name <- valley30cm_by_mukey$SHR7name[match(model_shp_0.5mAD30$mukey, valley30cm_by_mukey$mukey)]
model_shp_1.0mAD50$SHR7name <- valley30cm_by_mukey$SHR7name[match(model_shp_1.0mAD50$mukey, valley30cm_by_mukey$mukey)]
model_shp_2.0mAD50$SHR7name <- valley30cm_by_mukey$SHR7name[match(model_shp_2.0mAD50$mukey, valley30cm_by_mukey$mukey)]

tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='almond.mature'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='almond.mature'], sum)
sum(tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='almond.mature'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='almond.mature'], sum)) #1122779 acres total of almonds in shr study
tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='walnut.mature'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='walnut.mature'], sum)
sum(tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='walnut.mature'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='walnut.mature'], sum))
tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='alfalfa.CV'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='alfalfa.CV'], sum)
sum(tapply(model_shp_0.5mAD30$Acres[model_shp_0.5mAD30$crpnm=='alfalfa.CV'], model_shp_0.5mAD30$SHR7name[model_shp_0.5mAD30$crpnm=='alfalfa.CV'], sum))

GWresults_bySHR7 <- function(shp, gw_scenario) {
  alfalfa <- shp[which(shp$crpnm=='alfalfa.CV' & !is.na(shp$SHR7name)),]
  almonds <- shp[which(shp$crpnm=='almond.mature' & !is.na(shp$SHR7name)),]
  walnuts <- shp[which(shp$crpnm=='walnut.mature' & !is.na(shp$SHR7name)),]
  shapefile(alfalfa, filename=file.path(dataDir, 'GW shapefiles', 'alfalfa', paste0('alfalfa', gw_scenario, '.shp')))
  shapefile(almonds, filename=file.path(dataDir, 'GW shapefiles', 'almonds', paste0('almonds', gw_scenario, '.shp')))
  shapefile(walnuts, filename=file.path(dataDir, 'GW shapefiles', 'walnuts', paste0('walnuts', gw_scenario, '.shp')))
}
GWresults_bySHR7(model_shp_0.5mAD30, '0.5mAD30')
GWresults_bySHR7(model_shp_1.0mAD50, '1.0mAD50')
GWresults_bySHR7(model_shp_2.0mAD50, '2.0mAD50')

almonds_2.0mAD50_SHR7 <- shapefile()