ecoDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/ecoregions'
ssurgoDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/ssurgo_data'
list.files(ssurgoDir)
library(aqp)
#demo(aqp)
#demo(slope_effect_hz_thickness)
library(raster)

#read in map unit (mu) tabular data
mu_data <- read.csv(file.path(ssurgoDir, 'ca_mapunit_data.csv'), stringsAsFactors = FALSE)
dim(mu_data) #18865 map units across CA
mu_data_Fresno <- mu_data[mu_data$areasymbol %in% c('ca651', 'ca653', 'ca654'),]
dim(mu_data_Fresno) #812 map units

#read in map unit spatial data
list.files(file.path(ssurgoDir, 'ca_mapunits'))
mu_shp <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits.shp'))
unique(mu_shp$areasymbol)
mu_shp_fresno <- mu_shp[mu_shp$areasymbol %in% c('ca651', 'ca653', 'ca654'), ]
length(unique(mu_shp_fresno$mukey)) #812 is the trick
shapefile(mu_shp_fresno, file.path(ssurgoDir, 'ca_mapunits', 'fresno_mapunits.shp'))
fresno_area <- aggregate(mu_shp_fresno)
shapefile(fresno_area, file.path(ssurgoDir, 'ca_mapunits', 'fresno_mu_area.shp'))
plot(fresno_area)
plot(mu_shp_fresno)

#read in component (comp) data
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE)
dim(comp_data) #91193
length(unique(comp_data$mukey)) #18,861 
length(unique(comp_data$cokey)) #91,193 unique components
comp_data_Fresno <- comp_data[comp_data$mukey %in% mu_data_Fresno$mukey,]
dim(comp_data_Fresno) #3563 rows
length(unique(comp_data_Fresno$mukey)) #812 map units match above
length(unique(comp_data_Fresno$cokey)) #3563 unique components

#read in horizon data
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
horizon_data_Fresno <- horizon_data[horizon_data$cokey %in% comp_data_Fresno$cokey,]
dim(horizon_data_Fresno) #3993 horizons
unique(horizon_data_Fresno$hzname)

#read in ecoregions
list.files(ecoDir)
eco_L3 <- shapefile(file.path(ecoDir, 'ca_eco_l3.shp'))
plot(eco_L3)
