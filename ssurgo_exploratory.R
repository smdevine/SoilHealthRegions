#left off at line 171
library(raster)
library(aqp)
#demo(aqp)
#demo(slope_effect_hz_thickness)
mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
ecoDir <- file.path(mainDir, 'soil health/ecoregions') #need to rename directory to postdoc to match laptop
ssurgoDir <- file.path(mainDir, 'soil health/ssurgo_data')
cropsDir <- file.path(mainDir, 'soil health/crops')
summaryDir <- file.path(mainDir, 'soil health/summaries/fresno_area_trial')
cropsCRS <- crs('+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
list.files(summaryDir)
list.files(ssurgoDir)

#define some functions
min_modified <- function(x) {
  if(all(is.na(x))) {
    return(NA)
  }
  else {min(x, na.rm = TRUE)}
}
max_modified <- function(x) {
  if(all(is.na(x))) {
    return(NA)
  }
  else {max(x, na.rm = TRUE)}
}
concat_names <- function(x, decat=FALSE) {
  if (decat) {x <- unlist(strsplit(x, split = '-'))}
  if(all(is.na(x))) {
    NA
  } else if (length(unique(x[!is.na(x)]))==1) {
        unique(x[!is.na(x)])
    } else {
        paste(unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))], collapse = '-')
    }
}


#read in map unit (mu) tabular data
mu_data <- read.csv(file.path(ssurgoDir, 'ca_mapunit_data.csv'), stringsAsFactors = FALSE)
#dim(mu_data) #18865 map units across CA
mu_data_Fresno <- mu_data[mu_data$areasymbol %in% c('ca651', 'ca653', 'ca654'),]
#dim(mu_data_Fresno) #812 map units
sum(grepl('association', mu_data_Fresno$muname)) #40 associations
sum(grepl('complex', mu_data_Fresno$muname)) #49 complexes


#read in map unit spatial data
list.files(file.path(ssurgoDir, 'ca_mapunits'))
#mu_shp <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits.shp'))
#unique(mu_shp$areasymbol)
#mu_shp_fresno <- mu_shp[mu_shp$areasymbol %in% c('ca651', 'ca653', 'ca654'), ]
#length(unique(mu_shp_fresno$mukey)) #812 is the trick
#shapefile(mu_shp_fresno, file.path(ssurgoDir, 'ca_mapunits', 'fresno_mapunits.shp'))
#fresno_area <- aggregate(mu_shp_fresno)
#shapefile(fresno_area, file.path(ssurgoDir, 'ca_mapunits', 'fresno_mu_area.shp'))
#plot(fresno_area)
#plot(mu_shp_fresno)
fresno_area <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'fresno_only', 'fresno_mu_area.shp'))
mu_shp_fresno <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'fresno_only', 'fresno_mapunits.shp'))
fresno_mu_aea <- spTransform(mu_shp_fresno, cropsCRS)
area651 <- mu_shp_fresno[mu_shp_fresno$areasymbol=='ca651',]
area651 <- aggregate(area651)
area651_aea <- spTransform(area651, cropsCRS)
shapefile(area651_aea, file.path(ssurgoDir, 'ca_mapunits', 'fresno_only', 'area651_aea.shp'))
area653 <- mu_shp_fresno[mu_shp_fresno$areasymbol=='ca653',]
area653 <- aggregate(area653)
area653_aea <- spTransform(area653, cropsCRS)
shapefile(area653_aea, file.path(ssurgoDir, 'ca_mapunits', 'fresno_only', 'area653_aea.shp'))
area654 <- mu_shp_fresno[mu_shp_fresno$areasymbol=='ca654',]
area654 <- aggregate(area654)
area654_aea <- spTransform(area654, cropsCRS)
shapefile(area654_aea, file.path(ssurgoDir, 'ca_mapunits', 'fresno_only', 'area654_aea.shp'))

#read in component (comp) data
list.files(file.path(ssurgoDir, 'component_data'))
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' ')) #treat blanks or a space as a NA
colnames(comp_data)
dim(comp_data) #91193
length(unique(comp_data$mukey)) #18,861 
length(unique(comp_data$cokey)) #91,193 unique cokeys
comp_data_Fresno <- comp_data[comp_data$mukey %in% mu_data_Fresno$mukey,]
dim(comp_data_Fresno) #3563 rows
length(unique(comp_data_Fresno$mukey)) #812 map units match above
length(unique(comp_data_Fresno$cokey)) #3563 unique components
unique(comp_data_Fresno$taxorder) #6 soil orders here
length(unique(comp_data_Fresno$compname)) #235 unique component names
unique(comp_data_Fresno$taxgrtgroup)
sum(!is.na(comp_data_Fresno$castorieindex))
sum(!is.na(comp_data_Fresno$castorieindex) & comp_data_Fresno$majcompflag=='Yes')
sum(!is.na(comp_data_Fresno$castorieindex) & comp_data_Fresno$majcompflag=='No ')
comp_data_Fresno[comp_data_Fresno$majcompflag=='No ' & !is.na(comp_data_Fresno$castorieindex), ]
summary(comp_data_Fresno$comppct_r[comp_data_Fresno$majcompflag=='Yes'])
sum(comp_data_Fresno$comppct_r[comp_data_Fresno$majcompflag=='Yes'] < 15) #2 instances of <15% comppct_r flagged as majcomps
comp_data_Fresno[comp_data_Fresno$majcompflag=='Yes' & comp_data_Fresno$comppct_r < 15,]
sum(comp_data_Fresno$comppct_r[comp_data_Fresno$majcompflag=='No '] >= 15) #82 are not flagged as majcomp but most don't have data
comp_data_Fresno[comp_data_Fresno$majcompflag=='No ' & comp_data_Fresno$comppct_r>=15,]
sum(comp_data_Fresno$majcompflag=='Yes' & comp_data_Fresno$comppct_r < 15)
sum(comp_data_Fresno$majcompflag=='No ' & comp_data_Fresno$comppct_r>=15 & !is.na(comp_data_Fresno$castorieindex)) #hay tres
#fix this
comp_data_Fresno$majcompflag[comp_data_Fresno$majcompflag=='No ' & comp_data_Fresno$comppct_r>=15 & !is.na(comp_data_Fresno$castorieindex)] <- 'Yes'
comp_data_Fresno$majcompflag[comp_data_Fresno$majcompflag=='Yes' & comp_data_Fresno$comppct_r < 15] <- 'No '

restrictions <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_restrictions.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
#restrictions$mukey <- comp_data$mukey[match(restrictions$cokey, comp_data$cokey)]
head(restrictions)
dim(restrictions) #some cokeys have more than one restrictive layer
length(unique(restrictions$cokey))
restrictions_Fresno <- restrictions[restrictions$cokey %in% comp_data_Fresno$cokey, ]
restrictions_Fresno$majcompflag <- comp_data_Fresno$majcompflag[match(restrictions_Fresno$cokey, comp_data_Fresno$cokey)]

parentmat <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_parentmat.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
parentmat_Fresno <- parentmat[parentmat$cokey %in% comp_data_Fresno$cokey, ]
head(parentmat_Fresno)
unique(parentmat_Fresno$pmkind) #7 different kinds
table(parentmat$pmkind)
parentmat_Fresno$majcompflag <- comp_data_Fresno$majcompflag[match(parentmat_Fresno$cokey, comp_data_Fresno$cokey)]
table(parentmat_Fresno$majcompflag) #lots of minor component data
dim(parentmat_Fresno)
length(unique(parentmat_Fresno$cokey)) #lots of cokeys with more than one parent material, thus...
parentmat_by_cokey <- data.frame(cokey = row.names(tapply(parentmat_Fresno$pmkind[parentmat_Fresno$majcompflag=='Yes'], parentmat_Fresno$cokey[parentmat_Fresno$majcompflag=='Yes'],  function(x) unique(x))), pmkinds = as.character(tapply(parentmat_Fresno$pmkind[parentmat_Fresno$majcompflag=='Yes'], parentmat_Fresno$cokey[parentmat_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
unique(parentmat_by_cokey$pmkinds)

productivity <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_forprod.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
head(productivity)
productivity_Fresno[productivity_Fresno$cokey==16607839,]
productivity_Fresno <- productivity[productivity$cokey %in% comp_data_Fresno$cokey, ]
dim(productivity_Fresno)#132
length(unique(productivity_Fresno$cokey)) #lots of duplicates
productivity_Fresno <- productivity_Fresno[!duplicated(productivity_Fresno$cokey),]
productivity_Fresno #very limited dataset

hillslope_pos <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_hillslopeprof.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
hillslope_pos_Fresno <- hillslope_pos[hillslope_pos$cokey %in% comp_data_Fresno$cokey, ]
dim(hillslope_pos_Fresno)
length(unique(hillslope_pos_Fresno$cokey))
unique(hillslope_pos_Fresno$hillslopeprof)
hillslope_pos_Fresno$majcompflag <- comp_data_Fresno$majcompflag[match(hillslope_pos_Fresno$cokey, comp_data_Fresno$cokey)]
hillslope_pos_Fresno$mukey <- comp_data_Fresno$mukey[match(hillslope_pos_Fresno$cokey, comp_data_Fresno$cokey)]
table(hillslope_pos_Fresno$majcompflag)
table(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='Yes']) / sum(hillslope_pos_Fresno$majcompflag=='Yes')
table(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='No ']) / sum(hillslope_pos_Fresno$majcompflag=='No ')

#this data seems problematic
#for major components
hillslope_pos_by_cokey <- data.frame(cokey = row.names(tapply(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='Yes'], hillslope_pos_Fresno$cokey[hillslope_pos_Fresno$majcompflag=='Yes'], function(x) unique(x))), hillpos= as.character(tapply(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='Yes'], hillslope_pos_Fresno$cokey[hillslope_pos_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
unique(hillslope_pos_by_cokey$hillpos)
hillslope_pos_by_mukey <- data.frame(mukey = row.names(tapply(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='Yes'], hillslope_pos_Fresno$mukey[hillslope_pos_Fresno$majcompflag=='Yes'], function(x) unique(x))), hillpos= as.character(tapply(hillslope_pos_Fresno$hillslopeprof[hillslope_pos_Fresno$majcompflag=='Yes'], hillslope_pos_Fresno$mukey[hillslope_pos_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
unique(hillslope_pos_by_mukey$hillpos)
table(hillslope_pos_by_mukey$hillpos)
sum(is.na(hillslope_pos_by_mukey$hillpos)) #only 1 NA
hillslope_pos_by_mukey$mukey[hillslope_pos_by_mukey$hillpos=='Footslope-Summit']
hillslope_pos_Fresno[hillslope_pos_Fresno$mukey=='464368',]
comp_data_Fresno[comp_data_Fresno$mukey=='464368', ]
mu_data_Fresno[mu_data_Fresno$mukey=='464368', ]
horizon_data_Fresno[horizon_data_Fresno$cokey=='16626076', ]

#make this from Fresno only dataset, even though I don't keep tagging it with Fresno
reskinds_by_cokey <- data.frame(cokey = row.names(tapply(restrictions_Fresno$reskind[restrictions_Fresno$majcompflag=='Yes'], restrictions_Fresno$cokey[restrictions_Fresno$majcompflag=='Yes'],  function(x) unique(x))), reskinds = as.character(tapply(restrictions_Fresno$reskind[restrictions_Fresno$majcompflag=='Yes'], restrictions_Fresno$cokey[restrictions_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE) #function(x) unique(x)))) #function(x) if(all(is.na(x))) {NA} else if (length(unique(x[!is.na(x)]))==1) {unique(x[!is.na(x)])} else{paste(unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))], collapse = '-')})
dim(reskinds_by_cokey) #658 majcomp cokeys with restrictions
length(unique(reskinds_by_cokey$cokey))
lapply(reskinds_by_cokey, class)
unique(reskinds_by_cokey$reskinds) #17 different reskind combos for Fresno
reskinds_by_cokey$mukey <- comp_data_Fresno$mukey[match(reskinds_by_cokey$cokey, comp_data_Fresno$cokey)]
dim(reskinds_by_cokey)
length(unique(reskinds_by_cokey$mukey)) #469 mukeys
reskinds_by_cokey$comp_pct <- comp_data_Fresno$comppct_r[match(reskinds_by_cokey$cokey, comp_data_Fresno$cokey)]
reskinds_by_cokey$majcompflag <- comp_data_Fresno$majcompflag[match(reskinds_by_cokey$cokey, comp_data_Fresno$cokey)]
summary(as.factor(reskinds_by_cokey$majcompflag)) #minor components were left out in the original creation of the table above
summary(as.factor(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) length(x)))) #up to 3, makes sense given that there are at most 3 major components in a map-unit for this AOI
#reskinds needs to be unconcatenated first before finding unique
reskinds_by_mukey <- data.frame(mukey = row.names(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) unique(x))), reskinds = as.character(tapply(reskinds_by_cokey$reskind, reskinds_by_cokey$mukey, concat_names, decat=TRUE)), stringsAsFactors = FALSE)
unique(reskinds_by_mukey$reskinds) #now this appears to be ok, using the decat=TRUE

majcomps_no_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$majcompflag, comp_data_Fresno$mukey, function(x) sum(x=='Yes'))), majcomp_no = as.numeric(tapply(comp_data_Fresno$majcompflag, comp_data_Fresno$mukey, function(x) sum(x=='Yes'))), stringsAsFactors = FALSE)
unique(majcomps_no_by_mukey$majcomp_no)

majcompnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$compname[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], concat_names)), majcompnames=as.character(tapply(comp_data_Fresno$compname[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
unique(majcompnames_by_mukey$majcompnames) #277 unique major component name combos

compnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$compname, comp_data_Fresno$mukey, function(x) unique(x))), compnames = as.character(tapply(comp_data_Fresno$compname, comp_data_Fresno$mukey, concat_names)), stringsAsFactors = FALSE)
unique(compnames_by_mukey$compnames) #422 different map-unit combos accounting for both major and minor components; decat=TRUE has no effect
dim(compnames_by_mukey) #but 812 mukeys

majcokeys_by_mukey <- tapply(comp_data_Fresno$cokey[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], function(x) x)

majcomp_taxorders_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$taxorder[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], concat_names)), taxorders=as.character(tapply(comp_data_Fresno$taxorder[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)

domcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$comppct_r, comp_data_Fresno$mukey, function(x) max(x, na.rm=TRUE))), docomppct=as.numeric(tapply(comp_data_Fresno$comppct_r, comp_data_Fresno$mukey, function(x) max(x, na.rm=TRUE))), stringsAsFactors = FALSE)
summary(domcomp_pct_by_mukey$docomppct)

storie_rng_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$castorieindex[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], function(x) max_modified(x)- min_modified(x))), storierng=as.numeric(tapply(comp_data_Fresno$castorieindex[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], function(x) max_modified(x)- min_modified(x))), stringsAsFactors=FALSE)
hist(storie_rng_by_mukey$storierng)

storie_mean_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_Fresno$castorieindex[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], function(x) mean(x, na.rm=TRUE))), storiemn=as.numeric(tapply(comp_data_Fresno$castorieindex[comp_data_Fresno$majcompflag=='Yes'], comp_data_Fresno$mukey[comp_data_Fresno$majcompflag=='Yes'], function(x) mean(x, na.rm=TRUE))), stringsAsFactors=FALSE)
hist(storie_mean_by_mukey$storiemn)
sum(is.na(storie_mean_by_mukey$storiemn)) #30 are NA

#read in horizon data
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
horizon_data_Fresno <- horizon_data[horizon_data$cokey %in% comp_data_Fresno$cokey,]
dim(horizon_data_Fresno) #3993 horizons
unique(horizon_data_Fresno$hzname)
horizon_data_Fresno$majcompflag <- comp_data_Fresno$majcompflag[match(horizon_data_Fresno$cokey, comp_data_Fresno$cokey)]
table(horizon_data_Fresno$majcompflag) #493 minor component have horizon data
horizon_data_Fresno[horizon_data_Fresno$majcompflag=='No ', ]

#read in ecoregions
list.files(ecoDir)
eco_L3 <- shapefile(file.path(ecoDir, 'ca_eco_l3.shp'))
plot(eco_L3)
eco_L3_wgs84 <- spTransform(eco_L3, crs(fresno_area))
eco_L3_Fresno <- crop(eco_L3_wgs84, fresno_area)
plot(eco_L3_Fresno)
eco_L3_Fresno$US_L3NAME

eco_L4 <- shapefile(file.path(ecoDir, 'ca_eco_l4.shp'))
eco_L4_wgs84 <- spTransform(eco_L4, crs(fresno_area))
eco_L4_Fresno <- crop(eco_L4_wgs84, fresno_area)
plot(eco_L4_Fresno)
unique(eco_L4_Fresno$US_L4NAME)

#read in crops
list.files(cropsDir)
crops <- shapefile(file.path(cropsDir, 'i15_Crop_Mapping_2014_Final_LandIQonAtlas.shp'))
fresno_area_aea <- spTransform(fresno_area, crs(crops))
shapefile(fresno_area_aea, file.path(ssurgoDir, 'ca_mapunits/fresno_only/ca651_653_654_aea.shp'))
area(fresno_area_aea) / 10000 * 2.47105 #3,303,267 acres

#add some variable to this subset of map units
length(unique(fresno_mu_aea$mukey)) #812 mukeys
fresno_mu_aea$area_ac <- area(fresno_mu_aea) / 10000 * 2.47105 #acres calc
sum(fresno_mu_aea$area_ac) #3303267 matches above
fresno_mu_aea$majcomps_no <- majcomps_no_by_mukey$majcomp_no[match(fresno_mu_aea$mukey, majcomps_no_by_mukey$mukey)]
round(tapply(fresno_mu_aea$area_ac, fresno_mu_aea$majcomps_no, sum), 0)
#0       1        2        3 
#1742 2,223,661  668,190  409,675  (67.5% of survey area has 1 major component per map unit, 20.0% has 2 major components, 12.4% has 3 major components)
fresno_mu_aea$taxorders <- majcomp_taxorders_by_mukey$taxorders[match(fresno_mu_aea$mukey, majcomp_taxorders_by_mukey$mukey)]
unique(fresno_mu_aea$taxorders)
taxorders_area <- data.frame(taxorders=row.names(tapply(fresno_mu_aea$area_ac, fresno_mu_aea$taxorders, sum)), acres=as.numeric(tapply(fresno_mu_aea$area_ac, fresno_mu_aea$taxorders, sum)), stringsAsFactors = FALSE)
sum(taxorders_area$acres) #3259727 less than above because NA class is dropped (see next 2 calcs)
sum(fresno_mu_aea$area_ac) - sum(taxorders_area$acres) #43540
sum(fresno_mu_aea$area_ac[is.na(fresno_mu_aea$taxorders)]) #43540.03
taxorders_area[nrow(taxorders_area)+1, 'taxorders'] <- 'Undefined'
taxorders_area[nrow(taxorders_area), 'acres'] <- sum(fresno_mu_aea$area_ac[is.na(fresno_mu_aea$taxorders)]) 
taxorders_area <- taxorders_area[order(taxorders_area$acres, decreasing=TRUE), ]
taxorders_area
#write.csv(taxorders_area, file.path(summaryDir, 'taxorders_area.csv'), row.names=FALSE)

fresno_mu_aea$domcomp_pct <- domcomp_pct_by_mukey$docomppct[match(fresno_mu_aea$mukey, majcomps_no_by_mukey$mukey)]

#add more naming information
fresno_mu_aea$muname <- mu_data_Fresno$muname[match(fresno_mu_aea$mukey, mu_data_Fresno$mukey)]
fresno_mu_aea$complex <- ifelse(grepl('complex', fresno_mu_aea$muname), 'Yes', 'No') #374 yes
fresno_mu_aea$association <- ifelse(grepl('association', fresno_mu_aea$muname), 'Yes', 'No') #272 yes

#add Storie index range and mean (relative to major components only)
fresno_mu_aea$storiemn <- storie_mean_by_mukey$storiemn[match(fresno_mu_aea$mukey, storie_mean_by_mukey$mukey)]
fresno_mu_aea$storierng <- storie_rng_by_mukey$storierng[match(fresno_mu_aea$mukey, storie_mean_by_mukey$mukey)]


#all compnames by mukey 
#sum(grepl('-', unique(comp_data_Fresno$compname)))
#this counts rock outcrop when it's a minor component
fresno_mu_aea$Rock_OC <- ifelse(grepl('Rock outcrop', compnames_by_mukey$compnames[match(fresno_mu_aea$mukey, compnames_by_mukey$mukey)]), 'Yes', 'No')
summary(as.factor(fresno_mu_aea$Rock_OC)) #1953 polygons have rock outcrop


fresno_mu_aea$Lithic <- ifelse(grepl('Lithic bedrock', reskinds_by_mukey$reskinds[match(fresno_mu_aea$mukey, reskinds_by_mukey$mukey)]) | fresno_mu_aea$Rock_OC=='Yes', 'Yes', 'No') #2247 were 'yes' after accounting for mukeys with more than one cokey with restrictions and those with rock outrcrop but no major components with a subsoil lithic contact
summary(as.factor(fresno_mu_aea$Lithic))

fresno_mu_aea$Paralithic <- ifelse(grepl('Paralithic bedrock', reskinds_by_mukey$reskinds[match(fresno_mu_aea$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
summary(as.factor(fresno_mu_aea$Paralithic)) #5283 are yes

#what about rock outcrop? this below needs to be refined to account for percentage of rock outcrop also
fresno_lithic_comppct <- data.frame(mukey=row.names(tapply(comp_data_Fresno$comppct_r[comp_data_Fresno$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)]], comp_data_Fresno$mukey[comp_data_Fresno$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)]], sum)), compct_sum = as.numeric(tapply(comp_data_Fresno$comppct_r[comp_data_Fresno$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)]], comp_data_Fresno$mukey[comp_data_Fresno$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)]], sum)), stringsAsFactors = FALSE)
summary(fresno_lithic_comppct$compct_sum)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#19.00   20.00   30.00   44.85   85.00   95.00
sum(fresno_lithic_comppct$compct_sum==95) #2 are 95
fresno_lithic_comppct[fresno_lithic_comppct$compct_sum==95,]
comp_data_Fresno[comp_data_Fresno$mukey==464390,]
reskinds_by_cokey[reskinds_by_cokey$mukey==464390,]
rock_OC_cokeys <- comp_data_Fresno$cokey[comp_data_Fresno$compname=='Rock outcrop']
sum(rock_OC_cokeys %in% reskinds_by_cokey$cokey) #51 rock outcrop cokeys are in reskind table
reskinds_by_cokey[reskinds_by_cokey$cokey %in% rock_OC_cokeys, ] #one mistakenly has 'paralithic bedrock' as reskind
comp_data_Fresno[comp_data_Fresno$cokey==16607858,]

sum(grepl('Rock outcrop', compnames_by_mukey$compnames)) #143 here
sum(grepl('Rock outcrop', majcompnames_by_mukey$majcompnames)) #68 here
fresno_rockOC_comppct <- data.frame(mukey=row.names(tapply(comp_data_Fresno$comppct_r[comp_data_Fresno$compname=='Rock outcrop'], comp_data_Fresno$mukey[comp_data_Fresno$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), compct_sum = as.numeric(tapply(comp_data_Fresno$comppct_r[comp_data_Fresno$compname=='Rock outcrop'], comp_data_Fresno$mukey[comp_data_Fresno$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), stringsAsFactors = FALSE)
head(fresno_rockOC_comppct)
comp_data_Fresno[comp_data_Fresno$mukey==463312, ]
summary(fresno_rockOC_comppct$compct_sum)
test <- fresno_lithic_comppct$compct_sum[match(fresno_mu_aea$mukey[fresno_mu_aea$Lithic=='Yes'], fresno_lithic_comppct$mukey)]
test2 <- fresno_rockOC_comppct$compct_sum[match(fresno_mu_aea$mukey[fresno_mu_aea$Lithic=='Yes'], fresno_lithic_comppct$mukey)]
length(fresno_mu_aea$Lthc_pct[fresno_mu_aea$Lithic=='Yes'])
length(test)
length(test2)
fresno_mu_aea$Lthc_pct <- 0
fresno_mu_aea$Lthc_pct[fresno_mu_aea$Lithic=='Yes'] <- apply(cbind(fresno_lithic_comppct$compct_sum[match(fresno_mu_aea$mukey[fresno_mu_aea$Lithic=='Yes'], fresno_lithic_comppct$mukey)], fresno_rockOC_comppct$compct_sum[match(fresno_mu_aea$mukey[fresno_mu_aea$Lithic=='Yes'], fresno_lithic_comppct$mukey)]), 1, sum, na.rm = TRUE)
summary(fresno_mu_aea$Lthc_pct)
sum(fresno_mu_aea$Lthc_pct > 100)


#what about mukeys with more than one major component?
fresno_mu_aea$Lthc_cm[fresno_mu_aea$Lithic=='Yes'] <- restrictions_Fresno$resdept_r[restrictions_Fresno$reskind=="Lithic bedrock"][match(fresno_mu_aea$mukey[fresno_mu_aea$Lithic=='Yes'], restrictions_Fresno$mukey[restrictions_Fresno$reskind=="Lithic bedrock"])]


sum(grepl('\\bL|lithic\\b', reskinds_by_cokey$reskinds)) #569 in Fresno only
#write to file
shapefile(fresno_mu_aea, file.path(ssurgoDir, 'ca_mapunits/fresno_only/ca651_653_654mu_aea.shp'), overwrite=TRUE)

#check it!
fresno_mu_aea[fresno_mu_aea$mukey==467090, ]
comp_data_Fresno[comp_data_Fresno$mukey==467090,]

fresno_mu_aea[337,]





crops_fresno <- crop(crops, fresno_area_aea)
shapefile(crops_fresno, file.path(cropsDir, 'fresno_only', 'crops_fresno.shp'))
unique(crops_fresno$Crop2014)
sum(crops_fresno$Acres) #1,718.419 acres, so 52% of total area
crop_acreage_fresno <- data.frame(crop=row.names(tapply(crops_fresno$Acres, crops_fresno$Crop2014, function(x) sum(x))), acres=as.numeric(tapply(crops_fresno$Acres, crops_fresno$Crop2014, function(x) sum(x))))
crop_acreage_fresno <- crop_acreage_fresno[order(crop_acreage_fresno$acres, decreasing = TRUE), ]
write.csv(crop_acreage_fresno, file.path(summaryDir, 'crop_acreage_DWR14.csv'), row.names = FALSE)


#some test code
paste(unlist(majcomp_taxorders_by_mukey[which(row.names(majcomp_taxorders_by_mukey)=='467090')], use.names = FALSE), collapse = '-')
majcomp_taxorders_by_mukey

test[337]
