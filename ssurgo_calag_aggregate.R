#currently, only 0-30 cm aggregation done, but 0-10 and 0-100 cm are available if needed
#TO-DOs
#caculate minimum reskind type [DONE]
#check misc_res_by_cokey if applying to new datasets
#check effect of excluding minor components on AWC calc
#check how resdep is calculated in multi-component mukeys when some components have restrictions and some don't
#add tabulation of majcomps with data
library(raster)
library(aqp)
# library(extrafont)
# library(extrafontdb)
#font_import() only needs to be done once
# loadfonts()
#demo(aqp)
#demo(slope_effect_hz_thickness)
laptop <- TRUE
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
ecoDir <- file.path(mainDir, 'soil health/ecoregions') #need to rename directory to postdoc to match laptop
ssurgoDir <- file.path(mainDir, 'soil health/ssurgo_data')
cropsDir <- file.path(mainDir, 'soil health/crops')
summaryDir <- file.path(mainDir, 'soil health/summaries/valley_final')
# FiguresDir <- file.path(mainDir, 'soil health/Figures/ag_land')
cropsCRS <- crs('+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
# res_plots <- 150
CAcountiesDir <- 'C:/Users/smdevine/Desktop/SpatialData/CA_counties/government_units'
list.files(summaryDir)
list.files(ssurgoDir)

#define an assumption for classification exercise
assumed_depth <- 200

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
#this used to count up component percentages by reskind
reskind_comppct <- function(reskind, comp_df, reskind_df) {
  data.frame(mukey=row.names(tapply(comp_df$comppct_r[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], comp_df$mukey[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], sum)), compct_sum = as.numeric(tapply(comp_df$comppct_r[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], comp_df$mukey[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], sum)), stringsAsFactors = FALSE)
}

#functions to work with ProfileApply
wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=TRUE)
  m
}

kgOrgC_sum <- function(x, slice_it=FALSE, depth, rm.NAs=TRUE, om_to_c=1.72) {
  if (slice_it) {
    x <- horizons(x)[1:depth, ]
    depths(x) <- cokey ~ hzdept_r + hzdepb_r
  }
  thick <- x$hzdepb_r - x$hzdept_r
  sum((thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100), na.rm = rm.NAs)
}

awc_sum <- function(x, rm.NAs=TRUE) {
  thick <- x$hzdepb_r - x$hzdept_r
  sum(thick * x$awc_r, na.rm = rm.NAs)
}
textural.class.calc <- function(sand, silt, clay, QC_param=1) {
  ifelse(is.na(sand) | is.na(silt) | is.na(clay), NA,
    ifelse(sand + silt + clay > (100 + QC_param) | sand + silt + clay < (100 - QC_param), paste('proportions do not sum to 100+-', QC_param), 
      ifelse(silt + 1.5 * clay < 15, 'sand',
        ifelse(silt + 1.5 * clay >= 15 & silt + 2 * clay < 30, 'loamy sand',
          ifelse((clay >= 7 & clay < 20 & sand > 52 & silt + 2 * clay >= 30) | (clay < 7 & silt < 50 & silt + 2 * clay >= 30), 'sandy loam',
            ifelse(clay >= 7 & clay < 27 & silt >=28 & silt < 50 & sand <= 52, 'loam',
              ifelse((silt >= 50 & clay >= 12 & clay < 27) | (silt >=50 & silt < 80 & clay < 12), 'silt loam',
                ifelse(silt >= 80 & clay < 12, 'silt',
                  ifelse(clay >= 20 & clay < 35 & silt < 28 & sand > 45, 'sandy clay loam',
                    ifelse(clay >= 27 & clay < 40 & sand > 20 & sand <= 45, 'clay loam',
                      ifelse(clay >= 27 & clay < 40 & sand <= 20, 'silty clay loam',
                        ifelse(clay >= 35 & sand > 45, 'sandy clay',
                          ifelse(clay >= 40 & silt >= 40, 'silty clay',
                            ifelse(clay >= 40 & sand <= 45 & silt < 40, 'clay','undefined textural class'))))))))))))))
}

horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = c('claytotal_r', 'silttotal_r', 'sandtotal_r', 'om_r', 'cec7_r', 'dbthirdbar_r', 'fragvol_r_sum', 'kwfact', 'ec_r', 'ph1to1h2o_r', 'sar_r', 'caco3_r', 'gypsum_r', 'lep_r', 'ksat_r'), varnames = c('clay', 'silt', 'sand', 'om', 'cec', 'bd', 'frags', 'kwf', 'ec', 'pH', 'sar', 'caco3', 'gyp', 'lep', 'ksat')) { #lep is linear extensibility
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
  s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum)
  columnames <- c(columnames, paste0('kgOrg.m2_', depth, 'cm'), paste0('awc_', depth, 'cm')) 
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  s$compname <- comp_df$compname[match(s$cokey, comp_df$cokey)]
  s$mukey <- comp_df$mukey[match(s$cokey, comp_df$cokey)]
  s$comppct <- comp_df$comppct_r[match(s$cokey, comp_df$cokey)]
  s <- s[,c('mukey', 'cokey', 'compname', 'comppct', columnames)]
  s
}
MUaggregate <- function(df1, varname) {
  sapply(split(x=df1, f=df1$mukey), FUN=function(x) {if(sum(!is.na(x[[varname]]))==0) {NA} 
    else{sum(x$comppct[!is.na(x[[varname]])] * x[[varname]][!is.na(x[[varname]])] / sum(x$comppct[!is.na(x[[varname]])]))}
  })
}
MUAggregate_wrapper <- function(df1, varnames) {
  x <- sapply(varnames, FUN=MUaggregate, df1=df1)
  as.data.frame(cbind(mukey=as.integer(row.names(x)), x))
}
#read in CA shapefile
# CA_counties <- shapefile(file.path(CAcountiesDir, 'county_nrcs_a_ca.shp'))
# CA <- aggregate(CA_counties)
# CA_aea <- spTransform(CA, cropsCRS)

#read in map unit (mu) tabular data
mu_data <- read.csv(file.path(ssurgoDir, 'ca_mapunit_data.csv'), stringsAsFactors = FALSE)
# dim(mu_data) #18865 map units across CA
# sum(grepl('association', mu_data$muname)) #2051 associations
# sum(grepl('complex', mu_data$muname)) #5526 complexes

#read in map unit spatial data
# list.files(file.path(ssurgoDir, 'ca_mapunits'))
mu_shp <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits.shp'))
# unique(mu_shp$areasymbol)

#read in MLRA shapefile
# list.files(file.path(ssurgoDir, 'mlra'))
# mlra_shp <- shapefile(file.path(ssurgoDir, 'mlra', 'mlra_a_ca.shp')) #59 features
# plot(mlra_shp)
# valley_mlra_shp <- mlra_shp[mlra_shp$MLRA_NAME=='Sacramento and San Joaquin Valleys' | mlra_shp$MLRA_NAME=='Central California Coastal Valleys', ]
# mlra_aea_shp <- spTransform(mlra_shp, cropsCRS)
# shapefile(mlra_aea_shp, file.path(ssurgoDir, 'mlra', 'mlra_ca_aea.shp'))
# mlra_aea_shp <- shapefile(file.path(ssurgoDir, 'mlra', 'mlra_ca_aea.shp')) #Central California Coastal Valleys, Sacramento and San Joaquin Valleys, 
# plot(mlra_aea_shp)
# unique(mlra_aea_shp$MLRA_NAME)

#plot area of interest via MLRA shapefile
# tiff(file = file.path(FiguresDir, 'ca_mlra.tif'), pointsize = 10, family = 'Times New Roman', width = 4.5, height = 4.5, units = 'in', res = res_plots)
# par(mar=c(0.05, 0.05, 0.05, 0.05))
# plot(CA_aea)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Sacramento and San Joaquin Valleys', ], col='green', add=TRUE)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Central California Coastal Valleys',], col='green3', add=TRUE)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Central California Coast Range',], col='gold4', add=TRUE)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Sierra Nevada Foothills',], col='gold', add=TRUE)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Coastal Redwood Belt',], col='indianred4', add=TRUE)
# plot(mlra_aea_shp[mlra_aea_shp$MLRA_NAME=='Siskiyou-Trinity Area',], col='cadetblue', add=TRUE)
# legend('topright', legend=c('Sacramento and San Joaquin Valleys', 'Central California Coastal Valleys', 'Central California Coast Range', 'Sierra Nevada Foothills', 'Coastal Redwood Belt', 'Siskiyou-Trinity Area'), pch=22, col='black', pt.bg=c('green', 'green3', 'gold4', 'gold', 'indianred4', 'cadetblue'))
# dev.off()

#crop map unit shapfiles to mlra valley extents
# valley_mu_shp <- crop(mu_shp, valley_mlra_shp)
# plot(valley_mu_shp)
# sum(area(valley_mu_shp)) / 10000 * 2.47105 # 13,873,110 acres
# length(unique(valley_mu_shp$mukey)) #5,043 map-unit kinds
# shapefile(valley_mu_shp, file.path(ssurgoDir, 'ca_mapunits', 'valley_only', 'valley_mapunits.shp'))
valley_mu_shp <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'valley_only', 'valley_mapunits.shp'))
# mukeys_valley <- as.integer(unique(valley_mu_shp$mukey)) #5043
# write.csv(mukeys_valley, file.path(summaryDir, 'mukeys_unique_valley.csv'), row.names = FALSE)
mukeys_valley <- read.csv(file.path(summaryDir, 'mukeys_unique_valley.csv'), stringsAsFactors=FALSE)
mu_data_valley <- mu_data[mu_data$mukey %in% mukeys_valley$x, ]

#read in component (comp) data
#valley exploratory version replaces comp_data_Fresno with comp_data_valley
list.files(file.path(ssurgoDir, 'component_data'))
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' ')) #treat blanks or a space as a NA
colnames(comp_data)
# dim(comp_data) #91193
# length(unique(comp_data$mukey)) #18,861 
# length(unique(comp_data$cokey)) #91,193 unique cokeys
comp_data_valley <- comp_data[comp_data$mukey %in% mu_data_valley$mukey,]
if(sum(is.na(comp_data_valley$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data_valley$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(is.na(comp_data_valley$comppct_r))
comp_data_valley <- comp_data_valley[!is.na(comp_data_valley$comppct_r),]
dim(comp_data_valley) #23493 rows
length(unique(comp_data_valley$mukey)) #5043 map units match above
length(unique(comp_data_valley$cokey)) #23493 unique components
unique(comp_data_valley$taxorder) #9 soil orders here and NA
length(unique(comp_data_valley$compname)) #1153 unique component names
unique(comp_data_valley$taxgrtgroup)
sum(!is.na(comp_data_valley$castorieindex))
sum(!is.na(comp_data_valley$castorieindex) & comp_data_valley$majcompflag=='Yes')
sum(!is.na(comp_data_valley$castorieindex) & comp_data_valley$majcompflag=='No ')
comp_data_valley[comp_data_valley$majcompflag=='No ' & !is.na(comp_data_valley$castorieindex), ]
summary(comp_data_valley$comppct_r[comp_data_valley$majcompflag=='Yes'])
sum(comp_data_valley$comppct_r[comp_data_valley$majcompflag=='Yes'] < 15) #10 instance of <15% comppct_r flagged as majcomps
comp_data_valley[comp_data_valley$majcompflag=='Yes' & comp_data_valley$comppct_r < 15,]
sum(comp_data_valley$comppct_r[comp_data_valley$majcompflag=='No '] >= 15) #110 not flagged as majcomp in these valleys with > 15%
comp_data_valley[comp_data_valley$majcompflag=='No ' & comp_data_valley$comppct_r>=15, ]
sum(comp_data_valley$majcompflag=='No ' & comp_data_valley$comppct_r>=15 & !is.na(comp_data_valley$castorieindex)) #7

#fix a few majcompflag errors
comp_data_valley$majcompflag[comp_data_valley$majcompflag=='No ' & comp_data_valley$comppct_r>=15 & !is.na(comp_data_valley$castorieindex)] <- 'Yes'
comp_data_valley$majcompflag[comp_data_valley$majcompflag=='Yes' & comp_data_valley$comppct_r < 15] <- 'No '

restrictions <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_restrictions.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
# head(restrictions)
# dim(restrictions) #some cokeys have more than one restrictive layer
# length(unique(restrictions$cokey))
restrictions_valley <- restrictions[restrictions$cokey %in% comp_data_valley$cokey, ]
restrictions_valley$majcompflag <- comp_data_valley$majcompflag[match(restrictions_valley$cokey, comp_data_valley$cokey)]
# unique(restrictions_valley$reskind)
restrictions_valley$mukey <- comp_data_valley$mukey[match(restrictions_valley$cokey, comp_data_valley$cokey)]
sum(duplicated(restrictions_valley$cokey)) #368
restrictions_valley$comppct <- comp_data_valley$comppct_r[match(restrictions_valley$cokey, comp_data_valley$cokey)]

#depths to majcomp restrictions by reskind type
#which is necessary given that reskind has a NA
unique(restrictions_valley$reskind)
table(restrictions_valley$reskind[restrictions_valley$majcompflag=='Yes'])
lithic_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Lithic bedrock' & restrictions_valley$majcompflag=='Yes'), ]
head(lithic_by_cokey)
# dim(lithic_by_cokey)
# sum(duplicated(lithic_by_cokey$cokey))
lithic_by_mukey <- MUAggregate_wrapper(df1 = lithic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(lithic_by_mukey)

paralithic_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Paralithic bedrock' & restrictions_valley$majcompflag=='Yes'), ]
# dim(paralithic_by_cokey)
# sum(duplicated(paralithic_by_cokey$cokey))
paralithic_by_mukey <- MUAggregate_wrapper(df1 = paralithic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(paralithic_by_mukey)

duripan_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Duripan' & restrictions_valley$majcompflag=='Yes'), ]
# dim(duripan_by_cokey)
# sum(duplicated(duripan_by_cokey$cokey))
duripan_by_mukey <- MUAggregate_wrapper(df1 = duripan_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(duripan_by_mukey)

ATC_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Abrupt textural change' & restrictions_valley$majcompflag=='Yes'), ]
# dim(ATC_by_cokey)
# sum(duplicated(ATC_by_cokey$cokey))
ATC_by_mukey <- MUAggregate_wrapper(df1 = ATC_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(ATC_by_mukey)

natric_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Natric' & restrictions_valley$majcompflag=='Yes'), ]
# dim(natric_by_cokey)
# sum(duplicated(natric_by_cokey$cokey))
natric_by_mukey <- MUAggregate_wrapper(df1 = natric_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
natric_by_mukey

salic_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Salic' & restrictions_valley$majcompflag=='Yes'), ]
salic_by_mukey <- MUAggregate_wrapper(df1 = salic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
salic_by_mukey

SCTS_by_cokey <- restrictions_valley[which(restrictions_valley$reskind=='Strongly contrasting textural stratification' & restrictions_valley$majcompflag=='Yes'), ]
SCTS_by_mukey <- MUAggregate_wrapper(df1 = SCTS_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
SCTS_by_mukey

#now needs to be fixed so that if there's more than one restriction for a given cokey, the minimum is chosen before averaging with other cokeys
misc_res_by_cokey <- restrictions_valley[which(restrictions_valley$reskind %in% c('Densic material', 'Cemented horizon', 'Petrocalcic') & restrictions_valley$majcompflag=='Yes'), ]
# dim(misc_res_by_cokey) #
sum(duplicated(misc_res_by_cokey$cokey))
sum(duplicated(misc_res_by_cokey$mukey)) #5 will be averaged
misc_res_by_cokey[duplicated(misc_res_by_cokey$mukey),]
restrictions_valley[restrictions_valley$mukey %in% misc_res_by_cokey$mukey[duplicated(misc_res_by_cokey$mukey)],] #467123 has different depths for cokey 16620622
misc_res_by_mukey <- MUAggregate_wrapper(df1 = misc_res_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
misc_res_by_mukey[misc_res_by_mukey$mukey==461133, ]

#parent material info
parentmat <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_parentmat.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
parentmat_valley <- parentmat[parentmat$cokey %in% comp_data_valley$cokey, ]
# head(parentmat_valley)
# unique(parentmat_valley$pmkind) #11 different kinds
# table(parentmat_valley$pmkind)
parentmat_valley$majcompflag <- comp_data_valley$majcompflag[match(parentmat_valley$cokey, comp_data_valley$cokey)]
# table(parentmat_valley$majcompflag) #lots of minor component data
# dim(parentmat_valley)
# length(unique(parentmat_valley$cokey)) #lots of cokeys with more than one parent material, thus...
parentmat_by_cokey <- data.frame(cokey = row.names(tapply(parentmat_valley$pmkind[parentmat_valley$majcompflag=='Yes'], parentmat_valley$cokey[parentmat_valley$majcompflag=='Yes'],  function(x) unique(x))), pmkinds = as.character(tapply(parentmat_valley$pmkind[parentmat_valley$majcompflag=='Yes'], parentmat_valley$cokey[parentmat_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
# unique(parentmat_by_cokey$pmkinds)

productivity <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_forprod.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
# head(productivity)
productivity_valley <- productivity[productivity$cokey %in% comp_data_valley$cokey, ]
# dim(productivity_valley)
# length(unique(productivity_valley$cokey)) #some duplicates
productivity_valley <- productivity_valley[!duplicated(productivity_valley$cokey),]
# dim(productivity_valley)
# productivity_valley #very limited dataset

hillslope_pos <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_hillslopeprof.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))
hillslope_pos_valley <- hillslope_pos[hillslope_pos$cokey %in% comp_data_valley$cokey, ]
dim(hillslope_pos_valley)
# length(unique(hillslope_pos_valley$cokey))
# unique(hillslope_pos_valley$hillslopeprof)
hillslope_pos_valley$majcompflag <- comp_data_valley$majcompflag[match(hillslope_pos_valley$cokey, comp_data_valley$cokey)]
hillslope_pos_valley$mukey <- comp_data_valley$mukey[match(hillslope_pos_valley$cokey, comp_data_valley$cokey)]
table(hillslope_pos_valley$majcompflag)
table(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='Yes']) / sum(hillslope_pos_valley$majcompflag=='Yes')
table(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='No ']) / sum(hillslope_pos_valley$majcompflag=='No ')

#this data seems problematic
#for major components
hillslope_pos_by_cokey <- data.frame(cokey = row.names(tapply(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='Yes'], hillslope_pos_valley$cokey[hillslope_pos_valley$majcompflag=='Yes'], function(x) unique(x))), hillpos= as.character(tapply(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='Yes'], hillslope_pos_valley$cokey[hillslope_pos_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
# unique(hillslope_pos_by_cokey$hillpos)
hillslope_pos_by_mukey <- data.frame(mukey = row.names(tapply(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='Yes'], hillslope_pos_valley$mukey[hillslope_pos_valley$majcompflag=='Yes'], function(x) unique(x))), hillpos= as.character(tapply(hillslope_pos_valley$hillslopeprof[hillslope_pos_valley$majcompflag=='Yes'], hillslope_pos_valley$mukey[hillslope_pos_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
unique(hillslope_pos_by_mukey$hillpos)
table(hillslope_pos_by_mukey$hillpos)
sum(is.na(hillslope_pos_by_mukey$hillpos)) #only 6 NA
hillslope_pos_by_mukey$mukey[hillslope_pos_by_mukey$hillpos=='Footslope-Summit']
hillslope_pos_valley[hillslope_pos_valley$mukey==464368,]
comp_data_valley[comp_data_valley$mukey==464368, ]

#valley only dataset
reskinds_by_cokey <- data.frame(cokey = row.names(tapply(restrictions_valley$reskind[restrictions_valley$majcompflag=='Yes'], restrictions_valley$cokey[restrictions_valley$majcompflag=='Yes'],  function(x) unique(x))), reskinds = as.character(tapply(restrictions_valley$reskind[restrictions_valley$majcompflag=='Yes'], restrictions_valley$cokey[restrictions_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE) #function(x) unique(x)))) #function(x) if(all(is.na(x))) {NA} else if (length(unique(x[!is.na(x)]))==1) {unique(x[!is.na(x)])} else{paste(unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))], collapse = '-')})
# dim(reskinds_by_cokey) #3139 majcomp cokeys with restrictions
# length(unique(reskinds_by_cokey$cokey))
# lapply(reskinds_by_cokey, class)
# unique(reskinds_by_cokey$reskinds) #28 different reskind combos for valley
reskinds_by_cokey$mukey <- comp_data_valley$mukey[match(reskinds_by_cokey$cokey, comp_data_valley$cokey)]
length(unique(reskinds_by_cokey$mukey)) #2411 mukeys
reskinds_by_cokey$comp_pct <- comp_data_valley$comppct_r[match(reskinds_by_cokey$cokey, comp_data_valley$cokey)]
reskinds_by_cokey$majcompflag <- comp_data_valley$majcompflag[match(reskinds_by_cokey$cokey, comp_data_valley$cokey)]
#get rid of rock outcrop cokeys from this
reskinds_by_cokey$compname <- comp_data_valley$compname[match(reskinds_by_cokey$cokey, comp_data_valley$cokey)]
sum(reskinds_by_cokey$compname=='Rock outcrop') #hay 86
reskinds_by_cokey <- reskinds_by_cokey[reskinds_by_cokey$compname != 'Rock outcrop', ]
summary(as.factor(reskinds_by_cokey$majcompflag)) #minor components were left out in the original creation of the table above
summary(as.factor(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) length(x)))) #up to 4 per mukey, makes sense given that there are at most 3 major components in a map-unit for this AOI
#reskinds needs to be unconcatenated first before finding unique

reskinds_by_mukey <- data.frame(mukey = row.names(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) unique(x))), reskinds = as.character(tapply(reskinds_by_cokey$reskind, reskinds_by_cokey$mukey, concat_names, decat=TRUE)), stringsAsFactors = FALSE)
# unique(reskinds_by_mukey$reskinds) #now this appears to be ok, using the decat=TRUE
# table(reskinds_by_mukey$reskinds)

majcomps_no_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$majcompflag, comp_data_valley$mukey, function(x) sum(x=='Yes'))), majcomp_no = as.numeric(tapply(comp_data_valley$majcompflag, comp_data_valley$mukey, function(x) sum(x=='Yes'))), stringsAsFactors = FALSE)
unique(majcomps_no_by_mukey$majcomp_no)
sum(majcomps_no_by_mukey$majcomp_no==0) #2
comp_data_valley[comp_data_valley$mukey %in% majcomps_no_by_mukey$mukey[majcomps_no_by_mukey$majcomp_no==0], ] #its the ole' town map unit and Water

majcompnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$compname[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], concat_names)), majcompnames=as.character(tapply(comp_data_valley$compname[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
length(unique(majcompnames_by_mukey$majcompnames)) #1293 unique major component name combos

compnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$compname, comp_data_valley$mukey, function(x) unique(x))), compnames = as.character(tapply(comp_data_valley$compname, comp_data_valley$mukey, concat_names)), stringsAsFactors = FALSE)
length(unique(compnames_by_mukey$compnames)) #2960 different map-unit combos accounting for both major and minor components; decat=TRUE has no effect
dim(compnames_by_mukey) #but 5043 mukeys

majcokeys_by_mukey <- tapply(comp_data_valley$cokey[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) x)

majcomp_taxorders_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$taxorder[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], concat_names)), taxorders=as.character(tapply(comp_data_valley$taxorder[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
# unique(majcomp_taxorders_by_mukey$taxorders) #34 different

domcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$comppct_r, comp_data_valley$mukey, function(x) max(x, na.rm=TRUE))), docomppct=as.numeric(tapply(comp_data_valley$comppct_r, comp_data_valley$mukey, function(x) max(x, na.rm=TRUE))), stringsAsFactors = FALSE)
summary(domcomp_pct_by_mukey$docomppct)

majcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$comppct_r[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) sum(x))), majcomppct=as.numeric(tapply(comp_data_valley$comppct_r[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) sum(x))), stringsAsFactors = FALSE)
summary(majcomp_pct_by_mukey$majcomppct)

storie_rng_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$castorieindex[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) max_modified(x)- min_modified(x))), storierng=as.numeric(tapply(comp_data_valley$castorieindex[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) max_modified(x)- min_modified(x))), stringsAsFactors=FALSE)
hist(storie_rng_by_mukey$storierng)
summary(storie_rng_by_mukey$storierng)

storie_mean_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_valley$castorieindex[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) mean(x, na.rm=TRUE))), storiemn=as.numeric(tapply(comp_data_valley$castorieindex[comp_data_valley$majcompflag=='Yes'], comp_data_valley$mukey[comp_data_valley$majcompflag=='Yes'], function(x) mean(x, na.rm=TRUE))), stringsAsFactors=FALSE)
hist(storie_mean_by_mukey$storiemn)
sum(is.na(storie_mean_by_mukey$storiemn)) #922 are NA

#read in horizon data
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
horizon_data_valley <- horizon_data[horizon_data$cokey %in% comp_data_valley$cokey,]
# dim(horizon_data_valley) #25685 horizons
# unique(horizon_data_valley$hzname)
horizon_data_valley$majcompflag <- comp_data_valley$majcompflag[match(horizon_data_valley$cokey, comp_data_valley$cokey)]
horizon_data_valley$mukey <- comp_data_valley$mukey[match(horizon_data_valley$cokey, comp_data_valley$cokey)]
length(unique(horizon_data_valley$mukey)) #4900 mukeys with horizon data
table(horizon_data_valley$majcompflag) #4372 minor components have horizon data
head(horizon_data_valley[horizon_data_valley$majcompflag=='No ', ], 100)

#convert to Soil Profile collection class with no minor comps per filtering above
horizons_valley_majcomps <- horizon_data_valley[horizon_data_valley$majcompflag=='Yes', ]
# dim(horizons_valley_majcomps) #21313 horizons
# length(unique(horizons_valley_majcomps$cokey)) #6030 profiles (e.g. unique cokeys)
# length(unique(horizons_valley_majcomps$mukey)) #4897 after getting rid of minor components

#upgrade data.frame to a SoilProfileCollection object
#see https://r-forge.r-project.org/scm/viewvc.php/*checkout*/docs/aqp/aqp-intro.html?root=aqp
depths(horizons_valley_majcomps) <- cokey ~ hzdept_r + hzdepb_r
class(horizons_valley_majcomps)
print(horizons_valley_majcomps)
# reskinds_by_cokey_valley <- reskinds_by_cokey
# colnames(reskinds_by_cokey_valley)[2] <- 'kind'
# diagnostic_hz(horizons_valley_majcomps) <- reskinds_by_cokey_valley
# diagnostic_hz(horizons_valley_majcomps)
horizons_valley_majcomps$soil_depth <- profileApply(horizons_valley_majcomps, FUN = estimateSoilDepth, name='hzname', top='hzdept_r', bottom='hzdepb_r')

# apply custom functions and save results as a site-level attribute
# horizons_valley_majcomps$clay_wtd.mean_profile <- profileApply(horizons_valley_majcomps, FUN=wtd.mean, y='claytotal_r')
# horizons_valley_majcomps$kgOrg.m2_profile <- profileApply(horizons_valley_majcomps, FUN = kgOrgC_sum)

#slice into 1-cm increments to deepest extent as one approach
# sliced_horizons_valley <- slice(horizons_valley_majcomps, 0:222 ~ .)
# print(sliced_horizons_valley)
# horizons(sliced_horizons_valley)$awc_r[horizons(sliced_horizons_valley)$cokey==16607249]
# 
# sliced_horizons_valley$kgOrg.m2_30cm <- profileApply(sliced_horizons_valley, FUN = kgOrgC_sum, slice_it=TRUE, depth=30)
# hist(sliced_horizons_valley$kgOrg.m2_30cm)
# print(sliced_horizons_valley)
# depth_of_interest <- horizons(horizons_valley_majcomps)[horizons(horizons_valley_majcomps)$cokey==16607250,][1,]
# (depth_of_interest$om_r / 1.72) * depth_of_interest$dbthirdbar_r * (30 / 10) * (1 - depth_of_interest$fragvol_r_sum/100) #5.02
# print(sliced_horizons_valley)

#test with slab
# slab_results <- slab(horizons_valley_majcomps, fm = cokey ~ claytotal_r, slab.structure = c(0,10), slab.fun = mean) #+ silttotal_r + sandtotal_r + ksat_r + kwfact + ph1to1h2o_r + sar_r + caco3_r + gypsum_r + lep_r



#100 cm dataset
#Error: bad horizonation in IDs:16387801, 16597708 when running 100 cm aggregation
# comp_valley_100cm <- horizon_to_comp(horizon_SPC = horizons_valley_majcomps, depth = 100, comp_df = comp_data_valley)
# head(comp_valley_100cm)

#30 cm dataset
test_hz_aggregation <- horizon_to_comp(horizon_SPC = horizons_valley_majcomps[4], depth = 30, comp_df = comp_data_valley)
test_hz_aggregation
horizon_data_valley[horizon_data_valley$cokey==16264388,]
#om test: 1.791667 
2*25/30+0.75*5/30
#awc test: 4.75 cm
0.16*25 + 0.15*5

comp_valley_30cm <- horizon_to_comp(horizon_SPC = horizons_valley_majcomps, depth = 30, comp_df = comp_data_valley)
# write.csv(comp_valley_30cm, file.path(summaryDir, 'comp_valley_30cm_10.25.19.csv'), row.names = FALSE)
head(comp_valley_30cm)
dim(comp_valley_30cm) #6030
colnames(comp_valley_30cm)
lapply(comp_valley_30cm[ ,2:ncol(comp_valley_30cm)], summary)
lapply(comp_valley_30cm[,2:ncol(comp_valley_30cm)], function(x) unique(comp_valley_30cm$compname[is.na(x)]))
problematic_compnames <- unique(unlist(lapply(comp_valley_30cm[,c('compname', 'clay_30cm', 'om_30cm', 'cec_30cm', 'ksat_30cm', 'awc_30cm')], function(x) unique(comp_valley_30cm$compname[x==0]))))
problematic_compnames <- data.frame(compname=problematic_compnames[order(problematic_compnames)], stringsAsFactors = FALSE)
write.csv(problematic_compnames, file.path(summaryDir, 'problematic_compnames_10.25.19.csv'), row.names = FALSE)
unique(comp_valley_30cm$compname[which(comp_valley_30cm$om_30cm==0)])
zero_om_cokeys <- unique(comp_valley_30cm$cokey[which(comp_valley_30cm$om_30cm==0)])
horizon_data_zero_om <- horizon_data_valley[horizon_data_valley$cokey %in% zero_om_cokeys,]
horizon_data_zero_om$compname <- comp_valley_30cm$compname[match(horizon_data_zero_om$cokey, comp_valley_30cm$cokey)]
#write.csv(horizon_data_zero_om, file.path(summaryDir, 'horizon_data_zero_om_10.25.19.csv'), row.names=FALSE)

#convert zeroes for some variables (OM, CEC, Ks, AWC) to NA
#om
sum(comp_valley_30cm$om_30cm==0, na.rm = TRUE)
sum(is.na(comp_valley_30cm$om_30cm))

comp_valley_30cm$om_30cm[comp_valley_30cm$om_30cm==0] <- NA
#cec
sum(comp_valley_30cm$cec_30cm==0, na.rm = TRUE)
sum(is.na(comp_valley_30cm$cec_30cm))
comp_valley_30cm$cec_30cm[comp_valley_30cm$cec_30cm==0] <- NA
sum(is.na(comp_valley_30cm$cec_30cm))

#Ks
sum(comp_valley_30cm$ksat_30cm==0, na.rm = TRUE)
sum(is.na(comp_valley_30cm$ksat_30cm))
comp_valley_30cm$ksat_30cm[comp_valley_30cm$ksat_30cm==0] <- NA
sum(is.na(comp_valley_30cm$ksat_30cm))

#awc
sum(comp_valley_30cm$awc_30cm==0, na.rm = TRUE)
sum(is.na(comp_valley_30cm$awc_30cm))
comp_valley_30cm$awc_30cm[comp_valley_30cm$awc_30cm==0] <- NA
sum(is.na(comp_valley_30cm$awc_30cm))

length(unique(comp_valley_30cm$mukey)) #4897
length(unique(comp_data_valley$mukey)) #5043
length(unique(comp_valley_30cm$compname)) #899
length(unique(comp_data_valley$compname)) #1153
unique(comp_valley_30cm$compname)
#these variables are used for clustering later: c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm') #sar_30cm was removed
sum(is.na(comp_valley_30cm$om_30cm)) #210 comps have NA
sum(comp_valley_30cm$om_30cm==0, na.rm = TRUE) #66 comps impacted by om=0

#metadata about om data availability
comppct_by_mukey_om_data <- data.frame(mukey=row.names(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$om_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$om_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$om_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$om_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_om_data$comppct_tot)
sum(comppct_by_mukey_om_data$comppct_tot < 70) #290

#metadata about awc data availability
comppct_by_mukey_awc_data <- data.frame(mukey=row.names(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$awc_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$awc_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$awc_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$awc_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_awc_data$comppct_tot)
sum(comppct_by_mukey_awc_data$comppct_tot < 70) #286

#Ks
comppct_by_mukey_ksat_data <- data.frame(mukey=row.names(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$ksat_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$ksat_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$ksat_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$ksat_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_ksat_data$comppct_tot)
sum(comppct_by_mukey_ksat_data$comppct_tot < 70) #266

#cec
comppct_by_mukey_cec_data <- data.frame(mukey=row.names(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$cec_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$cec_30cm)], sum)), comppct_tot=as.numeric(tapply(comp_valley_30cm$comppct[!is.na(comp_valley_30cm$cec_30cm)], comp_valley_30cm$mukey[!is.na(comp_valley_30cm$cec_30cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_cec_data$comppct_tot)
sum(comppct_by_mukey_cec_data$comppct_tot < 70) #307

#10 cm dataset
# comp_valley_10cm <- horizon_to_comp(horizon_SPC = horizons_valley_majcomps, depth = 10, comp_df = comp_data_valley)
# head(comp_valley_10cm)
# hist(log(comp_valley_10cm$lep_10cm))
# hist(log(comp_valley_10cm$om_10cm))

#read in ecoregions
# list.files(ecoDir)
# eco_L3 <- shapefile(file.path(ecoDir, 'ca_eco_l3.shp'))
# plot(eco_L3)
# eco_L3_wgs84 <- spTransform(eco_L3, crs(fresno_area))
# eco_L3_Fresno <- crop(eco_L3_wgs84, fresno_area)
# plot(eco_L3_Fresno)
# eco_L3_Fresno$US_L3NAME
# 
# eco_L4 <- shapefile(file.path(ecoDir, 'ca_eco_l4.shp'))
# eco_L4_wgs84 <- spTransform(eco_L4, crs(fresno_area))
# eco_L4_Fresno <- crop(eco_L4_wgs84, fresno_area)
# plot(eco_L4_Fresno)
# unique(eco_L4_Fresno$US_L4NAME)

#read in crops
list.files(cropsDir)
# crops <- shapefile(file.path(cropsDir, 'i15_Crop_Mapping_2014_Final_LandIQonAtlas.shp'))

#add some variable to this subset of map units
# length(unique(valley_mu_shp$mukey)) #786 mukeys
valley_mu_shp$area_ac <- area(valley_mu_shp) / 10000 * 2.47105 #acres calc
# sum(valley_mu_shp$area_ac) #13873110 matches above
valley_mu_shp$mjcps_no <- majcomps_no_by_mukey$majcomp_no[match(valley_mu_shp$mukey, majcomps_no_by_mukey$mukey)]
# round(tapply(valley_mu_shp$area_ac, valley_mu_shp$majcomps_no, sum), 0)
#  0       1       2       3 
#1742 2,151,747  324,437  30,075  (85.8% of survey area has 1 major component per map unit, 12.9% has 2 major components, 1.2% has 3 major components)
valley_mu_shp$txorders <- majcomp_taxorders_by_mukey$taxorders[match(valley_mu_shp$mukey, majcomp_taxorders_by_mukey$mukey)]
# unique(valley_mu_shp$taxorders)
# taxorders_area <- data.frame(taxorders=row.names(tapply(valley_mu_shp$area_ac, valley_mu_shp$txorders, sum)), acres=as.numeric(tapply(valley_mu_shp$area_ac, valley_mu_shp$txorders, sum)), stringsAsFactors = FALSE)
# sum(taxorders_area$acres) #13497820 (375,290 less than above because NA class is dropped (see next 2 calcs)
# sum(valley_mu_shp$area_ac) - sum(taxorders_area$acres) #375290.2
# sum(valley_mu_shp$area_ac[is.na(valley_mu_shp$txorders)]) #375290.2
# taxorders_area[nrow(taxorders_area)+1, 'taxorders'] <- 'Undefined'
# taxorders_area[nrow(taxorders_area), 'acres'] <- sum(valley_mu_shp$area_ac[is.na(valley_mu_shp$txorders)])
# taxorders_area <- taxorders_area[order(taxorders_area$acres, decreasing=TRUE), ]
# taxorders_area
# write.csv(taxorders_area, file.path(summaryDir, 'taxorders_area_valley.csv'), row.names=FALSE)
valley_mu_shp$dmcmp_pct <- domcomp_pct_by_mukey$docomppct[match(valley_mu_shp$mukey, domcomp_pct_by_mukey$mukey)]
summary(valley_mu_shp$dmcmp_pct)
valley_mu_shp$mjcmp_pct <- majcomp_pct_by_mukey$majcomppct[match(valley_mu_shp$mukey, majcomp_pct_by_mukey$mukey)]
summary(valley_mu_shp$mjcmp_pct) #114 NAs
summary(valley_mu_shp$domcomp_pct / valley_mu_shp$majcomp_pct)

valley_mu_shp$compct_om <- comppct_by_mukey_om_data$comppct_tot[match(valley_mu_shp$mukey, comppct_by_mukey_om_data$mukey)]
summary(valley_mu_shp$compct_om) #5910 NAs
hist(valley_mu_shp$compct_om)
sum(valley_mu_shp$area_ac[is.na(valley_mu_shp$compct_om)]) #so, 303,581 acres have no organic matter data
sum(valley_mu_shp$area_ac[valley_mu_shp$compct_om < 70], na.rm = TRUE) #another 411,696 acres have some data but only for less than 70% of map unit
sum(valley_mu_shp$area_ac[valley_mu_shp$compct_om < 80], na.rm = TRUE) #or, 856,650 acres have some data but only for less than 80% of map unit
sum(valley_mu_shp$area_ac[valley_mu_shp$compct_om < 85], na.rm = TRUE) #1,930,013


#add more naming information
valley_mu_shp$muname <- mu_data_valley$muname[match(valley_mu_shp$mukey, mu_data_valley$mukey)]
unique(valley_mu_shp$muname[is.na(valley_mu_shp$compct_om)])
unique(valley_mu_shp$muname[valley_mu_shp$compct_om < 70])

valley_mu_shp$mjcmpnms <- majcompnames_by_mukey$majcompnames[match(valley_mu_shp$mukey, majcompnames_by_mukey$mukey)]
valley_mu_shp$complex <- ifelse(grepl('complex', valley_mu_shp$muname), 'Yes', 'No')
# table(valley_mu_shp$complex)
valley_mu_shp$associan <- ifelse(grepl('association', valley_mu_shp$muname), 'Yes', 'No') #272 yes
# table(valley_mu_shp$association)

#add Storie index range and mean (relative to major components only)
valley_mu_shp$storiemn <- storie_mean_by_mukey$storiemn[match(valley_mu_shp$mukey, storie_mean_by_mukey$mukey)]
# summary(valley_mu_shp$storiemn)
valley_mu_shp$storiern <- storie_rng_by_mukey$storierng[match(valley_mu_shp$mukey, storie_rng_by_mukey$mukey)]
# summary(valley_mu_shp$storierng)

#add concatenated restrictive info
valley_mu_shp$restrict <- reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]
table(valley_mu_shp$restrict)
sum(is.na(valley_mu_shp$restrict))
valley_mu_shp$restrict[is.na(valley_mu_shp$restrict)] <- 'None'
sum(valley_mu_shp$area_ac[valley_mu_shp$restrict=='None']) #8,978,326 map unit acres have no components with restrictions

valley_mu_shp$Rock_OC <- ifelse(grepl('Rock outcrop', compnames_by_mukey$compnames[match(valley_mu_shp$mukey, compnames_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Rock_OC)
sum(valley_mu_shp$area_ac[valley_mu_shp$Rock_OC=='Yes']) #531301.6 acres have at least some presence of rock oc

#not including rock OC in this one
# unique(reskinds_by_mukey$reskinds)
valley_mu_shp$Lithic <- ifelse(grepl('Lithic bedrock', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No') #632 were 'yes' after accounting for mukeys with more than one cokey with restrictions | valley_mu_shp$Rock_OC=='Yes' add to conditional
table(valley_mu_shp$Lithic)

valley_mu_shp$Paralith <- ifelse(grepl('Paralithic bedrock', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Paralith)

valley_mu_shp$Duripan <- ifelse(grepl('Duripan', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Duripan)

valley_mu_shp$ATC <- ifelse(grepl('Abrupt textural change', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No') #ATC=abrupt textural change
table(valley_mu_shp$ATC)

valley_mu_shp$Natric <- ifelse(grepl('Natric', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Natric)

valley_mu_shp$Salic <- ifelse(grepl('Salic', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Salic)

valley_mu_shp$SCTS <- ifelse(grepl('Strongly contrasting textural stratification', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$SCTS) 

valley_mu_shp$Misc_Res <- ifelse(grepl('Densic material|Cemented horizon|Petrocalcic', reskinds_by_mukey$reskinds[match(valley_mu_shp$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(valley_mu_shp$Misc_Res)
#another way to get a count
sum(table(valley_mu_shp$mukey[valley_mu_shp$mukey %in% misc_res_by_cokey$mukey])) 


#add awc info
valley_mu_shp$aws050wta <- mu_data_valley$aws050wta[match(valley_mu_shp$mukey, mu_data_valley$mukey)]
valley_mu_shp$aws100wta <- mu_data_valley$aws0100wta[match(valley_mu_shp$mukey, mu_data_valley$mukey)]
valley_mu_shp$aws150wta <- mu_data_valley$aws0150wta[match(valley_mu_shp$mukey, mu_data_valley$mukey)]

#not including rock OC in this definition
#however some rock OC components are listed as a lithic bedrock restrictive horizon so this is a complicating factor--removed above?
valley_lithic_comppct <- data.frame(mukey=row.names(tapply(comp_data_valley$comppct_r[comp_data_valley$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_valley$compname != 'Rock outcrop'], comp_data_valley$mukey[comp_data_valley$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_valley$compname != 'Rock outcrop'], sum)), compct_sum = as.numeric(tapply(comp_data_valley$comppct_r[comp_data_valley$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_valley$compname != 'Rock outcrop'], comp_data_valley$mukey[comp_data_valley$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_valley$compname != 'Rock outcrop'], sum)), stringsAsFactors = FALSE)
# head(valley_lithic_comppct)
# summary(valley_lithic_comppct$compct_sum)
# sum(valley_lithic_comppct$compct_sum==85) #20 are 85
# valley_lithic_comppct[valley_lithic_comppct$compct_sum==85,]
# comp_data_valley[comp_data_valley$mukey==463379,]
# reskinds_by_cokey[reskinds_by_cokey$mukey==463379,]
# rock_OC_cokeys <- comp_data_valley$cokey[comp_data_valley$compname=='Rock outcrop']
# length(rock_OC_cokeys)
# sum(rock_OC_cokeys %in% reskinds_by_cokey$cokey) #rock outcrop cokeys have been removed from reskind table
# sum(comp_data_valley$majcompflag[comp_data_valley$cokey %in% rock_OC_cokeys]=='Yes') #37 are major components

#rock outcrop component pct
# sum(grepl('Rock outcrop', compnames_by_mukey$compnames)) #70 here
# sum(grepl('Rock outcrop', majcompnames_by_mukey$majcompnames)) #37 here
valley_rockOC_comppct <- data.frame(mukey=row.names(tapply(comp_data_valley$comppct_r[comp_data_valley$compname=='Rock outcrop'], comp_data_valley$mukey[comp_data_valley$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), compct_sum = as.numeric(tapply(comp_data_valley$comppct_r[comp_data_valley$compname=='Rock outcrop'], comp_data_valley$mukey[comp_data_valley$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), stringsAsFactors = FALSE)
# dim(valley_rockOC_comppct)
# head(valley_rockOC_comppct)
# summary(valley_rockOC_comppct$compct_sum)
# comp_data_valley[comp_data_valley$mukey %in% valley_rockOC_comppct$mukey, ]

#paralithic component pct
valley_paralithic_comppct <- reskind_comppct(reskind = 'Paralithic', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
# summary(valley_paralithic_comppct$compct_sum)
# sum(valley_paralithic_comppct$compct_sum==100) #1 has 100% paralithic
# valley_paralithic_comppct$mukey[valley_paralithic_comppct$compct_sum==100]
# comp_data_valley[comp_data_valley$mukey==463500,]
# restrictions_valley[restrictions_valley$cokey %in% c(16608021, 16608022), ] #16608021 is "Rock land"

#duripan component pct
valley_duripan_comppct <- reskind_comppct(reskind = 'Duripan', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
# dim(valley_duripan_comppct)
# summary(valley_duripan_comppct$compct_sum)
#check it
# head(valley_duripan_comppct, 20)
# comp_data_valley[comp_data_valley$mukey==463412,] #has 70% Duripan
# restrictions_valley[restrictions_valley$cokey %in% comp_data_valley$cokey[comp_data_valley$mukey==463412],]
#-or-
# reskinds_by_cokey[reskinds_by_cokey$cokey %in% comp_data_valley$cokey[comp_data_valley$mukey==463412],]
#check one more
# valley_duripan_comppct$compct_sum[valley_duripan_comppct$mukey==464442] #has 90% Duripan
# comp_data_valley[comp_data_valley$mukey==464442,]
# reskinds_by_cokey[reskinds_by_cokey$cokey %in% comp_data_valley$cokey[comp_data_valley$mukey==464442],]

#abrupt textural contrast (ATC) comppct
valley_ATC_comppct <- reskind_comppct(reskind = 'Abrupt textural change', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
# dim(valley_ATC_comppct)
# summary(valley_ATC_comppct$compct_sum)

#natric comppct
valley_Natric_comppct <- reskind_comppct(reskind = 'Natric', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
dim(valley_Natric_comppct)
summary(valley_Natric_comppct$compct_sum)

##Salic comppct
valley_Salic_comppct <- reskind_comppct(reskind = 'Salic', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
dim(valley_Salic_comppct)
valley_Salic_comppct

#SCTS comppct
valley_SCTS_comppct <- reskind_comppct(reskind = 'Strongly contrasting textural stratification', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
valley_SCTS_comppct

valley_Misc_comppct <- reskind_comppct(reskind = 'Densic material|Cemented horizon|Petrocalcic', comp_df = comp_data_valley, reskind_df = reskinds_by_cokey)
valley_Misc_comppct

#add reskind comppct to mapunit
valley_mu_shp$Lthc_pct <- 0
valley_mu_shp$Lthc_pct[valley_mu_shp$Lithic=='Yes'] <- valley_lithic_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Lithic=='Yes'], valley_lithic_comppct$mukey)]
summary(valley_mu_shp$Lthc_pct)
valley_mu_shp$RckOC_pct <- 0
valley_mu_shp$RckOC_pct[valley_mu_shp$Rock_OC=='Yes'] <- valley_rockOC_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Rock_OC=='Yes'], valley_rockOC_comppct$mukey)]
summary(valley_mu_shp$RckOC_pct)
summary(rowSums(as.data.frame(valley_mu_shp[c('Lthc_pct', 'RckOC_pct')])))
valley_mu_shp$Plth_pct <- 0
valley_mu_shp$Plth_pct[valley_mu_shp$Paralith=='Yes'] <- valley_paralithic_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Paralith=='Yes'], valley_paralithic_comppct$mukey)]
summary(valley_mu_shp$Plth_pct)
valley_mu_shp$Drpn_pct <- 0
valley_mu_shp$Drpn_pct[valley_mu_shp$Duripan=='Yes'] <- valley_duripan_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Duripan=='Yes'], valley_duripan_comppct$mukey)]
summary(valley_mu_shp$Drpn_pct)
valley_mu_shp$ATC_pct <- 0
valley_mu_shp$ATC_pct[valley_mu_shp$ATC=='Yes'] <- valley_ATC_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$ATC=='Yes'], valley_ATC_comppct$mukey)]
summary(valley_mu_shp$ATC_pct)
valley_mu_shp$Natr_pct <- 0
valley_mu_shp$Natr_pct[valley_mu_shp$Natric=='Yes'] <- valley_Natric_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Natric=='Yes'], valley_Natric_comppct$mukey)]
summary(valley_mu_shp$Natr_pct)
valley_mu_shp$Salc_pct <- 0
valley_mu_shp$Salc_pct[valley_mu_shp$Salic=='Yes'] <- valley_Salic_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Salic=='Yes'], valley_Salic_comppct$mukey)]
valley_mu_shp$SCTS_pct <- 0
valley_mu_shp$SCTS_pct[valley_mu_shp$SCTS=='Yes'] <- valley_SCTS_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$SCTS=='Yes'], valley_SCTS_comppct$mukey)]
valley_mu_shp$MRes_pct <- 0
valley_mu_shp$MRes_pct[valley_mu_shp$Misc_Res=='Yes'] <- valley_Misc_comppct$compct_sum[match(valley_mu_shp$mukey[valley_mu_shp$Misc_Res=='Yes'], valley_Misc_comppct$mukey)]

#add restriction depth info
valley_mu_shp$Depth_aqp <- horizons_valley_majcomps$soil_depth[match(valley_mu_shp$mukey, horizons_valley_majcomps$mukey)]
valley_mu_shp$Lthc_dep <- assumed_depth
valley_mu_shp$Lthc_dep[valley_mu_shp$Lithic=='Yes'] <- lithic_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Lithic=='Yes'], lithic_by_mukey$mukey)]
valley_mu_shp$Plth_dep <- assumed_depth
valley_mu_shp$Plth_dep[valley_mu_shp$Paralith=='Yes'] <- paralithic_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Paralith=='Yes'], paralithic_by_mukey$mukey)]
valley_mu_shp$Drpn_dep <- assumed_depth
valley_mu_shp$Drpn_dep[valley_mu_shp$Duripan=='Yes'] <- duripan_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Duripan=='Yes'], duripan_by_mukey$mukey)]
valley_mu_shp$ATC_dep <- assumed_depth
valley_mu_shp$ATC_dep[valley_mu_shp$ATC=='Yes'] <- ATC_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$ATC=='Yes'], ATC_by_mukey$mukey)]
valley_mu_shp$Natr_dep <- assumed_depth
valley_mu_shp$Natr_dep[valley_mu_shp$Natric=='Yes'] <- natric_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Natric=='Yes'], natric_by_mukey$mukey)]
valley_mu_shp$Salc_dep <- assumed_depth
valley_mu_shp$Salc_dep[valley_mu_shp$Salic=='Yes'] <- salic_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Salic=='Yes'], salic_by_mukey$mukey)]
valley_mu_shp$SCTS_dep <- assumed_depth
valley_mu_shp$SCTS_dep[valley_mu_shp$SCTS=='Yes'] <- SCTS_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$SCTS=='Yes'], SCTS_by_mukey$mukey)]
valley_mu_shp$MRes_dep <- assumed_depth
valley_mu_shp$MRes_dep[valley_mu_shp$Misc_Res=='Yes'] <- misc_res_by_mukey$resdept_r[match(valley_mu_shp$mukey[valley_mu_shp$Misc_Res=='Yes'], misc_res_by_mukey$mukey)]
#minimum of all reskind depths
valley_mu_shp$MnRs_dep <- apply(as.data.frame(valley_mu_shp)[ ,c('Lthc_dep', 'Plth_dep', 'Drpn_dep', 'ATC_dep', 'Natr_dep', 'Salc_dep', 'SCTS_dep', 'MRes_dep')], 1, min)

#now add horizon aggregated data
#from previous comp level aggregation work
# lapply(comp_valley_10cm, class)
# colnames(comp_valley_10cm)[5:ncol(comp_valley_10cm)]
# valley_10cm_muagg <- MUAggregate_wrapper(df1=comp_valley_10cm, varnames = colnames(comp_valley_10cm)[5:ncol(comp_valley_10cm)])
# head(valley_10cm_muagg)
# dim(valley_10cm_muagg)
# lapply(valley_10cm_muagg, class)
valley_30cm_muagg <- MUAggregate_wrapper(df1=comp_valley_30cm, varnames = colnames(comp_valley_30cm)[5:ncol(comp_valley_30cm)])
# valley_100cm_muagg <- MUAggregate_wrapper(df1=comp_valley_100cm, varnames = colnames(comp_valley_100cm)[5:ncol(comp_valley_100cm)])

# names(valley_mu_shp)
# valley_mu_shp_10cm <- merge(valley_mu_shp, valley_10cm_muagg, by = 'mukey')
valley_mu_shp_30cm <- merge(valley_mu_shp, valley_30cm_muagg, by = 'mukey')
# valley_mu_shp_100cm <- merge(valley_mu_shp, valley_100cm_muagg, by = 'mukey')

#write 30cm to shapefile for cluster analysis
shapefile(valley_mu_shp_30cm, file.path(summaryDir, 'shapefiles with data', 'valley_30cm.shp'), overwrite=TRUE)

#write 10 cm to csv
# valley_10cm <- as.data.frame(valley_mu_shp_10cm)
# colnames(valley_10cm)
# write.csv(valley_10cm, file.path(summaryDir, 'valley_10cm_test.csv'), row.names = FALSE)
# lapply(valley_10cm, function(x) sum(is.na(x)))

#write 30 cm to csv
valley_30cm <- as.data.frame(valley_mu_shp_30cm)
# colnames(valley_30cm)
write.csv(valley_30cm, file.path(summaryDir, 'valley_30cm_data.csv'), row.names = FALSE)
# lapply(valley_30cm, function(x) sum(is.na(x)))

#write 100 cm to csv
# valley_100cm <- as.data.frame(valley_mu_shp_100cm)
# colnames(valley_100cm)
# write.csv(valley_100cm, file.path(summaryDir, 'valley_100cm_test.csv'), row.names = FALSE)
# plot(valley_100cm$aws100wta, valley_100cm$awc_100cm)
# sum(valley_100cm$awc_100cm  - valley_100cm$aws100wta > 0.2, na.rm = TRUE) #163 polygons have slightly different AWC values; most likely due to the fact that I did not include minor components
# unique(valley_100cm$muname[which(valley_100cm$awc_100cm  - valley_100cm$aws100wta > 0.2)])

#previous beginning of crops investigation
crops_fresno <- crop(crops, fresno_area_aea)
#shapefile(crops_fresno, file.path(cropsDir, 'fresno_only', 'crops_fresno.shp'))
unique(crops_fresno$Crop2014)
sum(crops_fresno$Acres) #1,718.419 acres, so 52% of total area
crop_acreage_fresno <- data.frame(crop=row.names(tapply(crops_fresno$Acres, crops_fresno$Crop2014, function(x) sum(x))), acres=as.numeric(tapply(crops_fresno$Acres, crops_fresno$Crop2014, function(x) sum(x))))
crop_acreage_fresno <- crop_acreage_fresno[order(crop_acreage_fresno$acres, decreasing = TRUE), ]
#write.csv(crop_acreage_fresno, file.path(summaryDir, 'crop_acreage_DWR14.csv'), row.names = FALSE)