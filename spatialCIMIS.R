#see "Daily reference evapotranspiration for California using satellite imagery and weather station measurement interpolation" by Hart et al. 2009
library(raster)
library(rgdal)
library(XML)
library(httr)
#NOTE: 1/1/2003-2/19/2003, 3/7/2003-3/9/2003, 4/9/2003-4/10/2003, 6/15/2003, 6/17/2003, 6/19/2003-6/22/2003, 6/24/2003, and a number more dates before 9/30/2003 are missing. 2004 water year to the present is good to go.
californiaDir <- 'C:/Users/smdevine/Desktop/SpatialData/CA_counties/government_units'
cellsofinterestDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/CIMIS'
#below is specific to downloading data
ca_ta <- showP4(showWKT("+init=epsg:3310"))
spatialCIMISdir <- 'D:/Dissertation/Allowable_Depletion/SpatialCIMIS'
#'C:/Users/smdevine/Desktop/Allowable_Depletion/SpatialCIMIS'
startyear <- '2018' #was 2003
endyear <- '2019'
startdate <- strptime(paste0("03/09/", startyear), '%m/%d/%Y') #was 6/14
enddate <- strptime(paste0("12/31/", endyear), '%m/%d/%Y') #was 9/12
datesequence <- seq.Date(from=as.Date(startdate), to=as.Date(enddate), by='day')
base_url <- 'http://cimis.casil.ucdavis.edu/cimis/'
varnames <- c('ETo.asc.gz', 'Tn.asc.gz', 'Tx.asc.gz', 'Tdew.asc.gz', 'U2.asc.gz')
dirnames <- c('ETo', 'Tmin', 'Tmax', 'Tdew', 'U2')
for (j in seq_along(varnames)) {
  varofinterest <- varnames[j]
  varDir <- dirnames[j]
  for (i in 1:length(datesequence)) {
    day <- format.Date(datesequence[i], '%d')
    mnth <- format.Date(datesequence[i], '%m')
    yr <- format.Date(datesequence[i], '%Y')
    if (file.exists(file.path(spatialCIMISdir, varDir)) == FALSE) {
      dir.create(file.path(spatialCIMISdir, varDir))
    }
    if (file.exists(file.path(spatialCIMISdir, varDir, yr)) == FALSE) {
      dir.create(file.path(spatialCIMISdir, varDir, yr))
    }
    url_download <- paste0(base_url, yr, '/', mnth, '/', day, '/', varofinterest)
    setwd(file.path(spatialCIMISdir, varDir, yr))
    error_catch <- try(download.file(url_download, destfile = varofinterest, quiet = TRUE, mode = 'wb'))
    if (class(error_catch)=='try-error') {
      print(error_catch)
      next
    }
    cimis_table <- read.table(varofinterest, skip=6, na.strings="-9999")
  #this two step operation is much faster than the "for loop" solution
    cimis_table <- t(cimis_table)
    cimis_values <- as.vector(cimis_table)
    cimis_raster <- raster(ncol=510, nrow=560, xmx=610000, xmn=-410000, ymn=-660000, ymx=460000, crs=ca_ta, vals=cimis_values)
    writeRaster(cimis_raster, filename = paste0(varDir, yr, mnth, day, '.tif'), format='GTiff', overwrite=TRUE) #each Gtiff will be 1092 KB
    file.remove(varofinterest)
  }
} #NOTE: ALL DATA MISSING FOR 2/23/18

#plot some files to check
yr <- 2017
setwd(file.path(spatialCIMISdir, yr))
raster_fnames <- list.files(pattern = glob2rx('*.tif'))
for(i in 1:100) {
  cimis_raster <- raster(raster_fnames[i])
  plot(cimis_raster, main=raster_fnames[i])
}
setwd(californiaDir)
list.files()
california <- shapefile("california_CA_TA.shp")
plot(california, add=T)

#calculate relative humidity
if (file.exists(file.path(spatialCIMISdir, 'minRH')) == FALSE) {
  dir.create(file.path(spatialCIMISdir, 'minRH'))
}
# startyear <- '2003'
# endyear <- '2017'
# startdate <- strptime(paste0("10/01/", startyear), '%m/%d/%Y')
# enddate <- strptime(paste0("06/13/", endyear), '%m/%d/%Y')
# datesequence <- seq.Date(from=as.Date(startdate), to=as.Date(enddate), by='day')
#used datesequence created 
for (i in 1:length(datesequence)) {
  day <- format.Date(datesequence[i], '%d')
  mnth <- format.Date(datesequence[i], '%m')
  yr <- format.Date(datesequence[i], '%Y')
  #setwd(file.path(spatialCIMISdir, 'Tmin', yr))
  #Tmin <- raster(paste0('Tmin', yr, mnth, day, '.tif'))
  if (file.exists(file.path(spatialCIMISdir, 'Tmax', yr, paste0('Tmax', yr, mnth, day, '.tif')))==FALSE) {
    print(paste0('No file for ', datesequence[i], '.'))
    next
  } else {
    Tmax <- raster(file.path(spatialCIMISdir, 'Tmax', yr, paste0('Tmax', yr, mnth, day, '.tif')))
    Tdew <- raster(file.path(spatialCIMISdir, 'Tdew', yr, paste0('Tdew', yr, mnth, day, '.tif')))
    Ea <- calc(Tdew, fun = function(x){0.6108*exp(17.27*x/(x+237.3))})
    Es <- calc(Tmax, fun= function(x){0.6108*exp(17.27*x/(x+237.3))})
  #Es <- overlay(x=Tmin, y=Tmax, fun=function(x,y){0.6108/2*(exp(17.27*x/(x+237.3))+exp(17.27*y/(y+237.3)))}) #this is if one desires the avg daily RH
    minRH <- overlay(x=Ea, y=Es, fun=function(x,y){100*x/y}) #with the above calculations for Ea and Es, this is an approximation of the daily minimum RH
    if (file.exists(file.path(spatialCIMISdir, 'minRH', yr)) == FALSE) {
      dir.create(file.path(spatialCIMISdir, 'minRH', yr))
    }
    writeRaster(minRH, filename = file.path(spatialCIMISdir, 'minRH', yr, paste0('minRH', yr, mnth, day, '.tif')), format='GTiff', overwrite=TRUE)
    print(i)
  }
}

#get SpatialCIMIS into matrix for cells of interest
#file naming convention for U2 rasters forgot '_' between 'U2' and date; hence the if/else statement below
SpCIMISExtract <- function(varname, startyear, endyear, startdate, enddate) {
  cellsofinterest <- read.csv(file.path(cellsofinterestDir, "CIMIS_cells_unique.csv"), stringsAsFactors = FALSE)
  cellsofinterest <- cellsofinterest[order(cellsofinterest$CIMIS_cells), ]
  print(dim(cellsofinterest))
  cellsofinterest_names <- paste0('cell_', as.character(cellsofinterest))
  startdate <- strptime(paste0(startdate, startyear), '%m/%d/%Y')
  enddate <- strptime(paste0(enddate, endyear), '%m/%d/%Y')
  datesequence <- seq.Date(from=as.Date(startdate), to=as.Date(enddate), by='day')
  cimis_data <- as.data.frame(matrix(nrow=length(datesequence), ncol=(length(cellsofinterest)+5)))
  colnames(cimis_data) <- c('dates', 'month', 'day', 'year', 'DOY', cellsofinterest_names)
  cimis_data$dates <- format.Date(datesequence, '%m_%d_%Y')
  cimis_data$month <- as.integer(format.Date(datesequence, '%m'))
  cimis_data$day <- as.integer(format.Date(datesequence, '%d'))
  cimis_data$year <- as.integer(format.Date(datesequence, '%Y'))
  cimis_data$DOY <- as.integer(format.Date(datesequence, '%j'))
  for (i in 1:length(datesequence)) {
    day <- format.Date(datesequence[i], '%d')
    mnth <- format.Date(datesequence[i], '%m')
    yr <- format.Date(datesequence[i], '%Y')
    if (file.exists(file.path(spatialCIMISdir, varname, yr, paste0(varname, yr, mnth, day, '.tif')))==FALSE) {
      print(paste0('No file for ', datesequence[i], '.'))
      next
    } else {
        spCIMIS <- raster(file.path(spatialCIMISdir, varname, yr, paste0(varname, yr, mnth, day, '.tif')))
      }
    cimis_data[i, 6:ncol(cimis_data)] <- extract(spCIMIS, cellsofinterest)
    print(i)
  }
#write.csv(cimis_data, paste0('SpatialCIMIS_', varname, '_data.csv'), row.names = F)
#cimis_data <- read.csv('SpatialCIMIS_data.csv') #this has 86,600,954 cells
  cimis_data2 <- cbind(cimis_data[ ,1:5], round(cimis_data[ ,6:ncol(cimis_data)], 3)) #cimis_data['cell_number'] preserves data.frame class
  write.csv(cimis_data2, file.path(cellsofinterestDir, paste0('SpatialCIMIS.', varname, 'update.rounded.csv')), row.names=FALSE)
}
# SpCIMISExtract('U2', '2003', '2018', '10/01/', '03/08/')
# SpCIMISExtract('minRH', '2003', '2018', '10/01/', '03/08/')
SpCIMISExtract('ETo', '2004', '2019', '01/01/', '12/31/')
list.files(cellsofinterestDir)
ETo_df <- read.csv(file.path(cellsofinterestDir, 'SpatialCIMIS.EToupdate.rounded.csv'), stringsAsFactors = FALSE)
dim(ETo_df)

#QC check
ETo <- ETo_df
dim(ETo)
sum(is.na(ETo[,6:ncol(ETo)])) #358,851 NAs as of 12/31/19; cell_148533 no longer a cell of interest, so all but 76 are from 2/23/18 missing data
sum(ETo[,6:ncol(ETo)] < 0, na.rm = TRUE) #618 negative numbers; this is 0.002% negative
sum(ETo[,6:ncol(ETo)] == 0, na.rm = TRUE) #61,383 observations equal to 0
max(ETo[,6:ncol(ETo)], na.rm = TRUE) #10.702 mm/day max
#strategy is to fill in 2/23/18 first; use mean of all previous days, excluding negative numbers
check2.23.data <- unlist(ETo[which(ETo$month==2 & ETo$day==23 & ETo$year!=2018),6:ncol(ETo)])
sum(check2.23.data < 0) #0 are negative
sum(is.na(check2.23.data)) #0 are NA, so don't have to worry about that for 2/23/18 correction
rm(check2.23.data)
test <- apply(ETo[which(ETo$month==2 & ETo$day==23 & ETo$year!=2018),6:ncol(ETo)], 2, function(x) {mean(x[x>0])})
summary(test)
head(test)
length(test)
rm(test)
ETo[ETo$dates=='02_23_2018', 6:ncol(ETo)] <- apply(ETo[which(ETo$month==2 & ETo$day==23 & ETo$year!=2018),6:ncol(ETo)], 2, function(x) {mean(x[x>0])}) #reduces NA total to 76
sum(is.na(ETo[,6:ncol(ETo)])) #346,037
#get rid of 2019 data since most of Dec 2019 still missing
ETo <- ETo[-which(ETo$year==2019),]
sum(is.na(ETo[,6:ncol(ETo)])) #59 remaining
na_presence_ETo <- lapply(ETo[,6:ncol(ETo)], function(x) {sum(is.na(x))})
na_presence_ETo <- na_presence_ETo[na_presence_ETo > 0]
na_presence_ETo #7 cells have NAs
for (i in seq_along(na_presence_ETo)) {
  print(na_presence_ETo[i])
  print(ETo$dates[is.na(ETo[names(na_presence_ETo)[i]])])
} #they are all different dates
for (i in seq_along(na_presence_ETo)) {
  dates <- ETo$dates[is.na(ETo[names(na_presence_ETo)[i]])]
  for (j in seq_along(dates)) {
    day <- ETo$day[ETo$dates==dates[j]]
    month <- ETo$month[ETo$dates==dates[j]]
    gap.fill.data <- ETo[which(ETo$month==month & ETo$day==day & ETo$dates != '02_23_2018'), names(na_presence_ETo)[i]]
    ETo[ETo$dates==dates[j], names(na_presence_ETo)[i]] <- mean(gap.fill.data[gap.fill.data > 0], na.rm = TRUE) #remove any negative values from averages
  }
}
sum(is.na(ETo[,6:ncol(ETo)]))==0 #all NAs removed 2004-2018
#now apply a negative value correction procedure
neg_presence_ETo <- lapply(ETo[,6:ncol(ETo)], function(x) {sum(x < 0)})
neg_presence_ETo <- neg_presence_ETo[neg_presence_ETo > 0]
neg_presence_ETo
length(neg_presence_ETo) #455 cells with negative ETo values
for (i in seq_along(neg_presence_ETo)) {
  print(neg_presence_ETo[i])
  #print(ETo$dates[is.na(ETo[names(neg_presence_ETo)[i]])])
} #they are all different dates
for (i in seq_along(neg_presence_ETo)) {
  dates <- ETo$dates[ETo[names(neg_presence_ETo)[i]] < 0]
  for (j in seq_along(dates)) {
    day <- ETo$day[ETo$dates==dates[j]]
    month <- ETo$month[ETo$dates==dates[j]]
    gap.fill.data <- ETo[which(ETo$month==month & ETo$day==day & ETo$dates != '02_23_2018'), names(neg_presence_ETo)[i]]
    ETo[ETo$dates==dates[j], names(neg_presence_ETo)[i]] <- mean(gap.fill.data[gap.fill.data > 0], na.rm = TRUE) #remove any negative values or zeroes before calculating averages to gap fill
  }
}

sum(ETo[,6:ncol(ETo)] < 0, na.rm = TRUE)==0 #good
sum(ETo[,6:ncol(ETo)] == 0, na.rm = TRUE) #still have 1777

#now, correct the zeroes
zero_presence_ETo <- lapply(ETo[,6:ncol(ETo)], function(x) {sum(x == 0)})
zero_presence_ETo <- zero_presence_ETo[zero_presence_ETo > 0] #identify cells that have at least 1 count of a zero value for ETo across all days
zero_presence_ETo
length(zero_presence_ETo) #27 cells with ETo values of zero 1 day
for (i in seq_along(zero_presence_ETo)) {
  print(zero_presence_ETo[i])
  print(ETo$dates[ETo[names(zero_presence_ETo)[i]]==0])
} #they are all different dates
for (i in seq_along(zero_presence_ETo)) {
  dates <- ETo$dates[ETo[names(zero_presence_ETo)[i]] == 0]
  for (j in seq_along(dates)) {
    day <- ETo$day[ETo$dates==dates[j]]
    month <- ETo$month[ETo$dates==dates[j]]
    gap.fill.data <- ETo[which(ETo$month==month & ETo$day==day & ETo$dates != '02_23_2018'), names(zero_presence_ETo)[i]]
    ETo[ETo$dates==dates[j], names(zero_presence_ETo)[i]] <- mean(gap.fill.data[gap.fill.data > 0], na.rm = TRUE) #remove any negative values or zeroes before calculating averages to gap fill
  }
}
sum(is.na(ETo[,6:ncol(ETo)])) #still good
sum(ETo[,6:ncol(ETo)] < 0, na.rm = TRUE) #still good
sum(ETo[,6:ncol(ETo)] == 0, na.rm = TRUE) #good now
max(ETo[,6:ncol(ETo)]) #max 10.702 mm/day;
hist(as.numeric(lapply(ETo[,6:ncol(ETo)], mean))) #centered over 3.8-4.0 mm/day
max(as.numeric(lapply(ETo[,6:ncol(ETo)], mean))) #max 6.3 mm/day across all days for one cell
min(as.numeric(lapply(ETo[,6:ncol(ETo)], mean))) #min 2.3 mm/day across all days for one cell
write.csv(ETo, file.path(cellsofinterestDir, 'ETo_daily_2004_2018_QCpass.csv'), row.names = FALSE)