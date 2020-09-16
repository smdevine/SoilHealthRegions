workDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture'
soilProps <- read.csv(file.path(workDir, 'soil_properties_SSURGO_by_SHR.csv'), stringsAsFactors = FALSE)
head(soilProps)
concat_func <- function(x,y)  {paste0(x, " (", y, ")")}
concat_func(soilProps$Sand, soilProps$Sand.IQR)
finalTable <- cbind(soilProps[1], Sand=concat_func(soilProps$Sand, soilProps$Sand.IQR), Silt=concat_func(soilProps$Silt, soilProps$Silt.IQR), Clay=concat_func(soilProps$Clay, soilProps$Clay.IQR), OM=concat_func(soilProps$OM, soilProps$OM.IQR), CEC=concat_func(soilProps$CEC, soilProps$CEC.IQR), pH=concat_func(soilProps$pH, soilProps$pH.IQR), EC=concat_func(soilProps$EC, soilProps$EC.IQR), LEP=concat_func(soilProps$LEP, soilProps$LEP.IQR), Ksat=concat_func(soilProps$Ksat, soilProps$Ksat.IQR), Storie=concat_func(soilProps$Storie, soilProps$Storie.IQR), deparse.level = 1, stringsAsFactors=FALSE)
write.csv(finalTable, file = file.path(workDir, 'ssurgo_properties_final_bySHR.csv'), row.names = FALSE)
