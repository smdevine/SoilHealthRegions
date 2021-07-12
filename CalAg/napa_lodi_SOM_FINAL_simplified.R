SOC_to_SOM <- 1.72
laptop <- TRUE
account_for_compost <- FALSE
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
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/Figures' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
  CalAgDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/napa lodi analysis'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}

clus_7_names <- c('6. Fine salt-affected', '3. Low OM with restrictive horizons', '4. High OM with restrictive horizons', '1. Coarse with no restrictions', '2. Loamy with no restrictions', '7. Shrink-swell', '5. Coarse-loamy salt-affected')

soil_data <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', "")) #colClasses = columnClass,
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))

#CDFA C validation for 7-region model
#compare 3 management types in cluster 6 across each depth segment
pts_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
pts_30cm$cc_type <- soil_data$perennial.Covercrop.vs.annual.covercrop[match(pts_30cm$Concatenate, soil_data$Concatenate)]
pts_30cm$SHR7name <- clus_7_names[pts_30cm$cluster_7]
pts_10cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_10cm.csv'), stringsAsFactors = FALSE)
pts_10cm$SHR7name <- clus_7_names[pts_10cm$cluster_7]
pts_50cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_50cm.csv'), stringsAsFactors = FALSE)
pts_50cm$SHR7name <- clus_7_names[pts_50cm$cluster_7]

# columnClass <- c(Concatenate='character')
unique(pts_30cm$SHR7name)

napa_lodi_data <- data.frame(ID=pts_30cm$Concatenate[!is.na(pts_30cm$SHR7name)], tillage=pts_30cm$tillage[!is.na(pts_30cm$SHR7name)], irrigation=pts_30cm$irrigated_vs_dryfarm[!is.na(pts_30cm$SHR7name)], compost=pts_30cm$compost_added[!is.na(pts_30cm$SHR7name)], cc_type=pts_30cm$cc_type[!is.na(pts_30cm$SHR7name)], stringsAsFactors = FALSE)
head(napa_lodi_data)
dim(napa_lodi_data)


soil_data_shrink_swell <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="7. Shrink-swell")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="7. Shrink-swell")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="7. Shrink-swell")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="7. Shrink-swell")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="7. Shrink-swell")], stringsAsFactors = FALSE)
head(soil_data_shrink_swell)
dim(soil_data_shrink_swell)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_shrink_swell$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_shrink_swell$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_shrink_swell$soil_C_0_10cm <- apply(soil_data_shrink_swell[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_shrink_swell$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_shrink_swell$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_shrink_swell$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_shrink_swell$soil_C_0_30cm <- apply(soil_data_shrink_swell[,c('soil_C_0_5cm', 'soil_C_5_10cm', 'soil_C_10_30cm')], 1, weighted.mean, w=c(5,5,20)/30)
soil_data_shrink_swell$soil_C_30_100cm <- apply(soil_data_shrink_swell[,c('soil_C_30_50cm', 'soil_C_50_100cm')], 1, weighted.mean, w=c(20,50)/70)
soil_data_shrink_swell <- soil_data_shrink_swell[order(soil_data_shrink_swell$tillage, soil_data_shrink_swell$irrigation, soil_data_shrink_swell$cc_type, soil_data_shrink_swell$compost), ]
soil_data_shrink_swell$SHR <- 7
colnames(soil_data_shrink_swell)
# write.csv(soil_data_shrink_swell, file.path(CalAgDir, 'soil_data_shrink_swell.csv'), row.names = FALSE)
shrink_swell <- data.frame(SOM_means=sapply(soil_data_shrink_swell[,12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_shrink_swell[,12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_shrink_swell[,12:13], function(x) {length(x[!is.na(x)])}), SOM_q25=sapply(soil_data_shrink_swell[,12:13]*SOC_to_SOM, quantile, probs=0.25, na.rm=TRUE), SOM_medians=sapply(soil_data_shrink_swell[,12:13]*SOC_to_SOM, median, na.rm=TRUE), SOM_q75=sapply(soil_data_shrink_swell[,12:13]*SOC_to_SOM, quantile, probs=0.75, na.rm=TRUE))
shrink_swell
# write.csv(shrink_swell, file.path(CalAgDir, 'shrink_swell_soilC.csv'), row.names = FALSE)


soil_data_loamy_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="2. Loamy with no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="2. Loamy with no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="2. Loamy with no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="2. Loamy with no restrictions")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="2. Loamy with no restrictions")], stringsAsFactors = FALSE)
head(soil_data_loamy_w_no_res)
dim(soil_data_loamy_w_no_res)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_loamy_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_no_res$soil_C_0_10cm <- apply(soil_data_loamy_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_loamy_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_loamy_w_no_res$soil_C_0_30cm <- apply(soil_data_loamy_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm', 'soil_C_10_30cm')], 1, weighted.mean, w=c(5,5,20)/30)
soil_data_loamy_w_no_res$soil_C_30_100cm <- apply(soil_data_loamy_w_no_res[,c('soil_C_30_50cm', 'soil_C_50_100cm')], 1, weighted.mean, w=c(20,50)/70)

soil_data_loamy_w_no_res <- soil_data_loamy_w_no_res[order(soil_data_loamy_w_no_res$tillage, soil_data_loamy_w_no_res$irrigation, soil_data_loamy_w_no_res$cc_type, soil_data_loamy_w_no_res$compost), ]
soil_data_loamy_w_no_res$SHR <- 2
colnames(soil_data_loamy_w_no_res)
# write.csv(soil_data_loamy_w_no_res, file.path(CalAgDir, 'soil_data_loamy_w_no_res.csv'), row.names = FALSE)

loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[,12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[,12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[,12:13], function(x) {length(x[!is.na(x)])}), SOM_q25=sapply(soil_data_loamy_w_no_res[,12:13]*SOC_to_SOM, quantile, probs=0.25, na.rm=TRUE), SOM_medians=sapply(soil_data_loamy_w_no_res[,12:13]*SOC_to_SOM, median, na.rm=TRUE), SOM_q75=sapply(soil_data_loamy_w_no_res[,12:13]*SOC_to_SOM, quantile, probs=0.75, na.rm=TRUE))
loamy_w_no_res
# write.csv(loamy_w_no_res, file.path(CalAgDir, 'loamy_w_no_res_soilC.csv'), row.names = FALSE)

till_irr_ann_loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])})) #2a in plot
till_irr_ann_loamy_w_no_res

notill_irr_ann_loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])})) #2b in plot
notill_irr_ann_loamy_w_no_res

notill_irr_per_loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])}))
notill_irr_per_loamy_w_no_res

till_dry_ann_loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])})) #2c in plot
till_dry_ann_loamy_w_no_res

notill_dry_ann_loamy_w_no_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])}))
notill_dry_ann_loamy_w_no_res

soil_data_coarse_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="1. Coarse with no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="1. Coarse with no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="1. Coarse with no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="1. Coarse with no restrictions")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="1. Coarse with no restrictions")], stringsAsFactors = FALSE)
soil_data_coarse_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_w_no_res$soil_C_0_10cm <- apply(soil_data_coarse_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_coarse_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_coarse_w_no_res$soil_C_0_30cm <- apply(soil_data_coarse_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm', 'soil_C_10_30cm')], 1, weighted.mean, w=c(5,5,20)/30)
soil_data_coarse_w_no_res$soil_C_30_100cm <- apply(soil_data_coarse_w_no_res[,c('soil_C_30_50cm', 'soil_C_50_100cm')], 1, weighted.mean, w=c(20,50)/70)

soil_data_coarse_w_no_res <- soil_data_coarse_w_no_res[order(soil_data_coarse_w_no_res$tillage, soil_data_coarse_w_no_res$irrigation, soil_data_coarse_w_no_res$cc_type, soil_data_coarse_w_no_res$compost), ]
soil_data_coarse_w_no_res$SHR <- 1
# write.csv(soil_data_coarse_w_no_res, file.path(CalAgDir, 'soil_data_coarse_w_no_res.csv'), row.names = FALSE)

coarse_w_no_res <- data.frame(SOM_means=sapply(soil_data_coarse_w_no_res[,12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_coarse_w_no_res[,12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_w_no_res[,12:13], function(x) {length(x[!is.na(x)])}), SOM_q25=sapply(soil_data_coarse_w_no_res[,12:13]*SOC_to_SOM, quantile, probs=0.25, na.rm=TRUE), SOM_medians=sapply(soil_data_coarse_w_no_res[,12:13]*SOC_to_SOM, median, na.rm=TRUE), SOM_q75=sapply(soil_data_coarse_w_no_res[,12:13]*SOC_to_SOM, quantile, probs=0.75, na.rm=TRUE))
coarse_w_no_res
# write.csv(coarse_w_no_res, file.path(CalAgDir, 'coarse_w_no_res_soilC.csv'), row.names = FALSE)

till_irr_coarse_w_no_res <- data.frame(SOM_means=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])}))
till_irr_coarse_w_no_res #excludes 3 dryfarm coarse w/no res in Napa

soil_data_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="4. High OM with restrictive horizons")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="4. High OM with restrictive horizons")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="4. High OM with restrictive horizons")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="4. High OM with restrictive horizons")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="4. High OM with restrictive horizons")], stringsAsFactors = FALSE)
soil_data_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_res$soil_C_0_5cm[soil_data_loamy_w_res$soil_C_0_5cm==0] <- NA #convert 0 value to NA
soil_data_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_res$soil_C_0_10cm <- apply(soil_data_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_loamy_w_res$soil_C_0_10cm[soil_data_loamy_w_res$ID==853] <- soil_data_loamy_w_res$soil_C_5_10cm[soil_data_loamy_w_res$ID==853]
soil_data_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_loamy_w_res$soil_C_0_30cm <- apply(soil_data_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm', 'soil_C_10_30cm')], 1, weighted.mean, w=c(5,5,20)/30)
soil_data_loamy_w_res$soil_C_30_100cm <- apply(soil_data_loamy_w_res[,c('soil_C_30_50cm', 'soil_C_50_100cm')], 1, weighted.mean, w=c(20,50)/70)

soil_data_loamy_w_res <- soil_data_loamy_w_res[order(soil_data_loamy_w_res$tillage, soil_data_loamy_w_res$irrigation, soil_data_loamy_w_res$cc_type, soil_data_loamy_w_res$compost), ]
soil_data_loamy_w_res$SHR <- 4
# write.csv(soil_data_loamy_w_res, file.path(CalAgDir, 'soil_data_loamy_w_res.csv'), row.names = FALSE)
colnames(soil_data_loamy_w_res)
loamy_w_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_res[,12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_res[,12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[,12:13], function(x) {length(x[!is.na(x)])}), SOM_q25=sapply(soil_data_loamy_w_res[,12:13]*SOC_to_SOM, quantile, probs=0.25, na.rm=TRUE), SOM_medians=sapply(soil_data_loamy_w_res[,12:13]*SOC_to_SOM, median, na.rm=TRUE), SOM_q75=sapply(soil_data_loamy_w_res[,12:13]*SOC_to_SOM, quantile, probs=0.75, na.rm=TRUE))
loamy_w_res
# write.csv(loamy_w_res, file.path(CalAgDir, 'loamy_w_res_soilC.csv'), row.names = FALSE)

notill_irr_ann_loamy_w_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])})) #4 in plot
notill_irr_ann_loamy_w_res

notill_irr_per_loamy_w_res <- data.frame(SOM_means=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),12:13], function(x) {length(x[!is.na(x)])}))
notill_irr_per_loamy_w_res

soil_data_coarse_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="3. Low OM with restrictive horizons")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="3. Low OM with restrictive horizons")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="3. Low OM with restrictive horizons")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="3. Low OM with restrictive horizons")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="3. Low OM with restrictive horizons")], stringsAsFactors = FALSE)

soil_data_coarse_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_loamy_w_res$soil_C_0_10cm <- apply(soil_data_coarse_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_coarse_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm[soil_data_coarse_loamy_w_res$soil_C_50_100cm==0] <- NA
soil_data_coarse_loamy_w_res$soil_C_0_30cm <- apply(soil_data_coarse_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm', 'soil_C_10_30cm')], 1, weighted.mean, w=c(5,5,20)/30)
soil_data_coarse_loamy_w_res$soil_C_30_100cm <- apply(soil_data_coarse_loamy_w_res[,c('soil_C_30_50cm', 'soil_C_50_100cm')], 1, weighted.mean, w=c(20,50)/70)

soil_data_coarse_loamy_w_res <- soil_data_coarse_loamy_w_res[order(soil_data_coarse_loamy_w_res$tillage, soil_data_coarse_loamy_w_res$irrigation, soil_data_coarse_loamy_w_res$cc_type, soil_data_coarse_loamy_w_res$compost), ]
soil_data_coarse_loamy_w_res$SHR <- 3
# write.csv(soil_data_coarse_loamy_w_res, file.path(CalAgDir, 'soil_data_coarse_loamy_w_res.csv'), row.names = FALSE)
coarse_loamy_w_res <- data.frame(SOM_means=sapply(soil_data_coarse_loamy_w_res[,12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_coarse_loamy_w_res[,12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_loamy_w_res[,12:13], function(x) {length(x[!is.na(x)])}), SOM_q25=sapply(soil_data_coarse_loamy_w_res[,12:13]*SOC_to_SOM, quantile, probs=0.25, na.rm=TRUE), SOM_medians=sapply(soil_data_coarse_loamy_w_res[,12:13]*SOC_to_SOM, median, na.rm=TRUE), SOM_q75=sapply(soil_data_coarse_loamy_w_res[,12:13]*SOC_to_SOM, quantile, probs=0.75, na.rm=TRUE))
coarse_loamy_w_res
# write.csv(coarse_loamy_w_res, file.path(CalAgDir, 'coarse_loamy_w_res_soilC.csv'), row.names = FALSE)

till_irr_ann_coarse_loamy_w_res <- data.frame(SOM_means=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated' & soil_data_coarse_loamy_w_res$cc_type!='perennial'),12:13]*SOC_to_SOM, mean, na.rm=TRUE), SOM_se=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'& soil_data_coarse_loamy_w_res$cc_type!='perennial'), 12:13]*SOC_to_SOM, function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'& soil_data_coarse_loamy_w_res$cc_type!='perennial'), 12:13], function(x) {length(x[!is.na(x)])}))
till_irr_ann_coarse_loamy_w_res

soil_data_allSHR <- rbind(soil_data_coarse_w_no_res, soil_data_loamy_w_no_res, soil_data_coarse_loamy_w_res, soil_data_loamy_w_res, soil_data_shrink_swell)
soil_data_allSHR$SHRcolor <- ifelse(soil_data_allSHR$SHR==1, 'lightgoldenrod', ifelse(soil_data_allSHR$SHR==2, 'tan4', ifelse(soil_data_allSHR$SHR==3, 'olivedrab3', ifelse(soil_data_allSHR$SHR==4, 'firebrick3', ifelse(soil_data_allSHR$SHR==7, 'violetred', 'black')))))
table(soil_data_allSHR$SHR)
soil_data_allSHR$Area <- soil_data$Area[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$Vineyard_Mgmt <- soil_data$Vineyard.Management[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data$years.in.vineyard[soil_data$years.in.vineyard=='unknown'] <- NA
soil_data$years.in.vineyard <- as.integer(soil_data$years.in.vineyard)
soil_data_allSHR$yrs_in_vineyard <- soil_data$years.in.vineyard[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$slope <- soil_data$Cluster.SLOPE[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$great_group <- soil_data$Cluster.GREATGROUP[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$soil_series <- soil_data$Soil.Series.SSURGO[match(soil_data_allSHR$ID, soil_data$Concatenate)]

# write.csv(soil_data_allSHR, file.path(CalAgDir, 'soil_data_allSHR.csv'), row.names = FALSE)
dim(soil_data_allSHR)
table(soil_data_allSHR$SHR)
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$tillage[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$irrigation[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$compost[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$cc_type[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$Area[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$soil_series[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$great_group[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$Vineyard_Mgmt[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {summary(soil_data_allSHR$slope[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {summary(soil_data_allSHR$yrs_in_vineyard[soil_data_allSHR$SHR==x])})
soil_data_allSHR[soil_data_allSHR$SHR==1,]
soil_data_allSHR[soil_data_allSHR$SHR==2,]

#see confidence interval formula for unknown mean and standard deviation at http://www.stat.yale.edu/Courses/1997-98/101/confint.htm
plot_da_error <- function(barnum, barcenters, df) {
  segments(x0=barcenters[barnum,], y0=df$means[barnum] + qt(p=0.975, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), x1=barcenters[barnum,], y1=df$means[barnum] + qt(p=0.025, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), lwd = 1.2)
}

#SOC concentration barplots
SOMmeans_SHRmgmt <- rbind(coarse_w_no_res$SOM_means, till_irr_ann_loamy_w_no_res$SOM_means, notill_irr_ann_loamy_w_no_res$SOM_means, notill_irr_per_loamy_w_no_res$SOM_means, till_dry_ann_loamy_w_no_res$SOM_means, notill_dry_ann_loamy_w_no_res$SOM_means, coarse_loamy_w_res$SOM_means, loamy_w_res$SOM_means, shrink_swell$SOM_means)
row.names(SOMmeans_SHRmgmt) <- c('1', '2a', '2b', '2c', '2d', '2e', '3', '4', '7')
colnames(SOMmeans_SHRmgmt) <- c('0-10 cm', '10-30 cm', '30-50 cm', '50-100 cm')
SOMmeans_SHRmgmt
SOM_SE_SHRmgmt <- rbind(coarse_w_no_res$SOM_se, till_irr_ann_loamy_w_no_res$SOM_se, notill_irr_ann_loamy_w_no_res$SOM_se, notill_irr_per_loamy_w_no_res$SOM_se, till_dry_ann_loamy_w_no_res$SOM_se, notill_dry_ann_loamy_w_no_res$SOM_se, coarse_loamy_w_res$SOM_se, loamy_w_res$SOM_se, shrink_swell$SOM_se)
SOM_SE_SHRmgmt
SOM_n_SHRmgmt <- rbind(coarse_w_no_res$n, till_irr_ann_loamy_w_no_res$n, notill_irr_ann_loamy_w_no_res$n, notill_irr_per_loamy_w_no_res$n, till_dry_ann_loamy_w_no_res$n, notill_dry_ann_loamy_w_no_res$n, coarse_loamy_w_res$n, loamy_w_res$n, shrink_swell$n)
SOM_n_SHRmgmt

tiff(file = file.path(FiguresDir, 'Napa_Lodi_soilSOM_barchart_two_depths_FINAL.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 4.25, units = 'in', res=800, compression = 'lzw')
par(mar=c(2.4, 3.5, 1.25, 0.5))
SOMmeans_bardims <- barplot(SOMmeans_SHRmgmt, beside=TRUE, legend.text=c('1. Coarse with no restrictions (n=14)', '2a. Loamy with no restrictions: tilled (n=18)', '2b. Loamy with no restrictions: no-till (n=5)', '2c. Loamy with no restrictions: no-till & perennial (n=13)', '2d. Loamy with no restrictions: tilled & dryfarm (n=12)', '2e. Loamy with no restrictions: no-till & dryfarm (n=6)', '3. Low OM with restrictive horizons (n=18)', '4. High OM with restrictive horizons (n=13)', '7. Shrink-swell (n=4)'), col=c('lightgoldenrod', 'tan4', 'tan4', 'tan4',  'tan4', 'tan4', 'olivedrab3', 'firebrick3', 'violetred'), ylab='', ylim=c(0,5.5), args.legend=list(x='topright', cex=1, pt.cex=1, bty='n', ncol=1, inset=c(0.02,-0.07), x.intersp=0.4), cex.names=1)
mtext(text = row.names(SOMmeans_SHRmgmt), side = 1, line=0.1, at=SOMmeans_bardims, cex=0.9)
mtext('Soil organic matter (%)', side=2, line=2.25)
segments(x0=SOMmeans_bardims[1,], y0=SOMmeans_SHRmgmt[1,] + qt(0.975,df=SOM_n_SHRmgmt[1,]) * SOM_SE_SHRmgmt[1,], x1=SOMmeans_bardims[1,], y1=SOMmeans_SHRmgmt[1,] + qt(0.025,df=SOM_n_SHRmgmt[1,]) * SOM_SE_SHRmgmt[1,], lwd = 1)
segments(x0=SOMmeans_bardims[2,], y0=SOMmeans_SHRmgmt[2,] + qt(0.975,df=SOM_n_SHRmgmt[2,]) * SOM_SE_SHRmgmt[2,], x1=SOMmeans_bardims[2,], y1=SOMmeans_SHRmgmt[2,] + qt(0.025,df=SOM_n_SHRmgmt[2,]) * SOM_SE_SHRmgmt[2,], lwd = 1)
segments(x0=SOMmeans_bardims[3,], y0=SOMmeans_SHRmgmt[3,] + qt(0.975,df=SOM_n_SHRmgmt[3,]) * SOM_SE_SHRmgmt[3,], x1=SOMmeans_bardims[3,], y1=SOMmeans_SHRmgmt[3,] + qt(0.025,df=SOM_n_SHRmgmt[3,]) * SOM_SE_SHRmgmt[3,], lwd = 1)
segments(x0=SOMmeans_bardims[4,], y0=SOMmeans_SHRmgmt[4,] + qt(0.975,df=SOM_n_SHRmgmt[4,]) * SOM_SE_SHRmgmt[4,], x1=SOMmeans_bardims[4,], y1=SOMmeans_SHRmgmt[4,] + qt(0.025,df=SOM_n_SHRmgmt[4,]) * SOM_SE_SHRmgmt[4,], lwd = 1)
segments(x0=SOMmeans_bardims[5,], y0=SOMmeans_SHRmgmt[5,] + qt(0.975,df=SOM_n_SHRmgmt[5,]) * SOM_SE_SHRmgmt[5,], x1=SOMmeans_bardims[5,], y1=SOMmeans_SHRmgmt[5,] + qt(0.025,df=SOM_n_SHRmgmt[5,]) * SOM_SE_SHRmgmt[5,], lwd = 1)
segments(x0=SOMmeans_bardims[6,], y0=SOMmeans_SHRmgmt[6,] + qt(0.975,df=SOM_n_SHRmgmt[6,]) * SOM_SE_SHRmgmt[6,], x1=SOMmeans_bardims[6,], y1=SOMmeans_SHRmgmt[6,] + qt(0.025,df=SOM_n_SHRmgmt[6,]) * SOM_SE_SHRmgmt[6,], lwd = 1)
segments(x0=SOMmeans_bardims[7,], y0=SOMmeans_SHRmgmt[7,] + qt(0.975,df=SOM_n_SHRmgmt[7,]) * SOM_SE_SHRmgmt[7,], x1=SOMmeans_bardims[7,], y1=SOMmeans_SHRmgmt[7,] + qt(0.025,df=SOM_n_SHRmgmt[7,]) * SOM_SE_SHRmgmt[7,], lwd = 1)
segments(x0=SOMmeans_bardims[8,], y0=SOMmeans_SHRmgmt[8,] + qt(0.975,df=SOM_n_SHRmgmt[8,]) * SOM_SE_SHRmgmt[8,], x1=SOMmeans_bardims[8,], y1=SOMmeans_SHRmgmt[8,] + qt(0.025,df=SOM_n_SHRmgmt[8,]) * SOM_SE_SHRmgmt[8,], lwd = 1)
segments(x0=SOMmeans_bardims[9,], y0=SOMmeans_SHRmgmt[9,] + qt(0.975,df=SOM_n_SHRmgmt[9,]) * SOM_SE_SHRmgmt[9,], x1=SOMmeans_bardims[9,], y1=SOMmeans_SHRmgmt[9,] + qt(0.025,df=SOM_n_SHRmgmt[9,]) * SOM_SE_SHRmgmt[9,], lwd = 1)
dev.off()

#make pdf or eps version
setEPS()
postscript(file = file.path(FiguresDir, 'EPS versions', 'Figure_3_FINAL.eps'), pointsize = 12, family = 'Times New Roman', width = 9, height = 4.25)
par(mar=c(2.4, 3.5, 1.25, 0.5))
SOMmeans_bardims <- barplot(SOMmeans_SHRmgmt, beside=TRUE, legend.text=c('1. Coarse with no restrictions (n=14)', '2a. Loamy with no restrictions: tilled (n=18)', '2b. Loamy with no restrictions: no-till (n=5)', '2c. Loamy with no restrictions: no-till & perennial (n=13)', '2d. Loamy with no restrictions: tilled & dryfarm (n=12)', '2e. Loamy with no restrictions: no-till & dryfarm (n=6)', '3. Low OM with restrictive horizons (n=18)', '4. High OM with restrictive horizons (n=13)', '7. Shrink-swell (n=4)'), col=c('lightgoldenrod', 'tan4', 'tan4', 'tan4',  'tan4', 'tan4', 'olivedrab3', 'firebrick3', 'violetred'), ylab='', ylim=c(0,5.5), args.legend=list(x=11.3, y=5.5, cex=1, pt.cex=1, bty='n', ncol=1, inset=c(-2,-2)), cex.names=1)
mtext(text = row.names(SOMmeans_SHRmgmt), side = 1, line=0.1, at=SOMmeans_bardims, cex=0.9)
mtext('Soil organic matter (%)', side=2, line=2.25)
segments(x0=SOMmeans_bardims[1,], y0=SOMmeans_SHRmgmt[1,] + qt(0.975,df=SOM_n_SHRmgmt[1,]) * SOM_SE_SHRmgmt[1,], x1=SOMmeans_bardims[1,], y1=SOMmeans_SHRmgmt[1,] + qt(0.025,df=SOM_n_SHRmgmt[1,]) * SOM_SE_SHRmgmt[1,], lwd = 1)
segments(x0=SOMmeans_bardims[2,], y0=SOMmeans_SHRmgmt[2,] + qt(0.975,df=SOM_n_SHRmgmt[2,]) * SOM_SE_SHRmgmt[2,], x1=SOMmeans_bardims[2,], y1=SOMmeans_SHRmgmt[2,] + qt(0.025,df=SOM_n_SHRmgmt[2,]) * SOM_SE_SHRmgmt[2,], lwd = 1)
segments(x0=SOMmeans_bardims[3,], y0=SOMmeans_SHRmgmt[3,] + qt(0.975,df=SOM_n_SHRmgmt[3,]) * SOM_SE_SHRmgmt[3,], x1=SOMmeans_bardims[3,], y1=SOMmeans_SHRmgmt[3,] + qt(0.025,df=SOM_n_SHRmgmt[3,]) * SOM_SE_SHRmgmt[3,], lwd = 1)
segments(x0=SOMmeans_bardims[4,], y0=SOMmeans_SHRmgmt[4,] + qt(0.975,df=SOM_n_SHRmgmt[4,]) * SOM_SE_SHRmgmt[4,], x1=SOMmeans_bardims[4,], y1=SOMmeans_SHRmgmt[4,] + qt(0.025,df=SOM_n_SHRmgmt[4,]) * SOM_SE_SHRmgmt[4,], lwd = 1)
segments(x0=SOMmeans_bardims[5,], y0=SOMmeans_SHRmgmt[5,] + qt(0.975,df=SOM_n_SHRmgmt[5,]) * SOM_SE_SHRmgmt[5,], x1=SOMmeans_bardims[5,], y1=SOMmeans_SHRmgmt[5,] + qt(0.025,df=SOM_n_SHRmgmt[5,]) * SOM_SE_SHRmgmt[5,], lwd = 1)
segments(x0=SOMmeans_bardims[6,], y0=SOMmeans_SHRmgmt[6,] + qt(0.975,df=SOM_n_SHRmgmt[6,]) * SOM_SE_SHRmgmt[6,], x1=SOMmeans_bardims[6,], y1=SOMmeans_SHRmgmt[6,] + qt(0.025,df=SOM_n_SHRmgmt[6,]) * SOM_SE_SHRmgmt[6,], lwd = 1)
segments(x0=SOMmeans_bardims[7,], y0=SOMmeans_SHRmgmt[7,] + qt(0.975,df=SOM_n_SHRmgmt[7,]) * SOM_SE_SHRmgmt[7,], x1=SOMmeans_bardims[7,], y1=SOMmeans_SHRmgmt[7,] + qt(0.025,df=SOM_n_SHRmgmt[7,]) * SOM_SE_SHRmgmt[7,], lwd = 1)
segments(x0=SOMmeans_bardims[8,], y0=SOMmeans_SHRmgmt[8,] + qt(0.975,df=SOM_n_SHRmgmt[8,]) * SOM_SE_SHRmgmt[8,], x1=SOMmeans_bardims[8,], y1=SOMmeans_SHRmgmt[8,] + qt(0.025,df=SOM_n_SHRmgmt[8,]) * SOM_SE_SHRmgmt[8,], lwd = 1)
segments(x0=SOMmeans_bardims[9,], y0=SOMmeans_SHRmgmt[9,] + qt(0.975,df=SOM_n_SHRmgmt[9,]) * SOM_SE_SHRmgmt[9,], x1=SOMmeans_bardims[9,], y1=SOMmeans_SHRmgmt[9,] + qt(0.025,df=SOM_n_SHRmgmt[9,]) * SOM_SE_SHRmgmt[9,], lwd = 1)
dev.off()
