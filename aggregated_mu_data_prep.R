#rule is 80% of a map unit must have data with no NAs (with NA ECs and no other NAs converted to 0) --or-- data coverage for at the dominant component (given by % with data with no NAs with same EC rule)
#rock outcrop with Ksat > 0 converted to NA so as not to affect map-unit aggregated data
#remaining question is whether or not CEC should be used in cluster analysis; would add 131,608 acres to analysis if excluded
laptop <- TRUE
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
ssurgoDir <- file.path(mainDir, 'soil health/ssurgo_data')
min_modified <- function(x) {
  if(all(is.na(x))) {
    return(NA)
  }
  else {min(x, na.rm = TRUE)}
}
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
valley30cm <- read.csv(file.path(dataDir, 'valley_30cm_data_2.7.20.csv'), stringsAsFactors = FALSE) #'valley_30cm_data_10.29.19.csv'
acres_by_mukey <- aggregate(area_ac ~ mukey, data = valley30cm, sum)
valley30cm_by_mukey <- valley30cm[!duplicated(valley30cm$mukey), ]
valley30cm_by_mukey$area_ac <- acres_by_mukey$area_ac[match(valley30cm_by_mukey$mukey, acres_by_mukey$mukey)]

#new checks based on meeting minimum data requirement for at least the dominant component percentage
c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm', 'area_ac', 'muname', 'mjcmpnms', 'compct_om', 'compct_cec', 'compct_ksat', 'compct_awc', 'compct_clay', 'compct_bd', 'compct_pH', 'compct_lep', 'dmcmp_pct', 'mukey', 'complex', 'associan') %in% colnames(valley30cm_by_mukey)
analysis_dataset <- valley30cm_by_mukey[ ,c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm', 'area_ac', 'muname', 'mjcmpnms', 'compct_om', 'compct_cec', 'compct_ksat', 'compct_awc', 'compct_clay', 'compct_bd', 'compct_pH', 'compct_lep', 'dmcmp_pct', 'mukey', 'complex', 'associan')] #don't include compct_ec because of rule below
analysis_dataset$count_NAs <- apply(analysis_dataset[,1:10], 1, function(x) sum(is.na(x)))
tapply(analysis_dataset$area_ac, analysis_dataset$count_NAs, sum) #0 NAs are 13,134,775 acres
analysis_dataset$ec_30cm[is.na(analysis_dataset$ec_30cm) & analysis_dataset$count_NAs==1] <- 0
analysis_dataset$count_NAs <- apply(analysis_dataset[,1:10], 1, function(x) sum(is.na(x)))
tapply(analysis_dataset$area_ac, analysis_dataset$count_NAs, sum) #with EC correction, now 13,229,694

colnames(analysis_dataset)
analysis_dataset$min_dat_cov <- apply(analysis_dataset[,which(colnames(analysis_dataset)=='compct_om'):which(colnames(analysis_dataset)=='compct_lep')], 1, min_modified)
analysis_dataset$dat_80pct <- ifelse(analysis_dataset$min_dat_cov >= 80, "Yes", "No")
summary(as.factor(analysis_dataset$dat_80pct)) #211 NAs if EC data is ignored
sum(analysis_dataset$area_ac[(analysis_dataset$dat_80pct=='Yes' & analysis_dataset$count_NAs==0)], na.rm = TRUE) #12,307,706 acres highest quality data
sum(analysis_dataset$area_ac[(analysis_dataset$dat_80pct=='Yes' & analysis_dataset$count_NAs==0) | (analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct & analysis_dataset$count_NAs==0)], na.rm = TRUE) #13,034,096 with additional rule
sum(analysis_dataset$area_ac[analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct & analysis_dataset$count_NAs==0 & analysis_dataset$min_dat_cov < 80], na.rm = TRUE) #because 726,390 meet the lower QC criteria
analysis_dataset_final <- analysis_dataset[which((analysis_dataset$dat_80pct=='Yes' & analysis_dataset$count_NAs==0) | (analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct & analysis_dataset$count_NAs==0)), ]
sum(analysis_dataset_final$area_ac) #check 13,034,096 acres

#investigate map units with only 1 NA using unedited version
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$clay_30cm)] & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct) #0
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$om_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #1064.5
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$cec_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #131607.7
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$bd_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #0
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$pH_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #8413
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$lep_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #0
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$ksat_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #0
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$awc_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct]) #11518
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==1 & is.na(analysis_dataset$ksat_30cm)]) #0
sum(analysis_dataset$area_ac[analysis_dataset$count_NAs==2 & is.na(analysis_dataset$ksat_30cm)]) #0
sum(analysis_dataset$area_ac[is.na(analysis_dataset$ksat_30cm)])
investigate_CEC <- analysis_dataset[which(analysis_dataset$count_NAs==1 & is.na(analysis_dataset$cec_30cm) & analysis_dataset$min_dat_cov >= analysis_dataset$dmcmp_pct), c('mukey', 'area_ac', 'muname', 'mjcmpnms')]

write.csv(investigate_CEC, file.path(dataDir, 'mapunits_with_NA_CEC.csv'), row.names = FALSE)

compkeys_invetigate_CEC <- comp_data$cokey[comp_data$majcompflag=='Yes' & comp_data$mukey %in% investigate_CEC$mukey]
investigate_CEC_horizons <- horizon_data[horizon_data$cokey %in% compkeys_invetigate_CEC, ]
investigate_CEC_horizons$mukey <- comp_data$mukey[match(investigate_CEC_horizons$cokey, comp_data$cokey)]
investigate_CEC_horizons <- investigate_CEC_horizons[,c(ncol(investigate_CEC_horizons), 1:(ncol(investigate_CEC_horizons)-1))]
write.csv(investigate_CEC_horizons, file.path(dataDir, 'horizons_with_NA_CEC.csv'), row.names = FALSE)

#add up complex and associations


#see scratch_checks_v2.R for previous checks
valley30cm_by_mukey[valley30cm_by_mukey$mukey==462578,]
comp_data[comp_data$mukey==462578,]
horizon_data[horizon_data$cokey==16997325,] #50%
horizon_data[horizon_data$cokey==16997326,] #30%
horizon_data[horizon_data$cokey==16997327,] #15%
horizon_data[horizon_data$cokey==16997328,] #minor component
#Ks check [should be 7.297105]
9*0.5+7.21*0.5 #8.105 for 16997326
9*0.5/0.95 + 8.105*0.3/0.95 + 0.005*0.15/0.95

#OM check [should be 1.109375]
1.5*0.5+0.25*0.5 #adding 1 extra cm in calculation
1.25*0.5/0.8 +  0.875*0.3/0.8

#awc check[should be 4.35]
0.13*30*0.5/0.8 + 0.17*30*0.3/0.8

#simple check
head(valley30cm_by_mukey[valley30cm_by_mukey$mjcps_no==1, ])
analysis_dataset_final[analysis_dataset_final$mukey==461027, ]
comp_data[comp_data$mukey==461027, ]
horizon_data[horizon_data$cokey==16777370, ]
#om calc
3 #considers 1 extra cm for some reason
#awc calc
0.21*30 #6.3 is answer
#simple check take 2
analysis_dataset_final[analysis_dataset_final$mukey==461063, ]
comp_data[comp_data$mukey==461063, ]
horizon_data[horizon_data$cokey==16776831, ]
#awc_30cm=4 (0.13 from 0-25; 0.15 from 25-94)
0.13*25+0.15*5

#rock OC check
mjcmp_rockOCs <- valley30cm_by_mukey$mukey[grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms)]
mjcmp_rockOCs

sum(mjcmp_rockOCs %in% analysis_dataset_final$mukey)
sum(analysis_dataset_final$area_ac[analysis_dataset_final$mukey %in% mjcmp_rockOCs]) #48,784 acres
analysis_rockOC_mukeys <- mjcmp_rockOCs[mjcmp_rockOCs %in% analysis_dataset_final$mukey]
analysis_rock_comp_data <- comp_data[comp_data$mukey %in% analysis_rockOC_mukeys & comp_data$compname== 'Rock outcrop', ]
horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
horizon_data_rockOC <- horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
dim(horizon_data_rockOC)
colnames(horizon_data_rockOC)
lapply(horizon_data_rockOC[ ,5:25], summary)
horizon_data_rockOC[,c('hzdept_r', 'hzdepb_r', 'ksat_r')]
rockOC_cokey_prob <- horizon_data_rockOC$cokey[which(horizon_data_rockOC$ksat_r > 0)]
rockOC_mukey_prob <- analysis_rock_comp_data$mukey[match(rockOC_cokey_prob, analysis_rock_comp_data$cokey)]
sum(analysis_dataset_final$area_ac[analysis_dataset_final$mukey %in% rockOC_mukey_prob]) #3897 acres affected by crappy rock OC Ksat data
analysis_dataset_final[analysis_dataset_final$mukey %in% rockOC_mukey_prob, ]
analysis_dataset_final[analysis_dataset_final$mukey==459994,]
comp_data[comp_data$mukey==459994,]
horizon_data[horizon_data$cokey==16453073,] #Gaviota
horizon_data[horizon_data$cokey==16453074,] #Rock outcrop with ksat of 71!  No longer affects final data


#write revised  mu aggregated data to file
valley30cm_by_mukey_final <- valley30cm_by_mukey[valley30cm_by_mukey$mukey %in% analysis_dataset_final$mukey, ]
sum(valley30cm_by_mukey_final$area_ac) #13034096 acres
sum(is.na(valley30cm_by_mukey_final$ec_30cm))
valley30cm_by_mukey_final$ec_30cm[is.na(valley30cm_by_mukey_final$ec_30cm)] <- 0
write.csv(valley30cm_by_mukey_final, file.path(dataDir, 'for cluster analysis', 'valley30cm_by_mukey_final_v2.csv'), row.names = FALSE) #only change in v2 is updated kg SOC m^-2 data
