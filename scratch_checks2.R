laptop <- TRUE
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
}
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
ssurgoDir <- file.path(mainDir, 'soil health/ssurgo_data')
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)

valley30cm <- read.csv(file.path(dataDir, 'valley_30cm_data.csv'), stringsAsFactors = FALSE)
acres_by_mukey <- aggregate(area_ac ~ mukey, data = valley30cm, sum)
valley30cm_by_mukey <- valley30cm[!duplicated(valley30cm$mukey), ]
valley30cm_by_mukey$area_ac <- acres_by_mukey$area_ac[match(valley30cm_by_mukey$mukey, acres_by_mukey$mukey)]
unique(valley30cm_by_mukey)
#a few checks before moving on to v2
sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$compct_om >= 80 & valley30cm_by_mukey$compct_om < 85)]) #1,073,363 acres
valley30cm_by_mukey$flag <- ifelse(valley30cm_by_mukey$compct_om >= 80, 1, 2)
sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$flag==1)]) #12,712,879
sum(valley30cm_by_mukey$area_ac) #13,873,110

# valley30cm_by_mukey <- valley30cm_by_mukey[which(valley30cm_by_mukey$compct_om >= 85), ]
# sum(is.na(valley30cm_by_mukey$aws050wta))
# sum(is.na(valley30cm_by_mukey$om_30cm))
# summary(valley30cm_by_mukey$compct_om)
# summary(valley30cm_by_mukey$mjcmp_pct)
analysis_preview <- valley30cm_by_mukey[which(valley30cm_by_mukey$flag==1), c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm', 'area_ac', 'muname', 'mjcmpnms', 'compct_om', 'dmcmp_pct', 'mukey', 'complex', 'associan')]
sum(analysis_preview$area_ac[analysis_preview$complex=='Yes'])
sum(analysis_preview$area_ac[analysis_preview$associan=='Yes'])
non_analysis_preview <- valley30cm_by_mukey[which(valley30cm_by_mukey$flag==2), c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm', 'area_ac', 'muname', 'mjcmpnms', 'compct_om', 'dmcmp_pct', 'mukey', 'complex', 'associan')]
sum(non_analysis_preview$area_ac)

analysis_preview$count_NAs <- apply(analysis_preview[,1:10], 1, function(x) sum(is.na(x)))
tapply(analysis_preview$area_ac, analysis_preview$count_NAs, sum)
apply(analysis_preview[,1:10], 2, function(x) sum(is.na(x)))
analysis_preview$ec_30cm[is.na(analysis_preview$ec_30cm) & analysis_preview$count_NAs==1] <- 0
analysis_preview$count_NAs <- apply(analysis_preview[,1:10], 1, function(x) sum(is.na(x)))
sum(analysis_preview$area_ac) #12,712,879 total acres
sum(na.omit(analysis_preview)$area_ac) #12,455,906 without EC or with corrected NA ECs;
#12,363,398 acres in v2 with EC and no SAR
#11,290,964 acres were assigned a soil health diagnostic indicators cluster class in v1
write.csv(na.omit(analysis_preview), file.path(dataDir, 'analysis_preview_10.24.19.csv'), row.names = FALSE)
analysis_preview <- na.omit(analysis_preview)

majcompnames <- unique(unlist(strsplit(analysis_preview$mjcmpnms, '-')))
majcompnames <- majcompnames[order(majcompnames)]
majcompnames <- data.frame(component_name=majcompnames, stringsAsFactors = FALSE)
write.csv(majcompnames, file.path(dataDir, 'majcompnames_unique_analysis_previes.csv'), row.names = FALSE)

problematic_compnames <- read.csv(file.path(dataDir, 'problematic_compnames_10.25.19.csv'), stringsAsFactors = FALSE)
analysis_prob_mukeys <- unlist(lapply(problematic_compnames$compname[2:nrow(problematic_compnames)], function(x) analysis_preview$mukey[grepl(x, analysis_preview$mjcmpnms)]))
sum(analysis_preview$area_ac[analysis_preview$mukey%in% analysis_prob_mukeys]) #418,271.5 acres potentially problematic
summary(analysis_preview$compct_om[analysis_preview$mukey%in% analysis_prob_mukeys])

analysis_rockOC_mukeys <- analysis_preview$mukey[grepl('Rock outcrop', analysis_preview$mjcmpnms)]
analysis_rockOC_mukeys
analysis_rock_comp_data <- comp_data[comp_data$mukey %in% analysis_rockOC_mukeys & comp_data$compname== 'Rock outcrop', ]
horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]

non_analysis_preview$count_NAs <- apply(non_analysis_preview[,1:10], 1, function(x) sum(is.na(x)))
tapply(non_analysis_preview$area_ac, non_analysis_preview$count_NAs, sum) #794,356.55 with no NAs; 796,308.86 after correcting NA EC with no other NAs
non_analysis_preview$ec_30cm[is.na(non_analysis_preview$ec_30cm) & non_analysis_preview$count_NAs==1] <- 0
apply(non_analysis_preview[,1:10], 2, function(x) sum(is.na(x)))
unique(non_analysis_preview$muname[non_analysis_preview$count_NAs==0])
unique(non_analysis_preview$mjcmpnms[non_analysis_preview$count_NAs==0])
hist(non_analysis_preview$compct_om[non_analysis_preview$count_NAs==0])
hist(non_analysis_preview$dmcmp_pct[non_analysis_preview$count_NAs==0])
sum(non_analysis_preview$area_ac[non_analysis_preview$count_NAs==0][non_analysis_preview$compct_om[non_analysis_preview$count_NAs==0] >= non_analysis_preview$dmcmp_pct[non_analysis_preview$count_NAs==0]]) #639,712 acres before correcting EC NAs with no other NA data; 640,556 after correcting EC NAs with no other NA data
imperfect_keep_preview <- non_analysis_preview[non_analysis_preview$count_NAs==0 & non_analysis_preview$compct_om >= non_analysis_preview$dmcmp_pct, ]
sum(imperfect_keep_preview$area_ac) #640,556 acres
sum(imperfect_keep_preview$area_ac[imperfect_keep_preview$complex=='Yes']) #300,081 acres are a complex
sum(imperfect_keep_preview$area_ac[imperfect_keep_preview$associan=='Yes']) #22825.48

# write.csv(imperfect_keep_preview, file.path(dataDir, 'imperfect_SSURGO_preview_10.22.19.csv'))


analysis_preview <- read.csv(file.path(dataDir, 'analysis_preview_10.23.19.csv'))
colnames(analysis_preview)
sum(grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms))
mjcmp_rockOCs <- valley30cm_by_mukey$mukey[grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms)]
mjcmp_rockOCs

sum(mjcmp_rockOCs %in% analysis_preview$mukey)
analysis_rockOC_mukeys <- mjcmp_rockOCs[mjcmp_rockOCs %in% analysis_preview$mukey]
analysis_rock_comp_data <- comp_data[comp_data$mukey %in% analysis_rockOC_mukeys & comp_data$compname== 'Rock outcrop', ]
horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
horizon_data_rockOC <- horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
colnames(horizon_data_rockOC)
lapply(horizon_data_rockOC[ ,5:25], summary)


valley30cm_by_mukey[valley30cm_by_mukey$mukey==462578,]
comp_data[comp_data$mukey==462578,]
horizon_data[horizon_data$cokey==16997325,]
horizon_data[horizon_data$cokey==16997326,]
horizon_data[horizon_data$cokey==16997327,]
horizon_data[horizon_data$cokey==16997328,] #minor component
#Ks check
9*0.5+7.21*0.5
9*0.5/0.95 + 8.105*0.3/0.95 + 0.005*10/30*0.15/0.95

#OM check
1.5*15/31+0.25*16/31 #adding 1 extra cm in calculation
1.25*0.5/0.8 +  0.8548387*0.3/0.8

#awc check[getting 3.785263]
0.13*31*0.5/0.95 + 0.17*31*0.3/0.95 #4.35 but getting 3.75 because rock outcrop is counting as a 0!

#simple check
head(valley30cm_by_mukey[valley30cm_by_mukey$mjcps_no==1, ])
analysis_preview[analysis_preview$mukey==461027, ]
comp_data[comp_data$mukey==461027, ]
horizon_data[horizon_data$cokey==16777370, ]
#om calc
3*30/31 + 0.5*1/31 #considers 1 extra cm for some reason
#awc calc
0.21*30 + 0.16*1 #6.46 is answer
#simple check take 2
analysis_preview[analysis_preview$mukey==461063, ]
comp_data[comp_data$mukey==461063, ]
horizon_data[horizon_data$cokey==16776831, ]
#awc_30cm=4.15 (0.13 from 0-25; 0.15 from 25-94)
0.13*25+0.15*6
