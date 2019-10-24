laptop <- TRUE
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  
}
analysis_preview <- read.csv(file.path(dataDir, 'analysis_preview_10.23.19.csv'))
dim(valley30cm_by_mukey)
sum(grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms))
mjcmp_rockOCs <- valley30cm_by_mukey$mukey[grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms)]
mjcmp_rockOCs

sum(mjcmp_rockOCs %in% analysis_preview$mukey)
analysis_rockOC_mukeys <- mjcmp_rockOCs[mjcmp_rockOCs %in% analysis_preview$mukey]
analysis_rock_comp_data <- comp_data[comp_data$mukey %in% analysis_rockOC_mukeys & comp_data$compname== 'Rock outcrop', ]
horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
valley30cm_by_mukey[valley30cm_by_mukey$mukey %in% analysis_rockOC_mukeys, ]
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

#simple check
head(valley30cm_by_mukey[valley30cm_by_mukey$mjcps_no==1, ])
analysis_preview[analysis_preview$mukey==461027, ]
comp_data[comp_data$mukey==461027, ]
horizon_data[horizon_data$cokey==16777370, ]
3*30/31 + 0.5*1/31 #considers 1 extra cm for some reason


comp_data_rockOC <- comp_data[comp_data$majcompflag=='Yes' & comp_data$compname== 'Rock outcrop',]
horizon_data_rockOC <- horizon_data[horizon_data$cokey %in% comp_data_rockOC$cokey, ]
colnames(horizon_data_rockOC)
lapply(horizon_data_rockOC[ ,5:25], summary)
lapply()
