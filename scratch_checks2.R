
dim(valley30cm_by_mukey)
sum(grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms))
mjcmp_rockOCs <- valley30cm_by_mukey$mukey[grepl('Rock outcrop', valley30cm_by_mukey$mjcmpnms)]
mjcmp_rockOCs

sum(mjcmp_rockOCs %in% analysis_preview$mukey)
analysis_rockOC_mukeys <- mjcmp_rockOCs[mjcmp_rockOCs %in% analysis_preview$mukey]
analysis_rock_comp_data <- comp_data[comp_data$mukey %in% analysis_rockOC_mukeys & comp_data$compname== 'Rock outcrop', ]
horizon_data[horizon_data$cokey %in% analysis_rock_comp_data$cokey, ]
valley30cm_by_mukey[valley30cm_by_mukey$mukey %in% analysis_rockOC_mukeys, ]

comp_data[comp_data$mukey==462578,]
horizon_data[horizon_data$cokey==16997325,]
horizon_data[horizon_data$cokey==16997326,]
horizon_data[horizon_data$cokey==16997327,]

comp_data_rockOC <- comp_data[comp_data$majcompflag=='Yes' & comp_data$compname== 'Rock outcrop',]
horizon_data_rockOC <- horizon_data[horizon_data$cokey %in% comp_data_rockOC$cokey, ]
colnames(horizon_data_rockOC)
lapply(horizon_data_rockOC[ ,5:25], summary)
lapply()
