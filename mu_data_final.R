laptop <- TRUE
library(vioplot)
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
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}
mar_settings <- c(4, 4.5, 1, 1)
om_to_oc <- 1.72
clus_7_names <- c('3. Coarse w/pans', '6. Fine saline-sodic', '5. Coarse saline-sodic', '1. Coarse w/no restrictions', '7. Fine shrink-swell', '2. Loamy w/no restrictions', '4. Loamy w/pans')
#produced in ssurgo_calag_cluster_v2.R
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'v2 results', 'valley30cm_by_mukey_cluster_v2.csv'), stringsAsFactors = FALSE)
colnames(valley30cm_by_mukey)
sum(valley30cm_by_mukey$area_ac)
unique(unlist(strsplit(unique(valley30cm_by_mukey$txorders), '-'))) #these orders are represented: "Alfisols" "Inceptisols" "Mollisols"  "Entisols"    "Vertisols" Ultisols" "Aridisols" "Andisols"
sum(valley30cm_by_mukey$area_ac[is.na(valley30cm_by_mukey$txorders)]) #only 3391.2 acres missing soil order

comp_data <- read.csv(file.path(mainDir, 'soil health/ssurgo_data/component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' ')) #treat blanks or a space as a NA)
colnames(comp_data)

#fix a few majcompflag errors as in original data prep
comp_data$majcompflag[comp_data$majcompflag=='No ' & comp_data$comppct_r>=15 & !is.na(comp_data$castorieindex)] <- 'Yes'
comp_data$majcompflag[comp_data$majcompflag=='Yes' & comp_data$comppct_r < 15] <- 'No '
comp_data[comp_data$majcompflag=='Yes' & is.na(comp_data$taxorder) & comp_data$mukey %in% valley30cm_by_mukey$mukey,]

SoilOrderArea_calc <- function(df_mu, df_comp, SoilOrder) {
  SoilOrder_by_mukey <- data.frame(comppct=tapply(df_comp$comppct_r[df_comp$taxorder==SoilOrder & df_comp$majcompflag=='Yes'], df_comp$mukey[df_comp$taxorder==SoilOrder & df_comp$majcompflag=='Yes'], sum))
  df_mu[[paste0(SoilOrder, '_pct')]] <- SoilOrder_by_mukey$comppct[match(df_mu$mukey, row.names(SoilOrder_by_mukey))]
  df_mu[[paste0(SoilOrder, '_pct')]][is.na(df_mu[[paste0(SoilOrder, '_pct')]])] <- 0
  df_mu[[paste0(SoilOrder, '_ac')]] <- df_mu$area_ac * (df_mu[[paste0(SoilOrder, '_pct')]] / df_mu$mjcmp_pct)
  result <- tapply(df_mu[[paste0(SoilOrder, '_ac')]], clus_7_names[match(df_mu$cluster_7, 1:7)], sum)
  print(result)
  df_mu
}
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Ultisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Vertisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Mollisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Inceptisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Entisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Aridisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Andisols')
valley30cm_by_mukey <- SoilOrderArea_calc(valley30cm_by_mukey, comp_data, 'Alfisols')

SoilOrder_summary <- as.data.frame(lapply(c('Alfisols', 'Andisols', 'Aridisols', 'Entisols', 'Inceptisols', 'Mollisols', 'Ultisols', 'Vertisols'), function(x) {
  tapply(valley30cm_by_mukey[[paste0(x, '_ac')]], clus_7_names[match(valley30cm_by_mukey$cluster_7, 1:7)], sum)
  }
), row.names = clus_7_names[order(clus_7_names)], col.names = c('Alfisols', 'Andisols', 'Aridisols', 'Entisols', 'Inceptisols', 'Mollisols', 'Ultisols', 'Vertisols'))
SoilOrder_summary
sum(valley30cm_by_mukey$area_ac) - sum(SoilOrder_summary) #95225.92 acres off
SoilOrder_summary$TOTAL <- tapply(valley30cm_by_mukey$area_ac, clus_7_names[match(valley30cm_by_mukey$cluster_7, 1:7)], sum)

write.csv(SoilOrder_summary, file.path(dataDir, 'soil survey facts', 'SoilOrders_by_SHR7.csv'), row.names=TRUE)


test <- valley30cm_by_mukey$area_ac - apply(valley30cm_by_mukey[,grepl('s_ac', colnames(valley30cm_by_mukey))], 1, sum)
which(test > 1)
valley30cm_by_mukey[5,]
comp_data[comp_data$mukey==461103,]

valley30cm_by_mukey[71,]
comp_data[comp_data$mukey==459435,]

valley30cm_by_mukey[80,]
comp_data[comp_data$mukey==459458,]

#calculate area
length(unique(valley30cm_by_mukey$muname)) #3942 unique map unit names
colnames(valley30cm_by_mukey)[grepl('sols_pct', colnames(valley30cm_by_mukey))]
test <- apply(valley30cm_by_mukey[ ,colnames(valley30cm_by_mukey)[grepl('sols_pct', colnames(valley30cm_by_mukey))]], 1, function(x) which.max(x))
table(test)
valley30cm_by_mukey$dom_orders <- apply(valley30cm_by_mukey[ ,colnames(valley30cm_by_mukey)[grepl('sols_pct', colnames(valley30cm_by_mukey))]], 1, function(x) sum(max(x) == x))
sum(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$dom_orders==1]) #12,857,758

