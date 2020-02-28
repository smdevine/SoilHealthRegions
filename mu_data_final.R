laptop <- FALSE
library(vioplot)
library(raster)
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
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'))
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
valley30cm_by_mukey$dom_orders <- apply(valley30cm_by_mukey[ ,colnames(valley30cm_by_mukey)[grepl('sols_pct', colnames(valley30cm_by_mukey))]], 1, function(x) sum(max(x) == x))
sum(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$dom_orders==1]) #12,857,758
table(valley30cm_by_mukey$dom_orders)

valley30cm_by_mukey$dom_order <- apply(valley30cm_by_mukey[ ,colnames(valley30cm_by_mukey)[grepl('sols_pct', colnames(valley30cm_by_mukey))]], 1, function(x) c('Ultisols', 'Vertisols', 'Mollisols',  'Inceptisols', 'Entisols', 'Aridisols', 'Andisols', 'Alfisols')[which.max(x)])
table(valley30cm_by_mukey$dom_order)
valley30cm_by_mukey$dom_order[valley30cm_by_mukey$dom_orders>1] <- NA
tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$dom_order, function(x) round(sum(x), 0))
tapply(valley30cm_by_mukey$om_30cm, valley30cm_by_mukey$dom_order, function(x) mean(x, na.rm = TRUE))

#TO-DO
#dominant soil order summary
DomSoilOrder_summary <- as.data.frame(lapply(c('Alfisols', 'Andisols', 'Aridisols', 'Entisols', 'Inceptisols', 'Mollisols', 'Ultisols', 'Vertisols'), function(x) {
  tapply(valley30cm_by_mukey[[paste0(x, '_ac')]], clus_7_names[match(valley30cm_by_mukey$cluster_7, 1:7)], sum)
}
), row.names = clus_7_names[order(clus_7_names)], col.names = c('Alfisols', 'Andisols', 'Aridisols', 'Entisols', 'Inceptisols', 'Mollisols', 'Ultisols', 'Vertisols'))
SoilOrder_summary
sum(valley30cm_by_mukey$area_ac) - sum(SoilOrder_summary) #95225.92 acres off
SoilOrder_summary$TOTAL <- tapply(valley30cm_by_mukey$area_ac, clus_7_names[match(valley30cm_by_mukey$cluster_7, 1:7)], sum)

write.csv(SoilOrder_summary, file.path(dataDir, 'soil survey facts', 'DominantSoilOrders_by_SHR7.csv'), row.names=TRUE)

#make a violin plot by soil order
#order: 
clus_7_colors <- c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3')
color_by_soil_order <- c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'black', 'lightblue1', 'deepskyblue', 'violetred')
col2rgb(color_by_soil_order)
soil_order_logic <- c('Entisols', 'Mollisols', 'Alfisols', 'Ultisols', 'Andisols', 'Inceptisols', 'Aridisols', 'Vertisols')
cbind(color_by_soil_order, soil_order_logic, t(col2rgb(color_by_soil_order)))
vioplot_mod_SoilOrder_validation <- function(df, varname, ylim_vioplot, area_fact, labnames, ylab, fname, mar) {
  tiff(file = file.path(FiguresDir, 'v2', 'soil order', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][which(df$dom_order=='Entisols')], times=round(df$area_ac[which(df$dom_order=='Entisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Mollisols')], times=round(df$area_ac[which(df$dom_order=='Mollisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Alfisols')], times=round(df$area_ac[which(df$dom_order=='Alfisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Ultisols')], times=round(df$area_ac[which(df$dom_order=='Ultisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Andisols')], times=round(df$area_ac[which(df$dom_order=='Andisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Inceptisols')], times=round(df$area_ac[which(df$dom_order=='Inceptisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Aridisols')], times=round(df$area_ac[which(df$dom_order=='Aridisols')]/area_fact, 0)), rep(df[[varname]][which(df$dom_order=='Vertisols')], times=round(df$area_ac[which(df$dom_order=='Vertisols')]/area_fact, 0)), col = color_by_soil_order, rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab) #col
  mtext('Soil orders (USDA-NRCS soil taxonomy)', side = 1, line = 2.25)
  # points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.1, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.7, pch=1, col='black')
  # kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  # kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  # points(x=kssl_means$xdim_vioplot-0.1, y=kssl_means$mean, pch=8, cex=0.9, col='orange')
  # points(x=cdfa_pts$xdim_vioplot7[!is.na(cdfa_pts[[cdfa_varname]])]+0.1, y=cdfa_pts[[cdfa_varname]][!is.na(cdfa_pts[[cdfa_varname]])], cex=0.7, pch=4, col='black')
  # cdfa_means <- data.frame(mean=tapply(cdfa_pts[[cdfa_varname]], cdfa_pts$cluster_7, mean, na.rm=TRUE))
  # cdfa_means$xdim_vioplot <- match(row.names(cdfa_means), plot_order)
  # points(x=cdfa_means$xdim_vioplot+0.1, y=cdfa_means$mean, pch=8, cex=0.9, col='darkblue')
  # text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5)
  # if(legend_plot) {
    # legend(x=legendloc, legend=c('SSURGO violin plots, area-weighted', 'KSSL data', 'KSSL mean', 'CDFA data', 'CDFA mean'), pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), cex=legend_cex)
  # }
  dev.off()
}
vioplot_mod_SoilOrder_validation(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(-0.2,12), area_fact = 10, ylab='Organic matter (%)', fname='SoilOrder_om_vioplots_validation.tif', mar=c(3.5, 4.25, 1, 1))
#additional function arguments kssl_df, kssl_varname, cdfa_pts, cdfa_varname, legend_plot, legendloc, legend_cex, sig_labels
summary(valley30cm_by_mukey$om_30cm[valley30cm_by_mukey$dom_order=='Aridisols'])
valley30cm_by_mukey[which(valley30cm_by_mukey$om_30cm>10 & valley30cm_by_mukey$dom_order=='Aridisols'),]


#add Order info to shapefile
length(unique(valley30cm_by_mukey$mukey))
length(unique(valley_mu_shp_30cm$mukey))
valley_mu_shp_30cm <- valley_mu_shp_30cm[valley_mu_shp_30cm$mukey %in% valley30cm_by_mukey$mukey, ]
length(unique(valley_mu_shp_30cm$mukey))
sum(valley_mu_shp_30cm$area_ac)
sum(valley30cm_by_mukey$area_ac)
unique(valley30cm_by_mukey$dom_order)
valley_mu_shp_30cm$dom_order <- valley30cm_by_mukey$dom_order[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
tapply(valley_mu_shp_30cm$area_ac, valley_mu_shp_30cm$dom_order, function(x) round(sum(x), 0))
tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$dom_order, function(x) round(sum(x), 0))
shapefile(valley_mu_shp_30cm, file.path(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster_SoilOrder.shp')))

