laptop <- TRUE
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
unique(unlist(strsplit(valley30cm_by_mukey$mjcmpnms, '-')))
unique_majcomps_by_SHR <- tapply(valley30cm_by_mukey$mjcmpnms, valley30cm_by_mukey$cluster_7, function(x) unique(unlist(strsplit(x, '-'))))
names(unique_majcomps_by_SHR) <- clus_7_names
unique_majcomps_by_SHR <- unique_majcomps_by_SHR[order(clus_7_names)]
sapply(unique_majcomps_by_SHR, length)
mean(sapply(unique_majcomps_by_SHR, length))

valley30cm_by_mukey$clus7_name <- clus_7_names[valley30cm_by_mukey$cluster_7]

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
sumDomOrderArea <- function(x, y) {
  data.frame(Alfisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Alfisols')]), Andisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Andisols')]), Aridisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Aridisols')]), Entisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Entisols')]), Inceptisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Inceptisols')]), Mollisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Mollisols')]), Ultisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Ultisols')]), Vertisols=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_order=='Vertisols')]), CoDominantOrders=sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$cluster_7==x & valley30cm_by_mukey$dom_orders>1)]), row.names = y)
}
sumDomOrderArea(x=1, y=clus_7_names[1])
DomSoilOrder_summary <- rbind(sumDomOrderArea(1, clus_7_names[1]), sumDomOrderArea(2, clus_7_names[2]), sumDomOrderArea(3, clus_7_names[3]), sumDomOrderArea(4, clus_7_names[4]), sumDomOrderArea(5, clus_7_names[5]), sumDomOrderArea(6, clus_7_names[6]), sumDomOrderArea(7, clus_7_names[7]))[order(clus_7_names),]
DomSoilOrder_summary$TOTAL <- tapply(valley30cm_by_mukey$area_ac, clus_7_names[match(valley30cm_by_mukey$cluster_7, 1:7)], sum)

write.csv(DomSoilOrder_summary, file.path(dataDir, 'soil survey facts', 'DominantSoilOrders_by_SHR7.csv'), row.names=TRUE)
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==1], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==1], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==2], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==2], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==3], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==3], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==4], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==4], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==5], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==5], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==6], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==6], function(x) as.integer(sum(x)))
tapply(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_7==7], valley30cm_by_mukey$dom_order[valley30cm_by_mukey$cluster_7==7], function(x) as.integer(sum(x)))


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


#get overall area summary by name
colnames(comp_data)
compname_area <- comp_data[comp_data$mukey %in% valley30cm_by_mukey$mukey & comp_data$majcompflag=='Yes', c('mukey', 'compname', 'comppct_r', 'taxorder')]
compname_area$mu_area_ac <- valley30cm_by_mukey$area_ac[match(compname_area$mukey, valley30cm_by_mukey$mukey)] 
compname_area$mjcmp_pct <- valley30cm_by_mukey$mjcmp_pct[match(compname_area$mukey, valley30cm_by_mukey$mukey)]
compname_area$comp_area_ha <- compname_area$mu_area_ac * (compname_area$comppct_r / compname_area$mjcmp_pct) / 2.47105
area_by_compname <- data.frame(compname=row.names(tapply(compname_area$comp_area_ha, compname_area$compname, sum)), area_ha=round(tapply(compname_area$comp_area_ha, compname_area$compname, sum), 1))
sum(area_by_compname$area_ha) #5274719
sum(valley30cm_by_mukey$area_ac/2.47105) #5274720
dim(area_by_compname)
length(unique(area_by_compname$compname))
area_by_compname$order <- compname_area$taxorder[match(area_by_compname$compname, compname_area$compname)]
table(area_by_compname$order)
sum(is.na(area_by_compname$order)) #13 NA
area_by_compname$compname[is.na(area_by_compname$order)]
unique(comp_data$taxorder[comp_data$compname=='Talus'])
unique(comp_data$taxorder[comp_data$compname=='Riverwash'])
write.csv(area_by_compname, file.path(dataDir, 'soil survey facts', 'area_taxorder_by_compname.csv'), row.names = FALSE)

#compname calc
#df_comp will be valley30cm_by_mukey
#df_comp2 will be comp_data
CompnameArea_calc <- function(df_comp, df_comp2, SHR_name, fname) {
  print(length(unique_majcomps_by_SHR[[match(SHR_name, clus_7_names[order(clus_7_names)])]]))
  mukeys <- df_comp$mukey[df_comp$clus7_name==SHR_name]
  df_comp3 <- df_comp2[df_comp2$mukey %in% mukeys & df_comp2$majcomp=='Yes',] #majcomp flags fixed above
  df_comp3$mu_area_ac <- df_comp$area_ac[match(df_comp3$mukey, df_comp$mukey)]
  df_comp3$mjcmp_pct <- df_comp$mjcmp_pct[match(df_comp3$mukey, df_comp$mukey)]
  df_comp3$comp_area_ha <- df_comp3$mu_area_ac * (df_comp3$comppct_r / df_comp3$mjcmp_pct) / 2.47105
  print(length(unique(df_comp3$compname)))
  print(as.integer(sum(df_comp$area_ac[df_comp$clus7_name==SHR_name]) / 2.47105))
  area_by_MajCompnames <- data.frame(compname=row.names(tapply(df_comp3$comp_area_ha, df_comp3$compname, sum)), area_ha=round(tapply(df_comp3$comp_area_ha, df_comp3$compname, sum), 1), row.names = seq_along(unique_majcomps_by_SHR[[match(SHR_name, clus_7_names[order(clus_7_names)])]]))
  print(as.integer(sum(area_by_MajCompnames$area_ha)))
  write.csv(area_by_MajCompnames, file.path(dataDir, 'soil survey facts', fname), row.names=FALSE)
  area_by_MajCompnames
}
clus_7_names[order(clus_7_names)]
coarse_w_no_res <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "1. Coarse w/no restrictions", fname='coarse_w_no_res_MajCompnames_area.csv')
loamy_w_no_res <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "2. Loamy w/no restrictions", fname='loamy_w_no_res_MajCompnames_area.csv')
coarse_w_pans <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "3. Coarse w/pans", fname='coarse_w_pans_MajCompnames_area.csv')
loamy_w_pans <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "4. Loamy w/pans", fname='loamy_w_pans_MajCompnames_area.csv')
coarse_saline_sodic <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "5. Coarse saline-sodic", fname='coarse_saline_sodic_MajCompnames_area.csv')
fine_saline_sodic <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "6. Fine saline-sodic", fname='fine_saline_sodic_MajCompnames_area.csv')
fine_shrink_swell <- CompnameArea_calc(df_comp = valley30cm_by_mukey, df_comp2 = comp_data, SHR_name = "7. Fine shrink-swell", fname='fine_shrink_swell_MajCompnames_area.csv')

majcompnames_area_summary <- list(coarse_w_no_res, loamy_w_no_res, coarse_w_pans, loamy_w_pans, coarse_saline_sodic, fine_saline_sodic, fine_shrink_swell)
names(majcompnames_area_summary) <- clus_7_names[order(clus_7_names)]
top30_compnames_by_SHR <- do.call(cbind, lapply(majcompnames_area_summary, function(x) {head(x[order(x$area_ha, decreasing = TRUE), ], 30)}))
write.csv(top30_compnames_by_SHR, file.path(dataDir, 'soil survey facts', 'top30_compnames_by_SHR.csv'), row.names = FALSE)
lapply(majcompnames_area_summary, function(x) {'Hanford' %in% x$compname}) #in coarse w/pans and coarse w/restrictions
lapply(majcompnames_area_summary, function(x) {'Panoche' %in% x$compname}) #in coarse w/pans and coarse w/restrictions
lapply(majcompnames_area_summary, function(x) {'San Joaquin' %in% x$compname})
lapply(majcompnames_area_summary, function(x) {'Flamen' %in% x$compname})
lapply(majcompnames_area_summary, function(x) {'Fresno' %in% x$compname})
lapply(majcompnames_area_summary, function(x) {'Calflax' %in% x$compname})
lapply(majcompnames_area_summary, function(x) {'Capay' %in% x$compname})
valley30cm_by_mukey$mukey[grepl('Capay', valley30cm_by_mukey$mjcmpnms) & valley30cm_by_mukey$clus7_name=='2. Loamy w/no restrictions']
comp_data[comp_data$mukey %in% valley30cm_by_mukey$mukey[grepl('Capay', valley30cm_by_mukey$mjcmpnms) & valley30cm_by_mukey$clus7_name=='2. Loamy w/no restrictions'], ]
comp_data[comp_data$mukey %in% valley30cm_by_mukey$mukey[grepl('Capay', valley30cm_by_mukey$mjcmpnms)], ]
valley30cm_by_mukey[grepl('Capay', valley30cm_by_mukey$mjcmpnms) & valley30cm_by_mukey$clus7_name=='2. Loamy w/no restrictions', ]

valley30cm_by_mukey[grepl('Capay', valley30cm_by_mukey$mjcmpnms), c(56:73, 103)]

sum(fine_shrink_swell$compname %in% fine_saline_sodic$compname) #24
fine_conflicts_majcompnames <- cbind(fine_shrink_swell[fine_shrink_swell$compname %in% fine_saline_sodic$compname,], fine_saline_sodic[fine_saline_sodic$compname %in% fine_shrink_swell$compname, ])
colnames(fine_conflicts_majcompnames) <- c('fine_shrink_swell_compname', 'area_ha', 'fine_saline_sodic_compname', 'area_ha')
write.csv(fine_conflicts_majcompnames, file.path(dataDir, 'soil survey facts', 'fine_conficts_majcompnames.csv'), row.names = FALSE)

sum(coarse_saline_sodic$compname %in% fine_saline_sodic$compname) #29
sum(coarse_w_no_res$compname %in% loamy_w_no_res$compname) #63
sum(coarse_w_pans$compname %in% loamy_w_pans$compname) #86

#find how many components make up a certain percentage of a SHR
compname_n_percentage <- function(compnames, fname, perc) {
  compnames <- compnames[order(compnames$area_ha, decreasing = TRUE),]
  print(nrow(compnames))
  compnames$cumperc <- cumsum(compnames$area_ha) / sum(compnames$area_ha)
  compnames$tot_comp_area <- area_by_compname$area_ha[match(compnames$compname, area_by_compname$compname)]
  compnames$SHR_purity_perc <- round(100* compnames$area_ha / compnames$tot_comp_area, 3)
  compnames$taxorder <- area_by_compname$order[match(compnames$compname, area_by_compname$compname)]
  # abs(compnames$cumperc-perc)
  if (perc==1) {
    compnames_trimmed <- compnames
  } else {compnames_trimmed <- compnames[1:which.min(abs(compnames$cumperc-perc)), ]}
  print(nrow(compnames_trimmed))
  write.csv(compnames_trimmed, file.path(dataDir, 'soil survey facts', fname, paste0(fname, '_top', perc*100, '_areal_percentage.csv')), row.names = FALSE)
}
compname_n_percentage(coarse_w_no_res, 'Coarse_w_no_res', 1)
compname_n_percentage(loamy_w_no_res, 'Loamy_w_no_res', 1)
compname_n_percentage(coarse_w_pans, 'Coarse_w_res', 1)
compname_n_percentage(loamy_w_pans, 'Loamy_w_res', 1)
compname_n_percentage(coarse_saline_sodic, 'Coarse_saline_sodic', 1)
compname_n_percentage(fine_saline_sodic, 'Fine_saline_sodic', 1)
compname_n_percentage(fine_shrink_swell, 'Fine_shrink_swell', 1)

compname_n_percentage(coarse_w_no_res, 'Coarse_w_no_res', 0.75)
compname_n_percentage(loamy_w_no_res, 'Loamy_w_no_res', 0.75)
compname_n_percentage(coarse_w_pans, 'Coarse_w_res', 0.75)
compname_n_percentage(loamy_w_pans, 'Loamy_w_res', 0.75)
compname_n_percentage(coarse_saline_sodic, 'Coarse_saline_sodic', 0.75)
compname_n_percentage(fine_saline_sodic, 'Fine_saline_sodic', 0.75)
compname_n_percentage(fine_shrink_swell, 'Fine_shrink_swell', 0.75)

compname_n_percentage(coarse_w_no_res, 'Coarse_w_no_res', 0.5)
compname_n_percentage(loamy_w_no_res, 'Loamy_w_no_res', 0.5)
compname_n_percentage(coarse_w_pans, 'Coarse_w_res', 0.5)
compname_n_percentage(loamy_w_pans, 'Loamy_w_res', 0.5)
compname_n_percentage(coarse_saline_sodic, 'Coarse_saline_sodic', 0.5)
compname_n_percentage(fine_saline_sodic, 'Fine_saline_sodic', 0.5)
compname_n_percentage(fine_shrink_swell, 'Fine_shrink_swell', 0.5)
