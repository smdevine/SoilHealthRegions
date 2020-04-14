library(raster)
library(aqp)
library(vioplot)
library(extrafont)
library(extrafontdb)
loadfonts(device = 'win')
laptop <- FALSE
if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
  workDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  workDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}
list.files(workDir)
columnClass <- c(Concatenate='character')
soil_data <- read.csv(file.path(workDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, colClasses = columnClass, na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', ""))
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))
soil_data_pts <- data.frame(ID=unique(soil_data$Concatenate),stringsAsFactors = FALSE)
soil_data_pts$Lat_WGS84 <- soil_data$GPS.N[match(soil_data_pts$ID, soil_data$Concatenate)]
soil_data_pts$Lon_WGS84 <- soil_data$GPS.W[match(soil_data_pts$ID, soil_data$Concatenate)]
soil_data_pts <- SpatialPointsDataFrame(coords = soil_data_pts[,c('Lon_WGS84', 'Lat_WGS84')], data = soil_data_pts['ID'], proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
# plot(soil_data_pts)
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster.shp'))
crs(valley_mu_shp_30cm)
pts_ssurgo_extract_30cm <- extract(valley_mu_shp_30cm, soil_data_pts, df=TRUE)
table(pts_ssurgo_extract_30cm$cluster_7)
sum(is.na(pts_ssurgo_extract_30cm$cluster_7))
pts_ssurgo_extract_30cm$Concatenate <- soil_data_pts$ID[pts_ssurgo_extract_30cm$point.ID] #point.ID is the row number
pts_ssurgo_extract_30cm[which(pts_ssurgo_extract_30cm$cluster_9==2),]
pts_ssurgo_extract_30cm[which(pts_ssurgo_extract_30cm$cluster_9==4),]

#aggregate horizon data to 0-30 cm
colnames(soil_data)
horizon_data <- soil_data[,c(4:5,15:18,20:24,27)]
colnames(horizon_data)
colnames(horizon_data)[2:12] <- c('Depth', 'gravel', 'clay', 'silt', 'sand', 'totN', 'totC', 'ph_h2O', 'DOC', 'DON', 'db')
lapply(horizon_data, class)
lapply(horizon_data[,3:12], summary)
horizon_data$totC[which(horizon_data$totC<0.01)] <- NA
horizon_data$totN[which(horizon_data$totN<0.005)] <- NA
unique(horizon_data$Depth)
horizon_data$hzn_top <- ifelse(horizon_data$Depth=='0_5', 0, ifelse(horizon_data$Depth=='5_10', 5, ifelse(horizon_data$Depth=='10_30', 10, ifelse(horizon_data$Depth=='30_50', 30, ifelse(horizon_data$Depth=='50_100', 50, NA)))))
horizon_data$hzn_bot <- ifelse(horizon_data$Depth=='0_5', 5, ifelse(horizon_data$Depth=='5_10', 10, ifelse(horizon_data$Depth=='10_30', 30, ifelse(horizon_data$Depth=='30_50', 50, ifelse(horizon_data$Depth=='50_100', 100, NA)))))

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}
kgOrgC_sum_v2 <- function(x, depth, rm.NAs=FALSE) {
  thick <- x$hzn_bot - x$hzn_top
  sum((thick / 10) * x$totC * x$db * (1 - x$gravel / 100), na.rm = rm.NAs)
}
# awc_sum_v2 <- function(x, rm.NAs=FALSE) {
#   thick <- x$hzn_bot - x$hzn_top
#   sum(thick * x$whc, na.rm = rm.NAs)
# }
# horizon_SPC <- kssl_horizons_subset[112:113]
# depth <- 30
# vars_of_interest <- c('clay')
# varnames <- 'clay'
horizon_to_comp_v3 <- function(horizon_SPC, depth, vars_of_interest = c('clay', 'silt', 'sand', 'totC', 'totN', 'db', 'gravel', 'ph_h2O',  'DOC', 'DON'), varnames = c('clay', 'silt', 'sand', 'totC', 'totN', 'bd', 'frags', 'pH_H2O', 'DOC', 'DON')) {
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
  horizons(sliced_SPC) <- horizons(sliced_SPC)[order(as.integer(sliced_SPC$Concatenate), sliced_SPC$hzn_top),] #because integer codes are coerced to character by slice with the sliced data.frame re-ordered contrary to order of site-level data 
  stopifnot(unique(sliced_SPC$Concatenate)==site(sliced_SPC)$Concatenate)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean_v2, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum_v2)
  #s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum_v2, rm.NAs=FALSE) #was TRUE for SSURGO
  #columnames <- c(columnames, paste0('awc_', depth, 'cm')) 
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  s
}

depths(horizon_data) <- Concatenate ~ hzn_top + hzn_bot
horizon_data
class(horizon_data)
pts_30cm <- horizon_to_comp_v3(horizon_SPC = horizon_data, depth = 30)
colnames(pts_30cm)
lapply(pts_30cm[,2:12], summary)
dim(pts_30cm)

pts_10cm <- horizon_to_comp_v3(horizon_SPC = horizon_data, depth = 10)
colnames(pts_10cm)
lapply(pts_10cm[,2:12], summary)
dim(pts_10cm)

pts_50cm <- horizon_to_comp_v3(horizon_SPC = horizon_data, depth = 50)
colnames(pts_50cm)
lapply(pts_50cm[,2:12], summary)
dim(pts_50cm)

pts_100cm <- horizon_to_comp_v3(horizon_SPC = horizon_data, depth = 100)
colnames(pts_100cm)
lapply(pts_100cm[,2:12], summary)
#check the aggregation
pts_30cm[126,]
#clay 30 cm: 51, 49.9, 51.15
51*5/30+49.9*5/30+51.15*20/30
#totC 30cm: 2.737411672, 2.671205503, 2.004525189
2.737411672*5/30+2.671205503*5/30+2.004525189*20/30
pts_30cm[15,]
#clay 30 cm: 9.183016977, 9.562007334, 9.034136208
9.183016977*5/30+9.562007334*5/30+9.034136208*20/30
#totC 30cm: 1.765431644, 1.497917658, 0.889017189
1.765431644*5/30+1.497917658*5/30+0.889017189*20/30
pts_30cm[61, ]
#clay 30 cm: 21.1, 21.1, 21.95
21.1*5/30+21.1*5/30+21.95*20/30
#totC 30cm: 3.552454827, 3.04512963, 2.922880435
3.552454827*5/30+3.04512963*5/30+2.922880435*20/30

#add cluster info
pts_30cm$cluster_2 <- pts_ssurgo_extract_30cm$cluster_2[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_3 <- pts_ssurgo_extract_30cm$cluster_3[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_4 <- pts_ssurgo_extract_30cm$cluster_4[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_5 <- pts_ssurgo_extract_30cm$cluster_5[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_6 <- pts_ssurgo_extract_30cm$cluster_6[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_7 <- pts_ssurgo_extract_30cm$cluster_7[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_8 <- pts_ssurgo_extract_30cm$cluster_8[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_9 <- pts_ssurgo_extract_30cm$cluster_9[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_10 <- pts_ssurgo_extract_30cm$cluster_10[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_11 <- pts_ssurgo_extract_30cm$cluster_11[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_30cm$cluster_12 <- pts_ssurgo_extract_30cm$cluster_12[match(pts_30cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]


colnames(pts_30cm)
tapply(pts_30cm$totC_30cm, pts_30cm$cluster_9, summary)
tapply(pts_30cm$pH_H2O_30cm, pts_30cm$cluster_9, summary)
tapply(pts_30cm$clay_30cm, pts_30cm$cluster_9, summary)
table(pts_30cm$cluster_9[!is.na(pts_30cm$clay_30cm)])
tapply(pts_30cm$bd_30cm, pts_30cm$cluster_9, summary)
tapply(pts_30cm$bd_30cm, pts_30cm$cluster_9, summary)
write.csv(pts_ssurgo_extract_30cm, file.path(workDir, 'FINAL', 'CDFA_pts_ssurgo_30cm_extract_FINAL.csv'), row.names = FALSE)

# pts_ssurgo_extract_30cm <- read.csv(file.path(workDir, 'CDFA_pts_ssurgo_30cm_extract.csv'), stringsAsFactors = FALSE) 
# pts_30cm <- read.csv(file.path(workDir, 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)

#add cluster info to 10 cm soil depth aggregation
pts_10cm$cluster_2 <- pts_ssurgo_extract_30cm$cluster_2[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_3 <- pts_ssurgo_extract_30cm$cluster_3[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_4 <- pts_ssurgo_extract_30cm$cluster_4[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_5 <- pts_ssurgo_extract_30cm$cluster_5[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_6 <- pts_ssurgo_extract_30cm$cluster_6[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_7 <- pts_ssurgo_extract_30cm$cluster_7[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_8 <- pts_ssurgo_extract_30cm$cluster_8[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_9 <- pts_ssurgo_extract_30cm$cluster_9[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_10 <- pts_ssurgo_extract_30cm$cluster_10[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_11 <- pts_ssurgo_extract_30cm$cluster_11[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_10cm$cluster_12 <- pts_ssurgo_extract_30cm$cluster_12[match(pts_10cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]

#add cluster info to 50 cm soil depth aggregation
pts_50cm$cluster_2 <- pts_ssurgo_extract_30cm$cluster_2[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_3 <- pts_ssurgo_extract_30cm$cluster_3[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_4 <- pts_ssurgo_extract_30cm$cluster_4[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_5 <- pts_ssurgo_extract_30cm$cluster_5[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_6 <- pts_ssurgo_extract_30cm$cluster_6[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_7 <- pts_ssurgo_extract_30cm$cluster_7[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_8 <- pts_ssurgo_extract_30cm$cluster_8[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_9 <- pts_ssurgo_extract_30cm$cluster_9[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_10 <- pts_ssurgo_extract_30cm$cluster_10[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_11 <- pts_ssurgo_extract_30cm$cluster_11[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]
pts_50cm$cluster_12 <- pts_ssurgo_extract_30cm$cluster_12[match(pts_50cm$Concatenate, pts_ssurgo_extract_30cm$Concatenate)]

#add management info
pts_30cm$tillage <- soil_data$till.or.no.till[match(pts_30cm$Concatenate, soil_data$Concatenate)]
unique(soil_data$compost.added.)
soil_data$compost.added.[which(soil_data$compost.added.=="")] <- NA
pts_30cm$compost_added <- soil_data$compost.added.[match(pts_30cm$Concatenate, soil_data$Concatenate)]
unique(soil_data$irrigated.vs.dryfarm)
pts_30cm$irrigated_vs_dryfarm <- soil_data$irrigated.vs.dryfarm[match(pts_30cm$Concatenate, soil_data$Concatenate)]
unique(soil_data$organic.vs.conventional.vs.biodynamic)
soil_data$organic.vs.conventional.vs.biodynamic[which(soil_data$organic.vs.conventional.vs.biodynamic=='organic')] <- 'Organic'
pts_30cm$management_type <- soil_data$organic.vs.conventional.vs.biodynamic[match(pts_30cm$Concatenate, soil_data$Concatenate)]
write.csv(pts_30cm, file.path(workDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), row.names = FALSE)

#10 cm management info
pts_10cm$tillage <- soil_data$till.or.no.till[match(pts_10cm$Concatenate, soil_data$Concatenate)]
pts_10cm$compost_added <- soil_data$compost.added.[match(pts_10cm$Concatenate, soil_data$Concatenate)]
pts_10cm$irrigated_vs_dryfarm <- soil_data$irrigated.vs.dryfarm[match(pts_10cm$Concatenate, soil_data$Concatenate)]
pts_10cm$management_type <- soil_data$organic.vs.conventional.vs.biodynamic[match(pts_10cm$Concatenate, soil_data$Concatenate)]
write.csv(pts_10cm, file.path(workDir, 'FINAL', 'CDFA_samples_cluster_10cm.csv'), row.names = FALSE)

#50 cm management info
pts_50cm$tillage <- soil_data$till.or.no.till[match(pts_50cm$Concatenate, soil_data$Concatenate)]
pts_50cm$compost_added <- soil_data$compost.added.[match(pts_50cm$Concatenate, soil_data$Concatenate)]
pts_50cm$irrigated_vs_dryfarm <- soil_data$irrigated.vs.dryfarm[match(pts_50cm$Concatenate, soil_data$Concatenate)]
pts_50cm$management_type <- soil_data$organic.vs.conventional.vs.biodynamic[match(pts_50cm$Concatenate, soil_data$Concatenate)]
write.csv(pts_50cm, file.path(workDir, 'FINAL', 'CDFA_samples_cluster_50cm.csv'), row.names = FALSE)

#make vioplots
order_lgnd_9 <- c(3,6,5,9,4,2) #reflecting the clusters that are missing
# vector_include <- table(pts_30cm$cluster_9[!is.na(pts_30cm$clay_30cm)]) > 4
# plot_order2 <- (1:9)[order_lgnd_9]
# ((1:9)[order_lgnd_9])[as.logical(vector_include)]
# order_lgnd_9[as.logical(vector_include)[match(as.character(order_lgnd_9), names(vector_include))]]
# plot_order2[as.logical(vector_include)[match(as.character(plot_order2), names(vector_include))]]
vioplot_soil_clus9 <- function(df, varname, ylim_vioplot, plot_order, input_include, ylab, fname, mar, group_names) {
  plot_order2 <- (1:9)[plot_order]
  if(is.null(names(input_include))) {
    names(input_include) <- plot_order
  }
  color_include <- plot_order[as.logical(input_include)[match(as.character(plot_order), names(input_include))]]
  print(color_include)
  plot_list <- list(df[[varname]][df$cluster_9==plot_order2[1]], df[[varname]][df$cluster_9==plot_order2[2]], df[[varname]][df$cluster_9==plot_order2[3]], df[[varname]][df$cluster_9==plot_order2[4]], df[[varname]][df$cluster_9==plot_order2[5]], df[[varname]][df$cluster_9==plot_order2[6]])[as.logical(input_include)[match(as.character(plot_order2), names(input_include))]]
  tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(plot_list, col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[color_include], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab, names = group_names)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
colnames(pts_30cm)
#order_lgnd was defined for radarchart

vioplot_soil_clus9(pts_30cm, 'clay_30cm', ylim_vioplot = c(0.5,55), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$clay_30cm)]) > 4, ylab='0-30 cm clay (%)', fname='class9_clay_0_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
summary(pts_30cm$silt_30cm)
table(pts_30cm$cluster_9[!is.na(pts_30cm$silt_30cm)])
vioplot_soil_clus9(pts_30cm, 'silt_30cm', ylim_vioplot = c(0,70), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$silt_30cm)]) > 4, ylab='0-30 cm silt (%)', fname='class9_silt_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
summary(pts_30cm$sand_30cm)
vioplot_soil_clus9(pts_30cm, 'sand_30cm', ylim_vioplot = c(0,75), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$sand_30cm)]) > 4, ylab='0-30 cm sand (%)', fname='class9_sand_0_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
summary(pts_30cm$frags_30cm)
vioplot_soil_clus9(pts_30cm, 'frags_30cm', ylim_vioplot = c(0,75), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$frags_30cm)]) > 4, ylab='0-30 cm gravel (%)', fname='class9_gravel_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
vioplot_soil_clus9(pts_30cm, 'totC_30cm', ylim_vioplot = c(0.1,4), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$totC_30cm)]) > 4,  ylab='0-30 cm soil organic carbon (%)', fname='class9_OC_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
vioplot_soil_clus9(pts_30cm, 'pH_H2O_30cm', ylim_vioplot = c(5.5,8), plot_order = order_lgnd_9, input_include = c(rep(TRUE, 4), rep(FALSE, 2)), ylab='0-30 cm soil pH', fname='class9_pH_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
summary(pts_30cm$bd_30cm)
table(pts_30cm$cluster_9[!is.na(pts_30cm$bd_30cm)])
tapply(pts_30cm$bd_30cm, pts_30cm$cluster_9, summary)
vioplot_soil_clus9(pts_30cm, 'bd_30cm', ylim_vioplot = c(0.9,1.75), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$bd_30cm)]) > 4,  ylab=expression('0-30 cm bulk density (g soil cm'^-3*'soil)'), fname='class9_BD_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names=c(2,3,4,5))
vioplot_soil_clus9(pts_30cm, 'totN_30cm', ylim_vioplot = c(0,0.24), plot_order = order_lgnd_9, input_include = table(pts_30cm$cluster_9[!is.na(pts_30cm$totN_30cm)]) > 4,  ylab='0-30 cm soil nitrogen, total (%)', fname='class9_totN_30cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))

table(pts_30cm$cluster_9[!is.na(pts_30cm$DOC_30cm)])
table(pts_30cm$cluster_9[!is.na(pts_30cm$totN_30cm)])

#vioplots for 10 cm
summary(pts_10cm$DOC_10cm)
vioplot_soil_clus9(pts_10cm, 'DOC_10cm', ylim_vioplot = c(7,170), plot_order = order_lgnd_9, input_include = table(pts_10cm$cluster_9[!is.na(pts_10cm$DOC_10cm)]) > 4,  ylab='Dissolved organic carbon (ppm) 0-10 cm soil', fname='class9_DOC_10cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
vioplot_soil_clus9(pts_10cm, 'DOC_10cm', ylim_vioplot = c(7,170), plot_order = order_lgnd_9, input_include = table(pts_10cm$cluster_9[!is.na(pts_10cm$DOC_10cm)]) > 4,  ylab='0-10 cm dissolved organic carbon (ppm)', fname='class9_DOC_10cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
vioplot_soil_clus9(pts_10cm, 'totC_10cm', ylim_vioplot = c(0.1,4), plot_order = order_lgnd_9, input_include = table(pts_10cm$cluster_9[!is.na(pts_10cm$totC_10cm)]) > 4,  ylab='0-10 cm soil organic carbon (%)', fname='class9_OC_10cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))
vioplot_soil_clus9(pts_10cm, 'clay_10cm', ylim_vioplot = c(0.5,55), plot_order = order_lgnd_9, input_include = table(pts_10cm$cluster_9[!is.na(pts_10cm$clay_10cm)]) > 4, ylab='0-10 cm clay (%)', fname='class9_clay_0_10cm_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(2,3,4,5))

#tabulate management reps overall and by cluster
#overall with %C and cluster info
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm))])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm))]))
sum(!is.na(pts_30cm$totC_30cm))
sum(!is.na(pts_30cm$totC_30cm)&is.na(pts_30cm$cluster_9))
#tilled samples with %C data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till')]))
#no till samples with %C data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='no till')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='no till')]))

#conv. samples with %C data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional')]))
#org. samples with %C data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic')]))
#biodynamic samples with %C data count by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic')]))
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & is.na(pts_30cm$management_type))]))

#irrigated samples
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated')]))

#dryland samples
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm')]))

#compost added samples
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$compost_added=='yes')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$compost_added=='yes')]))

#no compost added samples
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$compost_added=='none')])
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$compost_added=='none')]))
sum(table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & is.na(pts_30cm$compost_added))]))

#summary of groups by irrigation x tillage x compost
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$tillage=='till' & pts_30cm$compost_added=='yes')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$tillage=='no till' & pts_30cm$compost_added=='yes')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$tillage=='till' & pts_30cm$compost_added=='none')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$tillage=='no till' & pts_30cm$compost_added=='none')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$tillage=='till' & pts_30cm$compost_added=='yes')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$tillage=='no till' & pts_30cm$compost_added=='yes')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$tillage=='till' & pts_30cm$compost_added=='none')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$tillage=='no till' & pts_30cm$compost_added=='none')])

#total by mgmt type with C and cluster data
sum(!is.na(pts_30cm$totC_30cm) & !is.na(pts_30cm$management_type) & !is.na(pts_30cm$cluster_9))
#org. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic')])
#conv. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional')])
#biodynamic samples with C
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic')])
#org. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic')])
#conv. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional')])
#biodynamic samples with C
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic')])

#org. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic' & pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic' & pts_30cm$tillage=='no till')])

#conv. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional'& pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional'& pts_30cm$tillage=='no till')])

#biodynamic samples with C
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic'& pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic'& pts_30cm$tillage=='no till')])

#org. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic' & pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Organic' & pts_30cm$tillage=='no till')])

#conv. samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional'& pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='conventional'& pts_30cm$tillage=='no till')])

#biodynamic samples with C
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic'& pts_30cm$tillage=='till')])
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$totC_30cm) & pts_30cm$management_type=='Biodynamic'& pts_30cm$tillage=='no till')])

#tilled samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$kgOrg.m2_30cm) & pts_30cm$tillage=='till')])
#no till samples with C content data counts by cluster 9 ID
table(pts_30cm$cluster_9[which(!is.na(pts_30cm$kgOrg.m2_30cm) & pts_30cm$tillage=='no till')])

#30 cm
#(a) Compare Till x Irrigated x Compost across these regions: cluster 3 (n=9) vs. cluster 5 (n=9) vs. cluster 6 (n=14)
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==3, na.rm=TRUE) #9
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==5, na.rm=TRUE)#9
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6, na.rm=TRUE) #14
vioplot(pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==3)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==5)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], names=c(3,5,6))

#(b)	Compare tillage in cluster 6 (controlling for irrigation and compost): till x irrigated x compost (n=14); no-till x irrigated x compost (n=13); dryfarm x tillage x compost (n=12)
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6, na.rm=TRUE) #14
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6, na.rm=TRUE)#13
sum(!is.na(pts_30cm$totC_30cm) & pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6, na.rm=TRUE) #12
vioplot(pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_9==6)], names=c('T-Irr', 'NT-Irr', 'T-DF'))
#across all clusters
vioplot(pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes')], pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes')], pts_30cm$totC_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes')], names=c('T-Irr', 'NT-Irr', 'T-DF'))


#(c)	Compare cluster 6 (n=13) vs. cluster 9 (n=8): controlling for tillage (no-till), irrigation and compost): 
vioplot(pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==9)])
t.test(pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], pts_30cm$totC_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==9)])

#compare 3 management types in cluster 6 across each depth segment
soil_data_cluster6_of_9 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_9==6)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_9==6)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_9==6)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_9==6)], stringsAsFactors = FALSE)
head(soil_data_cluster6_of_9)
dim(soil_data_cluster6_of_9)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_cluster6_of_9$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster6_of_9$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster6_of_9$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster6_of_9$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster6_of_9$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster6_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
colnames(soil_data_cluster6_of_9)
lapply(soil_data_cluster6_of_9[,5:9], summary)
lapply(soil_data_cluster6_of_9[,5:9], function(x) tapply(x, soil_data_cluster6_of_9$tillage, mean, na.rm=TRUE))
lapply(soil_data_cluster6_of_9[,5:9], function(x) tapply(x, soil_data_cluster6_of_9$irrigation, mean, na.rm=TRUE))
lapply(soil_data_cluster6_of_9[,5:9], function(x) tapply(x, soil_data_cluster6_of_9$compost, mean, na.rm=TRUE))


till_irr_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='till' & soil_data_cluster6_of_9$irrigation=='irrigated' & soil_data_cluster6_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='till' & soil_data_cluster6_of_9$irrigation=='irrigated' & soil_data_cluster6_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))
notill_irr_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='no till' & soil_data_cluster6_of_9$irrigation=='irrigated' & soil_data_cluster6_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='no till' & soil_data_cluster6_of_9$irrigation=='irrigated' & soil_data_cluster6_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))
till_dry_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='till' & soil_data_cluster6_of_9$irrigation=='dryfarm' & soil_data_cluster6_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_9[which(soil_data_cluster6_of_9$tillage=='till' & soil_data_cluster6_of_9$irrigation=='dryfarm' & soil_data_cluster6_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))


soil_data_cluster3_of_9 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_9==3)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_9==3)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_9==3)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_9==3)], stringsAsFactors = FALSE)
soil_data_cluster3_of_9$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster3_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster3_of_9$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster3_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster3_of_9$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster3_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster3_of_9$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster3_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster3_of_9$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster3_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
till_irr_clus3 <- data.frame(C_means=sapply(soil_data_cluster3_of_9[which(soil_data_cluster3_of_9$tillage=='till' & soil_data_cluster3_of_9$irrigation=='irrigated' & soil_data_cluster3_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster3_of_9[which(soil_data_cluster3_of_9$tillage=='till' & soil_data_cluster3_of_9$irrigation=='irrigated' & soil_data_cluster3_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))

soil_data_cluster9_of_9 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_9==9)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_9==9)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_9==9)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_9==9)], stringsAsFactors = FALSE)
soil_data_cluster9_of_9$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster9_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster9_of_9$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster9_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster9_of_9$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster9_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster9_of_9$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster9_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster9_of_9$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster9_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
notill_irr_clus9 <- data.frame(C_means=sapply(soil_data_cluster9_of_9[which(soil_data_cluster9_of_9$tillage=='no till' & soil_data_cluster9_of_9$irrigation=='irrigated' & soil_data_cluster9_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster9_of_9[which(soil_data_cluster9_of_9$tillage=='no till' & soil_data_cluster9_of_9$irrigation=='irrigated' & soil_data_cluster9_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))

soil_data_cluster5_of_9 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_9==5)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_9==5)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_9==5)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_9==5)], stringsAsFactors = FALSE)
soil_data_cluster5_of_9$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster5_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster5_of_9$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster5_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster5_of_9$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster5_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster5_of_9$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster5_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster5_of_9$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster5_of_9$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
till_irr_clus5 <- data.frame(C_means=sapply(soil_data_cluster5_of_9[which(soil_data_cluster5_of_9$tillage=='till' & soil_data_cluster5_of_9$irrigation=='irrigated' & soil_data_cluster5_of_9$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster5_of_9[which(soil_data_cluster5_of_9$tillage=='till' & soil_data_cluster5_of_9$irrigation=='irrigated' & soil_data_cluster5_of_9$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))


depths_cm <- c(2.5, 7.5, 20, 40, 75)
tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', 'CDFA_soilC_profile.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(till_irr_clus6$C_means, -depths_cm, type='b', xlim=c(0.1,2.6), pch=24, bg='tan4', col='black', cex=1, xlab='', ylab='', yaxt='n')
axis(side = 2, at=seq(from=-100, to=0, by=20), labels=-seq(from=-100, to=0, by=20))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(till_dry_clus6$C_means, -depths_cm, type='b', pch=24, bg='tan4', col='black', lty=2)
points(notill_irr_clus6$C_means, -depths_cm, type='b', pch=23, bg='tan4', col='black')
points(till_irr_clus3$C_means, -depths_cm, type = 'b', pch=24, bg='tan2', col='black')
points(till_irr_clus5$C_means, -depths_cm, type = 'b', pch=24, bg='gold', col='black')
points(notill_irr_clus9$C_means, -depths_cm, type='b', pch=23, bg='firebrick3', col='black')
legend('bottomright', legend=c('2. Loams w/ no res. & mod OM-low SS', '3a. Loams w/ no res. & mod OM-mod SS', '3b. Loams w/ no res. & mod OM-mod SS: No-till', '3c. Loams w/ no res. & mod OM-mod SS: Dry farm', '4. Loams w/ res. & low OM', '5. Loams w/ res. & mod OM: No-till'), pch=c(24,24,23,24,24,23), lty=c(1,1,1,2,1,1), pt.bg=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), col='black')
dev.off()

#make barplot of 30cm SOC kg m^-2
plot_da_error <- function(barnum, barcenters, df) {
  segments(x0=barcenters[barnum,], y0=df$means[barnum] + qnorm(0.975) * df$sd[barnum] / sqrt(df$n[barnum]), x1=barcenters[barnum,], y1=df$means[barnum] + qnorm(0.025) * df$sd[barnum] / sqrt(df$n[barnum]), lwd = 1.2)
}
kgOrgC_30cm_data <- data.frame(means=c(mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==3)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_9==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==5)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==9)], na.rm=TRUE)), sd=c(sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==3)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_9==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==5)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==9)], na.rm=TRUE)), n=c(sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==3)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_9==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==5)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_9==9)]))), row.names = c('2. Loams w/ no res. & mod OM-low SS', '3a. Loams w/ no res. & mod OM-mod SS', '3b. Loams w/ no res. & mod OM-mod SS: No-till', '3c. Loams w/ no res. & mod OM-mod SS: Dry farm', '4. Loams w/ res. & low OM', '5. Loams w/ res. & mod OM: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', 'CDFA_soilC_kg_0_30cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_30cm <- barplot(height=kgOrgC_30cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('2', '3a', '3b', '3c', '4', '5'), ylim=c(0,8.5))
mtext(text=expression('0-30 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(2, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(3, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(4, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(5, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(6, barcenters_30cm, kgOrgC_30cm_data)
dev.off()

#plot 50 cm C contents
kgOrgC_50cm_data <- data.frame(means=c(mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==3)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_9==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==5)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==9)], na.rm=TRUE)), sd=c(sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==3)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_9==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==5)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==9)], na.rm=TRUE)), n=c(sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==3)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_9==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==5)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_9==9)]))), row.names = c('2. Loams w/ no res. & mod OM-low SS', '3a. Loams w/ no res. & mod OM-mod SS', '3b. Loams w/ no res. & mod OM-mod SS: No-till', '3c. Loams w/ no res. & mod OM-mod SS: Dry farm', '4. Loams w/ res. & low OM', '5. Loams w/ res. & mod OM: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', 'CDFA_soilC_kg_0_50cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_50cm <- barplot(height=kgOrgC_50cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('2', '3a', '3b', '3c', '4', '5'), ylim=c(0,13))
mtext(text=expression('0-50 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(2, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(3, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(4, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(5, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(6, barcenters_50cm, kgOrgC_50cm_data)
dev.off()

#plot 50 cm C contents
kgOrgC_10cm_data <- data.frame(means=c(mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==3)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_9==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==5)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==9)], na.rm=TRUE)), sd=c(sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==3)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_9==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==5)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==9)], na.rm=TRUE)), n=c(sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==3)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_9==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==5)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_9==9)]))), row.names = c('2. Loams w/ no res. & mod OM-low SS', '3a. Loams w/ no res. & mod OM-mod SS', '3b. Loams w/ no res. & mod OM-mod SS: No-till', '3c. Loams w/ no res. & mod OM-mod SS: Dry farm', '4. Loams w/ res. & low OM', '5. Loams w/ res. & mod OM: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', 'CDFA_soilC_kg_0_10cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_10cm <- barplot(height=kgOrgC_10cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('2', '3a', '3b', '3c', '4', '5'), ylim=c(0,3.9))
mtext(text=expression('0-10 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(2, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(3, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(4, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(5, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(6, barcenters_10cm, kgOrgC_10cm_data)
dev.off()

#%C by tillage and all samples
vioplot(pts_30cm$totC_30cm[pts_30cm$tillage=='till'], pts_30cm$totC_30cm[pts_30cm$tillage=='no till'], names = c('Tillage', 'No-till'))

#%C by tillage for cluster 6 of 9
vioplot(pts_30cm$totC_30cm[pts_30cm$tillage=='till'& pts_30cm$cluster_9==6], pts_30cm$totC_30cm[pts_30cm$tillage=='no till'& pts_30cm$cluster_9==6], names = c('Tillage', 'No-till'))
#C content by tillage and all samples
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$tillage=='till'], pts_30cm$kgOrg.m2_30cm[pts_30cm$tillage=='no till'], names = c('Tillage', 'No-till'))

#C content by tillage for cluster 6 of 9
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$tillage=='till'& pts_30cm$cluster_9==6], pts_30cm$kgOrg.m2_30cm[pts_30cm$tillage=='no till'& pts_30cm$cluster_9==6], names = c('Tillage', 'No-till'))

#%C by org. vs. conv. and all samples
vioplot(pts_30cm$totC_30cm[pts_30cm$management_type=='conventional'], pts_30cm$totC_30cm[pts_30cm$management_type=='Organic'], names = c('Conv.', 'Organic'))

#%C by conv. & org. for cluster 6 of 9
vioplot(pts_30cm$totC_30cm[pts_30cm$management_type=='conventional' & pts_30cm$cluster_9==6], pts_30cm$totC_30cm[pts_30cm$management_type=='Organic' & pts_30cm$cluster_9==6], names = c('Conv.', 'Organic'))
#C content by conv. & org. for all samples
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional'], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic'], names = c('Conv.', 'Organic'))

#C content by org. vs. conv. for cluster 6 of 9
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional' & pts_30cm$cluster_9==6], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic' & pts_30cm$cluster_9==6], names = c('Conv.', 'Organic'))
#C content by tillage x mgmt type for all samples
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='no till'], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='till'], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='no till'], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='till'], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#C content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='no till' & pts_30cm$cluster_9==6], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='till' & pts_30cm$cluster_9==6], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='no till' & pts_30cm$cluster_9==6], pts_30cm$kgOrg.m2_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='till' & pts_30cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#gravel content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_30cm$frags_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='no till' & pts_30cm$cluster_9==6], pts_30cm$frags_30cm[pts_30cm$management_type=='conventional' & pts_30cm$tillage=='till' & pts_30cm$cluster_9==6], pts_30cm$frags_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='no till' & pts_30cm$cluster_9==6], pts_30cm$frags_30cm[pts_30cm$management_type=='Organic' & pts_30cm$tillage=='till' & pts_30cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))

#10 cm
#%C by tillage and all samples
vioplot(pts_10cm$totC_10cm[pts_10cm$tillage=='till'], pts_10cm$totC_10cm[pts_10cm$tillage=='no till'], names = c('Tillage', 'No-till'))
#tilled samples with %C data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$totC_10cm) & pts_10cm$tillage=='till')])
#no till samples with %C data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$totC_10cm) & pts_10cm$tillage=='no till')])
#%C by tillage for cluster 6 of 9
vioplot(pts_10cm$totC_10cm[pts_10cm$tillage=='till'& pts_10cm$cluster_9==6], pts_10cm$totC_10cm[pts_10cm$tillage=='no till'& pts_10cm$cluster_9==6], names = c('Tillage', 'No-till'))
#C content by tillage and all samples
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$tillage=='till'], pts_10cm$kgOrg.m2_10cm[pts_10cm$tillage=='no till'], names = c('Tillage', 'No-till'))
#tilled samples with C content data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$kgOrg.m2_10cm) & pts_10cm$tillage=='till')])
#no till samples with C content data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$kgOrg.m2_10cm) & pts_10cm$tillage=='no till')])
#C content by tillage for cluster 6 of 9
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$tillage=='till'& pts_10cm$cluster_9==6], pts_10cm$kgOrg.m2_10cm[pts_10cm$tillage=='no till'& pts_10cm$cluster_9==6], names = c('Tillage', 'No-till'))

#%C by org. vs. conv. and all samples
vioplot(pts_10cm$totC_10cm[pts_10cm$management_type=='conventional'], pts_10cm$totC_10cm[pts_10cm$management_type=='Organic'], names = c('Conv.', 'Organic'))
#conv. samples with %C data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$totC_10cm) & pts_10cm$management_type=='conventional')])
#org. samples with %C data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$totC_10cm) & pts_10cm$management_type=='Organic')])
#%C by conv. & org. for cluster 6 of 9
vioplot(pts_10cm$totC_10cm[pts_10cm$management_type=='conventional' & pts_10cm$cluster_9==6], pts_10cm$totC_10cm[pts_10cm$management_type=='Organic' & pts_10cm$cluster_9==6], names = c('Conv.', 'Organic'))
#C content by conv. & org. for all samples
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional'], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic'], names = c('Conv.', 'Organic'))
#conv. samples with C content data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$kgOrg.m2_10cm) & pts_10cm$management_type=='conventional')])
#org. samples with C content data counts by cluster 9 ID
table(pts_10cm$cluster_9[which(!is.na(pts_10cm$kgOrg.m2_10cm) & pts_10cm$management_type=='Organic')])
#C content by org. vs. conv. for cluster 6 of 9
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional' & pts_10cm$cluster_9==6], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic' & pts_10cm$cluster_9==6], names = c('Conv.', 'Organic'))
#C content by tillage x mgmt type for all samples
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='no till'], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='till'], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='no till'], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='till'], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#C content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='no till' & pts_10cm$cluster_9==6], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='till' & pts_10cm$cluster_9==6], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='no till' & pts_10cm$cluster_9==6], pts_10cm$kgOrg.m2_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='till' & pts_10cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#gravel content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_10cm$frags_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='no till' & pts_10cm$cluster_9==6], pts_10cm$frags_10cm[pts_10cm$management_type=='conventional' & pts_10cm$tillage=='till' & pts_10cm$cluster_9==6], pts_10cm$frags_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='no till' & pts_10cm$cluster_9==6], pts_10cm$frags_10cm[pts_10cm$management_type=='Organic' & pts_10cm$tillage=='till' & pts_10cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))

#50 cm data
#%C by tillage and all samples
vioplot(pts_50cm$totC_50cm[pts_50cm$tillage=='till'], pts_50cm$totC_50cm[pts_50cm$tillage=='no till'], names = c('Tillage', 'No-till'))
#tilled samples with %C data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$totC_50cm) & pts_50cm$tillage=='till')])
#no till samples with %C data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$totC_50cm) & pts_50cm$tillage=='no till')])
#%C by tillage for cluster 6 of 9
vioplot(pts_50cm$totC_50cm[pts_50cm$tillage=='till'& pts_50cm$cluster_9==6], pts_50cm$totC_50cm[pts_50cm$tillage=='no till'& pts_50cm$cluster_9==6], names = c('Tillage', 'No-till'))
#C content by tillage and all samples
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$tillage=='till'], pts_50cm$kgOrg.m2_50cm[pts_50cm$tillage=='no till'], names = c('Tillage', 'No-till'))
#tilled samples with C content data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$kgOrg.m2_50cm) & pts_50cm$tillage=='till')])
#no till samples with C content data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$kgOrg.m2_50cm) & pts_50cm$tillage=='no till')])
#C content by tillage for cluster 6 of 9
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$tillage=='till'& pts_50cm$cluster_9==6], pts_50cm$kgOrg.m2_50cm[pts_50cm$tillage=='no till'& pts_50cm$cluster_9==6], names = c('Tillage', 'No-till'))
#gravel content by tillage for cluster 6 of 9
vioplot(pts_50cm$frags_50cm[pts_50cm$tillage=='till'& pts_50cm$cluster_9==6], pts_50cm$frags_50cm[pts_50cm$tillage=='no till'& pts_50cm$cluster_9==6], names = c('Tillage', 'No-till'))
#clay by tillage for cluster 6 of 9
vioplot(pts_50cm$clay_50cm[pts_50cm$tillage=='till'& pts_50cm$cluster_9==6], pts_50cm$clay_50cm[pts_50cm$tillage=='no till'& pts_50cm$cluster_9==6], names = c('Tillage', 'No-till'))

#%C by org. vs. conv. and all samples
vioplot(pts_50cm$totC_50cm[pts_50cm$management_type=='conventional'], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic'], names = c('Conv.', 'Organic'))
#conv. samples with %C data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$totC_50cm) & pts_50cm$management_type=='conventional')])
#org. samples with %C data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$totC_50cm) & pts_50cm$management_type=='Organic')])
#%C by conv. & org. for cluster 6 of 9
vioplot(pts_50cm$totC_50cm[pts_50cm$management_type=='conventional' & pts_50cm$cluster_9==6], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic' & pts_50cm$cluster_9==6], names = c('Conv.', 'Organic'))
#%C by tillage x mgmt type for all samples
vioplot(pts_50cm$totC_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='no till'], pts_50cm$totC_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='till'], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='no till'], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='till'], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#%C by tillage x mgmt type for cluster 6 of 9
vioplot(pts_50cm$totC_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$totC_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$totC_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))

#C content by conv. & org. for all samples
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional'], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic'], names = c('Conv.', 'Organic'))
#conv. samples with C content data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$kgOrg.m2_50cm) & pts_50cm$management_type=='conventional')])
#org. samples with C content data counts by cluster 9 ID
table(pts_50cm$cluster_9[which(!is.na(pts_50cm$kgOrg.m2_50cm) & pts_50cm$management_type=='Organic')])
#C content by org. vs. conv. for cluster 6 of 9
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional' & pts_50cm$cluster_9==6], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic' & pts_50cm$cluster_9==6], names = c('Conv.', 'Organic'))
#C content by tillage x mgmt type for all samples
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='no till'], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='till'], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='no till'], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='till'], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#C content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$kgOrg.m2_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))
#gravel content by tillage x mgmt type for cluster 6 of 9
vioplot(pts_50cm$frags_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$frags_50cm[pts_50cm$management_type=='conventional' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], pts_50cm$frags_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='no till' & pts_50cm$cluster_9==6], pts_50cm$frags_50cm[pts_50cm$management_type=='Organic' & pts_50cm$tillage=='till' & pts_50cm$cluster_9==6], names = c('Conv. \nno-till', 'Conv. \ntillage', 'Organic \nno-till', 'Organic \ntillage'))


order_lgnd_9 <- c(1,3,6,5,9,4,7,8,2)
points_by_cluster <- data.frame(name=1:9, cluster_count=0)
cluster_count <- as.data.frame(table(pts_30cm$cluster_9[!is.na(pts_30cm$totC_30cm)]))
points_by_cluster$cluster_count[points_by_cluster$name%in%cluster_count$Var1] <- cluster_count$Freq
points_by_cluster

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', 'CDFA_samples_haveCdata_cluster9.tif'), pointsize = 11, family = 'Times New Roman', width = 6.5, height = 5, units = 'in', res=800, compression = 'lzw')
par(mar=c(8, 4, 1, 1))
barplot(points_by_cluster$cluster_count[order_lgnd_9], col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[order_lgnd_9], ylab = 'Number of sample points w/ 0-30 cm C data', legend.text=c('1. Sandy soils', '9. Shrink-swell clays', '2. Loams w/ no res. & mod OM-low SS', '6. Loams w/ res. & high OM', '4. Loams w/ res. & low OM', '3. Loams w/ no res. & mod OM-mod SS', '7. Saline-sodic loams', '8. Saline-sodic clays', '5. Loams w/ res. & mod OM')[order_lgnd_9], cex.axis = 1, cex.names = 1, cex.lab = 1, args.legend = list(x=11, y=-5, cex=1, ncol=2))
dev.off()

#cluster bar plot by organic, conventional, tilled, and no-till


#look at effects by land management
table(pts_30cm$tillage[which(pts_30cm$cluster_9==6)])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$tillage=='till')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$tillage=='no till')])

table(pts_30cm$tillage[which(pts_30cm$cluster_9==5)])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==5 & pts_30cm$tillage=='till')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==5 & pts_30cm$tillage=='no till')])


table(pts_30cm$compost_added[which(pts_30cm$cluster_9==6)])
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$compost_added=='yes')])
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$compost_added=='none')])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$compost_added=='yes')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$compost_added=='none')])

table(pts_30cm$compost_added[which(pts_30cm$cluster_9==5)])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==5 & pts_30cm$compost_added=='yes')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==5 & pts_30cm$compost_added=='none')])


table(pts_30cm$irrigated_vs_dryfarm[pts_30cm$cluster_9==6])
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$irrigated_vs_dryfarm=='irrigated')])
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$irrigated_vs_dryfarm=='dryfarm')])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$irrigated_vs_dryfarm=='irrigated')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$irrigated_vs_dryfarm=='dryfarm')])

table(pts_30cm$irrigated_vs_dryfarm[pts_30cm$cluster_9==6])

table(pts_30cm$management_type[pts_30cm$cluster_9==6])
unique(pts_30cm$management_type)
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$management_type=='conventional')])
summary(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$management_type=='Organic')])
t.test(pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$management_type=='conventional')], pts_30cm$totC_30cm[which(pts_30cm$cluster_9==6 & pts_30cm$management_type=='Organic')])

table(pts_30cm$management_type[pts_30cm$cluster_9==5])

#look at 10 cm effects
#look at management effects at 10 cm on total C
table(pts_10cm$cluster_9)

table(pts_10cm$tillage[which(pts_10cm$cluster_9==6)])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='till')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='no till')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='till')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='no till')])

table(pts_10cm$tillage[which(pts_10cm$cluster_9==5)])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='till')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='no till')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='till')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='no till')])


table(pts_10cm$compost_added[which(pts_10cm$cluster_9==6)])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='yes')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='none')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='yes')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='none')])

table(pts_10cm$compost_added[which(pts_10cm$cluster_9==5)])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='yes')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='none')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='yes')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='none')])

table(pts_10cm$irrigated_vs_dryfarm[pts_10cm$cluster_9==6])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='irrigated')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='dryfarm')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='irrigated')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='dryfarm')])

table(pts_10cm$irrigated_vs_dryfarm[pts_10cm$cluster_9==6])

table(pts_10cm$management_type[pts_10cm$cluster_9==6])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='conventional')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='Organic')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='conventional')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='Organic')])

table(pts_10cm$management_type[pts_10cm$cluster_9==5])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='conventional')])
summary(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='Organic')])
t.test(pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='conventional')], pts_10cm$totC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='Organic')])

#look at management effects on DOC through 10 cm
table(pts_10cm$tillage[which(pts_10cm$cluster_9==6)])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='till')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='no till')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='till')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$tillage=='no till')])

table(pts_10cm$tillage[which(pts_10cm$cluster_9==5)])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='till')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='no till')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='till')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$tillage=='no till')])


table(pts_10cm$compost_added[which(pts_10cm$cluster_9==6)])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='yes')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='none')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='yes')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$compost_added=='none')])

table(pts_10cm$compost_added[which(pts_10cm$cluster_9==5)])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='yes')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='none')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='yes')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$compost_added=='none')])

table(pts_10cm$irrigated_vs_dryfarm[pts_10cm$cluster_9==6])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='irrigated')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='dryfarm')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='irrigated')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$irrigated_vs_dryfarm=='dryfarm')])

t.test(pts_10cm$DOC_10cm[which(pts_10cm$management_type=='conventional')], pts_10cm$DOC_10cm[which(pts_10cm$management_type=='Organic')])
table(pts_30cm[,c('management_type', 'cluster_9')])
table(pts_10cm$management_type[pts_10cm$cluster_9==6])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='conventional')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='Organic')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='conventional')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==6 & pts_10cm$management_type=='Organic')])

table(pts_10cm$management_type[pts_10cm$cluster_9==5])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='conventional')])
summary(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='Organic')])
t.test(pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='conventional')], pts_10cm$DOC_10cm[which(pts_10cm$cluster_9==5 & pts_10cm$management_type=='Organic')])
