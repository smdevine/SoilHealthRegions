laptop <- TRUE
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
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}

cluster_fk <- read.csv(file.path(dataDir, 'v2 results', 'clusters2_to_12_df_for_radarchart.csv'), stringsAsFactors = FALSE)
cluster_7_df <- cluster_fk[cluster_fk$clusters==7,1:10]
radarchart_7 <- rbind(apply(cluster_7_df, 2, max), apply(cluster_7_df, 2, min), cluster_7_df)
clus_7_colors <- c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3')
order_lgnd_7 <- c(4,6,1,7,3,2,5)
clus_7_lines <- c(3,2,2,1,2,1,3)
clus_7_names <- c('3. Coarse with pans', '6. Fine saline-sodic', '5. Coarse saline-sodic', '1. Coarse with no restrictions', '7. Fine shrink-swell', '2. Loamy with no restrictions', '4. Loamy with pans')
#clus_7_names <- c('3. Very coarse w/pans', '6. Loamy saline-sodic', '5. Very coarse saline-sodic', '1. Very coarse w/no restrictions', '7. Loamy shrink-swell', '2. Coarse w/no resrictions', '4. Coarse w/pans') #this is v1
clus_7_names[order_lgnd_7]
fix_legend_columns <- function(x) {
  c(x[1:2], NA, x[3:4], NA, x[5:7])
}
fix_legend_columns(clus_7_names[order_lgnd_7])
tiff(file = file.path(FiguresDir, 'v2', '7 region', 'valley_7_classes_spider.tif'), family = 'Times New Roman', width = 9, height = 6, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(2.5, 0.1, 0.1, 0.1))
radarchart(radarchart_7[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = clus_7_lines, pcol = clus_7_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
legend(x=-1.6, y=-1.2, legend = fix_legend_columns(clus_7_names[order_lgnd_7]), col=fix_legend_columns(clus_7_colors[order_lgnd_7]), lty=fix_legend_columns(clus_7_lines[order_lgnd_7]), lwd = 3, ncol = 3, bty='n')
dev.off()

#CDFA C validation for 7-region model
#compare 3 management types in cluster 6 across each depth segment
pts_30cm <- read.csv(file.path(kerriDir, 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE) 
pts_10cm <- read.csv(file.path(kerriDir, 'CDFA_samples_cluster_10cm.csv'), stringsAsFactors = FALSE)
pts_50cm <- read.csv(file.path(kerriDir, 'CDFA_samples_cluster_50cm.csv'), stringsAsFactors = FALSE) 
table(pts_30cm$cluster_7)
# columnClass <- c(Concatenate='character')
soil_data <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', "")) #colClasses = columnClass,
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))
soil_data_cluster6_of_7 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_7==6)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_7==6)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_7==6)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_7==6)], stringsAsFactors = FALSE)
head(soil_data_cluster6_of_7)
dim(soil_data_cluster6_of_7)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_cluster6_of_7$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster6_of_7$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster6_of_7$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster6_of_7$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster6_of_7$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster6_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
colnames(soil_data_cluster6_of_7)
lapply(soil_data_cluster6_of_7[,5:9], summary)
lapply(soil_data_cluster6_of_7[,5:9], function(x) tapply(x, soil_data_cluster6_of_7$tillage, mean, na.rm=TRUE))
lapply(soil_data_cluster6_of_7[,5:9], function(x) tapply(x, soil_data_cluster6_of_7$irrigation, mean, na.rm=TRUE))
lapply(soil_data_cluster6_of_7[,5:9], function(x) tapply(x, soil_data_cluster6_of_7$compost, mean, na.rm=TRUE))


till_irr_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='till' & soil_data_cluster6_of_7$irrigation=='irrigated' & soil_data_cluster6_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='till' & soil_data_cluster6_of_7$irrigation=='irrigated' & soil_data_cluster6_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))
notill_irr_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='no till' & soil_data_cluster6_of_7$irrigation=='irrigated' & soil_data_cluster6_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='no till' & soil_data_cluster6_of_7$irrigation=='irrigated' & soil_data_cluster6_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))
till_dry_clus6 <- data.frame(C_means=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='till' & soil_data_cluster6_of_7$irrigation=='dryfarm' & soil_data_cluster6_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster6_of_7[which(soil_data_cluster6_of_7$tillage=='till' & soil_data_cluster6_of_7$irrigation=='dryfarm' & soil_data_cluster6_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))


soil_data_cluster4_of_7 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_7==4)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_7==4)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_7==4)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_7==4)], stringsAsFactors = FALSE)
soil_data_cluster4_of_7$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster4_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster4_of_7$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster4_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster4_of_7$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster4_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster4_of_7$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster4_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster4_of_7$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster4_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
till_irr_clus4 <- data.frame(C_means=sapply(soil_data_cluster4_of_7[which(soil_data_cluster4_of_7$tillage=='till' & soil_data_cluster4_of_7$irrigation=='irrigated' & soil_data_cluster4_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster4_of_7[which(soil_data_cluster4_of_7$tillage=='till' & soil_data_cluster4_of_7$irrigation=='irrigated' & soil_data_cluster4_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))

soil_data_cluster7_of_7 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_7==7)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_7==7)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_7==7)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_7==7)], stringsAsFactors = FALSE)
soil_data_cluster7_of_7$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster7_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster7_of_7$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster7_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster7_of_7$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster7_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster7_of_7$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster7_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster7_of_7$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster7_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
notill_irr_clus7 <- data.frame(C_means=sapply(soil_data_cluster7_of_7[which(soil_data_cluster7_of_7$tillage=='no till' & soil_data_cluster7_of_7$irrigation=='irrigated' & soil_data_cluster7_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster7_of_7[which(soil_data_cluster7_of_7$tillage=='no till' & soil_data_cluster7_of_7$irrigation=='irrigated' & soil_data_cluster7_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))

soil_data_cluster1_of_7 <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$cluster_7==1)], tillage=pts_30cm$tillage[which(pts_30cm$cluster_7==1)], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$cluster_7==1)], compost=pts_30cm$compost_added[which(pts_30cm$cluster_7==1)], stringsAsFactors = FALSE)
soil_data_cluster1_of_7$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_cluster1_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_cluster1_of_7$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_cluster1_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_cluster1_of_7$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_cluster1_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_cluster1_of_7$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_cluster1_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_cluster1_of_7$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_cluster1_of_7$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
till_irr_clus5 <- data.frame(C_means=sapply(soil_data_cluster1_of_7[which(soil_data_cluster1_of_7$tillage=='till' & soil_data_cluster1_of_7$irrigation=='irrigated' & soil_data_cluster1_of_7$compost=='yes'),5:9], mean), C_se=sapply(soil_data_cluster1_of_7[which(soil_data_cluster1_of_7$tillage=='till' & soil_data_cluster1_of_7$irrigation=='irrigated' & soil_data_cluster1_of_7$compost=='yes'),5:9], function(x) {sd(x)/sqrt(length(x[!is.na(x)]))}))


depths_cm <- c(2.5, 7.5, 20, 40, 75)
clus_7_names
tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_soilC_profile.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(till_irr_clus6$C_means, -depths_cm, type='b', xlim=c(0.1,2.6), pch=24, bg='tan4', col='black', cex=1, xlab='', ylab='', yaxt='n')
axis(side = 2, at=seq(from=-100, to=0, by=20), labels=-seq(from=-100, to=0, by=20))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(till_dry_clus6$C_means, -depths_cm, type='b', pch=24, bg='tan4', col='black', lty=2)
points(notill_irr_clus6$C_means, -depths_cm, type='b', pch=23, bg='tan4', col='black')
points(till_irr_clus4$C_means, -depths_cm, type = 'b', pch=24, bg='tan2', col='black')
points(till_irr_clus5$C_means, -depths_cm, type = 'b', pch=24, bg='gold', col='black')
points(notill_irr_clus7$C_means, -depths_cm, type='b', pch=23, bg='firebrick3', col='black')
legend('bottomright', legend=c('1. Coarse w/no restrictions', '2a. Loamy w/no restrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'), pch=c(24,24,23,24,24,23), lty=c(1,1,1,2,1,1), pt.bg=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), col='black')
dev.off()
#option 2 names: c('1. Very coarse w/no restrictions', '2a. Coarse w/no resrictions', '2b. Coarse w/no resrictions: No-till', '2c. Coarse w/no resrictions: Dry farm', '3. Very coarse w/pans', '4. Coarse w/pans: No-till')

#make barplot of 30cm SOC kg m^-2
plot_da_error <- function(barnum, barcenters, df) {
  segments(x0=barcenters[barnum,], y0=df$means[barnum] + qnorm(0.975) * df$sd[barnum] / sqrt(df$n[barnum]), x1=barcenters[barnum,], y1=df$means[barnum] + qnorm(0.025) * df$sd[barnum] / sqrt(df$n[barnum]), lwd = 1.2)
}
kgOrgC_30cm_data <- data.frame(means=c(mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==4)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_7==6)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==1)], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==7)], na.rm=TRUE)), sd=c(sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==4)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_7==6)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==1)], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==7)], na.rm=TRUE)), n=c(sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==4)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$cluster_7==6)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==1)])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$cluster_7==7)]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_soilC_kg_0_30cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_30cm <- barplot(height=kgOrgC_30cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,8.5))
mtext(text=expression('0-30 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(2, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(3, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(4, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(5, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(6, barcenters_30cm, kgOrgC_30cm_data)
dev.off()

#plot 50 cm C contents
kgOrgC_50cm_data <- data.frame(means=c(mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==4)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_7==6)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==1)], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==7)], na.rm=TRUE)), sd=c(sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==4)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_7==6)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==1)], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==7)], na.rm=TRUE)), n=c(sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==4)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_7==6)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==1)])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==7)]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_soilC_kg_0_50cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_50cm <- barplot(height=kgOrgC_50cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,13))
mtext(text=expression('0-50 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(2, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(3, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(4, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(5, barcenters_50cm, kgOrgC_50cm_data)
plot_da_error(6, barcenters_50cm, kgOrgC_50cm_data)
text(x=barcenters_50cm, y=1, labels = c('A', 'BC', 'BC', 'C', 'AB', 'BC'), adj=0.5)
dev.off()

#plot 10 cm C contents
kgOrgC_10cm_data <- data.frame(means=c(mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==4)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_7==6)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==1)], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==7)], na.rm=TRUE)), sd=c(sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==4)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_7==6)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==1)], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==7)], na.rm=TRUE)), n=c(sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==4)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_7==6)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==1)])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==7)]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_soilC_kg_0_10cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_10cm <- barplot(height=kgOrgC_10cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,3.9))
mtext(text=expression('0-10 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(2, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(3, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(4, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(5, barcenters_10cm, kgOrgC_10cm_data)
plot_da_error(6, barcenters_10cm, kgOrgC_10cm_data)
text(x=barcenters_10cm, y=0.2, labels = c('A', 'ABC', 'BC', 'C', 'AB', 'BC'), adj=0.5)
dev.off()


#analysis of variance
pts_50cm$aov_code <- ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==4, 'Coarse w/no restrictions', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6, 'Loamy w/no restrictions', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==6, 'Loamy w/no restrictions: no-till', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$cluster_7==6, 'Loamy w/no restrictions: dry farm', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==1, 'Coarse w/pans', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$cluster_7==7, 'Loamy w/pans: no-till', NA))))))

tapply(pts_50cm$kgOrg.m2_50cm, pts_50cm$aov_code, mean, na.rm=TRUE)
kgOrgC_50cm_data
summary(aov(kgOrg.m2_50cm ~ aov_code, data=pts_50cm))
TukeyHSD(aov(kgOrg.m2_50cm ~ aov_code, data=pts_50cm))

#10 cm test
pts_10cm$aov_code <- ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==4, 'Coarse w/no restrictions', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6, 'Loamy w/no restrictions', ifelse(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==6, 'Loamy w/no restrictions: no-till', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$cluster_7==6, 'Loamy w/no restrictions: dry farm', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==1, 'Coarse w/pans', ifelse(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$cluster_7==7, 'Loamy w/pans: no-till', NA))))))

tapply(pts_10cm$kgOrg.m2_10cm, pts_10cm$aov_code, mean, na.rm=TRUE)
kgOrgC_10cm_data
summary(aov(kgOrg.m2_10cm ~ aov_code, data=pts_10cm))
TukeyHSD(aov(kgOrg.m2_10cm ~ aov_code, data=pts_10cm))
