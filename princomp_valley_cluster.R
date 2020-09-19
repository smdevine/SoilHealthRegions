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
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
}
mar_settings <- c(4, 4.5, 1, 1)
ec_zero_rule <- 1.5
list.files(dataDir)
# valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm.shp'))
# names(valley_mu_shp_30cm)
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'for cluster analysis', 'valley30cm_by_mukey_final_v2.csv'), stringsAsFactors = FALSE) #re-run 9/17/20 for NbClust analysis


#create data.frame for cluster analysis
df_for_clustering <- valley30cm_by_mukey
colnames(df_for_clustering)
rownames(df_for_clustering) <- df_for_clustering$mukey
# summary(df_for_clustering)
df_for_clustering <- df_for_clustering[ ,c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm')] #'sar_30cm' no longer included; frags have 160 mukeys with NA
# "storiemn", "Lthc_dep", "Plth_dep", "Drpn_dep", "ATC_dep", "Natr_dep", "Salc_dep", "SCTS_dep", "MRes_dep" 'kwf_30cm'

mapply(function(x, y) hist(x, main=y), x=df_for_clustering, y=colnames(df_for_clustering))
lapply(df_for_clustering, class)
lapply(df_for_clustering, function(x) sum(is.na(x)))
lapply(df_for_clustering, function(x) sum(x==0, na.rm = TRUE)) #verify why some are showing min res depth of 0
lapply(df_for_clustering, function(x) sum(x<0, na.rm = TRUE))
lapply(df_for_clustering, function(x) sum(x>0 & x < 0.01, na.rm = TRUE))
lapply(df_for_clustering, function(x) sum(x>=0.01 & x < 0.1, na.rm = TRUE))
lapply(df_for_clustering, function(x) sum(x>=0.1 & x < 0.5, na.rm = TRUE))
lapply(df_for_clustering, function(x) summary(x))

log_transform <- function(x, df, c) {
  df[[x]] <- log(df[[x]] + c, 10)
  df
}
#modify EC data
df_for_clustering$ec_30cm[df_for_clustering$ec_30cm < ec_zero_rule] <- 0
#manual solutions for c
hist(df_for_clustering$clay_30cm)
test <- log((df_for_clustering$clay_30cm + 5), 10)
hist(test)
c_arg <- c(0, 5, 2, 5, 0)
var_list <- c('om_30cm', 'cec_30cm', 'ec_30cm', 'lep_30cm', 'ksat_30cm') #'sar_30cm'
for (i in seq_along(var_list)) {
  df_for_clustering <- log_transform(x=var_list[i], df = df_for_clustering, c=c_arg[i]) #above variables are log transformed
}
# par(mar=c(4, 4, 4, 4))
mapply(function(x, y) hist(x, main=y), x=df_for_clustering, y=colnames(df_for_clustering))
dim(df_for_clustering)
df_for_clustering_scaled <- as.data.frame(scale(df_for_clustering))

#princomp was used in soilC_spatial.R tests in random vs. stratified sampling exercise
princomp_results <- princomp(df_for_clustering_scaled)
dim(princomp_results$scores)
biplot(princomp_results)
summary(princomp_results)
plot(princomp_results)
loadings(princomp_results)
summary(princomp_results)
lapply(princomp_results)
lapply(as.data.frame(princomp_results$scores), function(x) sum((x-mean(x))^2))

kmeans_test <- function(y, z) {sapply(1:y, function(x) {
  result <- kmeans(z, centers = x, iter.max = 100, nstart = 25)
  round(100 * result$betweenss / result$totss, 1)})
}
results <- replicate(30, kmeans_test(12, princomp_results$scores[,1:5]))
plot(1:12, rowMeans(results), type='b', xlab='Number of clusters', ylab='Soil variability captured by clustering (%)', cex=0.8, cex.axis=1, cex.lab=1)

#col_index selects how many princomps to use
col_index <- 1:7
set.seed(11431030)
cluster_2_pc <- kmeans(princomp_results$scores[,col_index], centers=2, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_3_pc <- kmeans(princomp_results$scores[,col_index], centers=3, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_4_pc <- kmeans(princomp_results$scores[,col_index], centers=4, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_5_pc <- kmeans(princomp_results$scores[,col_index], centers=5, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_6_pc <- kmeans(princomp_results$scores[,col_index], centers=6, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_7_pc <- kmeans(princomp_results$scores[,col_index], centers=7, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_8_pc <- kmeans(princomp_results$scores[,col_index], centers=8, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_9_pc <- kmeans(princomp_results$scores[,col_index], centers=9, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_10_pc <- kmeans(princomp_results$scores[,col_index], centers=10, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_11_pc <- kmeans(princomp_results$scores[,col_index], centers=11, iter.max = 200, nstart = 50)
set.seed(11431030)
cluster_12_pc <- kmeans(princomp_results$scores[,col_index], centers=12, iter.max = 200, nstart = 50)

stdev_calc <- function(SS, n) {
  sqrt(SS / (n - 1))
}

cluster_fk_pc <- do.call(rbind, mapply(function(x, y) {
  z <- as.data.frame(x[['centers']])
  z$clusters <- y
  z$clus_cl <- 1:nrow(z)
  z$withinss <- x[['withinss']]
  z$n <- x[['size']]
  z$sdev <- stdev_calc(SS = z$withinss, n = z$n)
  z}, x = list(cluster_2_pc, cluster_3_pc, cluster_4_pc, cluster_5_pc, cluster_6_pc, cluster_7_pc, cluster_8_pc, cluster_9_pc, cluster_10_pc, cluster_11_pc, cluster_12_pc), y = 2:12, SIMPLIFY = FALSE))

cluster_fk_pc$avg_dist_within <- cluster_fk_pc$withinss / cluster_fk_pc$n

valley30cm_by_mukey_PC <- valley30cm_by_mukey
valley30cm_by_mukey_PC$cluster_2pc <- cluster_2_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_2_pc$cluster))]
valley30cm_by_mukey_PC$cluster_3pc <- cluster_3_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_3_pc$cluster))]
valley30cm_by_mukey_PC$cluster_4pc <- cluster_4_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_4_pc$cluster))]
valley30cm_by_mukey_PC$cluster_5pc <- cluster_5_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_5_pc$cluster))]
valley30cm_by_mukey_PC$cluster_6pc <- cluster_6_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_6_pc$cluster))]
valley30cm_by_mukey_PC$cluster_7pc <- cluster_7_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_7_pc$cluster))]
valley30cm_by_mukey_PC$cluster_8pc <- cluster_8_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_8_pc$cluster))]
valley30cm_by_mukey_PC$cluster_9pc <- cluster_9_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_9_pc$cluster))]
valley30cm_by_mukey_PC$cluster_10pc <- cluster_10_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_10_pc$cluster))]
valley30cm_by_mukey_PC$cluster_11pc <- cluster_11_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_11_pc$cluster))]
valley30cm_by_mukey_PC$cluster_12pc <- cluster_12_pc$cluster[match(valley30cm_by_mukey_PC$mukey, names(cluster_12_pc$cluster))]

cluster_pc_area_summary <- c(tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_2pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_3pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_4pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_5pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_6pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_7pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_8pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_9pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_10pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_11pc, sum), tapply(valley30cm_by_mukey_PC$area_ac, valley30cm_by_mukey_PC$cluster_12pc, sum))
cluster_fk_pc$area_ac <- cluster_pc_area_summary
cluster_fk_pc$area_pct <- 100 * cluster_fk_pc$area_ac / sum(valley30cm_by_mukey_PC$area_ac)
write.csv(cluster_fk_pc, file.path(dataDir, 'FINAL results', 'pc analysis', 'clusters_pc7_2_to_12_df_for_radarchart.csv'), row.names = FALSE)

cluster_no <- 7
for (i in 1:cluster_no) {
  print(i)
  print(summary(valley30cm_by_mukey_PC[[paste0('cluster_',cluster_no)]][valley30cm_by_mukey_PC[[paste0('cluster_', cluster_no, 'pc')]]==i]))
}
cluster_no <- 7
for (i in 1:cluster_no) {
  hist(valley30cm_by_mukey_PC[[paste0('cluster_',cluster_no)]][valley30cm_by_mukey_PC[[paste0('cluster_', cluster_no, 'pc')]]==i], main=i)
}
cluster_no <- 9
for (i in 1:cluster_no) {
  print(i)
  print(table(valley30cm_by_mukey_PC[[paste0('cluster_',cluster_no)]][valley30cm_by_mukey_PC[[paste0('cluster_', cluster_no, 'pc')]]==i]))
}

cluster_radar_plot <- function(cluster_no) {
  radarchart_df <- cluster_fk_pc[cluster_fk_pc$clusters==cluster_no,]
  radarchart_df <- rbind(apply(radarchart_df, 2, max), apply(radarchart_df, 2, min), radarchart_df) #have to create a data.frame where first row is max, second row is min, and then followed by the data for radarchart
  tiff(file = file.path(FiguresDir, 'FINAL', 'pc analysis', paste0('valley_', cluster_no, '_classes_spider.tif')), family = 'Times New Roman', width = 6.5, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
  par(mar=rep(0.1, 4))
  radarchart(radarchart_df[,1:10], vlabels=1:10, maxmin = TRUE, plwd = 3)
  legend(x = -1.25, y = -0.25, legend = 1:cluster_no, col=1:cluster_no, lty=1:cluster_no, lwd = 3)
  dev.off()
}
cluster_radar_plot(2)
cluster_radar_plot(3)
cluster_radar_plot(4)
cluster_radar_plot(5)
cluster_radar_plot(6)
cluster_radar_plot(7)
cluster_radar_plot(8)
cluster_radar_plot(9) #color starts to repeat
cluster_radar_plot(10)
cluster_radar_plot(11)
cluster_radar_plot(12)

#add silhouette distances
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_2pc, mean))
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_3pc, mean))
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_4pc, mean))
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_5pc, mean))
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_6pc, mean))
lapply(valley30cm_by_mukey_PC[,c(56:62, 65:66, 71, 73)], function(x) tapply(x, valley30cm_by_mukey_PC$cluster_7pc, mean))
#name order corrected 9/17/20
# clus_9_names <- c('4. Coarse w/pans', '5. Loamy w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions',  '7. Coarse saline-sodic', '2. Coarse w/no restrictions', '8. Fine saline-sodic', '6. Loamy w/pans (high OM)', '9. Fine shrink-swell') #order corrected 4/7/20
# clus_8_names <- c('6. Coarse saline-sodic', '7. Fine saline-sodic', '5. Loamy w/pans', '2. Coarse w/no restrictions', '8. Fine shrink-swell', '4. Coarse w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions') #order corrected 4/7/20
clus_7_names <- c('6. Fine salt-affected', '3. Low OM & restrictive horizons', '4. High OM & restrictive horizons', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy salt-affected')
clus_6_names <- c('4. High OM w/restrictive horizons', '1. Coarse w/no restrictions', '2. Loamy w/no restrictions', '6. Shrink-swell', '3. Low OM w/restrictive horizons', '5. Salt-affected') 
clus_5_names <- c('4. Salt affected', '5. Fine shrink-swell', '1. Coarse w/no restrictions', '2. Loamy w/no restrictions', '3. Soils w/restrictive horizons')
clus_4_names <- c( '4. Shrink-swell', '3. Salt-affected', '2. Loamy', '1. Coarse')
clus_3_names <- c('1. Loamy', '2. Coarse', '3. Shrink-swell')
clus_2_names <- c('1. Coarse', '2. Fine')
valley30cm_by_mukey_silhouette <- valley30cm_by_mukey_PC
add_silhouette_data <- function(clus_df, clus_no, clus_names) {
  result <- silhouette(clus_df$cluster, dist(princomp_results$scores[,col_index]))
  valley30cm_by_mukey_silhouette[[paste0('clus', clus_no, '_sil_width')]] <- result[,3]
  valley30cm_by_mukey_silhouette[[paste0('clus', clus_no, '_names')]] <- clus_names[valley30cm_by_mukey_silhouette[[paste0('cluster_', clus_no, 'pc')]]]
  valley30cm_by_mukey_silhouette
}
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_2_pc, '2', clus_2_names)
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_3_pc, '3', clus_3_names)
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_4_pc, '4', clus_4_names)
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_5_pc, '5', clus_5_names)
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_6_pc, '6', clus_6_names)
valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_7_pc, '7', clus_7_names)
# valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_8_pc, '8', clus_8_names)
# valley30cm_by_mukey_silhouette <- add_silhouette_data(cluster_9_pc, '9', clus_9_names)
clus7_silhouette <- silhouette(cluster_7$cluster, dist(df_for_clustering_scaled))
summary_sil7 <- summary(clus7_silhouette)
mean(summary_sil7$clus.avg.widths)
mean(clus7_silhouette[,3])
median(clus7_silhouette[,3])
plot(clus7_silhouette)

lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 85)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus2_names, mean))
lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 87)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus3_names, mean))
lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 89)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus4_names, mean))
lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 91)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus5_names, mean))
lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 93)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus6_names, mean))
lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 95)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus7_names, mean))
# lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 95)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus8_names, mean))
# lapply(valley30cm_by_mukey_silhouette[,c(37, 56:62, 65:66, 70:71, 73, 97)], function(x) tapply(x, valley30cm_by_mukey_silhouette$clus9_names, mean))
lapply(valley30cm_by_mukey_silhouette[,seq(85,99,2)], function(x) mean(x))
lapply(valley30cm_by_mukey_silhouette[,seq(85,99,2)], function(x) mean(x[valley30cm_by_mukey_silhouette$clus7_names!='4. High OM & restrictive horizons']))
lapply(valley30cm_by_mukey_silhouette[,seq(85,99,2)], function(x) mean(x[valley30cm_by_mukey_silhouette$clus7_names!='6. Loamy w/pans (high OM) ']))

#gap stat
gap_stats <- clusGap(princomp_results$scores[,col_index], FUN = kmeans, nstart = 50, K.max = 20, B = 100, iter.max = 200)

fviz_gap_stat(gap_stats, maxSE = list('firstSEmax', SE.factor=1))
fviz_gap_stat(gap_stats, maxSE = list('firstSEmax', SE.factor=3))
fviz_gap_stat(gap_stats, maxSE = list('Tibs2001SEmax'))

#test NbClust function
library(NbClust)
NbClust_result_euc <- NbClust(data=princomp_results$scores[,col_index], distance = 'euclidean', min.nc=2, max.nc = 20, method= 'kmeans', index = 'all')
write.csv(NbClust_result_euc$All.index, file.path(dataDir, 'FINAL results', 'pc analysis', 'NbClust_cluster_metrics.csv'), row.names = TRUE)
write.csv(NbClust_result_euc$Best.nc, file.path(dataDir, 'FINAL results', 'pc analysis', 'NbClust_optimal_cluster_no.csv'), row.names = TRUE)

#y nada cambia
best_cluster <- kmeansruns(data = princomp_results$scores, krange = 1:12, criterion = 'asw', iter.max = 200, runs = 100, nstart=50)
best_cluster
best_cluster$bestk #2!

gap_stats <- clusGap(princomp_results$scores, FUN = kmeans, nstart = 50, K.max = 12, B = 500, iter.max = 200)
fviz_gap_stat(gap_stats) #12!


#violin plots for 9 classes based 7 pc
order_lgnd_pc7 <- c(1,3,6,5,9,4,7,8,2)
vioplot_mod_clus9_pc <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:9)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', 'pc analysis', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_9pc==plot_order2[1]], times=round(df$area_ac[df$cluster_9pc==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[2]], times=round(df$area_ac[df$cluster_9pc==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[3]], times=round(df$area_ac[df$cluster_9pc==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[4]], times=round(df$area_ac[df$cluster_9pc==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[5]], times=round(df$area_ac[df$cluster_9pc==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[6]], times=round(df$area_ac[df$cluster_9pc==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[7]], times=round(df$area_ac[df$cluster_9pc==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[8]], times=round(df$area_ac[df$cluster_9pc==plot_order2[8]]/area_fact, 0)), rep(df[[varname]][df$cluster_9pc==plot_order2[9]], times=round(df$area_ac[df$cluster_9pc==plot_order2[9]]/area_fact, 0)), col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
#order_lgnd was defined for radarchart
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_pc7, area_fact = 10, ylab='Clay (%)', fname='class9_clay_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_pc7, area_fact = 10, ylab='Organic matter (%)', fname='class9_OM_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
valley30cm_by_mukey_PC$logom_30cm <- log(valley30cm_by_mukey_PC$om_30cm)
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class9_LogOM_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class9_CEC_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_pc7, area_fact = 10, ylab='soil pH', fname='class9_pH_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class9_AWC_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class9_EC_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class9_BD_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_pc7, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class9_MnRs_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus9(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
valley30cm_by_mukey_PC$logks_30cm <- log(valley30cm_by_mukey_PC$ksat_30cm)
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_pc7, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class9_logKs_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9_pc(valley30cm_by_mukey_PC, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_pc7, area_fact = 10, ylab='Linear extensibility (%)', fname='class9_lep_vioplots_pc7.tif', mar=c(3.5, 4.25, 1, 1))

var_list <- c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm')
for (i in 1:9) {
  print(i)
  print(summary(valley30cm_by_mukey$clay_30cm[valley30cm_by_mukey$cluster_9==i]))
  print(summary(valley30cm_by_mukey_PC$clay_30cm[valley30cm_by_mukey_PC$cluster_9pc==i]))
}
