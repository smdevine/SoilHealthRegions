#this modified to work with aggregated dataset from ssurgo_calag_aggregate.R
#TO-DO
#(1) run cluster analysis on Salinas only [DONE]
#(2) summarize classes for 4 and 5 
#(3) log transform om and ksat [DONE]
#(4) identify outliers
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
mar_settings <- c(3.5, 4.5, 1, 1)
ec_zero_rule <- 1.5
list.files(file.path(dataDir, 'for cluster analysis'))
# valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm.shp'))
# names(valley_mu_shp_30cm)
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'for cluster analysis', 'valley30cm_by_mukey_final_v2.csv'), stringsAsFactors = FALSE)


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
summary(df_for_clustering_scaled)
colnames(df_for_clustering_scaled)
mapply(function(x, y) hist(x, main=y), x=df_for_clustering_scaled, y=colnames(df_for_clustering_scaled))
lapply(df_for_clustering_scaled, function(x) summary(x))
lapply(df_for_clustering_scaled, function(x) sd(x))
lapply(na.omit(df_for_clustering_scaled), function(x) sum((x-mean(x))^2))
sum((df_for_clustering_scaled$MnRs_dep - mean(df_for_clustering_scaled$MnRs_dep))^2)

#look at relationship between scaled and untransformed data
plot(valley30cm_by_mukey$ec_30cm, df_for_clustering_scaled$ec_30cm)
plot(valley30cm_by_mukey$om_30cm, df_for_clustering_scaled$om_30cm)

#ec
df_exploration <- df_for_clustering
df_exploration$area_ac <- valley30cm_by_mukey$area_ac[match(rownames(df_exploration), valley30cm_by_mukey$mukey)]
sum(df_exploration$area_ac)
#classification scheme from https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044788.pdf
sum(df_exploration$area_ac[df_exploration$ec_30cm < 2]) / sum(df_exploration$area_ac) #76% considered 'non-saline'
sum(df_exploration$area_ac[df_exploration$ec_30cm >= 2 & df_exploration$ec_30cm < 4]) / sum(df_exploration$area_ac) #9.0% very slightly saline
sum(df_exploration$area_ac[df_exploration$ec_30cm >= 4 & df_exploration$ec_30cm < 8]) / sum(df_exploration$area_ac) #7.1% slightly saline
sum(df_exploration$area_ac[df_exploration$ec_30cm >= 8 & df_exploration$ec_30cm < 16]) / sum(df_exploration$area_ac) #6.3 moderately saline
sum(df_exploration$area_ac[df_exploration$ec_30cm >= 16]) / sum(df_exploration$area_ac) #1.2% strongly saline

kmeans_test <- function(y, z) {sapply(1:y, function(x) {
  result <- kmeans(z, centers = x, iter.max = 200, nstart = 50)
  round(100 * result$betweenss / result$totss, 1)})
}
results <- replicate(100, kmeans_test(20, df_for_clustering_scaled))
# dim(results)
# rowMeans(results)
# apply(results, 1, sd)
tiff(file = file.path(FiguresDir, 'v2', 'validation plots', 'kmeans_comparison_3.31.20.tif'), family = 'Times New Roman', width = 6.5, height = 3.5, pointsize = 11, units = 'in', res=800, compression='lzw')
par(mar=c(2, 4, 0.5, 0.5))
plot(1:20, rowMeans(results), type='b', xlab='', ylab='', cex=0.8, cex.axis=1, cex.lab=1, ylim=c(0,82), xaxt='n')
mtext(text = 'Number of soil health regions (cluster size in conceptual model)', side=1, line=0.75)
mtext(text = 'SSURGO variability captured by clusters (%)', side=2, line=2.5, at=35)
text(1:20, rowMeans(results), labels=as.character(1:20), pos=3, offset=0.5)
text(x=1, y=80, 'a', adj=c(0,0))
dev.off()

set.seed(11431030)
cluster_2 <- kmeans(df_for_clustering_scaled, centers=2, iter.max = 200, nstart = 50)
cluster_2
set.seed(11431030)
cluster_3 <- kmeans(df_for_clustering_scaled, centers=3, iter.max = 200, nstart = 50)
cluster_3
set.seed(11431030)
cluster_4 <- kmeans(df_for_clustering_scaled, centers=4, iter.max = 200, nstart = 50)
cluster_4
set.seed(11431030)
cluster_5 <- kmeans(df_for_clustering_scaled, centers=5, iter.max = 200, nstart = 50)
cluster_5
set.seed(11431030)
cluster_6 <- kmeans(df_for_clustering_scaled, centers=6, iter.max = 200, nstart = 50)
cluster_6
set.seed(11431030)
cluster_7 <- kmeans(df_for_clustering_scaled, centers=7, iter.max = 200, nstart = 50)
cluster_7
set.seed(11431030)
cluster_8 <- kmeans(df_for_clustering_scaled, centers=8, iter.max = 200, nstart = 50)
cluster_8
set.seed(11431030)
cluster_9 <- kmeans(df_for_clustering_scaled, centers=9, iter.max = 200, nstart = 50)
cluster_9
set.seed(11431030)
cluster_10 <- kmeans(df_for_clustering_scaled, centers=10, iter.max = 200, nstart = 50)
cluster_10
set.seed(11431030)
cluster_11 <- kmeans(df_for_clustering_scaled, centers=11, iter.max = 200, nstart = 50)
cluster_11
set.seed(11431030)
cluster_12 <- kmeans(df_for_clustering_scaled, centers=12, iter.max = 200, nstart = 50)
cluster_12

#bind them
stdev_calc <- function(SS, n) {
  sqrt(SS / (n - 1))
}
# x <- cluster_2
# y <- 2
cluster_fk <- do.call(rbind, mapply(function(x, y) {
  z <- as.data.frame(x[['centers']])
  z$clusters <- y
  z$clus_cl <- 1:nrow(z)
  z$withinss <- x[['withinss']]
  z$n <- x[['size']]
  z$sdev <- stdev_calc(SS = z$withinss, n = z$n)
  z}, x = list(cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, cluster_8, cluster_9, cluster_10, cluster_11, cluster_12), y = 2:12, SIMPLIFY = FALSE))

# var_quantiles <- as.data.frame(apply(na.omit(df_for_clustering_scaled)[,1:11], 2, function(x) quantile(x, probs=c(0.1, 0.25, 0.75, 0.9)))) #use whole dataset to establish breaks

# cluster_fk_labels <- as.data.frame(mapply(function(x, y, z) if(z=='MnRs_dep') {ifelse(x < y[1], 'v. shallow', ifelse(x < y[2], 'shallow', ifelse(x < y[3], 'mod', ifelse(x < y[4], 'deep', 'v. deep'))))} else {ifelse(x < y[1], 'v. low', ifelse(x < y[2], 'low', ifelse(x < y[3], 'mod', ifelse(x < y[4], 'high', 'v. high'))))}, x = cluster_fk[,1:11], y=var_quantiles, z=colnames(var_quantiles), SIMPLIFY = FALSE)) #sar, lep, ec, and res dep have problematic distributions
# colnames(cluster_fk)
# dim(cluster_fk)
cluster_fk$avg_dist_within <- cluster_fk$withinss / cluster_fk$n
# cbind(cluster_fk[cluster_fk$clusters==5, 13:16], cluster_fk_labels[(nrow(cluster_fk_labels)-4):nrow(cluster_fk_labels), ])

#make a radarchart plotting function
cluster_radar_plot <- function(cluster_no) {
  radarchart_df <- cluster_fk[cluster_fk$clusters==cluster_no,]
  radarchart_df <- rbind(apply(radarchart_df, 2, max), apply(radarchart_df, 2, min), radarchart_df) #have to create a data.frame where first row is max, second row is min, and then followed by the data for radarchart
  tiff(file = file.path(FiguresDir, 'v2', paste0('valley_', cluster_no, '_classes_spider.tif')), family = 'Times New Roman', width = 6.5, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
  par(mar=rep(0.1, 4))
  radarchart(radarchart_df[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'Saline', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
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

valley30cm_by_mukey$cluster_2 <- cluster_2$cluster[match(valley30cm_by_mukey$mukey, names(cluster_2$cluster))]
valley30cm_by_mukey$cluster_3 <- cluster_3$cluster[match(valley30cm_by_mukey$mukey, names(cluster_3$cluster))]
valley30cm_by_mukey$cluster_4 <- cluster_4$cluster[match(valley30cm_by_mukey$mukey, names(cluster_4$cluster))]
valley30cm_by_mukey$cluster_5 <- cluster_5$cluster[match(valley30cm_by_mukey$mukey, names(cluster_5$cluster))]
valley30cm_by_mukey$cluster_6 <- cluster_6$cluster[match(valley30cm_by_mukey$mukey, names(cluster_6$cluster))]
valley30cm_by_mukey$cluster_7 <- cluster_7$cluster[match(valley30cm_by_mukey$mukey, names(cluster_7$cluster))]
valley30cm_by_mukey$cluster_8 <- cluster_8$cluster[match(valley30cm_by_mukey$mukey, names(cluster_8$cluster))]
valley30cm_by_mukey$cluster_9 <- cluster_9$cluster[match(valley30cm_by_mukey$mukey, names(cluster_9$cluster))]
valley30cm_by_mukey$cluster_10 <- cluster_10$cluster[match(valley30cm_by_mukey$mukey, names(cluster_10$cluster))]
valley30cm_by_mukey$cluster_11 <- cluster_11$cluster[match(valley30cm_by_mukey$mukey, names(cluster_11$cluster))]
valley30cm_by_mukey$cluster_12 <- cluster_12$cluster[match(valley30cm_by_mukey$mukey, names(cluster_12$cluster))]
write.csv(valley30cm_by_mukey, file.path(dataDir, 'v2 results', 'valley30cm_by_mukey_cluster_v2.csv'), row.names = FALSE) #this produced a different order of cluster labels even with setting seed--perhaps because R version has changed since last run?

#read in cluster valley file by mukey with cluster info
valley30cm_by_mukey_orig <- read.csv(file.path(dataDir, 'v2 results', 'valley30cm_by_mukey_cluster.csv'), stringsAsFactors = FALSE)
all(valley30cm_by_mukey_orig$mukey==valley30cm_by_mukey$mukey) #it's in same order
head(cbind(valley30cm_by_mukey$cluster_7, valley30cm_by_mukey_orig$cluster_7), 20)
sum(3!=valley30cm_by_mukey$cluster_7[valley30cm_by_mukey_orig$cluster_7==7]) #14 are different
sum(5!=valley30cm_by_mukey$cluster_7[valley30cm_by_mukey_orig$cluster_7==6]) #2 are different
sum(2!=valley30cm_by_mukey$cluster_7[valley30cm_by_mukey_orig$cluster_7==1]) #3 are different
sum(4!=valley30cm_by_mukey$cluster_7[valley30cm_by_mukey_orig$cluster_7==4]) #0 are different

#conclusion: there was a change in the set.seed as a result of updating R version
colnames(valley30cm_by_mukey)[74:84]
colnames(valley30cm_by_mukey_orig)[74:84]
valley30cm_by_mukey[,74:84] <- valley30cm_by_mukey_orig[,74:84]
write.csv(valley30cm_by_mukey, file.path(dataDir, 'v2 results', 'valley30cm_by_mukey_cluster_v2.csv'), row.names = FALSE)

#add cluster info to shapefile
valley_mu_shp_30cm$cluster_2 <- valley30cm_by_mukey$cluster_2[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
tapply(valley_mu_shp_30cm$area_ac, valley_mu_shp_30cm$cluster_2, sum)
valley_mu_shp_30cm$cluster_3 <- valley30cm_by_mukey$cluster_3[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_4 <- valley30cm_by_mukey$cluster_4[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_5 <- valley30cm_by_mukey$cluster_5[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_6 <- valley30cm_by_mukey$cluster_6[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_7 <- valley30cm_by_mukey$cluster_7[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_8 <- valley30cm_by_mukey$cluster_8[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_9 <- valley30cm_by_mukey$cluster_9[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_10 <- valley30cm_by_mukey$cluster_10[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_11 <- valley30cm_by_mukey$cluster_11[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
valley_mu_shp_30cm$cluster_12 <- valley30cm_by_mukey$cluster_12[match(valley_mu_shp_30cm$mukey, valley30cm_by_mukey$mukey)]
sum(tapply(valley_mu_shp_30cm$area_ac, valley_mu_shp_30cm$cluster_10, sum))
sum(valley_mu_shp_30cm$area_ac[is.na(valley_mu_shp_30cm$cluster_10)])
# shapefile(valley_mu_shp_30cm, file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'), overwrite=TRUE)

cluster_area_summary <- c(tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_2, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_3, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_4, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_5, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_6, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_7, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_8, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_9, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_10, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_11, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_12, sum))
#area check
sum(cluster_area_summary[1:2])
sum(valley30cm_by_mukey$area_ac)
cluster_fk$area_ac <- cluster_area_summary
cluster_fk$area_pct <- 100 * cluster_fk$area_ac / sum(valley30cm_by_mukey$area_ac)
# write.csv(cluster_fk, file.path(dataDir, 'v2 results', 'clusters2_to_12_df_for_radarchart.csv'), row.names = FALSE)
# cluster_fk <- read.csv(file.path(dataDir, 'v2', 'clusters2_to_12_df_for_radarchart.csv'))

# colnames(cluster_fk)
cluster_9_df <- cluster_fk[cluster_fk$clusters==9,1:10]
radarchart_9 <- rbind(apply(cluster_9_df, 2, max), apply(cluster_9_df, 2, min), cluster_9_df)

#9-class cluster option
tiff(file = file.path(FiguresDir, 'v2', 'valley_9_classes_spider_11.8.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_9[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = c(1, 2, 1, 3, 3, 1, 2, 2, 3), pcol = c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3'), vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
order_lgnd <- c(1,3,6,5,9,4,7,8,2)
legend(x=-1.8, y=-1.2, legend = c('1. Sandy soils', '9. Shrink-swell clays', '2. Loams w/ no res. & mod OM-low SS', '6. Loams w/ res. & high OM', '4. Loams w/ res. & low OM', '3. Loams w/ no res. & mod OM-mod SS', '7. Saline-sodic loams', '8. Saline-sodic clays', '5. Loams w/ res. & mod OM')[order_lgnd], col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[order_lgnd], lty=c(1, 2, 1, 3, 3, 1, 2, 2, 3)[order_lgnd], lwd = 3, ncol = 3, bty='n')
dev.off()

# #rgb color codes
col2rgb(c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3'))

#revised 10 class option
cluster_10_df <- cluster_fk[cluster_fk$clusters==10,1:10]
radarchart_10 <- rbind(apply(cluster_10_df, 2, max), apply(cluster_10_df, 2, min), cluster_10_df)
clus_10_colors <- c('violetred', 'tan4', 'lightgoldenrod', 'black', 'sandybrown', 'deepskyblue', 'gold', 'firebrick3', 'tan2', 'lightblue1')
col2rgb(clus_10_colors)
order_lgnd_10 <- c(3,5,9,2,7,8,4,10,6,1)

tiff(file = file.path(FiguresDir, 'v2', 'valley_10_classes_spider_11.11.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_10[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = c(2,1,1,3,1,2,3,3,1,2), pcol = clus_10_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
order_lgnd_10 <- c(3,5,9,2,7,8,4,10,6,1)
legend(x=-1.8, y=-1.2, legend = c('10. Shrink-swell clays', '4. Loams w/no res & mod OM-mod SS', '1. Sandy soils', '7. Loams w/res & high OM', '2. Loams w/no res & low OM', '9. Saline-sodic clays', '5. Loams w/res & low OM', '6. Loams w/res & mod OM', '3. Loams w/no res & mod OM-low SS', '8. Saline-sodic loams')[order_lgnd_10], col=clus_10_colors[order_lgnd_10], lty=c(2,1,1,3,1,2,3,3,1,2)[order_lgnd_10], lwd = 3, ncol = 3, bty='n')
dev.off()

#revised 7-class option
cluster_7_df <- cluster_fk[cluster_fk$clusters==7,1:10]
radarchart_7 <- rbind(apply(cluster_7_df, 2, max), apply(cluster_7_df, 2, min), cluster_7_df)
clus_7_colors <- c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3')
order_lgnd_7 <- c(4,6,1,7,3,2,5)
clus_7_lines <- c(3,2,2,1,2,1,3)
tiff(file = file.path(FiguresDir, 'v2', 'valley_7_classes_spider_11.15.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_7[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = clus_7_lines, pcol = clus_7_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
legend(x=-1.6, y=-1.2, legend = c('3. Loams w/res & low OM', '6. Saline-sodic clays', '5. Saline-sodic loams', '1. Sandy loam soils', '7. Shrink-swell clays', '2. Loams w/no res', '4. Loams w/ res. & mod OM')[order_lgnd_7], col=clus_7_colors[order_lgnd_7], lty=clus_7_lines[order_lgnd_7], lwd = 3, ncol = 3, bty='n')
dev.off()
#Loams w/no res could have mod OM-mod SS added
col2rgb(c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3'))

#revised 8-class option
cluster_8_df <- cluster_fk[cluster_fk$clusters==8,1:10]
radarchart_8 <- rbind(apply(cluster_8_df, 2, max), apply(cluster_8_df, 2, min), cluster_8_df)
clus_8_colors <- c('firebrick3', 'lightgoldenrod', 'deepskyblue', 'gold', 'tan2', 'lightblue1', 'tan4', 'violetred')
order_lgnd_8 <- c(2,5,7,4,1,6,3,8)
clus_8_lines <- c(3,1,2,3,1,2,1,2)
tiff(file = file.path(FiguresDir, 'v2', 'valley_8_classes_spider_11.15.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_8[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = clus_8_lines, pcol = clus_8_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
legend(x=-1.6, y=-1.2, legend = c('5. Loams w/ res. & mod OM', '1. Sandy soils', '7. Saline-sodic clays', '4. Loams w/ res. & low OM', '2. Loams w/ no res. & mod. OM-low SS', '6. Saline-sodic loams', '3. Loams w/ no res. & mod OM-mod SS', '8. Shrink-swell clays')[order_lgnd_8], col=clus_8_colors[order_lgnd_8], lty=clus_8_lines[order_lgnd_8], lwd = 3, ncol = 3, bty='n')
dev.off()
#Loams w/no res could have mod OM-mod SS added
col2rgb(clus_8_colors)

#revised 5-class option
cluster_5_df <- cluster_fk[cluster_fk$clusters==5,1:10]
radarchart_5 <- rbind(apply(cluster_5_df, 2, max), apply(cluster_5_df, 2, min), cluster_5_df)
clus_5_colors <- c('tan4', 'lightgoldenrod', 'violetred', 'lightblue1', 'firebrick3')
order_lgnd_5 <- c(2,1,5,4,3)
clus_5_lines <- c(1,1,3,3,2)
tiff(file = file.path(FiguresDir, 'v2', 'valley_5_classes_spider_11.15.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_5[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = clus_5_lines, pcol = clus_5_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
legend(x=-1.2, y=-1.2, legend = c('2. Loams w/no res', '1. Sandy & sandy loam soils', '5. Shrink-swell clays', '4. Saline-sodic loams', '3. Loams w/res')[order_lgnd_5], col=clus_5_colors[order_lgnd_5], lty=clus_5_lines[order_lgnd_5], lwd = 3, ncol = 2, bty='n')
dev.off()
#Loams w/no res could have mod OM-mod SS added
col2rgb(clus_5_colors)

#revised 6-class option
cluster_6_df <- cluster_fk[cluster_fk$clusters==6,1:10]
radarchart_6 <- rbind(apply(cluster_6_df, 2, max), apply(cluster_6_df, 2, min), cluster_6_df)
clus_6_colors <- c('lightblue1', 'violetred', 'lightgoldenrod', 'firebrick3', 'gold', 'tan4')
order_lgnd_6 <- c(3,6,5,4,1,2)
clus_6_lines <- c(3,3,1,2,2,1)
tiff(file = file.path(FiguresDir, 'v2', 'valley_6_classes_spider_11.15.19.tif'), family = 'Times New Roman', width = 9, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_6[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'cec_30cm', 'ec_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = clus_6_lines, pcol = clus_6_colors, vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'CEC', 'Salinity', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
legend(x=-1.2, y=-1.2, legend = c('5. Saline-sodic loams', '6. Shrink-swell clays', '1. Sandy (sandy loam) soils', '4. Loams w/res & mod OM', '3. Loams w/res & low OM', '2. Loams w/no res')[order_lgnd_6], col=clus_6_colors[order_lgnd_6], lty=clus_6_lines[order_lgnd_6], lwd = 3, ncol = 2, bty='n')
dev.off()
#Loams w/no res could have mod OM-mod SS added
col2rgb(clus_6_colors)

#calc intra-cluster distances for 7-region model
compute_intra_clus_dist <- function(x, y) {
  dist(rbind(x, y), method = 'euclidean')
}
df_for_clustering_scaled_t <- as.data.frame(t(df_for_clustering_scaled))
intra_clus_dist <- mapply(function(a, b) {compute_intra_clus_dist(x=cluster_9$centers[a,], y=b)}, a=cluster_9$cluster, b=df_for_clustering_scaled_t) #only loops across columns if df_for_clustering_scaled_t is a data.frame
hist(intra_clus_dist)
summary(intra_clus_dist)
names(intra_clus_dist)
valley30cm_by_mukey$intra_clus9_dist <- intra_clus_dist
tapply(valley30cm_by_mukey$intra_clus9_dist, valley30cm_by_mukey$cluster_9, summary)
for (i in 1:9) {
  hist(valley30cm_by_mukey$intra_clus9_dist[valley30cm_by_mukey$cluster_9==i], main=hc_9_labels[i], xlab='intra-cluster euc. distance')
}

compute_intra_clus_dist(x=df_for_clustering_scaled[1,], y=cluster_9$centers[cluster_9$cluster[1],])
compute_intra_clus_dist(x=df_for_clustering_scaled_t[,1], y=cluster_9$centers[cluster_9$cluster[1],])
  
#see daisy in the cluster package with more possibilities in the case of mixed (continuous / categorical) variables.
dist_test_4 <- dist(cluster_4$centers, method = 'euclidean')
hc_4 <- hclust(dist_test_4)
plot(hc_4)

dist_test_5 <- dist(cluster_5$centers, method = 'euclidean')
hc_5 <- hclust(dist_test_5)
hc_5_labels <- c('2. Loams w/no res', '1. Sandy & sandy loam soils', '5. Shrink-swell clays', '4. Saline-sodic loams', '3. Loams w/res')
plot(hc_5, labels=hc_5_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class5.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_5, labels=hc_5_labels, axes=FALSE, ann=FALSE)
dev.off()

dist_test_6 <- dist(cluster_6$centers, method = 'euclidean')
hc_6 <- hclust(dist_test_6)
hc_6_labels <- c('5. Saline-sodic loams', '6. Shrink-swell clays', '1. Sandy (sandy loam) soils', '4. Loams w/res & mod OM', '3. Loams w/res & low OM', '2. Loams w/no res')
plot(hc_6, labels= hc_6_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class6.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_6, labels=hc_6_labels, axes=FALSE, ann=FALSE)
dev.off()

dist_test_7 <- dist(cluster_7$centers, method = 'euclidean')
hc_7 <- hclust(dist_test_7)
hc_7_labels <- c('3. Loams w/res & low OM', '6. Saline-sodic clays', '5. Saline-sodic loams', '1. Sandy loam soils', '7. Shrink-swell clays', '2. Loams w/no res', '4. Loams w/ res. & mod OM')
plot(hc_7, labels=hc_7_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class7.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_7, labels=hc_7_labels, axes=FALSE, ann=FALSE)
dev.off()

dist_test_8 <- dist(cluster_8$centers, method = 'euclidean')
hc_8 <- hclust(dist_test_8)
hc_8_labels <- c('Loams w/ res. & mod OM', 'Sandy soils', 'Saline-sodic clays', 'Loams w/ res. & low OM', 'Loams w/ no res. & mod. OM-low SS', 'Saline-sodic loams', 'Loams w/ no res. & mod OM-mod SS', 'Shrink-swell clays')
plot(hc_8, labels=hc_8_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class8.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_8, labels=hc_8_labels, axes=FALSE, ann=FALSE)
dev.off()

dist_test_9 <- dist(cluster_9$centers, method = 'euclidean')
hc_9 <- hclust(dist_test_9)
hc_9_labels <- c('Sandy soils (6.7%)', 'Shrink-swell clays (14.8%)', 'Loams w/ no res. & mod OM-low SS (21.5%)', 'Loams w/ res. & high OM (0.9%)', 'Loams w/ res. & low OM (15.6%)', 'Loams w/ no res. & mod OM-mod SS (18.4%)', 'Saline-sodic loams (8.8%)', 'Saline-sodic clays (6.9%)', 'Loams w/ res. & mod OM (6.3%)')
plot(hc_9, labels=hc_9_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class9.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_9, labels=hc_9_labels, axes=FALSE, ann=FALSE)
dev.off()

#see ?plot.hclust

dist_test_10 <- dist(cluster_10$centers, method = 'euclidean')
hc_10 <- hclust(dist_test_10)
hc_10_labels <- c('10. Shrink-swell clays', '4. Loams w/ no res. & mod OM-mod SS', '1. Sandy soils', '7. Loams w/ res. & high OM', '2. Loams w/ no res. & low OM', '9. Saline-sodic clays', '5. Loams w/ res. & low OM', '6. Loams w/ res. & mod OM', '3. Loams w/ no res. & mod OM-low SS', '8. Saline-sodic loams')
plot(hc_10, labels=hc_10_labels)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_class10.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.2,0.2,0.2,0.2))
plot(hc_10, labels=hc_10_labels, axes=FALSE, ann=FALSE)
dev.off()

dist_test_11 <- dist(cluster_11$centers, method = 'euclidean')
hc_11 <- hclust(dist_test_11)
hc_11_labels <- c('Shrink-swell clays', 'Loams w/ res. & low OM', 'Loams w/ res. & mod OM-low SS', 'Loams w/ no res. & low OM', 'Loams w/ res. & mod OM-mod SS', 'Loams w/ res. & high OM', 'Saline-sodic clays',  'Saline-sodic loams', 'Sands', 'Loams w/ no res. & mod OM-mod SS', 'Loams w/ no res. & mod OM-low SS')
plot(hc_11, labels=hc_11_labels)

#figure out correlation across cluster sizes 2-7
cluster_fk_2_7 <- rbind(cluster_2$centers, cluster_3$centers, cluster_4$centers, cluster_5$centers, cluster_6$centers, cluster_7$centers)
rownames(cluster_fk_2_7) <- paste0(c(rep('2_', 2), rep('3_', 3), rep('4_', 4), rep('5_', 5), rep('6_', 6), rep('7_', 7)), rownames(cluster_fk_2_7))
dist_test_2_7 <- dist(cluster_fk_2_7, method = 'euclidean')
hc_2_7 <- hclust(dist_test_2_7)
plot(hc_2_7, axes=FALSE, ann=FALSE, hang=0.1)
tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_classes2_to7.tif'), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0,0,0,0))
plot(hc_2_7, axes=FALSE, ann=FALSE)
dev.off()

#figure out correlation across cluster sizes 5-10
make_hclus_labels <- function(x, y) {paste0(x, y)}
cluster_fk_5_11 <- rbind(cluster_5$centers, cluster_6$centers, cluster_7$centers, cluster_8$centers, cluster_9$centers, cluster_10$centers, cluster_11$centers)
dist_test_5_11 <- dist(cluster_fk_5_11, method = 'euclidean')
hc_5_11 <- hclust(dist_test_5_11)
plot(hc_5_11, labels=c(make_hclus_labels(hc_5_labels, '-5'), make_hclus_labels(hc_6_labels, '-6'), make_hclus_labels(hc_7_labels, '-7'), make_hclus_labels(hc_8_labels, '-8'), make_hclus_labels(hc_9_labels, '-9'), make_hclus_labels(hc_10_labels, '-10'), make_hclus_labels(hc_11_labels, '-11')), axes=FALSE, ann=FALSE, hang=0.1)

tiff(file = file.path(FiguresDir, 'v2', 'cluster_dendrogram_classes5_to11.tif'), family = 'Times New Roman', width = 9, height = 7, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0,0,0,0))
plot(hc_5_11, labels=c(make_hclus_labels(hc_5_labels, '-5'), make_hclus_labels(hc_6_labels, '-6'), make_hclus_labels(hc_7_labels, '-7'), make_hclus_labels(hc_8_labels, '-8'), make_hclus_labels(hc_9_labels, '-9'), make_hclus_labels(hc_10_labels, '-10'), make_hclus_labels(hc_11_labels, '-11')), axes=FALSE, ann=FALSE)
dev.off()

#look at variables by cluster
#OM
tapply(valley30cm_by_mukey$om_30cm, valley30cm_by_mukey$cluster_7, summary)
class_9_om_summary <- tapply(valley30cm_by_mukey$om_30cm, valley30cm_by_mukey$cluster_9, summary)
names(class_9_om_summary) <- hc_9_labels
sum(grepl('Vertisols', valley30cm_by_mukey$txorders[valley30cm_by_mukey$cluster_7==5]), na.rm = TRUE) / sum(!is.na(valley30cm_by_mukey$txorders[valley30cm_by_mukey$cluster_7==5])) #68% labeled as Verts have Verts as a taxinomic order
sum(grepl('Vertisols', valley30cm_by_mukey$txorders[valley30cm_by_mukey$cluster_9==2]), na.rm = TRUE) / sum(!is.na(valley30cm_by_mukey$txorders[valley30cm_by_mukey$cluster_9==2])) #71.4%

#calculate intracluster distances


#another way to get optimal number of clusters
best_cluster <- kmeansruns(data = df_for_clustering_scaled, krange = 1:20, criterion = 'asw', iter.max = 200, runs = 100, nstart=50)
best_cluster
best_cluster$bestk #2!!!
best_cluster_ch <- kmeansruns(data = df_for_clustering_scaled, krange = 1:20, criterion = 'ch', iter.max = 200, runs = 100, nstart=50)
best_cluster_ch
best_cluster_ch$bestk

km.boot_2 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=2, runs=25) #runs is same as nstart
km.boot_2
km.boot_3 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=3, runs=25) #runs is same as nstart
km.boot_3
km.boot_4 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=4, runs=25) #runs is same as nstart
km.boot_4


#find optimum number of clusters based on gap statistic
gap_stats <- clusGap(df_for_clustering_scaled, FUN = kmeans, nstart = 50, K.max = 20, B = 500, iter.max = 200)
fviz_gap_stat(gap_stats)

#visualize stats by cluster using untransformed data

#compare key properties across classes
colnames(df_for_clustering_scaled)
axis_labels <- c("Depth to Restriction", 'Clay', 'Organic matter', 'CEC', 'Bulk density', 'Salinity', 'pH', 'Shrink-swell', 'Saturated conductivity', 'Available water capacity')
cor_matrix_pearson <- cor(df_for_clustering_scaled, method = 'pearson')
rownames(cor_matrix_pearson) <- axis_labels
colnames(cor_matrix_pearson) <- axis_labels
write.csv(cor_matrix_pearson, file.path(dataDir, 'v2 results', 'soil property correlations', 'pearson_scaled.csv'), row.names = TRUE)
cor_pearson_pvals <- cor.mtest(df_for_clustering_scaled, method='pearson')
rownames(cor_pearson_pvals$p) <- axis_labels
colnames(cor_pearson_pvals$p) <- axis_labels
write.csv(cor_pearson_pvals$p, file.path(dataDir, 'v2 results', 'soil property correlations', 'pearson_pval_scaled.csv'), row.names = TRUE)

cor_matrix_spearman <- cor(df_for_clustering_scaled, method = 'spearman')
rownames(cor_matrix_spearman) <- axis_labels
colnames(cor_matrix_spearman) <- axis_labels
write.csv(cor_matrix_spearman, file.path(dataDir, 'v2 results', 'soil property correlations', 'spearman_scaled.csv'), row.names = TRUE)
cor_spearman_pvals <- cor.mtest(cor_matrix_spearman)
rownames(cor_spearman_pvals$p) <- axis_labels
colnames(cor_spearman_pvals$p) <- axis_labels
write.csv(cor_spearman_pvals$p, file.path(dataDir, 'v2 results', 'soil property correlations', 'spearman_pval_scaled.csv'), row.names = TRUE)
cor.test(df_for_clustering_scaled$om_30cm, df_for_clustering_scaled$cec_30cm, method = 'pearson')


cex.before <- par("cex") #saves current cex setting for plotting
cex.reduced <- 0.5 #desired cex setting for plotting p-values
mag.factor <- 2.5 #fudge factor to increase size of axis and legend text

#produced
axis_labels <- c("Depth to\nRestriction", 'Clay', 'Organic\nmatter', 'CEC', 'Bulk \ndensity', 'Salinity', 'pH', 'Shrink-\n  swell', 'Saturated\nconductivity', 'Available water\n   capacity')
tiff(file = file.path(FiguresDir, 'v2', 'soil_property_corrplot_pvals.tif'), family = 'Times New Roman', width = 6.5, height = 6.5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0.1,0.1,0.1,0.1), cex = cex.reduced)
corrplot(cor_matrix, p.mat = cor_matrix_pvals$p, insig = "p-value", sig.level = -1, tl.cex = par("cex") * mag.factor, cl.cex = par("cex") * mag.factor, tl.col = 'black')
par(cex = cex.before)
dev.off()

#violin plots
#om
boxplot(om_30cm ~ cluster_5, data = valley30cm_by_mukey, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_6, data = valley30cm_by_mukey, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_7, data = valley30cm_by_mukey, ylim=c(0, 5))

#vioplot for cluster 5
c(1,2,3,4,5)[c(2,1,5,4,3)]
hc_5_labels[c(2,1,5,4,3)]

#5 class vioplot
clus_5_colors <- c('tan4', 'lightgoldenrod', 'violetred', 'lightblue1', 'firebrick3')
order_lgnd_5 <- c(2, 1, 5, 4, 3)
clus_5_colors[order_lgnd_5]
vioplot_mod_clus5 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:5)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_5==plot_order2[1]], times=round(df$area_ac[df$cluster_5==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[2]], times=round(df$area_ac[df$cluster_5==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[3]], times=round(df$area_ac[df$cluster_5==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[4]], times=round(df$area_ac[df$cluster_5==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[5]], times=round(df$area_ac[df$cluster_5==plot_order2[5]]/area_fact, 0)), col=clus_5_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}

vioplot_mod_clus5(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_5, area_fact = 10, ylab='Clay (%)', fname='class5_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_5, area_fact = 10, ylab='Organic matter (%)', fname='class5_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus5(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class5_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class5_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_5, area_fact = 10, ylab='soil pH', fname='class5_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class5_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class5_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class5_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_5, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class5_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus5(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
# valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus5(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_5, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class5_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus5(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_5, area_fact = 10, ylab='Linear extensibility (%)', fname='class5_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))

#6 class vioplot
clus_6_colors <- c('lightblue1', 'violetred', 'lightgoldenrod', 'firebrick3', 'gold', 'tan4')
order_lgnd_6 <- c(3,6,5,4,1,2)
clus_6_colors[order_lgnd_6]
vioplot_mod_clus6 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:6)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_6==plot_order2[1]], times=round(df$area_ac[df$cluster_6==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_6==plot_order2[2]], times=round(df$area_ac[df$cluster_6==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_6==plot_order2[3]], times=round(df$area_ac[df$cluster_6==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_6==plot_order2[4]], times=round(df$area_ac[df$cluster_6==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_6==plot_order2[5]], times=round(df$area_ac[df$cluster_6==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_6==plot_order2[6]], times=round(df$area_ac[df$cluster_6==plot_order2[6]]/area_fact, 0)), col=clus_6_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}

vioplot_mod_clus6(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_6, area_fact = 10, ylab='Clay (%)', fname='class6_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_6, area_fact = 10, ylab='Organic matter (%)', fname='class6_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus6(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class6_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class6_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_6, area_fact = 10, ylab='soil pH', fname='class6_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class6_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class6_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class6_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_6, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class6_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus6(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus6(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_6, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class6_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus6(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_6, area_fact = 10, ylab='Linear extensibility (%)', fname='class6_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))


#7 class vioplot
clus_7_colors <- c('gold', 'deepskyblue', 'lightblue1', 'lightgoldenrod', 'violetred', 'tan4', 'firebrick3')
order_lgnd_7 <- c(4,6,1,7,3,2,5)
vioplot_mod_clus7 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:7)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_7==plot_order2[1]], times=round(df$area_ac[df$cluster_7==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[2]], times=round(df$area_ac[df$cluster_7==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[3]], times=round(df$area_ac[df$cluster_7==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[4]], times=round(df$area_ac[df$cluster_7==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[5]], times=round(df$area_ac[df$cluster_7==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[6]], times=round(df$area_ac[df$cluster_7==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_7==plot_order2[7]], times=round(df$area_ac[df$cluster_7==plot_order2[7]]/area_fact, 0)), col=clus_7_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}

vioplot_mod_clus7(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_7, area_fact = 10, ylab='Clay (%)', fname='class7_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_7, area_fact = 10, ylab='Organic matter (%)', fname='class7_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus7(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class7_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class7_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_7, area_fact = 10, ylab='soil pH', fname='class7_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class7_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class7_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class7_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_7, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class7_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus7(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus7(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class7_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus7(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_7, area_fact = 10, ylab='Linear extensibility (%)', fname='class7_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))

#8 class vioplot
clus_8_colors <- c('firebrick3', 'lightgoldenrod', 'deepskyblue', 'gold', 'tan2', 'lightblue1', 'tan4', 'violetred')
order_lgnd_8 <- c(2,5,7,4,1,6,3,8)
vioplot_mod_clus8 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:8)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_8==plot_order2[1]], times=round(df$area_ac[df$cluster_8==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[2]], times=round(df$area_ac[df$cluster_8==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[3]], times=round(df$area_ac[df$cluster_8==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[4]], times=round(df$area_ac[df$cluster_8==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[5]], times=round(df$area_ac[df$cluster_8==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[6]], times=round(df$area_ac[df$cluster_8==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[7]], times=round(df$area_ac[df$cluster_8==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_8==plot_order2[8]], times=round(df$area_ac[df$cluster_8==plot_order2[8]]/area_fact, 0)), col=clus_8_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
vioplot_mod_clus8(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_8, area_fact = 10, ylab='Clay (%)', fname='class8_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_8, area_fact = 10, ylab='Organic matter (%)', fname='class8_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus8(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class8_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class8_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_8, area_fact = 10, ylab='soil pH', fname='class8_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class8_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class8_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class8_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_8, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class8_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus8(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus8(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_8, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class8_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus8(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_8, area_fact = 10, ylab='Linear extensibility (%)', fname='class8_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))

#9 class vioplot
vioplot_mod_clus9 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:9)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_9==plot_order2[1]], times=round(df$area_ac[df$cluster_9==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[2]], times=round(df$area_ac[df$cluster_9==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[3]], times=round(df$area_ac[df$cluster_9==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[4]], times=round(df$area_ac[df$cluster_9==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[5]], times=round(df$area_ac[df$cluster_9==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[6]], times=round(df$area_ac[df$cluster_9==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[7]], times=round(df$area_ac[df$cluster_9==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[8]], times=round(df$area_ac[df$cluster_9==plot_order2[8]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[9]], times=round(df$area_ac[df$cluster_9==plot_order2[9]]/area_fact, 0)), col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
#order_lgnd was defined for radarchart
vioplot_mod_clus9(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd, area_fact = 10, ylab='Clay (%)', fname='class9_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd, area_fact = 10, ylab='Organic matter (%)', fname='class9_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus9(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class9_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class9_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd, area_fact = 10, ylab='soil pH', fname='class9_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class9_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class9_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class9_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class9_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus9(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus9(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class9_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus9(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd, area_fact = 10, ylab='Linear extensibility (%)', fname='class9_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))


test <- vioplot(rep(valley30cm_by_mukey$om_30cm[valley30cm_by_mukey$cluster_9==1], times=round(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_9==1]/10, 0)))
vioplot(rep(valley30cm_by_mukey$om_30cm[valley30cm_by_mukey$cluster_9==4], times=round(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_9==4]/10, 0)))
quantile(rep(valley30cm_by_mukey$om_30cm[valley30cm_by_mukey$cluster_9==2], times=round(valley30cm_by_mukey$area_ac[valley30cm_by_mukey$cluster_9==2]/10, 0)), probs=0.1)


#10 class vioplot
hc_10_labels
vioplot_mod_clus10 <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar) {
  plot_order2 <- (1:10)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(rep(df[[varname]][df$cluster_10==plot_order2[1]], times=round(df$area_ac[df$cluster_10==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[2]], times=round(df$area_ac[df$cluster_10==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[3]], times=round(df$area_ac[df$cluster_10==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[4]], times=round(df$area_ac[df$cluster_10==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[5]], times=round(df$area_ac[df$cluster_10==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[6]], times=round(df$area_ac[df$cluster_10==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[7]], times=round(df$area_ac[df$cluster_10==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[8]], times=round(df$area_ac[df$cluster_10==plot_order2[8]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[9]], times=round(df$area_ac[df$cluster_10==plot_order2[9]]/area_fact, 0)), rep(df[[varname]][df$cluster_10==plot_order2[10]], times=round(df$area_ac[df$cluster_10==plot_order2[10]]/area_fact, 0)), col=clus_10_colors[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}

vioplot_mod_clus10(valley30cm_by_mukey, 'clay_30cm', ylim_vioplot = c(0.5,70), plot_order = order_lgnd_10, area_fact = 10, ylab='Clay (%)', fname='class10_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'om_30cm', ylim_vioplot = c(0.1,12), plot_order = order_lgnd_10, area_fact = 10, ylab='Organic matter (%)', fname='class10_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
valley30cm_by_mukey$logom_30cm <- log(valley30cm_by_mukey$om_30cm)
vioplot_mod_clus10(valley30cm_by_mukey, 'logom_30cm', ylim_vioplot = c(-3.15,3.61),  plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Organic matter (Log'[10]~'%)'), fname='class10_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'cec_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class10_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'pH_30cm', ylim_vioplot = c(5.2,10.1), plot_order = order_lgnd_10, area_fact = 10, ylab='soil pH', fname='class10_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class10_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'ec_30cm', ylim_vioplot = c(0,50), plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class10_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'bd_30cm', ylim_vioplot = c(0.75,1.8), plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class10_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'MnRs_dep', ylim_vioplot = c(0,200), plot_order = order_lgnd_10, area_fact = 10, ylab='Minimum depth to restrictive layer', fname='class10_MnRs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
# vioplot_mod_clus10(valley30cm_by_mukey, 'ksat_30cm', ylim_vioplot = c(0,350), plot_order = 1:9, area_fact = 10, labnames = NULL)
# valley30cm_by_mukey$logks_30cm <- log(valley30cm_by_mukey$ksat_30cm)
vioplot_mod_clus10(valley30cm_by_mukey, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_10, area_fact = 10, ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class10_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
vioplot_mod_clus10(valley30cm_by_mukey, 'lep_30cm', ylim_vioplot = c(0,17), plot_order = order_lgnd_10, area_fact = 10, ylab='Linear extensibility (%)', fname='class10_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1))



#boxplot 5 class
boxplot_mod_clus5 <- function(df, varname, ylim_boxplot, plot_order, area_fact, labnames) {
  plot_order2 <- c(1,2,3,4,5)[plot_order]
  par(mar=c(6, 2, 2, 2))
  boxplot(rep(df[[varname]][df$cluster_5==plot_order2[1]], times=round(df$area_ac[df$cluster_5==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[2]], times=round(df$area_ac[df$cluster_5==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[3]], times=round(df$area_ac[df$cluster_5==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[4]], times=round(df$area_ac[df$cluster_5==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[5]], times=round(df$area_ac[df$cluster_5==plot_order2[5]]/area_fact, 0)), ylim = ylim_boxplot, ann=FALSE)
  mtext(text=labnames[plot_order], side = 1, line=3, at=1:5)
}
boxplot_mod_clus5(valley30cm_by_mukey, 'om_30cm', ylim_boxplot = c(0,5), plot_order = c(2,1,5,4,3), area_fact = 10, labnames = c('Loams w/\nno res.\nmod OM-mod SS', 'Sandy\nsoils', 'Shrink-swell\nclays', 'Saline-sodic\nloams', 'Loams w/\nres &\nmod OM'))

boxplot_mod_clus6 <- function(df, varname, ylim_boxplot, plot_order, area_fact, labnames) {
  plot_order2 <- c(1,2,3,4,5)[plot_order]
  par(mar=c(6, 2, 2, 2))
  boxplot(rep(df[[varname]][df$cluster_5==plot_order2[1]], times=round(df$area_ac[df$cluster_5==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[2]], times=round(df$area_ac[df$cluster_5==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[3]], times=round(df$area_ac[df$cluster_5==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[4]], times=round(df$area_ac[df$cluster_5==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_5==plot_order2[5]], times=round(df$area_ac[df$cluster_5==plot_order2[5]]/area_fact, 0)), ylim = ylim_boxplot, ann=FALSE)
  mtext(text=labnames[plot_order], side = 1, line=3, at=1:5)
}
boxplot_mod_clus5(valley30cm_by_mukey, 'om_30cm', ylim_boxplot = c(0,5), plot_order = c(2,1,5,4,3), area_fact = 10, labnames = c('Loams w/\nno res.\nmod OM-mod SS', 'Sandy\nsoils', 'Shrink-swell\nclays', 'Saline-sodic\nloams', 'Loams w/\nres &\nmod OM'))

#boxplot 9 class
boxplot_mod_clus9 <- function(df, varname, ylim_boxplot, plot_order, area_fact, labnames) {
  plot_order2 <- c(1,2,3,4,5,6,7,8,9)[plot_order]
  boxplot(rep(df[[varname]][df$cluster_9==plot_order2[1]], times=round(df$area_ac[df$cluster_9==plot_order2[1]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[2]], times=round(df$area_ac[df$cluster_9==plot_order2[2]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[3]], times=round(df$area_ac[df$cluster_9==plot_order2[3]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[4]], times=round(df$area_ac[df$cluster_9==plot_order2[4]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[5]], times=round(df$area_ac[df$cluster_9==plot_order2[5]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[6]], times=round(df$area_ac[df$cluster_9==plot_order2[6]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[7]], times=round(df$area_ac[df$cluster_9==plot_order2[7]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[8]], times=round(df$area_ac[df$cluster_9==plot_order2[8]]/area_fact, 0)), rep(df[[varname]][df$cluster_9==plot_order2[9]], times=round(df$area_ac[df$cluster_9==plot_order2[9]]/area_fact, 0)), ylim = ylim_boxplot, ann=FALSE)
}
boxplot_mod_clus9(valley30cm_by_mukey, 'om_30cm', ylim_boxplot = c(0,8), plot_order = 1:9, area_fact = 10, labnames = 1:9)
