#this modified to work with aggregated dataset from ssurgo_calag_aggregate.R
#TO-DO
#(1) run cluster analysis on Salinas only [DONE]
#(2) summarize classes for 4 and 5 
#(3) log transform om and ksat [DONE]
#(4) identify outliers
laptop <- FALSE
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
# 
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  
}
mar_settings <- c(4, 4.5, 1, 1)
list.files(dataDir)
# valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm.shp'))
# names(valley_mu_shp_30cm)
valley30cm <- read.csv(file.path(dataDir, 'valley_30cm_data.csv'), stringsAsFactors = FALSE)
# dim(valley30cm)
# colnames(valley30cm)
# unique(valley30cm$mjcmpnms)
# sum(grepl('Rock outcrop', valley30cm$mjcmpnms)) #929
# unique(valley30cm$mjcmpnms[grepl('Rock outcrop', valley30cm$mjcmpnms)]) #72 unique majcomp combos
# unique(valley30cm$mukey[grepl('Rock outcrop', valley30cm$mjcmpnms)]) #148 unique mukeys
# sum(valley30cm$area_ac[grepl('Rock outcrop', valley30cm$mjcmpnms)]) #53,523 have rock OC as majcomp
acres_by_mukey <- aggregate(area_ac ~ mukey, data = valley30cm, sum)
# table(valley30cm$mukey)
# length(unique(valley30cm$mukey)) #5043
valley30cm_by_mukey <- valley30cm[!duplicated(valley30cm$mukey), ]
# dim(valley30cm_by_mukey)
# length(unique(valley30cm_by_mukey$mukey))
valley30cm_by_mukey$area_ac <- acres_by_mukey$area_ac[match(valley30cm_by_mukey$mukey, acres_by_mukey$mukey)]
# sum(valley30cm_by_mukey$area_ac) #13873110 matches
valley30cm_by_mukey$area_proportion <- valley30cm_by_mukey$area_ac / sum(valley30cm_by_mukey$area_ac)
# hist(valley30cm_by_mukey$area_proportion)
# valley30cm_by_mukey[which.max(valley30cm_by_mukey$area_proportion),] #Tulare clay wins the max extent prize!
# colnames(valley30cm_by_mukey)
# unique(valley30cm_by_mukey$muname[is.na(valley30cm_by_mukey$aws050wta)])
# sum(valley30cm_by_mukey$area_ac[is.na(valley30cm_by_mukey$aws050wta)]) #270,261
# valley30cm_by_mukey_v2 <- valley30cm_by_mukey[!is.na(valley30cm_by_mukey$aws050wta), ]
# sum(valley30cm_by_mukey_v2$area_ac)
# sum(valley30cm_by_mukey_v2$area_ac[grepl('Rock outcrop', valley30cm_by_mukey_v2$mjcmpnms)])
#apply QC using criteria that 85% of mapunit needs to have OM data, as is SSURGO standard
# sum(valley30cm_by_mukey$area_ac[which(valley30cm_by_mukey$compct_om >= 85)]) / sum(valley30cm_by_mukey$area_ac) #which is to handle NAs
valley30cm_by_mukey <- valley30cm_by_mukey[which(valley30cm_by_mukey$compct_om >= 85), ]
# sum(is.na(valley30cm_by_mukey$aws050wta))
# sum(is.na(valley30cm_by_mukey$om_30cm))
# summary(valley30cm_by_mukey$compct_om)
# summary(valley30cm_by_mukey$mjcmp_pct)
analysis_preview <- na.omit(valley30cm_by_mukey[ ,c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'sar_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm', 'area_ac')])
# sum(analysis_preview$area_ac) #11,290,964 acres will be assigned a soil health diagnostic indicators cluster class


#create data.frame for cluster analysis
df_for_clustering <- valley30cm_by_mukey
rownames(df_for_clustering) <- df_for_clustering$mukey
# summary(df_for_clustering)
df_for_clustering <- df_for_clustering[,c('MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'ec_30cm', 'pH_30cm', 'sar_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm')] #frags have 160 mukeys with NA
# "storiemn", "Lthc_dep", "Plth_dep", "Drpn_dep", "ATC_dep", "Natr_dep", "Salc_dep", "SCTS_dep", "MRes_dep" 'kwf_30cm'
sum()
mapply(function(x, y) hist(x, main=y), x=df_for_clustering, y=colnames(df_for_clustering))
lapply(df_for_clustering, class)
lapply(df_for_clustering, function(x) sum(x==0, na.rm = TRUE)) #verify why some are showing min res depth of 0
lapply(df_for_clustering, function(x) sum(x<0, na.rm = TRUE))
lapply(df_for_clustering, function(x) sum(x>0 & x < 0.01, na.rm = TRUE))
lapply(df_for_clustering, function(x) sum(x>0 & x < 0.1, na.rm = TRUE))
df_for_clustering[df_for_clustering$MnRs_dep==100,]
df_for_clustering[which(df_for_clustering$ec_30cm > 0 & df_for_clustering$ec_30cm < 0.01), ]
df_for_clustering[which(df_for_clustering$bd_30cm < 0.5), ] #doesn't matter since these excluded from na.omit
log_transform <- function(x, df) {
  df[[x]] <- log(df[[x]] + 0.01, 10)
  df
}
var_list <- c('om_30cm', 'cec_30cm', 'ec_30cm', 'sar_30cm', 'lep_30cm', 'ksat_30cm')
for (i in seq_along(var_list)) {
  df_for_clustering <- log_transform(x=var_list[i], df = df_for_clustering) #above variables are log transformed
}
mapply(function(x, y) hist(x, main=y), x=df_for_clustering, y=colnames(df_for_clustering))
dim(df_for_clustering)
df_for_clustering_scaled <- as.data.frame(scale(df_for_clustering))
summary(df_for_clustering_scaled)
colnames(df_for_clustering_scaled)
mapply(function(x, y) hist(x, main=y), x=df_for_clustering_scaled, y=colnames(df_for_clustering_scaled))


kmeans_test <- function(y, z) {sapply(1:y, function(x) {
  result <- kmeans(na.omit(z), centers = x, iter.max = 100, nstart = 25)
  round(100 * result$betweenss / result$totss, 1)})
}
results <- replicate(100, kmeans_test(20, df_for_clustering_scaled))
dim(results)
rowMeans(results)
apply(results, 1, sd)
tiff(file = file.path(FiguresDir, 'kmeans_comparison.tif'), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=mar_settings)
plot(1:20, rowMeans(results), type='b', xlab='Number of clusters', ylab='Soil variability captured by clustering (%)', cex=0.8, cex.axis=1, cex.lab=1)
#text(1:12, rowMeans(results), labels=as.character(1:12), pos=1, offset=0.5)
dev.off()

cluster_2 <- kmeans(na.omit(df_for_clustering_scaled), centers=2, iter.max = 200, nstart = 50)
cluster_2
cluster_3 <- kmeans(na.omit(df_for_clustering_scaled), centers=3, iter.max = 200, nstart = 50)
cluster_3
cluster_4 <- kmeans(na.omit(df_for_clustering_scaled), centers=4, iter.max = 200, nstart = 50)
cluster_4
cluster_5 <- kmeans(na.omit(df_for_clustering_scaled), centers=5, iter.max = 200, nstart = 50)
cluster_5
cluster_6 <- kmeans(na.omit(df_for_clustering_scaled), centers=6, iter.max = 200, nstart = 50)
cluster_6
cluster_7 <- kmeans(na.omit(df_for_clustering_scaled), centers=7, iter.max = 200, nstart = 50)
cluster_7
cluster_8 <- kmeans(na.omit(df_for_clustering_scaled), centers=8, iter.max = 200, nstart = 50)
cluster_8
cluster_9 <- kmeans(na.omit(df_for_clustering_scaled), centers=9, iter.max = 200, nstart = 50)
cluster_9
cluster_10 <- kmeans(na.omit(df_for_clustering_scaled), centers=10, iter.max = 200, nstart = 50)
cluster_10
cluster_11 <- kmeans(na.omit(df_for_clustering_scaled), centers=11, iter.max = 200, nstart = 50)
cluster_11
cluster_12 <- kmeans(na.omit(df_for_clustering_scaled), centers=12, iter.max = 200, nstart = 50)
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
colnames(cluster_fk)
dim(cluster_fk)
cluster_fk$avg_dist_within <- cluster_fk$withinss / cluster_fk$n
cluster_fk
# cbind(cluster_fk[cluster_fk$clusters==5, 13:16], cluster_fk_labels[(nrow(cluster_fk_labels)-4):nrow(cluster_fk_labels), ])

#make a radarchart plotting function
cluster_radar_plot <- function(cluster_no) {
  radarchart_df <- cluster_fk[cluster_fk$clusters==cluster_no,]
  radarchart_df <- rbind(apply(radarchart_df, 2, max), apply(radarchart_df, 2, min), radarchart_df) #have to create a data.frame where first row is max, second row is min, and then followed by the data for radarchart
  tiff(file = file.path(FiguresDir, paste0('valley_', cluster_no, '_classes_spider.tif')), family = 'Times New Roman', width = 6.5, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
  par(mar=rep(0.1, 4))
  radarchart(radarchart_df[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'ec_30cm', 'sar_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], vlabels=c('Clay', 'Organic\n matter', ' Shrink-\n  swell', 'Saline', 'Sodic',  "Depth to Restriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
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
# write.csv(valley30cm_by_mukey, file.path(dataDir, 'valley30cm_by_mukey_cluster.csv'), row.names = FALSE)

cluster_area_summary <- c(tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_2, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_3, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_4, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_5, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_6, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_7, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_8, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_9, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_10, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_11, sum), tapply(valley30cm_by_mukey$area_ac, valley30cm_by_mukey$cluster_12, sum))
cluster_fk$area_ac <- cluster_area_summary
cluster_fk$area_pct <- 100 * cluster_fk$area_ac / sum(analysis_preview$area_ac)
# write.csv(cluster_fk, file.path(FiguresDir, 'clusters2_to_12_df_for_radarchart.csv'), row.names = FALSE)
cluster_fk <- read.csv(file.path(FiguresDir, 'clusters2_to_12_df_for_radarchart.csv'))

cluster_10 <- cluster_fk[cluster_fk$clusters==10,1:11]
# cluster_10$cluster_name <- c('Loams w/ res.', 'Saline-sodic w/ high 2:1', 'Sands', 'Loams w/ no res.', 'Mod 2:1 w/ no res.', 'Saline-sodic w/ mod 2:1', 'High 2:1', 'Mod. 2:1 w/ res.', 'Saline loams', 'Saline-sodic w/ low 2:1')
radarchart_10 <- rbind(apply(cluster_10, 2, max), apply(cluster_10, 2, min), cluster_10)

tiff(file = file.path(FiguresDir, 'valley_10_classes_spider_FINAL.tif'), family = 'Times New Roman', width = 6.5, height = 6.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(xpd=TRUE, mar=c(4, 0.1, 0.1, 0.1))
radarchart(radarchart_10[,c('clay_30cm', 'om_30cm', 'lep_30cm', 'ec_30cm', 'sar_30cm', 'pH_30cm',  "MnRs_dep", 'ksat_30cm', 'awc_30cm')], plty = c(3, 1, 1, 1, 1, 2, 2, 3, 1, 1), pcol = c('tan3', 'lightblue4', 'lightgoldenrod', 'tan3', 'tan4', 'lightblue2', 'black', 'tan4', 'gold', 'lightblue1'), vlabels=c('Clay', 'Organic\nmatter', ' Shrink-\n  swell', 'Saline', 'Sodic', 'pH',  "Depth to\nRestriction", expression('K'['s']), 'Available water\n   capacity'), maxmin = TRUE, plwd = 3)
order_lgnd <- c(3,4,1,9,2,6,10,7,5,8)
legend(x=-1.4, y=-1.2, legend =c('Loams w/ res.', 'Saline-sodic w/ high 2:1', 'Sands', 'Loams w/ no res.', 'Mod 2:1 w/ no res.', 'Saline-sodic w/ mod 2:1', 'High 2:1', 'Mod. 2:1 w/ res.', 'Saline loams', 'Saline-sodic w/ low 2:1')[order_lgnd], col=c('tan3', 'lightblue4', 'lightgoldenrod', 'tan3', 'tan4', 'lightblue2', 'black', 'tan4', 'gold', 'lightblue1')[order_lgnd], lty=c(3, 1, 1, 1, 1, 2, 2, 3, 1, 1)[order_lgnd], lwd = 3, ncol = 3)
dev.off()

#rgb color codes
col2rgb(c('tan3', 'lightblue4', 'lightgoldenrod', 'tan3', 'tan4', 'lightblue2', 'black', 'tan4', 'gold', 'lightblue1')[order_lgnd])

#see daisy in the cluster package with more possibilities in the case of mixed (continuous / categorical) variables.
dist_test_4 <- dist(cluster_4$centers, method = 'euclidean')
hc_4 <- hclust(dist_test_4)
plot(hc_4)

dist_test_5 <- dist(cluster_5$centers, method = 'euclidean')
hc_5 <- hclust(dist_test_5)
plot(hc_5)

dist_test_6 <- dist(cluster_6$centers, method = 'euclidean')
hc_6 <- hclust(dist_test_6)
plot(hc_6)

dist_test_7 <- dist(cluster_7$centers, method = 'euclidean')
hc_7 <- hclust(dist_test_7)
plot(hc_7)

#test to figure out correlation across cluster sizes
cluster_fk_2_7 <- rbind(cluster_2$centers, cluster_3$centers, cluster_4$centers, cluster_5$centers, cluster_6$centers, cluster_7$centers)
rownames(cluster_fk_2_7) <- paste0(c(rep('2_', 2), rep('3_', 3), rep('4_', 4), rep('5_', 5), rep('6_', 6), rep('7_', 7)), rownames(cluster_fk_2_7))
dist_test_2_7 <- dist(cluster_fk_2_7, method = 'euclidean')
hc_2_7 <- hclust(dist_test_2_7)
plot(hc_2_7, axes=FALSE, ann=FALSE, hang=0.1)
tiff(file = file.path(FiguresDir, 'cluster_dendrogram_classes2_to7.tif'), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
par(mar=c(0,0,0,0))
plot(hc_2_7, axes=FALSE, ann=FALSE)
dev.off()

hcd <- as.dendrogram(hc_5)
par(mfrow=c(2,1))
plot(cut(hcd, h=100)$upper, 
     main="Upper tree of cut at h=100")
plot(cut(hcd, h=100)$lower[[2]], 
     main="Second branch of lower tree with cut at h=100")

rect.hclust(hc_5, k = 4, border = "red")
#test plotting
hc <- hclust(dist(USArrests), "ave")
plot(hc)
plot(hc, hang = -1)


#another way to get optimal number of clusters
best_cluster <- kmeansruns(data = na.omit(df_for_clustering_scaled), krange = 1:12, criterion = 'asw', iter.max = 100, runs = 100, nstart=25)
best_cluster
best_cluster$bestk

km.boot_2 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=2, runs=25) #runs is same as nstart
km.boot_2
km.boot_3 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=3, runs=25) #runs is same as nstart
km.boot_3
km.boot_4 <- clusterboot(na.omit(df_for_clustering_scaled), B=500, bootmethod="boot", clustermethod=kmeansCBI, krange=4, runs=25) #runs is same as nstart
km.boot_4

#kmeans check
kmeans_test_2 <- function(y, z) {sapply(1:y, function(x) {
  result <- kmeans(na.omit(z), centers = x, iter.max = 20, nstart = 25)
  round(sum(result$withinss), 0)})
}
results_2 <- replicate(100, kmeans_test_2(12, df_for_clustering_scaled))
plot(1:12, rowMeans(results_2), type='b')
text(1:12, rowMeans(results_2), labels=as.character(1:12), pos=1, offset=0.5)

#find optimum number of clusters based on gap statistic
gap_stats <- clusGap(na.omit(df_for_clustering_scaled), FUN = kmeans, nstart = 25, K.max = 12, B = 500, iter.max = 20)
fviz_gap_stat(gap_stats)

valley_mu_shp_30cm$cluster_2 <- cluster_2$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_2$cluster))]
valley_mu_shp_30cm$cluster_3 <- cluster_3$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_3$cluster))]
valley_mu_shp_30cm$cluster_4 <- cluster_4$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_4$cluster))]
valley_mu_shp_30cm$cluster_5 <- cluster_5$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_5$cluster))]
valley_mu_shp_30cm$cluster_6 <- cluster_6$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_6$cluster))]
valley_mu_shp_30cm$cluster_7 <- cluster_7$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_7$cluster))]
# valley_mu_shp_30cm$cluster_8 <- cluster_8$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_8$cluster))]
# valley_mu_shp_30cm$cluster_10 <- cluster_10$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_10$cluster))]
# valley_mu_shp_30cm$cluster_13 <- cluster_13$cluster[match(valley_mu_shp_30cm$mukey, names(cluster_13$cluster))]
shapefile(valley_mu_shp_30cm, file.path(dataDir, 'valley_30cm_cluster.shp'), overwrite=TRUE)
sum(valley_mu_shp_30cm$area_ac[!is.na(valley_mu_shp_30cm$cluster_4)]) / sum(valley_mu_shp_30cm$area_ac)

CustomBP <- function(x) {
  lower.x <- quantile(x, 0.1, na.rm = TRUE)
  q1.x <- quantile(x, 0.25, na.rm = TRUE)
  med.x <- median(x, na.rm = TRUE)
  q3.x <- quantile(x, 0.75, na.rm = TRUE)
  upper.x <- quantile(x, 0.9, na.rm = TRUE)
  c(lower.x, q1.x, med.x, q3.x, upper.x)
}
#compare key properties across classes
cor_matrix <- cor(na.omit(df_for_clustering))
write.csv(cor_matrix, file.path(dataDir, 'soil_property_correlations.csv'), row.names = FALSE)
cor_matrix_pvals <- cor.mtest(cor_matrix)
cex.before <- par("cex") #saves current cex setting for plotting
cex.reduced <- 0.5 #desired cex setting for plotting p-values
mag.factor <- 2.5 #fudge factor to increase size of axis and legend text
par(cex = cex.reduced)
corrplot(cor_matrix, p.mat = cor_matrix_pvals$p, insig = "p-value", sig.level = -1, tl.cex = par("cex") * mag.factor, cl.cex = par("cex") * mag.factor)
par(cex = cex.before)

#some plotting
plot(cluster_results$clay_30cm, log(cluster_results$ksat_30cm), col=cluster_results$cluster_4+1)
legend()
plot(cluster_results$clay_30cm, cluster_results$cec_30cm, col=cluster_results$cluster_4+1)
plot(cluster_results$clay_30cm, cluster_results$awc_30cm, col=cluster_results$cluster_4+1)
plot(cluster_results$clay_30cm, cluster_results$bd_30cm, col=cluster_results$cluster_4+1)
plot(cluster_results$clay_30cm, cluster_results$lep_30cm, col=cluster_results$cluster_4+1)

#om
boxplot(om_30cm ~ cluster_2, data = cluster_results, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_3, data = cluster_results, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_4, data = cluster_results, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_5, data = cluster_results, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_6, data = cluster_results, ylim=c(0, 5))
boxplot(om_30cm ~ cluster_7, data = cluster_results, ylim=c(0, 5))
#clay
boxplot(clay_30cm ~ cluster_2, data = cluster_results)
boxplot(clay_30cm ~ cluster_3, data = cluster_results)
boxplot(clay_30cm ~ cluster_4, data = cluster_results)
boxplot(clay_30cm ~ cluster_5, data = cluster_results)
boxplot(clay_30cm ~ cluster_6, data = cluster_results)
boxplot(clay_30cm ~ cluster_7, data = cluster_results)
#awc
boxplot(awc_30cm ~ cluster_2, data = cluster_results)
boxplot(awc_30cm ~ cluster_3, data = cluster_results)
boxplot(awc_30cm ~ cluster_4, data = cluster_results)
boxplot(awc_30cm ~ cluster_5, data = cluster_results)
boxplot(awc_30cm ~ cluster_6, data = cluster_results)
boxplot(awc_30cm ~ cluster_7, data = cluster_results)
#cec
boxplot(cec_30cm ~ cluster_2, data = cluster_results)
boxplot(cec_30cm ~ cluster_3, data = cluster_results)
boxplot(cec_30cm ~ cluster_4, data = cluster_results)
boxplot(cec_30cm ~ cluster_5, data = cluster_results)
boxplot(cec_30cm ~ cluster_6, data = cluster_results)
boxplot(cec_30cm ~ cluster_7, data = cluster_results)
#bd
boxplot(bd_30cm ~ cluster_2, data = cluster_results)
boxplot(bd_30cm ~ cluster_3, data = cluster_results)
boxplot(bd_30cm ~ cluster_4, data = cluster_results)
boxplot(bd_30cm ~ cluster_5, data = cluster_results)
boxplot(bd_30cm ~ cluster_6, data = cluster_results)
boxplot(bd_30cm ~ cluster_7, data = cluster_results)

summary(lm(clay_30cm ~ as.factor(cluster_4), data = cluster_results, weights = area_ac))
summary(lm(clay_30cm ~ as.factor(cluster_5), data = cluster_results))
summary(lm(clay_30cm ~ as.factor(cluster_6), data = cluster_results))
clay30cm_lm <- lm(clay_30cm ~ as.factor(cluster_4), data = cluster_results)
anova(clay30cm_lm)
clay30cm_aov <- aov(clay_30cm ~ as.factor(cluster_4), data = cluster_results)
TukeyHSD(clay30cm_aov)
OM30cm_aov <- aov(om_30cm ~ as.factor(cluster_6), data= cluster_results)
TukeyHSD(OM30cm_aov)
summary(lm(om_30cm ~ as.factor(cluster_4), data = cluster_results))
summary(lm(om_30cm ~ as.factor(cluster_5), data = cluster_results))
summary(lm(om_30cm ~ as.factor(cluster_6), data = cluster_results))



table(cluster_results$cluster_4)
table(cluster_results$cluster_5)
table(cluster_results$cluster_6)
tapply(cluster_results$area_ac, cluster_results$cluster_4, sum)
tapply(cluster_results$area_ac, cluster_results$cluster_5, sum)
tapply(cluster_results$area_ac, cluster_results$cluster_6, sum)
cluster_results$area_ac <- as.integer(round(cluster_results$area_ac, 0))
cluster_results_expanded <- cluster_results[rep(rownames(cluster_results), times = cluster_results$area_ac), ]
colnames(cluster_results_expanded)
dim(cluster_results_expanded)
cluster_results_expanded$area_ac <- 1
lapply(cluster_results_expanded[, c("storiemn", 'aws150wta', "Lthc_dep", "Plth_dep", "Drpn_dep", "ATC_dep", "Natr_dep", "Salc_dep", "SCTS_dep", "MRes_dep", 'MnRs_dep', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'kwf_30cm', 'ec_30cm', 'pH_30cm', 'sar_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm')], function(x) tapply(x, cluster_results_expanded$cluster_4, summary))
boxplot(cluster_results_expanded$clay_30cm ~ cluster_results_expanded$cluster_4)
boxplot(cluster_results_expanded$clay_30cm ~ cluster_results_expanded$cluster_5)
boxplot(cluster_results_expanded$clay_30cm ~ cluster_results_expanded$cluster_6)
boxplot(cluster_results_expanded$om_30cm ~ cluster_results_expanded$cluster_4)
boxplot(cluster_results_expanded$om_30cm ~ cluster_results_expanded$cluster_5)
boxplot(cluster_results_expanded$om_30cm ~ cluster_results_expanded$cluster_6)
boxplot(cluster_results_expanded$storiemn ~ cluster_results_expanded$cluster_4)
boxplot(cluster_results_expanded$storiemn ~ cluster_results_expanded$cluster_5)
boxplot(cluster_results_expanded)


MakeBPs <- function(cropname, years=2004:2016) {
  make.bp <- function(bp_name, yaxis_lab, varname, bxfill) {
    png(file.path(dataDir, 'plots', paste0('boxplot_', varname, '.png', sep = '')), family = 'Times New Roman', width = 6.5, height = 4.5, units = 'in', pointsize=11, res = 600)
    par(mai=c(0.9, 0.9, 0.2, 0.2))
    bxp(bp_name, outline = FALSE, boxfill=bxfill, las=2, ylab='', xlab='', cex.axis=1.2)
    mtext(text='Class', side=1, line=3.5, cex = 1.2)
    mtext(text=yaxis_lab, side=2, line=3.5, cex = 1.2)
    dev.off()
  }
  define.bp.stats <- function(bp_name, varname) { #this is to replace default boxplot stats with custom function
    for (i in 1:ncol(bp_name$stats)) { 
      df <- var_df[which(var_df$Model.Year==years[i]),]
      bp_name$stats[,i] <- boxplot(rep(df[[varname]], times=df$cellcounts30m2))
    }
  }
  loadfonts(quiet=TRUE, device='win')
  fnames <- list.files(path=file.path(resultsDir, cropname), pattern = glob2rx('*.csv'))
  for (j in seq_along(fnames)) {
    var_df <- read.csv(file.path(resultsDir, cropname, fnames[j]), stringsAsFactors = FALSE)
    scenario_name <- gsub('_FAO56results_points_rounded.csv', '', fnames[j])
    scenario_name <- paste0('scenario_', gsub(cropname, '', scenario_name))
    bp_GW <- boxplot(GW.ET.growing ~ Model.Year, data=var_df, plot=FALSE)
    bp_BW <- boxplot(Irr.app.total ~ Model.Year, data=var_df, plot=FALSE)
    bp_P <- boxplot(P.WY ~ Model.Year, data=var_df, plot=FALSE)
    bp_ETo <- boxplot(ETo.WY ~ Model.Year, data=var_df, plot=FALSE)
    define.bp.stats(bp_GW, 'GW.ET.growing')
    define.bp.stats(bp_BW, 'Irr.app.total')
    define.bp.stats(bp_P, 'P.WY')
    define.bp.stats(bp_ETo, 'ETo.WY')
    if (!dir.exists(file.path(resultsDir, cropname, 'figures'))) {
      dir.create(file.path(resultsDir, cropname, 'figures'))
    }
    if (!dir.exists(file.path(resultsDir, cropname, 'figures', scenario_name))) {
      dir.create(file.path(resultsDir, cropname, 'figures', scenario_name))
    }
    make.bp(bp_name = bp_GW, yaxis_lab = 'Growing season green water (mm)', varname = 'GW.ET.growing', bxfill = 'green')
    make.bp(bp_name = bp_BW, yaxis_lab = 'Irrigation [blue water] demand (mm)', varname = 'Irr.app.total', bxfill = 'lightblue')
    make.bp(bp_name = bp_P, yaxis_lab = 'Water year precipitation (mm)', varname = 'P.WY', bxfill = bxfill = 'blue')
    make.bp(bp_name = bp_ETo, yaxis_lab = 'Reference evapotranspiration (mm)', varname = 'ETo.WY', bxfill = 'orange')
  } #bp_name, yaxis_lab, varname, bxfill
}

test <- data.frame(A=c(1, 2, 3), B=c(4, 5, 6), times=c(2, 10, 20))
test[rep(row.names(test), test$times), 1:2]
apply(test, 1, function(x) rbind(rep(x, times=x[3])))
lapply(split(test, rownames(test)), )
