library(raster)
dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_trial'
list.files(dataDir)
valley_mu_aea_30cm <- shapefile(file.path(dataDir, 'valley_30cm.shp'))
valley30cm <- read.csv(file.path(dataDir, 'valley_30cm_test.csv'), stringsAsFactors = FALSE)
dim(valley30cm)
#acres_by_mukey <- tapply(valley30cm$area_ac, valley30cm$mukey, sum)
acres_by_mukey <- aggregate(area_ac ~ mukey, data = valley30cm, sum)
table(valley30cm$mukey)
length(unique(valley30cm$mukey)) #786
valley30cm_by_mukey <- valley30cm[!duplicated(valley30cm$mukey), ]
dim(valley30cm_by_mukey)
length(unique(valley30cm_by_mukey$mukey))
valley30cm_by_mukey$area_ac <- acres_by_mukey$area_ac[match(valley30cm_by_mukey$mukey, acres_by_mukey$mukey)]
sum(valley30cm_by_mukey$area_ac) #250800
valley30cm_by_mukey$area_proportion <- valley30cm_by_mukey$area_ac / sum(valley30cm_by_mukey$area_ac)
hist(valley30cm_by_mukey$area_proportion)
valley30cm_by_mukey[which.max(valley30cm_by_mukey$area_proportion),]
sum(valley30cm$Drpn_dep == 200, na.rm = TRUE)
valley30cm_by_mukey$Drpn_dep
colnames(valley30cm_by_mukey)
valley30cm_by_mukey$muname[is.na(valley30cm_by_mukey$aws050wta)]
df_for_clustering <- valley30cm_by_mukey[!is.na(valley30cm_by_mukey$aws050wta), ]
rownames(df_for_clustering) <- df_for_clustering$mukey
summary(df_for_clustering)
df_for_clustering <- df_for_clustering[,c('area_proportion', "Lthc_dep", "Plth_dep", "Drpn_dep", "ATC_dep", "Natr_dep", "Salc_dep", "SCTS_dep", "MRes_dep", 'storiemn', 'clay_30cm', 'om_30cm', 'cec_30cm', 'bd_30cm', 'kwf_30cm', 'ec_30cm', 'pH_30cm', 'sar_30cm', 'lep_30cm', 'ksat_30cm', 'awc_30cm')] #frags have 160 mukeys with NA
df_for_clustering_scaled <- scale(df_for_clustering[,2:ncol(df_for_clustering)])
summary(df_for_clustering_scaled) #awkward
cluster_3 <- kmeans(na.omit(df_for_clustering_scaled), 3)
cluster_3
cluster_4 <- kmeans(na.omit(df_for_clustering_scaled), 4)
cluster_4
cluster_5 <- kmeans(na.omit(df_for_clustering_scaled), 5)
cluster_5
cluster_6 <- kmeans(na.omit(df_for_clustering_scaled), 6)
cluster_6
cluster_7 <- kmeans(na.omit(df_for_clustering_scaled), 7)
cluster_7
cluster_8 <- kmeans(na.omit(df_for_clustering_scaled), 8)
cluster_8
cluster_9 <- kmeans(na.omit(df_for_clustering_scaled), 9)
cluster_9
cluster_10 <- kmeans(na.omit(df_for_clustering_scaled), 10)
cluster_10
cluster_11 <- kmeans(na.omit(df_for_clustering_scaled), 11)
cluster_11
cluster_12 <- kmeans(na.omit(df_for_clustering_scaled), 12)
cluster_12
cluster_13 <- kmeans(na.omit(df_for_clustering_scaled), 13)
cluster_13
results <- sapply(1:20, function(x) {
  result <- kmeans(na.omit(df_for_clustering_scaled), x)
 round(100 * result$betweenss / result$totss, 1)})
plot(1:20, results, type='b')
valley30cm_by_mukey$cluster_4 <- cluster_4$cluster[match(valley30cm_by_mukey$mukey, names(cluster_4$cluster))]
sum(valley30cm_by_mukey$area_ac[!is.na(valley30cm_by_mukey$cluster_4)]) #1535232
2377819 / 2508000 #94.8% got a cluster if we ignore fragments

valley_mu_aea_30cm$cluster_4 <- cluster_4$cluster[match(valley_mu_aea_30cm$mukey, names(cluster_4$cluster))]
valley_mu_aea_30cm$cluster_6 <- cluster_6$cluster[match(valley_mu_aea_30cm$mukey, names(cluster_6$cluster))]
valley_mu_aea_30cm$cluster_8 <- cluster_8$cluster[match(valley_mu_aea_30cm$mukey, names(cluster_8$cluster))]
valley_mu_aea_30cm$cluster_10 <- cluster_10$cluster[match(valley_mu_aea_30cm$mukey, names(cluster_10$cluster))]
valley_mu_aea_30cm$cluster_13 <- cluster_13$cluster[match(valley_mu_aea_30cm$mukey, names(cluster_13$cluster))]
shapefile(valley_mu_aea_30cm, file.path(dataDir, 'valley_30cm.shp'), overwrite=TRUE)
sum(valley_mu_aea_30cm$area_ac[!is.na(valley_mu_aea_30cm$cluster_4)]) / sum(valley_mu_aea_30cm$area_ac)
