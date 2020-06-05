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
  resultsDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final/FINAL results/AOV results'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}
om_to_oc <- 1.72
crit_pH <- 7.8
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)
dom_order_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', "soil survey facts", 'dom_order_by_mukey.csv'), stringsAsFactors = FALSE)

kssl_ssurgo_extract <- read.csv(file.path(ksslDir, 'kssl_pts_ssurgo_30cm_extract_FINAL.csv'), stringsAsFactors = FALSE)
kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_FINAL.csv'), stringsAsFactors = FALSE) #updated 2/11/20 to consider total C when org C not available but soil pH sufficiently low for calculating SOC content
kssl_points_30cm$labsampnum <- kssl_ssurgo_extract$labsampnum[match(kssl_points_30cm$pedon_key, kssl_ssurgo_extract$pedon_key)]
sum(grepl('UCD', kssl_points_30cm$labsampnum)) #54
colnames(kssl_points_30cm)
dim(kssl_points_30cm)
kssl_points_30cm$mukey <- kssl_ssurgo_extract$mukey[match(kssl_points_30cm$pedon_key, kssl_ssurgo_extract$pedon_key)]
kssl_points_30cm$dom_order <- dom_order_by_mukey$dom_order[match(kssl_points_30cm$mukey, dom_order_by_mukey$mukey)]
tapply(kssl_points_30cm$om_30cm, kssl_points_30cm$dom_order, summary)
table(kssl_points_30cm$dom_order)
sum(is.na(kssl_points_30cm$clay_30cm)) #55 are NA
sum(is.na(kssl_points_30cm$oc_30cm)) #88 are NA
sum(is.na(kssl_points_30cm$cluster_7))
sampnumbers_kgOrgC <- table(kssl_points_30cm$cluster_7[!is.na(kssl_points_30cm$kgOrg.m2_30cm)])
names(sampnumbers_kgOrgC) <- clus_7_names
sampnumbers_kgOrgC

sampnumbers_OC <- table(kssl_points_30cm$cluster_7[!is.na(kssl_points_30cm$oc_30cm)])
names(sampnumbers_OC) <- clus_7_names
sampnumbers_OC

#replace OC with totC when pH sufficiently low
kssl_points_30cm$oc_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)] <- kssl_points_30cm$c_tot_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)]
sum(is.na(kssl_points_30cm$oc)) #NAs reduced from 88 to 56

#replace estimated om with
kssl_points_30cm$om_30cm_v2 <- kssl_points_30cm$oc_30cm * om_to_oc
plot(kssl_points_30cm$om_30cm, kssl_points_30cm$om_30cm_v2)
sum(is.na(kssl_points_30cm$om_30cm_v2))
sum(is.na(kssl_points_30cm$om_30cm))
kssl_points_30cm$om_30cm <- kssl_points_30cm$om_30cm_v2
kssl_points_30cm$om_30cm_v2 <- NULL
kssl_points_30cm[139:150,]
valley30cm_by_mukey[valley30cm_by_mukey$mukey==460870,]
dom_order_by_mukey[dom_order_by_mukey$mukey==460870,]
sum(is.na(dom_order_by_mukey$dom_order)) #89
lapply(kssl_points_30cm, function(x) sum(!is.na(x)))

# kssl_points_100cm <- read.csv(file.path(ksslDir, 'kssl_cluster_100cm_NArm.csv'), stringsAsFactors = FALSE)
# colnames(kssl_points_100cm)
# table(kssl_points_100cm$cluster_7[!is.na(kssl_points_100cm$kgOrg.m2_100cm)])
# test <- table(kssl_points_100cm$cluster_7[!is.na(kssl_points_100cm$kgOrg.m2_100cm)])
# names(test) <- clus_7_names
# test

df_combined1 <-  kssl_points_30cm[ ,c('clay_30cm', 'om_30cm', 'bd_13b_30cm', 'pH_H2O_30cm', 'kgOrg.m2_30cm', colnames(kssl_points_30cm)[grepl('cluster_', colnames(kssl_points_30cm))], 'dom_order')]
dim(df_combined1) #369 rows
colnames(df_combined1)[3] <- 'bd_30cm'
df_combined1$source <- 'KSSL'

kerri_points_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
colnames(kerri_points_30cm)
kerri_points_30cm$om_30cm <- kerri_points_30cm$totC_30cm * om_to_oc
kerri_ssurgo_extract <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_pts_ssurgo_30cm_extract_FINAL.csv'), stringsAsFactors = FALSE)
colnames(kerri_ssurgo_extract)
kerri_points_30cm$mukey <- kerri_ssurgo_extract$mukey[match(kerri_points_30cm$Concatenate, kerri_ssurgo_extract$Concatenate)]
kerri_points_30cm$dom_order <- dom_order_by_mukey$dom_order[match(kerri_points_30cm$mukey, dom_order_by_mukey$mukey)]
tapply(kerri_points_30cm$om_30cm, kerri_points_30cm$dom_order, summary)
table(kerri_points_30cm$dom_order)

kerri_metadata <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE)
unique(kerri_metadata$Area)
length(unique(kerri_metadata$Concatenate))#127 unique points, 102 in study area

kerri_points_30cm$vineyard_region <- kerri_metadata$Area[match(kerri_points_30cm$Concatenate, kerri_metadata$Concatenate)]
kerri_points_30cm$vineyard_name <- kerri_metadata$Vineyard.Management[match(kerri_points_30cm$Concatenate, kerri_metadata$Concatenate)]
table(kerri_points_30cm$vineyard_region) #30 in Lodi, 97 in Napa
table(kerri_points_30cm$vineyard_region[!is.na(kerri_points_30cm$cluster_7)]) #27 in Lodi, 76 in Napa in area of interest
unique(kerri_points_30cm$vineyard_name[kerri_points_30cm$vineyard_region=='Napa']) #16 vineyards in Napa (paper says ninteen)
unique(kerri_points_30cm$vineyard_name[kerri_points_30cm$vineyard_region=='Lodi']) #5 vineyards (paper says nine)
tapply(kerri_points_30cm$Concatenate[kerri_points_30cm$vineyard_region=='Napa'], kerri_points_30cm$vineyard_name[kerri_points_30cm$vineyard_region=='Napa'], function(x) length(x))
tapply(kerri_points_30cm$Concatenate[kerri_points_30cm$vineyard_region=='Lodi'], kerri_points_30cm$vineyard_name[kerri_points_30cm$vineyard_region=='Lodi'], function(x) length(x))
kerri_points_30cm[kerri_points_30cm$vineyard_name=='Big Ranch Vineyard',]
tapply(kerri_points_30cm$cluster_7, kerri_points_30cm$vineyard_name, function(x) length(unique(x)))
tapply(kerri_points_30cm$cluster_7, kerri_points_30cm$vineyard_name, function(x) unique(x))
tapply(kerri_points_30cm$cluster_7, kerri_points_30cm$vineyard_region, function(x) unique(x))
tapply(kerri_points_30cm$totC_30cm, kerri_points_30cm$vineyard_region, function(x) length(unique(x)))
kerri_points_30cm$shr7_name <- cluster_

df_combined2 <- kerri_points_30cm[!is.na(kerri_points_30cm$cluster_7), c('clay_30cm', 'om_30cm', 'bd_30cm', 'pH_H2O_30cm', 'kgOrg.m2_30cm', colnames(kerri_points_30cm)[grepl('cluster_', colnames(kerri_points_30cm))], 'dom_order')]
dim(df_combined2) #103 rows
df_combined2$source <- 'CDFA'

df_combined <- rbind(df_combined1, df_combined2)
dim(df_combined) #472
table(df_combined$dom_order[!is.na(df_combined$om_30cm)])
df_combined$om_30cm[df_combined$dom_order=='Ultisols']
lapply(df_combined, function(x) sum(!is.na(x)))

clus_9_names <- c('4. Coarse w/pans', '5. Loamy w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions',  '7. Coarse saline-sodic', '2. Coarse w/no restrictions', '8. Fine saline-sodic', '6. Loamy w/pans (high OM)', '9. Fine shrink-swell') #order corrected 4/7/20

clus_8_names <- c('6. Coarse saline-sodic', '7. Fine saline-sodic', '5. Loamy w/pans', '2. Coarse w/no restrictions', '8. Fine shrink-swell', '4. Coarse w/pans', '3. Loamy w/no restrictions', '1. Very coarse w/no restrictions') #order corrected 4/7/20

clus_6_names <- c('3. Coarse w/pans', '4. Loamy w/pans', '5. Loamy saline-sodic', '1. Coarse w/no restrictions', '6. Fine shrink-swell', '2. Loamy w/no restrictions') #order corrected 4/7/20

clus_5_names <- c('3. Loamy w/pans', '4. Loamy saline-sodic', '2. Loamy w/no restrictions', '5. Fine shrink-swell', '1. Coarse w/no restrictions')

clus_4_names <- c('1. Coarse w/no restrictions', '2. Loamy w/pans', '3. Loamy saline-sodic', '4. Fine shrink-swell')

clus_3_names <- c('1. Loamy w/pans', '3. Fine shrink-swell', '2. Coarse saline sodic')

test <- TukeyHSD(aov(kssl_points_30cm$clay_30cm ~ as.factor(kssl_points_30cm$cluster_7)))
test[[1]][,4]

compare_region_means <- function(df, names, cluster_no, y, alpha) {
  x <- names[match(df[[paste0('cluster_', cluster_no)]], 1:cluster_no)]
  y_vector <- df[[y]]
  print('These points have data (by soil order)')
  print(table(x[!is.na(y_vector)]))
  print(paste('Total of', sum(table(x[!is.na(y_vector)])), 'points with data.'))
  print(paste(sum(TukeyHSD(aov(y_vector ~ x))[[1]][,4] <= alpha), 'out of ', length(TukeyHSD(aov(y_vector ~ x))[[1]][,4]), 'are significant.'))
  # print(paste0(y, '_cluster_', as.character(cluster_no), '_aov_results.txt'))
  if(!dir.exists(file.path(resultsDir, paste(as.character(cluster_no), 'region')))) {dir.create(file.path(resultsDir, paste(as.character(cluster_no), 'region')))}
  sink(file=file.path(resultsDir, paste(as.character(cluster_no), 'region'), paste0(y, '_cluster_', as.character(cluster_no), '_aov_results_TEST.txt')))
  print(paste('ANOVA and Tukey HSD results for', y, 'and', cluster_no, '-region model:'))
  print(paste(sum(TukeyHSD(aov(y_vector ~ x))[[1]][,4] <= alpha), 'out of ', length(TukeyHSD(aov(y_vector ~ x))[[1]][,4]), 'are significant.'))
  print(paste0(y, '_cluster_', as.character(cluster_no), '_aov_results.txt'))
  print(summary(aov(y_vector ~ x)))
  print(TukeyHSD(aov(y_vector ~ x), ordered=FALSE))
  sink()
  #plot(TukeyHSD(aov(y_vector ~ x), ordered=FALSE))
  # model.tables(aov(y_vector ~ x))
}

#9 region
compare_region_means(df=df_combined1, names=clus_9_names, cluster_no=9, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_9_names, cluster_no=9, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_9_names, cluster_no=9, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_9_names, cluster_no=9, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_9_names, cluster_no=9, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_9_names, cluster_no=9, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_9_names, cluster_no=9, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_9_names, cluster_no=9, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_9_names, cluster_no=9, y='kgOrg.m2_30cm', alpha = 0.05)


#8 region
compare_region_means(df=df_combined1, names=clus_8_names, cluster_no=8, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_8_names, cluster_no=8, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_8_names, cluster_no=8, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_8_names, cluster_no=8, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_8_names, cluster_no=8, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_8_names, cluster_no=8, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_8_names, cluster_no=8, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_8_names, cluster_no=8, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_8_names, cluster_no=8, y='kgOrg.m2_30cm', alpha = 0.05)

#7 region
compare_region_means(df=df_combined1, names=clus_7_names, cluster_no=7, y='clay_30cm', alpha = 0.05) #kssl only
compare_region_means(df=df_combined, names=clus_7_names, cluster_no=7, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_7_names, cluster_no=7, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_7_names, cluster_no=7, y='pH_H2O_30cm', alpha = 0.05) #verified that re-aggregation on 2/11/20 did not change soil pH results (only OM and kg SOC)
compare_region_means(df=kssl_points_30cm, names=clus_7_names, cluster_no=7, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_7_names, cluster_no=7, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_7_names, cluster_no=7, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_7_names, cluster_no=7, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_7_names, cluster_no=7, y='kgOrg.m2_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_100cm, names=clus_7_names, cluster_no=7, y='kgOrg.m2_100cm', alpha = 0.05)

meanclay_clus7 <- tapply(df_combined$clay_30cm, clus_7_names[match(df_combined$cluster_7, 1:7)], mean, na.rm=TRUE)
meanclay_clus7[order(meanclay_clus7)]

meanOM_clus7 <- tapply(df_combined$om_30cm, clus_7_names[match(df_combined$cluster_7, 1:7)], mean, na.rm=TRUE)
meanOM_clus7[order(meanOM_clus7)]

meankgOC_clus7 <- tapply(df_combined$kgOrg.m2_30cm, clus_7_names[match(df_combined$cluster_7, 1:7)], mean, na.rm=TRUE)
meankgOC_clus7[order(meankgOC_clus7)]

meanBD_clus7 <- tapply(df_combined$bd_30cm, clus_7_names[match(df_combined$cluster_7, 1:7)], mean, na.rm=TRUE)
meanBD_clus7[order(meanBD_clus7)]

meanpH_clus7 <- tapply(df_combined$pH_H2O_30cm, clus_7_names[match(df_combined$cluster_7, 1:7)], mean, na.rm=TRUE)
meanpH_clus7[order(meanpH_clus7)]

meanCEC_clus7 <- tapply(kssl_points_30cm$cec_7_30cm, clus_7_names[match(kssl_points_30cm$cluster_7, 1:7)], mean, na.rm=TRUE)
meanCEC_clus7[order(meanCEC_clus7)]

meanLEP_clus7 <- tapply(kssl_points_30cm$lep_30cm, clus_7_names[match(kssl_points_30cm$cluster_7, 1:7)], mean, na.rm=TRUE)
meanLEP_clus7[order(meanLEP_clus7)]

meanEC_clus7 <- tapply(kssl_points_30cm$ec_30cm, clus_7_names[match(kssl_points_30cm$cluster_7, 1:7)], mean, na.rm=TRUE)
meanEC_clus7[order(meanEC_clus7)]

meanAWC_clus7 <- tapply(kssl_points_30cm$awc_30cm, clus_7_names[match(kssl_points_30cm$cluster_7, 1:7)], mean, na.rm=TRUE)
meanAWC_clus7[order(meanAWC_clus7)]

meankgOC_clus7 <- tapply(kssl_points_30cm$kgOrg.m2_30cm, clus_7_names[match(kssl_points_30cm$cluster_7, 1:7)], mean, na.rm=TRUE)
meankgOC_clus7[order(meankgOC_clus7)]

#6 region
compare_region_means(df=df_combined1, names=clus_6_names, cluster_no=6, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_6_names, cluster_no=6, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_6_names, cluster_no=6, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_6_names, cluster_no=6, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_6_names, cluster_no=6, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_6_names, cluster_no=6, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_6_names, cluster_no=6, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_6_names, cluster_no=6, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_6_names, cluster_no=6, y='kgOrg.m2_30cm', alpha = 0.05)


meanclay_clus6 <- tapply(df_combined$clay_30cm, clus_6_names[match(df_combined$cluster_6, 1:6)], mean, na.rm=TRUE)
meanclay_clus6[order(meanclay_clus6)]

meanOM_clus6 <- tapply(df_combined$om_30cm, clus_6_names[match(df_combined$cluster_6, 1:6)], mean, na.rm=TRUE)
meanOM_clus6[order(meanOM_clus6)]

meanpH_clus6 <- tapply(df_combined$pH_H2O_30cm, clus_6_names[match(df_combined$cluster_6, 1:6)], mean, na.rm=TRUE)
meanpH_clus6[order(meanpH_clus6)]

meanCEC_clus6 <- tapply(kssl_points_30cm$cec_7_30cm, clus_6_names[match(kssl_points_30cm$cluster_6, 1:6)], mean, na.rm=TRUE)
meanCEC_clus6[order(meanCEC_clus6)]

meanLEP_clus6 <- tapply(kssl_points_30cm$lep_30cm, clus_6_names[match(kssl_points_30cm$cluster_6, 1:6)], mean, na.rm=TRUE)
meanLEP_clus6[order(meanLEP_clus6)]

meanEC_clus6 <- tapply(kssl_points_30cm$ec_30cm, clus_6_names[match(kssl_points_30cm$cluster_6, 1:6)], mean, na.rm=TRUE)
meanEC_clus6[order(meanEC_clus6)]

meanAWC_clus6 <- tapply(kssl_points_30cm$awc_30cm, clus_6_names[match(kssl_points_30cm$cluster_6, 1:6)], mean, na.rm=TRUE)
meanAWC_clus6[order(meanAWC_clus6)]

#5 region
compare_region_means(df=df_combined1, names=clus_5_names, cluster_no=5, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_5_names, cluster_no=5, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_5_names, cluster_no=5, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_5_names, cluster_no=5, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_5_names, cluster_no=5, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_5_names, cluster_no=5, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_5_names, cluster_no=5, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_5_names, cluster_no=5, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_5_names, cluster_no=5, y='kgOrg.m2_30cm', alpha = 0.05)


#4 region
compare_region_means(df=df_combined1, names=clus_4_names, cluster_no=4, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_4_names, cluster_no=4, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_4_names, cluster_no=4, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_4_names, cluster_no=4, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_4_names, cluster_no=4, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_4_names, cluster_no=4, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_4_names, cluster_no=4, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_4_names, cluster_no=4, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_4_names, cluster_no=4, y='kgOrg.m2_30cm', alpha = 0.05)

#3 region
compare_region_means(df=df_combined1, names=clus_3_names, cluster_no=3, y='clay_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_3_names, cluster_no=3, y='om_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_3_names, cluster_no=3, y='bd_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_3_names, cluster_no=3, y='pH_H2O_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_3_names, cluster_no=3, y='cec_7_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_3_names, cluster_no=3, y='lep_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_3_names, cluster_no=3, y='ec_30cm', alpha = 0.05)
compare_region_means(df=kssl_points_30cm, names=clus_3_names, cluster_no=3, y='awc_30cm', alpha = 0.05)
compare_region_means(df=df_combined, names=clus_3_names, cluster_no=3, y='kgOrg.m2_30cm', alpha = 0.05)

#make group contrast by cluster figure
contrast_summary <- data.frame(cluster_no=3:9, contrast_no=c(3,6,10,15,21,28,36), clay=c(3,4,7,8,14,17,17), om=c(2,4,6,9,11,13,17), bd=c(0,0,2,2,2,3,2), pH=c(2,6,9,12,16,18,18), cec=c(3,5,8,9,13,16,16), lep=c(2,3,6,7,11,11,11), ec=c(1,3,2,1,5,6,6), kgSOC=c(2,4,5,8,11,8,12)) #updated on 4/27/20 so that clay analysis only includes KSSL data

mean(sapply(contrast_summary[1,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[1], 2)})) #0.63 for 3
mean(sapply(contrast_summary[2,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[2], 2)})) #0.63 for 4
mean(sapply(contrast_summary[3,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[3], 2)})) #0.56 for 5
mean(sapply(contrast_summary[4,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[4], 2)})) #0.47 for 6
mean(sapply(contrast_summary[5,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[5], 2)})) #0.49 for 7
mean(sapply(contrast_summary[6,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[6], 2)})) #0.41 for 8
mean(sapply(contrast_summary[7,3:ncol(contrast_summary)], function(x) {round(x/contrast_summary$contrast_no[7], 2)})) #0.35

tiff(file = file.path(FiguresDir, 'FINAL', 'sig_group_contrasts_3to9.tif'), family = 'Times New Roman', width = 6.5, height = 5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(4,4,4.5,1))
plot(x=3:9, y=100*contrast_summary$clay/contrast_summary$contrast_no, xlab='', ylab='', type='b', col='grey', ylim=c(0,100), pch=21, bg='grey')
mtext(text='Number of soil health regions in conceptual model', side=1, line=2.5)
mtext(text='Significant group contrasts, validation data (%)', side=2, line=2.5)
lines(x=3:9, y=100*contrast_summary$om/contrast_summary$contrast_no, col='black', type='b', pch=21, bg='black')
lines(x=3:9, y=100*contrast_summary$bd/contrast_summary$contrast_no, col='brown1', type='b', pch=21, bg='brown1')
lines(x=3:9, y=100*contrast_summary$pH/contrast_summary$contrast_no, col='red3', type='b', pch=21, bg='red3')
lines(x=3:9, y=100*contrast_summary$cec/contrast_summary$contrast_no, col='blue', type='b', pch=21, bg='blue')
lines(x=3:9, y=100*contrast_summary$lep/contrast_summary$contrast_no, col='darkviolet', type='b', pch=21, bg='darkviolet')
lines(x=3:9, y=100*contrast_summary$ec/contrast_summary$contrast_no, col='darkorange', type='b', pch=21, bg='darkorange')
lines(x=3:9, y=100*contrast_summary$kgSOC/contrast_summary$contrast_no, col='chocolate4', type='b', pch=21, bg='chocolate4')
abline(h=5, lty=2)
legend(x=3.75,y=128.5, legend = c('pH', 'Clay', 'CEC', 'LEP', 'OM (%)', expression('OC (kg m'^-2*')'), 'EC', 'BD', 'Null'), col=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'black'), lty=c(rep(1, 8), 2), pch=c(rep(21, 8), NA), pt.bg=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'black', NA), ncol = 3, xpd=TRUE)
dev.off()

tiff(file = file.path(FiguresDir, 'FINAL', 'sig_group_contrasts_3to9_not_rel.tif'), family = 'Times New Roman', width = 6.5, height = 3.5, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(4,4,0.5,0.5))
plot(x=3:9, y=contrast_summary$clay, xlab='', ylab='', type='b', col='grey', ylim=c(0,40), pch=21, bg='grey')#, xaxt='n')
mtext(text='Number of soil health regions in conceptual model', side=1, line=2.5)
mtext(text='Number of significant group contrasts, validation data', side=2, line=2.5, at=15)
lines(x=3:9, y=contrast_summary$om, col='black', type='b', pch=21, bg='black')
lines(x=3:9, y=contrast_summary$bd, col='brown1', type='b', pch=21, bg='brown1')
lines(x=3:9, y=contrast_summary$pH, col='red3', type='b', pch=21, bg='red3')
lines(x=3:9, y=contrast_summary$cec, col='blue', type='b', pch=21, bg='blue')
lines(x=3:9, y=contrast_summary$lep, col='darkviolet', type='b', pch=21, bg='darkviolet')
lines(x=3:9, y=contrast_summary$ec, col='darkorange', type='b', pch=21, bg='darkorange')
lines(x=3:9, y=contrast_summary$kgSOC, col='chocolate4', type='b', pch=21, bg='chocolate4')
lines(x=3:9, y=contrast_summary$contrast_no, col='grey', type='l', lty=2)
text(x=9, y=39, 'b')
legend(x=2.8, y=40.5, legend = c('pH', 'Clay', 'CEC', 'LEP', 'OM (%)', expression('OC (kg m'^-2*')'), 'EC', 'BD', 'Total possible contrasts'), col=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'grey'), lty=c(rep(1, 8), 2), pch=c(rep(21, 8), NA), pt.bg=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'black', NA), ncol = 2, inset = 0.005)
dev.off()

# tiff(file = file.path(FiguresDir, 'v2', 'validation plots', 'sig_group_contrasts_3to9_no_legend.tif'), family = 'Times New Roman', width = 6.5, height = 4, pointsize = 11, units = 'in', res=800, compression = 'lzw')
# par(mar=c(4,4,1,1))
# plot(x=3:9, y=100*contrast_summary$clay/contrast_summary$contrast_no, xlab='', ylab='', type='b', col='grey', ylim=c(0,100), pch=21, bg='grey')
# mtext(text='Number of soil health regions in conceptual model', side=1, line=2.5)
# mtext(text='Significant group contrasts of total possible, validation data (%)', side=2, line=2.5, at=42)
# lines(x=3:9, y=100*contrast_summary$om/contrast_summary$contrast_no, col='black', type='b', pch=21, bg='black')
# lines(x=3:9, y=100*contrast_summary$bd/contrast_summary$contrast_no, col='brown1', type='b', pch=21, bg='brown1')
# lines(x=3:9, y=100*contrast_summary$pH/contrast_summary$contrast_no, col='red3', type='b', pch=21, bg='red3')
# lines(x=3:9, y=100*contrast_summary$cec/contrast_summary$contrast_no, col='blue', type='b', pch=21, bg='blue')
# lines(x=3:9, y=100*contrast_summary$lep/contrast_summary$contrast_no, col='darkviolet', type='b', pch=21, bg='darkviolet')
# lines(x=3:9, y=100*contrast_summary$ec/contrast_summary$contrast_no, col='darkorange', type='b', pch=21, bg='darkorange')
# lines(x=3:9, y=100*contrast_summary$kgSOC/contrast_summary$contrast_no, col='chocolate4', type='b', pch=21, bg='chocolate4')
# abline(h=5, lty=2)
# text(x=9, y=95, 'b')
# # legend(x=3.75,y=128.5, legend = c('pH', 'Clay', 'CEC', 'LEP', 'OM (%)', expression('OC (kg m'^-2*')'), 'EC', 'BD', 'Null'), col=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'black'), lty=c(rep(1, 8), 2), pch=c(rep(21, 8), NA), pt.bg=c('red3', 'grey', 'blue', 'darkviolet', 'black', 'chocolate4', 'darkorange', 'brown1', 'black', NA), ncol = 3, xpd=TRUE)
# dev.off()

#compare order means
TukeyHSD()
compare_SoilOrder_means <- function(df, y, alpha) {
  y_vector <- df[[y]]
  x <- df[['dom_order']]
  print('These points have data (by soil order)')
  print(table(x[!is.na(y_vector)]))
  print(paste('Total of', sum(table(x[!is.na(y_vector)])), 'points with data.'))
  print(paste(sum(TukeyHSD(aov(y_vector ~ x))[[1]][,4] <= alpha), 'out of ', length(TukeyHSD(aov(y_vector ~ x))[[1]][,4]), 'possible contrasts are significant.'))
  sink(file=file.path(resultsDir, 'SoilOrders', paste0(y, '_SoilOrder_aov_results.txt')))
  print(paste('ANOVA and Tukey HSD results for', y, 'and soil orders'))
  print(paste(sum(TukeyHSD(aov(y_vector ~ x))[[1]][,4] <= alpha), 'out of ', length(TukeyHSD(aov(y_vector ~ x))[[1]][,4]), 'are significant.'))
  print(paste0(y, '_SoilOrder_aov_results.txt'))
  print(summary(aov(y_vector ~ x)))
  print(TukeyHSD(aov(y_vector ~ x), ordered=FALSE))
  sink()
  #plot(TukeyHSD(aov(y_vector ~ x), ordered=FALSE))
  # model.tables(aov(y_vector ~ x))
}

compare_SoilOrder_means(df=df_combined1, y='clay_30cm', alpha = 0.05) #kssl only
compare_SoilOrder_means(df=df_combined, y='om_30cm', alpha = 0.05)
compare_SoilOrder_means(df=df_combined, y='bd_30cm', alpha = 0.05)
compare_SoilOrder_means(df=df_combined, y='pH_H2O_30cm', alpha = 0.05) #verified that re-aggregation on 2/11/20 did not change soil pH results (only OM and kg SOC)
compare_SoilOrder_means(df=kssl_points_30cm, y='cec_7_30cm', alpha = 0.05)
compare_SoilOrder_means(df=kssl_points_30cm, y='lep_30cm', alpha = 0.05)
compare_SoilOrder_means(df=kssl_points_30cm, y='ec_30cm', alpha = 0.05)
compare_SoilOrder_means(df=kssl_points_30cm, y='awc_30cm', alpha = 0.05)
compare_SoilOrder_means(df=df_combined, y='kgOrg.m2_30cm', alpha = 0.05)

pH_by_order <- tapply(df_combined$pH_H2O_30cm, df_combined$dom_order, mean, na.rm=TRUE)
pH_by_order[order(pH_by_order)]
om_by_order <- tapply(df_combined$om_30cm, df_combined$dom_order, mean, na.rm=TRUE)
om_by_order[order(om_by_order)]
clay_by_order <- tapply(df_combined$clay_30cm, df_combined$dom_order, mean, na.rm=TRUE)
clay_by_order[order(clay_by_order)]

#check treatment types by soil health region
dim(kerri_points_30cm)
sum(is.na(kerri_points_30cm$cluster_7))
napa_lodi_pts_30cm <- kerri_points_30cm[!is.na(kerri_points_30cm$cluster_7),]
napa_lodi_pts_30cm$SHR7name <- clus_7_names[napa_lodi_pts_30cm$cluster_7]
table(napa_lodi_pts_30cm$SHR7name)
table(napa_lodi_pts_30cm$tillage)
tapply(napa_lodi_pts_30cm$tillage, napa_lodi_pts_30cm$SHR7name, table)
tapply(napa_lodi_pts_30cm$tillage, napa_lodi_pts_30cm$vineyard_region, table)
table(napa_lodi_pts_30cm$irrigated_vs_dryfarm)
tapply(napa_lodi_pts_30cm$irrigated_vs_dryfarm, napa_lodi_pts_30cm$SHR7name, table)
tapply(napa_lodi_pts_30cm$irrigated_vs_dryfarm, napa_lodi_pts_30cm$vineyard_region, table)
table(napa_lodi_pts_30cm$compost_added)
tapply(napa_lodi_pts_30cm$compost_added, napa_lodi_pts_30cm$SHR7name, table)
tapply(napa_lodi_pts_30cm$compost_added, napa_lodi_pts_30cm$vineyard_region, table)
table(napa_lodi_pts_30cm$management_type)
tapply(napa_lodi_pts_30cm$management_type, napa_lodi_pts_30cm$SHR7name, table)
tapply(napa_lodi_pts_30cm$management_type, napa_lodi_pts_30cm$vineyard_region, table)
table(napa_lodi_pts_30cm$vineyard_region)
tapply(napa_lodi_pts_30cm$vineyard_region, napa_lodi_pts_30cm$SHR7name, table)
