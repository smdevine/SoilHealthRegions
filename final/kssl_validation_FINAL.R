#TO-DO: run soil horizon aggregation by excluding any 0-30 cm variables with missing data for a horizon
library(raster)
library(aqp)
library(vioplot)
library(extrafont)
library(extrafontdb)
loadfonts(device = 'win')
laptop <- TRUE

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
}
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')
list.files(ksslDir)
kssl_shp <- shapefile(file.path(ksslDir, 'shapefiles', 'kssl_CA.shp'))
crs(kssl_shp)
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'FINAL results', 'shapefiles with data', 'valley_30cm_cluster.shp'))
crs(valley_mu_shp_30cm)
kssl_ssurgo_30cm <- extract(valley_mu_shp_30cm, kssl_shp, df=TRUE)
# sum(kssl_ssurgo_30cm$cluster_10 == old_data$cluster_10, na.rm = TRUE) #made no difference whether NAD83 or WGS84 was specified
# old_data <- kssl_ssurgo_30cm
# head(kssl_ssurgo_30cm)
# dim(kssl_ssurgo_30cm)
# table(kssl_ssurgo_30cm$cluster_5)
# table(kssl_ssurgo_30cm$cluster_6)
# table(kssl_ssurgo_30cm$cluster_7)
# table(kssl_ssurgo_30cm$cluster_8)
# table(kssl_ssurgo_30cm$cluster_9)
# table(kssl_ssurgo_30cm$cluster_10)
kssl_ssurgo_30cm <- kssl_ssurgo_30cm[!is.na(kssl_ssurgo_30cm$cluster_7), ]
kssl_ssurgo_30cm$pedon_key <- kssl_shp$pedn_ky[kssl_ssurgo_30cm$point.ID] #point.ID is the row number
kssl_ssurgo_30cm$labsampnum <- kssl_shp$pdlbsmp[kssl_ssurgo_30cm$point.ID]
# colnames(kssl_ssurgo_30cm)

kssl_horizons <- read.csv(file.path(ksslDir, 'ca_kssl_horizons.csv'), stringsAsFactors = FALSE)
sum(is.na(kssl_horizons$oc) & !is.na(kssl_horizons$c_tot))
sum(is.na(kssl_horizons$oc) & !is.na(kssl_horizons$c_tot) & kssl_horizons$ph_h2o < 7.8, na.rm = TRUE)

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}

kgOrgC_sum_v2 <- function(x, depth, rm.NAs=FALSE, crit_pH=7.8) {
  thick <- x$hzn_bot - x$hzn_top
  if(any(is.na(x$oc)) & all(!is.na(x$ph_h2o))) {
    if (all(x$ph_h2o < crit_pH)) {
      sum((thick / 10) * x$c_tot * x$db_13b * (1 - x$frags / 100), na.rm = rm.NAs)
    } else {NA}
  } else {sum((thick / 10) * x$oc * x$db_13b * (1 - x$frags / 100), na.rm = rm.NAs)}
}

awc_sum_v2 <- function(x, rm.NAs=FALSE) {
  thick <- x$hzn_bot - x$hzn_top
  sum(thick * x$whc, na.rm = rm.NAs)
}
# horizon_SPC <- kssl_horizons_subset[112:113]
# depth <- 30
# vars_of_interest <- c('clay')
# varnames <- 'clay'
horizon_to_comp_v2 <- function(horizon_SPC, depth, vars_of_interest = c('clay', 'silt', 'sand', 'oc', 'c_tot', 'estimated_om', 'cec7', 'cec82', 'db_13b', 'db_od', 'frags', 'ec_12pre', 'ph_h2o', 'ph_cacl2', 'sar', 'caco3', 'gypl20', 'COLEws', 'whc', 'w3cld', 'w15l2', 'Ks'), varnames = c('clay', 'silt', 'sand', 'oc', 'c_tot', 'om', 'cec_7', 'cec_8.2', 'bd_13b', 'bd_od', 'frags', 'ec', 'pH_H2O', 'pH_CaCl2', 'sar', 'caco3', 'gyp', 'lep', 'awc', 'w3cld', 'w15l2', 'ksat')) { #lep is linear extensibility
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
  horizons(sliced_SPC) <- horizons(sliced_SPC)[order(as.integer(sliced_SPC$pedon_key), sliced_SPC$hzn_top),] #because integer codes are coerced to character by slice with the sliced data.frame re-ordered contrary to order of site-level data 
  stopifnot(unique(sliced_SPC$pedon_key)==site(sliced_SPC)$pedon_key)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean_v2, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum_v2)
  s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum_v2, rm.NAs=FALSE) #was TRUE for SSURGO
  #columnames <- c(columnames, paste0('awc_', depth, 'cm')) 
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  s
}
#Error: bad horizonation in IDs: 2676, 7183, 7324, 7330, 12668, 14777, 18190, 19977, 20883, 20886, 20887, 72040, 17783 
kssl_horizons[kssl_horizons$pedon_key==2676,] #has two complete profiles: exclude 40A22466, 40A22467, 40A22468
kssl_horizons[kssl_horizons$pedon_key==7183,] #two extra horizons 79P00471, 79P00472
kssl_horizons[kssl_horizons$pedon_key==7324,] #one extra horizon 79P01177
kssl_horizons[kssl_horizons$pedon_key==7330,] #four extra horizons: 79P01239, 79P01240, 79P01241, 79P01242
kssl_horizons[kssl_horizons$pedon_key==12668,] #one extra horizon: 85P05356
kssl_horizons[kssl_horizons$pedon_key==14777,] #one extra horizon: 88P01488
kssl_horizons[kssl_horizons$pedon_key==18190,] #extra horizons: 91P04438, 91P04439, 91P04440, 91P04441, 91P04442, 91P04443, 91P04444, 91P04445
kssl_horizons[kssl_horizons$pedon_key==19977,] #93P02017, 93P02018
kssl_horizons[kssl_horizons$pedon_key==20883,] #94P02018, 94P02019, 94P02020
kssl_horizons[kssl_horizons$pedon_key==20886,] #94P02028
kssl_horizons[kssl_horizons$pedon_key==20887,] #94P02038
kssl_horizons[kssl_horizons$pedon_key==72040,] #15N01738, 15N01739
kssl_horizons[kssl_horizons$pedon_key==17783,] #all problematic 91P02043, 91P02044, 91P02045
kssl_horizons[kssl_horizons$pedon_key==14162,]#87P02575,87P02577,87P02578,87P02583
kssl_horizons[kssl_horizons$pedon_key==20635,] #94P00526, 94P00527, 94P00528, 94P00529 do not make sense in terms of depths
kssl_horizons[kssl_horizons$pedon_key==10150,] #83P01308
kssl_horizons[kssl_horizons$pedon_key==14771,] #88P01429
kssl_horizons[kssl_horizons$pedon_key==2747,] #40A23004 but has gap between 91-110 cm
kssl_horizons[kssl_horizons$pedon_key==19974,] #93P01981, 93P01982
kssl_horizons[kssl_horizons$pedon_key==20636,] #missing 0-48 cm so remove entire pedon below
kssl_horizons[kssl_horizons$pedon_key==14773,] #88P01449
horizons_to_exclude <- c('91P02043', '91P02044', '91P02045', '15N01738', '15N01739', '94P02038', '94P02028', '94P02018', '94P02019', '94P02020', '93P02017', '93P02018',  '91P04438', '91P04439', '91P04440', '91P04441', '91P04442', '91P04443', '91P04444', '91P04445', '88P01488', '85P05356', '79P01239', '79P01240', '79P01241', '79P01242', '79P01177', '79P00471', '79P00472', '40A22466', '40A22467', '40A22468', '87P02575', '87P02577', '87P02578', '87P02583', '94P00526', '94P00527', '94P00528', '94P00529', '83P01308', '88P01429', '40A23004', '93P01981', '93P01982', '88P01449')
kssl_horizons_subset <- kssl_horizons[kssl_horizons$pedon_key %in% kssl_ssurgo_30cm$pedon_key, ]
# dim(kssl_horizons_subset)
kssl_horizons_subset <- kssl_horizons_subset[!(kssl_horizons_subset$labsampnum %in% horizons_to_exclude), ]
# dim(kssl_horizons_subset)
kssl_horizons_subset <- kssl_horizons_subset[kssl_horizons_subset$pedon_key!=20636, ]

depths(kssl_horizons_subset) <- pedon_key ~ hzn_top + hzn_bot
# class(kssl_horizons_subset)
# horizons(kssl_horizons_subset)[kssl_horizons_subset$pedon_key==73286,]
kssl_horizons_subset$soil_depth <- profileApply(kssl_horizons_subset, FUN = estimateSoilDepth, name='hzn_desgn', top='hzn_top', bottom='hzn_bot')
# which(kssl_horizons_subset$soil_depth==0)
# which(site(kssl_horizons_subset)==73286)
# kssl_horizons_subset[113]
# test_hz_aggregation_v1 <- horizon_to_comp_v2(horizon_SPC = kssl_horizons_subset[1:113], depth = 30)
# head(test_hz_aggregation_v1)
# kssl_horizons[kssl_horizons$pedon_key %in% c(165,166,204) & kssl_horizons$hzn_top <=30,]
kssl_points_30cm <- horizon_to_comp_v2(horizon_SPC = kssl_horizons_subset, depth = 30)
colnames(kssl_points_30cm)

kssl_points_30cm[kssl_points_30cm$pedon_key==165,]
kssl_horizons[kssl_horizons$pedon_key==165 & kssl_horizons$hzn_top <=30,]
#awc check: 3.3 cm
15*0.15+15*0.07
#clay check: 16.9%
12.5*0.5+21.3*0.5
#oc check:0.77%
1.13*0.5+0.41*0.5
#kg OC check
(15/10)*1.13*1.62+(15/10)*0.41*1.71
kssl_points_30cm[kssl_points_30cm$pedon_key==73286,]
kssl_horizons[kssl_horizons$pedon_key==73286 & kssl_horizons$hzn_top <=30,]
kssl_points_30cm[185,]
kssl_horizons[kssl_horizons$pedon_key==14771 & kssl_horizons$hzn_top <=30,]
#clay:11.943%
21.5*4/30+16.3*11/30+6.7*10/30+5.2*5/30
#oc:0.656%
1.87*4/30+0.71*11/30+0.37*10/30+0.14*5/30
sum(kssl_points_30cm$soil_depth==0) #12
kssl_points_30cm[kssl_points_30cm$pedon_key==7327, ]

kssl_points_30cm <- kssl_points_30cm[!(kssl_points_30cm$soil_depth==0), ]
kssl_points_30cm$cluster_2 <- kssl_ssurgo_30cm$cluster_2[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_3 <- kssl_ssurgo_30cm$cluster_3[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_4 <- kssl_ssurgo_30cm$cluster_4[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_5 <- kssl_ssurgo_30cm$cluster_5[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_6 <- kssl_ssurgo_30cm$cluster_6[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_7 <- kssl_ssurgo_30cm$cluster_7[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_8 <- kssl_ssurgo_30cm$cluster_8[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_9 <- kssl_ssurgo_30cm$cluster_9[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_10 <- kssl_ssurgo_30cm$cluster_10[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_11 <- kssl_ssurgo_30cm$cluster_11[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_30cm$cluster_12 <- kssl_ssurgo_30cm$cluster_12[match(kssl_points_30cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
sum(is.na(kssl_points_30cm$cluster_7))
dim(kssl_points_30cm) #369 kssl points
summary(kssl_points_30cm$kgOrg.m2_30cm)
sum(!is.na(kssl_points_30cm$oc_30cm))
sum(is.na(kssl_points_30cm$oc_30cm) & !is.na(kssl_points_30cm$c_tot_30cm))
sum(is.na(kssl_points_30cm$oc_30cm) & !is.na(kssl_points_30cm$c_tot_30cm) & kssl_points_30cm$pH_H2O_30cm <= 8.2, na.rm = TRUE) #39 points have total C but not organic C; 31 have pH <= 7.5; 32 have pH <= 7.8; 34 have pH <= 8.0

write.csv(kssl_points_30cm, file.path(ksslDir,  'kssl_cluster_30cm_FINAL.csv'), row.names = FALSE)

tapply(kssl_points_30cm$oc_30cm, kssl_points_30cm$cluster_9, summary)
tapply(kssl_points_30cm$pH_H2O_30cm, kssl_points_30cm$cluster_9, summary)
tapply(kssl_points_30cm$cec_7_30cm, kssl_points_30cm$cluster_9, summary)
sum(!is.na(kssl_points_30cm$kgOrg.m2_30cm)) #120 of 370 have 0-30 cm content data (was 106)
write.csv(kssl_ssurgo_30cm, file.path(ksslDir, 'kssl_pts_ssurgo_30cm_extract_FINAL.csv'), row.names = FALSE)


#read in 30 cm data
om_to_oc <- 1.72
crit_pH <- 7.8
clus_7_colors <- c('deepskyblue', 'gold', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1')
order_lgnd_7 <- c(4,5,2,3,7,1,6)
valley30cm_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', 'valley30cm_by_mukey_cluster_FINAL.csv'), stringsAsFactors = FALSE)
dom_order_by_mukey <- read.csv(file.path(dataDir, 'FINAL results', "soil survey facts", 'dom_order_by_mukey.csv'), stringsAsFactors = FALSE)

kssl_ssurgo_extract <- read.csv(file.path(ksslDir, 'kssl_pts_ssurgo_30cm_extract_FINAL.csv'), stringsAsFactors = FALSE)
kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_FINAL.csv'), stringsAsFactors = FALSE) #updated 2/11/20 to consider total C when org C not available but soil pH sufficiently low for calculating SOC content
kssl_points_30cm$xdim_vioplot7 <- match(kssl_points_30cm$cluster_7, order_lgnd_7)
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

df_combined1 <-  kssl_points_30cm[ ,c('clay_30cm', 'om_30cm', 'bd_13b_30cm', 'pH_H2O_30cm', 'kgOrg.m2_30cm', colnames(kssl_points_30cm)[grepl('cluster_', colnames(kssl_points_30cm))], 'dom_order', 'xdim_vioplot7')]
dim(df_combined1) #369 rows
colnames(df_combined1)[3] <- 'bd_30cm'
df_combined1$source <- 'KSSL'

kerri_points_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
kerri_points_30cm$xdim_vioplot7 <- match(kerri_points_30cm$cluster_7, order_lgnd_7)
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
kerri_points_30cm$cc_meta <- kerri_metadata$covercrop.vs.resident.vegetation[match(kerri_points_30cm$Concatenate, kerri_metadata$Concatenate)]
kerri_points_30cm$cc_type <- kerri_metadata$perennial.Covercrop.vs.annual.covercrop[match(kerri_points_30cm$Concatenate, kerri_metadata$Concatenate)]
table(kerri_points_30cm$cc_type)

tapply(kerri_points_30cm$compost_added, kerri_points_30cm$xdim_vioplot7, table)
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$compost_added[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$clay_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$compost_added[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==3] ~ kerri_points_30cm$compost_added[kerri_points_30cm$xdim_vioplot7==3])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==2] ~ kerri_points_30cm$compost_added[kerri_points_30cm$xdim_vioplot7==2])
plot(kerri_points_30cm$clay_30cm[kerri_points_30cm$xdim_vioplot7==2], kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==2])

plot(kerri_points_30cm$bd_30cm[kerri_points_30cm$xdim_vioplot7==2], kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==2])

tapply(kerri_points_30cm$cc_type, kerri_points_30cm$xdim_vioplot7, table)
boxplot(kerri_points_30cm$om_30cm ~ kerri_points_30cm$cc_type)
plot(ifelse(kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==7]=='perennial', 3, ifelse(kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==7]=='annual', 2, 1)), kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==7])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$clay_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==3] ~ kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==3])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==2] ~ kerri_points_30cm$cc_type[kerri_points_30cm$xdim_vioplot7==2])
kerri_points_30cm[kerri_points_30cm$cc_type=='perennial' & !is.na(kerri_points_30cm$cluster_7),]
table(kerri_points_30cm$compost_added[kerri_points_30cm$cc_type=='perennial' & kerri_points_30cm$xdim_vioplot7==2])

#till vs. no-till
tapply(kerri_points_30cm$tillage, kerri_points_30cm$xdim_vioplot7, table)
boxplot(kerri_points_30cm$om_30cm ~ kerri_points_30cm$tillage)
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$tillage[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$clay_30cm[kerri_points_30cm$xdim_vioplot7==4] ~ kerri_points_30cm$tillage[kerri_points_30cm$xdim_vioplot7==4])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==3] ~ kerri_points_30cm$tillage[kerri_points_30cm$xdim_vioplot7==3])
boxplot(kerri_points_30cm$om_30cm[kerri_points_30cm$xdim_vioplot7==2] ~ kerri_points_30cm$tillage[kerri_points_30cm$xdim_vioplot7==2])

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
kerri_points_30cm$SHM_rating <- sum(kerri_points_30cm$compost_added=='yes', kerri_points_30cm$tillage=='no till', kerri_points_30cm$irrigated_vs_dryfarm=='dryfarm',  
# kerri_points_30cm$shr7_name <- cluster_

df_combined2 <- kerri_points_30cm[!is.na(kerri_points_30cm$cluster_7), c('clay_30cm', 'om_30cm', 'bd_30cm', 'pH_H2O_30cm', 'kgOrg.m2_30cm', colnames(kerri_points_30cm)[grepl('cluster_', colnames(kerri_points_30cm))], 'dom_order', 'xdim_vioplot7')]
dim(df_combined2) #103 rows
df_combined2$source <- 'CDFA'

df_combined <- rbind(df_combined1, df_combined2)
dim(df_combined) #472
table(df_combined$dom_order[!is.na(df_combined$om_30cm)])
df_combined$om_30cm[df_combined$dom_order=='Ultisols']
lapply(df_combined, function(x) sum(!is.na(x)))

#create kssl point vioplots
vioplot_mod_clus7_KSSL_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, legend_plot, legendloc, legend_cex, sig_labels, fig_label, legend_text, fig_height, removeOutlier, OutlierThreshold) {
  plot_order2 <- (1:7)[plot_order]
  if (removeOutlier) {
    df <- df[which(df[[varname]] < OutlierThreshold), ]
  } else{
    df <- df[!is.na(df[[varname]]), ]
  }
  tiff(file = file.path(FiguresDir, 'FINAL', 'CalAg validation plots', fname), family = 'Times New Roman', width = 3.25, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[[varname]][df$cluster_7==plot_order2[1]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[2]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[3]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[4]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[5]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[6]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[7]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, col=clus_7_colors[plot_order], wex=1.2, cex=0.8, rectCol = 'gray', ylim = ylim_vioplot, ylab = NULL)
  mtext('Soil health region', side = 1, line = 2)
  mtext(ylab, side = 2, line = 2)
  points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.15, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.5, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.15, y=kssl_means$mean, pch=8, cex=0.5, col='orange')
  text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5, cex=0.9)
  legend('topright', fig_label, bty='n', inset=0.005)
  if(legend_plot) {
    legend(x=legendloc, legend=legend_text, pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), pt.cex=legend_cex, bty = 'n')
  }
  dev.off()
}
vioplot_mod_clus7_KSSL_validation(kssl_points_30cm, 'clay_30cm', ylim_vioplot = c(-2.5,80), plot_order = order_lgnd_7, area_fact = 10, ylab='Clay (%)', fname='class7_clay_vioplots_KSSL_pts_only_validation.tif', mar=c(0.02, 3.25, 0.25, 0.25), kssl_df = kssl_points_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = 'clay_30cm', sig_labels = c('A', 'B', 'A', 'AB', 'A', 'C', 'D'), fig_label = 'a', fig_height = 3, legend_text =  c('KSSL point data violin plots', 'KSSL point data', 'KSSL mean'))

vioplot_mod_clus7_KSSL_validation(kssl_points_30cm, 'kgOrg.m2_30cm', ylim_vioplot = c(-0.15,12), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('0-30 cm SOC (kg m'^-2*')'), fname='class7_kgSOC_0_30cm_vioplots_KSSL_pts_only_validation.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "kgOrg.m2_30cm", sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'b', fig_height = 3, removeOutlier = FALSE)

TukeyHSD(aov(kssl_points_30cm$oc_30cm[which(kssl_points_30cm$oc_30cm < 3.5)] * 1.72 ~ kssl_points_30cm$SHR7name[which(kssl_points_30cm$oc_30cm < 3.5)]), ordered=FALSE)
sum(TukeyHSD(aov(kssl_points_30cm$oc_30cm[which(kssl_points_30cm$oc_30cm < 3.5)] * 1.72 ~ kssl_points_30cm$SHR7name[which(kssl_points_30cm$oc_30cm < 3.5)]))[[1]][,4] <= 0.05) #9
vioplot_mod_clus7_KSSL_validation(kssl_points_30cm, 'oc_30cm', ylim_vioplot = c(0, 3.5), plot_order = order_lgnd_7, area_fact = 10, ylab='Organic matter (%)', fname='class7_OM_vioplots_KSSL_pts_only_validation_0_4.tif', mar=c(3.5, 4.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "oc_30cm", sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'b', fig_height = 3, removeOutlier=TRUE, OutlierThreshold=3.5)

#add Napa-Lodi data to vioplots
vioplot_mod_clus7_validation <- function(df, varname, ylim_vioplot, plot_order, area_fact, labnames, ylab, fname, mar, kssl_df, kssl_varname, cdfa_pts, cdfa_varname, legend_plot, legendloc, legend_cex, sig_labels, fig_label, legend_text, fig_height, removeOutlier, OutlierThreshold) {
  plot_order2 <- (1:7)[plot_order]
  if (removeOutlier) {
    df <- df[which(df[[varname]] < OutlierThreshold), ]
  } else{
    df <- df[!is.na(df[[varname]]), ]
  }
  tiff(file = file.path(FiguresDir, 'FINAL', 'CalAg validation plots', fname), family = 'Times New Roman', width = 3.25, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[[varname]][df$cluster_7==plot_order2[1]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[2]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[3]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[4]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[5]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[6]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, df[[varname]][df$cluster_7==plot_order2[7]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, col=clus_7_colors[plot_order], wex=1.2, cex=0.8, rectCol = 'gray', ylim = ylim_vioplot, ylab = NULL)
  mtext('Soil health region', side = 1, line = 2)
  mtext(ylab, side = 2, line = 2)
  points(x=kssl_df$xdim_vioplot7[!is.na(kssl_df[[kssl_varname]])]-0.15, y=kssl_df[[kssl_varname]][!is.na(kssl_df[[kssl_varname]])]*if(kssl_varname=='totC_30cm'){1.72} else{1}, cex=0.5, pch=1, col='black')
  kssl_means <- data.frame(mean=tapply(kssl_df[[kssl_varname]]*if(kssl_varname=='totC_30cm'){1.72} else{1}, kssl_df$cluster_7, mean, na.rm=TRUE))
  kssl_means$xdim_vioplot <- match(row.names(kssl_means), plot_order)
  points(x=kssl_means$xdim_vioplot-0.15, y=kssl_means$mean, pch=8, cex=0.5, col='orange')
  points(x=cdfa_pts$xdim_vioplot7[!is.na(cdfa_pts[[cdfa_varname]])]+0.15, y=cdfa_pts[[cdfa_varname]][!is.na(cdfa_pts[[cdfa_varname]])], cex=0.5, pch=4, col='black')
  cdfa_means <- data.frame(mean=tapply(cdfa_pts[[cdfa_varname]], cdfa_pts$cluster_7, mean, na.rm=TRUE))
  cdfa_means$xdim_vioplot <- match(row.names(cdfa_means), plot_order)
  points(x=cdfa_means$xdim_vioplot+0.15, y=cdfa_means$mean, pch=8, cex=0.5, col='darkblue')
  text(x=1:7, y=ylim_vioplot[1], labels = sig_labels, adj=0.5, cex=0.9)
  legend('topright', fig_label, bty='n', inset=0.005)
  if(legend_plot) {
    legend(x=legendloc, legend=c('Validation point data violin plots', 'KSSL points', 'KSSL average', 'Napa-Lodi points', 'Napa-Lodi average'), pch=c(NA,1,8,4,8), col = c(NA, 'black', 'orange', 'black', 'darkblue'), pt.cex=legend_cex, bty='n')
  }
  dev.off()
}
# vioplot_mod_clus7_validation(df_combined, 'clay_30cm', ylim_vioplot = c(-2.5,80), plot_order = order_lgnd_7, area_fact = 10, ylab='Clay (%)', fname='class7_clay_vioplots_validation_pts_all.tif', mar=c(3, 3, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=TRUE, legendloc='topleft', legend_cex = 0.5, legend_text='', kssl_varname = 'clay_30cm', cdfa_varname = 'clay_30cm', sig_labels = c('A', 'B', 'A', 'AB', 'A', 'C', 'D'), fig_label = 'a', removeOutlier = FALSE, fig_height = 3)

vioplot_mod_clus7_validation(df_combined, 'om_30cm', ylim_vioplot = c(0,4.4), plot_order = order_lgnd_7, area_fact = 10, ylab='Organic matter (%)', fname='class7_om_vioplots_validation_pts_all.tif', mar=c(3, 3, 1, 1), kssl_df = kssl_points_30cm, cdfa_pts=kerri_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.5, legend_text='', kssl_varname = 'om_30cm', cdfa_varname = 'om_30cm', sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'a', removeOutlier = TRUE, OutlierThreshold = 4.5, fig_height = 3.5)
hist(kssl_points_30cm$om_30cm)

vioplot_mod_clus7_validation(df_combined, 'kgOrg.m2_30cm', ylim_vioplot = c(-0.15,10.1), plot_order = order_lgnd_7, area_fact = 10, ylab=expression('0-30 cm SOC (kg m'^-2*')'), fname='class7_kgSOC_0_30cm_vioplots_validation_pts_all.tif', mar=c(3, 3.25, 1, 1), kssl_df = kssl_points_30cm, legend_plot=FALSE, legendloc='topleft', legend_cex = 0.9 , kssl_varname = "kgOrg.m2_30cm", cdfa_pts=kerri_points_30cm, cdfa_varname = 'kgOrg.m2_30cm', sig_labels = c('A', 'C', 'A', 'C', 'A', 'AB', 'BC'), fig_label = 'b', fig_height = 3.5, removeOutlier = FALSE, OutlierThreshold = 10)

#create KSSL shapefile that are within soil health regions
sum(!is.na(kssl_points_30cm$cluster_7))
kssl_points_30cm$pedon_key
kssl_SHR_shp <- kssl_shp
kssl_SHR_shp <- kssl_SHR_shp[kssl_SHR_shp$pedn_ky %in% kssl_points_30cm$pedon_key, ]
kssl_SHR_shp$SHR7name  <- kssl_points_30cm$SHR7name[match(kssl_SHR_shp$pedn_ky, kssl_points_30cm$pedon_key)]
# shapefile(kssl_SHR_shp, file.path(ksslDir, 'shapefiles', 'kssl_SHR7.shp'))
table(kssl_SHR_shp$SHR7name)
kssl_SHR_UCD_shp <- kssl_SHR_shp
kssl_SHR_UCD_shp <- kssl_SHR_UCD_shp[grepl('UCD', kssl_SHR_UCD_shp$pdlbsmp), ] #54
# shapefile(kssl_SHR_UCD_shp, file.path(ksslDir, 'shapefiles', 'kssl_UCD_SHR7.shp'))
table(kssl_SHR_UCD_shp$SHR7name)

#create 1 m soil property summary
kssl_points_100cm <- horizon_to_comp_v2(horizon_SPC = kssl_horizons_subset, depth = 100)
sum(kssl_points_100cm$soil_depth==0)
kssl_points_100cm <- kssl_points_100cm[!(kssl_points_100cm$soil_depth==0), ]
kssl_points_100cm$cluster_2 <- kssl_ssurgo_30cm$cluster_2[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_3 <- kssl_ssurgo_30cm$cluster_3[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_4 <- kssl_ssurgo_30cm$cluster_4[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_5 <- kssl_ssurgo_30cm$cluster_5[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_6 <- kssl_ssurgo_30cm$cluster_6[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_7 <- kssl_ssurgo_30cm$cluster_7[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_8 <- kssl_ssurgo_30cm$cluster_8[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_9 <- kssl_ssurgo_30cm$cluster_9[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_10 <- kssl_ssurgo_30cm$cluster_10[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_11 <- kssl_ssurgo_30cm$cluster_11[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
kssl_points_100cm$cluster_12 <- kssl_ssurgo_30cm$cluster_12[match(kssl_points_100cm$pedon_key, kssl_ssurgo_30cm$pedon_key)]
sum(is.na(kssl_points_100cm$cluster_7))
sum(!is.na(kssl_points_100cm$kgOrg.m2_100cm)) #only 84 of 369 have complete content data (was 79 before)
write.csv(kssl_points_100cm, file.path(ksslDir, 'kssl_cluster_100cm_NArm_v2.csv'), row.names = FALSE)

kssl_points_100cm <- read.csv(file.path(ksslDir, 'kssl_cluster_100cm_NArm_v2.csv'), stringsAsFactors = FALSE)

table(kssl_points_100cm$cluster_7)

clus_5_names <- c()
clus_6_names <- c()
clus_7_names <- c('3. Coarse w/pans', '6. Fine saline-sodic', '5. Coarse saline-sodic', '1. Coarse w/no restrictions', '7. Fine shrink-swell', '2. Loamy w/no restrictions', '4. Loamy w/pans')



compare_region_means <- function(df, names, cluster_no, y) {
  x <- names[match(df[[paste0('cluster_', cluster_no)]], 1:cluster_no)]
  y <- df[[y]]
  print(summary(aov(y ~ x)))
  print(TukeyHSD(aov(y ~ x), ordered=TRUE))
}
compare_region_means(kssl_points_30cm, clus_7_names, 7, 'clay_30cm')
compare_region_means(kssl_points_30cm, clus_7_names, 7, 'om_30cm')
compare_region_means(kssl_points_30cm, clus_7_names, 7, 'lep_30cm')
compare_region_means(kssl_points_30cm, clus_7_names, 7, 'lep_30cm')
#analysis of variance by cluster
summary(aov(clay_30cm ~ as.factor(cluster_2), kssl_points_30cm))
summary(aov(clay_30cm ~ as.factor(cluster_3), kssl_points_30cm))
summary(aov(clay_30cm ~ as.factor(cluster_4), kssl_points_30cm))
summary(aov(clay_30cm ~ as.factor(cluster_5), kssl_points_30cm))
TukeyHSD(aov(clay_30cm ~ as.factor(cluster_5), kssl_points_30cm), ordered=TRUE)
summary(aov(clay_30cm ~ cluster_6, kssl_points_30cm))
summary(aov(clay_30cm ~ cluster_7, kssl_points_30cm))
TukeyHSD(aov(clay_30cm ~ as.factor(cluster_7), kssl_points_30cm), ordered=TRUE)
summary(aov(clay_30cm ~ cluster_8, kssl_points_30cm))
summary(aov(clay_30cm ~ cluster_9, kssl_points_30cm))



order_lgnd_9 <- c(1,3,6,5,9,7,8,2)
vioplot_kssl_clus9 <- function(df, varname, ylim_vioplot, plot_order, labnames, ylab, fname, mar, group_names) {
  plot_order2 <- (1:9)[plot_order]
  tiff(file = file.path(FiguresDir, 'v2', 'kssl', fname), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[[varname]][df$cluster_9==plot_order2[1]], df[[varname]][df$cluster_9==plot_order2[2]], df[[varname]][df$cluster_9==plot_order2[3]], df[[varname]][df$cluster_9==plot_order2[4]], df[[varname]][df$cluster_9==plot_order2[5]], df[[varname]][df$cluster_9==plot_order2[6]], df[[varname]][df$cluster_9==plot_order2[7]], df[[varname]][df$cluster_9==plot_order2[8]], col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[plot_order], rectCol = 'gray', ylim = ylim_vioplot, ylab = ylab, names = group_names)
  mtext('Soil health region', side = 1, line = 2.25)
  dev.off()
}
colnames(kssl_points_30cm)
#order_lgnd was defined for radarchart
vioplot_kssl_clus9(kssl_points_30cm, 'clay_30cm', ylim_vioplot = c(0.5,75), plot_order = order_lgnd_9, ylab='Clay (%)', fname='class9_clay_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
vioplot_kssl_clus9(kssl_points_30cm, 'om_30cm', ylim_vioplot = c(0.1,9), plot_order = order_lgnd_9,  ylab='Organic matter (%)', fname='class9_OM_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
kssl_points_30cm$logom_30cm <- log(kssl_points_30cm$om_30cm)
vioplot_kssl_clus9(kssl_points_30cm, 'logom_30cm', ylim_vioplot = c(-2,2),  plot_order = order_lgnd_9,  ylab=expression('Organic matter (Log'[10]~'%)'), fname='class9_LogOM_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
vioplot_kssl_clus9(kssl_points_30cm, 'cec_7_30cm', ylim_vioplot = c(0.5,60), plot_order = order_lgnd_9,  ylab=expression('Cation exchange capacity (mEq 100g'^-1*')'), fname='class9_CEC_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
vioplot_kssl_clus9(kssl_points_30cm, 'pH_H2O_30cm', ylim_vioplot = c(4.9,10.1), plot_order = order_lgnd_9,  ylab='soil pH', fname='class9_pH_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
vioplot_kssl_clus9(kssl_points_30cm, 'awc_30cm', ylim_vioplot = c(0.55,7.1), plot_order = order_lgnd_9,  ylab=expression('Available water capacity (cm H'[2]*'O 30 cm'^-1~'soil)'), fname='class9_AWC_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
summary(kssl_points_30cm$ec_30cm)
vioplot_kssl_clus9(kssl_points_30cm, 'ec_30cm', ylim_vioplot = c(0,20.5), plot_order = order_lgnd_9,  ylab=expression('Electrical conductivity (mmhos cm'^-1*')'), fname='class9_EC_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
summary(kssl_points_30cm$bd_13b_30cm)
tapply(kssl_points_30cm$bd_13b_30cm, kssl_points_30cm$cluster_9, summary)
vioplot_kssl_clus9(kssl_points_30cm, 'bd_13b_30cm', ylim_vioplot = c(1.0,1.9), plot_order = order_lgnd_9,  ylab=expression('Bulk density (g soil cm'^-3*'soil)'), fname='class9_BD_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
vioplot_kssl_clus9(kssl_points_30cm, 'ksat_30cm', ylim_vioplot = c(0.2,2.8), plot_order = order_lgnd_9, ylab=expression('Hydraulic conductivity (units unknown)'), fname='class9_ksat_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))
#kssl_points_30cm$logks_30cm <- log(kssl_points_30cm$ksat_30cm)
#vioplot_kssl_clus9(kssl_points_30cm, 'logks_30cm', ylim_vioplot = c(-4.58,5.85), plot_order = order_lgnd_9,  ylab=expression('Saturated conductivity (log'[10]~mu*'m H'[2]*'O s'^-1*')'), fname='class9_logKs_vioplots.tif', mar=c(3.5, 4.25, 1, 1))
summary(kssl_points_30cm$lep_30cm)
kssl_points_30cm$lep_30cm <- kssl_points_30cm$lep_30cm*100
tapply(kssl_points_30cm$lep_30cm, kssl_points_30cm$cluster_9, function(x) sum(!is.na(x)))
vioplot_kssl_clus9(kssl_points_30cm, 'lep_30cm', ylim_vioplot = c(-0.1,17), plot_order = order_lgnd_9,  ylab='Linear extensibility (%)', fname='class9_lep_vioplots.tif', mar=c(3.5, 4.25, 1, 1), group_names = c(1:5,7:9))

order_lgnd_9 <- c(1,3,6,5,9,4,7,8,2)
tiff(file = file.path(FiguresDir, 'v2', 'kssl', 'kssl_samples_cluster9.tif'), pointsize = 11, family = 'Times New Roman', width = 6.5, height = 5, units = 'in', res=800, compression = 'lzw')
par(mar=c(8, 4, 1, 1))
barplot(as.numeric(table(kssl_points_30cm$cluster_9)[order_lgnd_9]), col=c('lightgoldenrod', 'violetred', 'tan2', 'black', 'gold', 'tan4', 'lightblue1', 'deepskyblue', 'firebrick3')[order_lgnd_9], ylab = 'Number of KSSL points', legend.text=c('1. Sandy soils', '9. Shrink-swell clays', '2. Loams w/ no res. & mod OM-low SS', '6. Loams w/ res. & high OM', '4. Loams w/ res. & low OM', '3. Loams w/ no res. & mod OM-mod SS', '7. Saline-sodic loams', '8. Saline-sodic clays', '5. Loams w/ res. & mod OM')[order_lgnd_9], cex.axis = 1, cex.names = 1, cex.lab = 1, args.legend = list(x=11, y=-5, cex=1, ncol=2))
dev.off()
