#TO-DO: run soil horizon aggregation by excluding any 0-30 cm variables with missing data for a horizon
library(raster)
library(aqp)
library(vioplot)
library(extrafont)
library(extrafontdb)
loadfonts(device = 'win')
laptop <- FALSE

if (laptop) {
  mainDir <- 'C:/Users/smdevine/Desktop/post doc'
} else { #on UCD desktop
  mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
}
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
}

list.files(ksslDir)
kssl_shp <- shapefile(file.path(ksslDir, 'shapefiles', 'kssl_CA.shp'))
crs(kssl_shp)
valley_mu_shp_30cm <- shapefile(file.path(dataDir, 'shapefiles with data', 'valley_30cm_cluster.shp'))
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
kssl_ssurgo_30cm <- kssl_ssurgo_30cm[!is.na(kssl_ssurgo_30cm$cluster_9), ]
kssl_ssurgo_30cm$pedon_key <- kssl_shp$pedn_ky[kssl_ssurgo_30cm$point.ID] #point.ID is the row number
kssl_ssurgo_30cm$labsampnum <- kssl_shp$pdlbsmp[kssl_ssurgo_30cm$point.ID]
# colnames(kssl_ssurgo_30cm)

kssl_horizons <- read.csv(file.path(ksslDir, 'ca_kssl_horizons.csv'), stringsAsFactors = FALSE)

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}

awc_sum_v2 <- function(x, rm.NAs=FALSE) {
  thick <- x$hzn_bot - x$hzn_top
  sum(thick * x$whc, na.rm = rm.NAs)
}
# horizon_SPC <- kssl_horizons_subset[112:113]
# depth <- 30
# vars_of_interest <- c('clay')
# varnames <- 'clay'
horizon_to_comp_v2 <- function(horizon_SPC, depth, vars_of_interest = c('clay', 'silt', 'sand', 'oc', 'estimated_om', 'cec7', 'cec82', 'db_13b', 'db_od', 'frags', 'ec_12pre', 'ph_h2o', 'ph_cacl2', 'sar', 'caco3', 'gypl20', 'COLEws', 'whc', 'w3cld', 'w15l2', 'Ks'), varnames = c('clay', 'silt', 'sand', 'oc', 'om', 'cec_7', 'cec_8.2', 'bd_13b', 'bd_od', 'frags', 'ec', 'pH_H2O', 'pH_CaCl2', 'sar', 'caco3', 'gyp', 'lep', 'awc', 'w3cld', 'w15l2', 'ksat')) { #lep is linear extensibility
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
  #s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
  s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum_v2, rm.NAs=FALSE) #was TRUE for SSURGO
  columnames <- c(columnames, paste0('awc_', depth, 'cm')) 
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
horizons_to_exclude <- c('91P02043', '91P02044', '91P02045', '15N01738', '15N01739', '94P02038', '94P02028', '94P02018', '94P02019', '94P02020', '93P02017', '93P02018',  '91P04438', '91P04439', '91P04440', '91P04441', '91P04442', '91P04443', '91P04444', '91P04445', '88P01488', '85P05356', '79P01239', '79P01240', '79P01241', '79P01242', '79P01177', '79P00471', '79P00472', '40A22466', '40A22467', '40A22468')
kssl_horizons_subset <- kssl_horizons[kssl_horizons$pedon_key %in% kssl_ssurgo_30cm$pedon_key, ]
# dim(kssl_horizons_subset)
kssl_horizons_subset <- kssl_horizons_subset[!(kssl_horizons_subset$labsampnum %in% horizons_to_exclude), ]
# dim(kssl_horizons_subset)
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
dim(kssl_points_30cm)
head(kssl_points_30cm)
kssl_points_30cm[kssl_points_30cm$pedon_key==165,]
kssl_horizons[kssl_horizons$pedon_key==165 & kssl_horizons$hzn_top <=30,]
#awc check: 3.3 cm
15*0.15+15*0.07
#clay check: 16.9%
12.5*0.5+21.3*0.5
#oc check:0.77%
1.13*0.5+0.41*0.5
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

write.csv(kssl_points_30cm, file.path(ksslDir, 'kssl_cluster_30cm_NArm.csv'), row.names = FALSE)

tapply(kssl_points_30cm$oc_30cm, kssl_points_30cm$cluster_9, summary)
tapply(kssl_points_30cm$pH_H2O_30cm, kssl_points_30cm$cluster_9, summary)
tapply(kssl_points_30cm$cec_7_30cm, kssl_points_30cm$cluster_9, summary)
write.csv(kssl_ssurgo_30cm, file.path(ksslDir, 'kssl_pts_ssurgo_30cm_extract.csv'), row.names = FALSE)

kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_NArm.csv'), stringsAsFactors = FALSE)

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