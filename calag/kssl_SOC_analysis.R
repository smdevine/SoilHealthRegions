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
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
}
clus_7_colors <- c('deepskyblue', 'gold', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1')
order_lgnd_7 <- c(4,5,2,3,7,1,6)
clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')
list.files(ksslDir)
kssl_horizons_SHR <- read.csv(file.path(ksslDir, 'kssl_horizons_SHRonly.csv'), stringsAsFactors = FALSE)
kssl_climate_SHR <- read.csv(file.path(ksslDir, 'kssl_climate_SHRpts.csv'), stringsAsFactors = FALSE)
sum(kssl_horizons_SHR$hzn_mid > 100)
kssl_horizons_SHR$SHR7color <- clus_7_colors[kssl_horizons_SHR$SHR7code]
sum(kssl_horizons_SHR$hzn_mid >= 200)
kssl_horizons_SHR <- kssl_horizons_SHR[kssl_horizons_SHR$hzn_mid < 200, ]
sum(kssl_horizons_SHR$oc_est==0, na.rm = TRUE)
kssl_horizons_SHR <- kssl_horizons_SHR[-which(kssl_horizons_SHR$oc_est==0),]
sum(kssl_horizons_SHR$hzn_mid==0)
kssl_horizons_SHR <- kssl_horizons_SHR[kssl_horizons_SHR$hzn_mid!=0,]
sum(kssl_horizons_SHR$oc_est > 4, na.rm = TRUE)
kssl_horizons_SHR <- kssl_horizons_SHR[which(kssl_horizons_SHR$oc_est < 4), ]


plot(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est)], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est)], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est)], xlim=c(0,4))
plot(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est)], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est)], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est)], ylim=c(-100,0), xlim=c(0,4))
plot(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="1. Coarse & no restrictions"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="1. Coarse & no restrictions"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="1. Coarse & no restrictions"], xlim=c(0,3), ylim=c(0,-50))
points(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="2. Loamy & no restrictions"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="2. Loamy & no restrictions"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="2. Loamy & no restrictions"], xlim=c(0,3), ylim=c(0,-50))
points(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="7. Shrink-swell"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="7. Shrink-swell"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="7. Shrink-swell"], xlim=c(0,3), ylim=c(0,-100))
plot(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="3. Coarse-loamy & restrictive layers"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="3. Coarse-loamy & restrictive layers"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="3. Coarse-loamy & restrictive layers"], xlim=c(0,3), ylim=c(0,-100))
points(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="4. Loamy & restrictive layers"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="4. Loamy & restrictive layers"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="4. Loamy & restrictive layers"], xlim=c(0,3), ylim=c(0,-100))
plot(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="5. Coarse-loamy saline-sodic"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="5. Coarse-loamy saline-sodic"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="5. Coarse-loamy saline-sodic"], xlim=c(0,3), ylim=c(0,-100))
points(kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="6. Fine saline-sodic"], -kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="6. Fine saline-sodic"], col=kssl_horizons_SHR$SHR7color[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="6. Fine saline-sodic"], xlim=c(0,3), ylim=c(0,-100))

#exponential fits
#nls(y ~ yf + (y0 - yf) * exp(-alpha * t)
#C(z) = Cb + (C0 - Cb) * exp(-K * z)
shrink_swell_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="7. Shrink-swell"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="7. Shrink-swell"])
head(shrink_swell_data)
shrink_swell_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = shrink_swell_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(shrink_swell_nls)
plot(fitted(shrink_swell_nls), shrink_swell_data$depth)
points(shrink_swell_data$SOC, shrink_swell_data$depth, col='red')

coarse_no_res_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="1. Coarse & no restrictions"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="1. Coarse & no restrictions"])
head(coarse_no_res_data)
coarse_no_res_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = coarse_no_res_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(coarse_no_res_nls)
plot(fitted(coarse_no_res_nls), coarse_no_res_data$depth)
points(coarse_no_res_data$SOC, coarse_no_res_data$depth, col='red')

loamy_no_res_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="2. Loamy & no restrictions"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="2. Loamy & no restrictions"])
head(loamy_no_res_data)
dim(loamy_no_res_data)
sum(loamy_no_res_data$depth > 100)
loamy_no_res_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = loamy_no_res_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(loamy_no_res_nls)
plot(fitted(loamy_no_res_nls), loamy_no_res_data$depth)
points(loamy_no_res_data$SOC, loamy_no_res_data$depth, col='red')

coarse_loamy_res_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="3. Coarse-loamy & restrictive layers"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="3. Coarse-loamy & restrictive layers"])
head(coarse_loamy_res_data)

coarse_loamy_res_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = coarse_loamy_res_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(coarse_loamy_res_nls)
plot(fitted(coarse_loamy_res_nls), coarse_loamy_res_data$depth)
points(coarse_loamy_res_data$SOC, coarse_loamy_res_data$depth, col='red')

loamy_res_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="4. Loamy & restrictive layers"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="4. Loamy & restrictive layers"])
head(loamy_res_data)
loamy_res_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = loamy_res_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(loamy_res_nls)
plot(fitted(loamy_res_nls), loamy_res_data$depth)
points(loamy_res_data$SOC, loamy_res_data$depth, col='red')

fine_saline_sodic_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="6. Fine saline-sodic"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="6. Fine saline-sodic"])
head(fine_saline_sodic_data)
fine_saline_sodic_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = fine_saline_sodic_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(fine_saline_sodic_nls)
plot(fitted(fine_saline_sodic_nls), fine_saline_sodic_data$depth)
points(fine_saline_sodic_data$SOC, fine_saline_sodic_data$depth, col='red')

coarse_loamy_saline_sodic_data <- data.frame(SOC=kssl_horizons_SHR$oc_est[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="5. Coarse-loamy saline-sodic"], depth=kssl_horizons_SHR$hzn_mid[!is.na(kssl_horizons_SHR$oc_est) & kssl_horizons_SHR$SHR7name=="5. Coarse-loamy saline-sodic"])
head(coarse_loamy_saline_sodic_data)
coarse_loamy_saline_sodic_nls <- nls(SOC ~ Cb + (Cs - Cb) * exp(-K * depth), data = coarse_loamy_saline_sodic_data, start = list(Cb = 0.01, Cs = 5, K = 0.05))
summary(coarse_loamy_saline_sodic_nls)
plot(fitted(coarse_loamy_saline_sodic_nls), coarse_loamy_saline_sodic_data$depth)
points(coarse_loamy_saline_sodic_data$SOC, coarse_loamy_saline_sodic_data$depth, col='red')

#compare coarse, loamy, & shrink-swell
plot(fitted(coarse_no_res_nls), coarse_no_res_data$depth, col='lightgoldenrod', ylab='Depth (cm)', xlab='Soil organic carbon fitted (%)', xlim=c(0,3), ylim=c(0,150))
points(fitted(loamy_no_res_nls), loamy_no_res_data$depth, col='tan4')
points(fitted(shrink_swell_nls), shrink_swell_data$depth, col='violetred')
points(fitted(coarse_loamy_res_nls), coarse_loamy_res_data$depth, col='gold')
points(fitted(loamy_res_nls), loamy_res_data$depth, col='firebrick3')
points(fitted(fine_saline_sodic_nls), fine_saline_sodic_data$depth, col='deepskyblue')
points(fitted(coarse_loamy_saline_sodic_nls), coarse_loamy_saline_sodic_data$depth, col='lightblue1')

#kssl 30 cm data
om_to_oc <- 1.72
crit_pH <- 7.8
list.files(ksslDir)
kssl_points_30cm <- read.csv(file.path(ksslDir, 'kssl_cluster_30cm_FINAL.csv'), stringsAsFactors = FALSE)
head(kssl_points_30cm)
sum(is.na(kssl_points_30cm$oc_30cm))
kssl_points_30cm$oc_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)] <- kssl_points_30cm$c_tot_30cm[which(is.na(kssl_points_30cm$oc_30cm) & kssl_points_30cm$pH_H2O_30cm < crit_pH)]
sum(is.na(kssl_points_30cm$oc_30cm))
sum(is.na(kssl_points_30cm$cluster_10))
head(kssl_climate_SHR)
kssl_points_30cm$annual.P <- kssl_climate_SHR$annual.P[match(kssl_points_30cm$pedon_key, kssl_climate_SHR$pedn_ky)]
kssl_points_30cm$annual.T <- kssl_climate_SHR$annual.T[match(kssl_points_30cm$pedon_key, kssl_climate_SHR$pedn_ky)]
kssl_points_30cm$cec_to_clay <- kssl_points_30cm$clay_30cm / kssl_points_30cm$cec_7_30cm
plot(kssl_points_30cm$cec_to_clay, kssl_points_30cm$oc_30cm)

summary(lm(oc_30cm ~ cec_to_clay + annual.T + annual.P, data = kssl_points_30cm))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm))
summary(lm(oc_30cm ~  SHR7name + annual.T + annual.P, data = kssl_points_30cm))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5,]))
summary(lm(oc_30cm ~  SHR7name + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5,]))
plot(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5,]))

summary(lm(oc_30cm ~ clay_30cm, data = kssl_points_30cm))
summary(lm(oc_30cm ~ annual.T, data = kssl_points_30cm))
summary(lm(oc_30cm ~ annual.P, data = kssl_points_30cm))
kssl_points_30cm$SHR7name <- clus_7_names[kssl_points_30cm$cluster_7]
kssl_points_30cm$oc_30cm_log <- log(kssl_points_30cm$oc_30cm)
summary(lm(oc_30cm_log ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm))
plot(lm(oc_30cm_log ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm))
summary(lm(kssl_points_30cm$oc_30cm[kssl_points_30cm$SHR7name==] ~ kssl_points_30cm$clay_30cm + kssl_points_30cm$annual.T + kssl_points_30cm$annual.P))

lapply(clus_7_names[order(clus_7_names)], function(x) {
  print(x)
  summary(lm(kssl_points_30cm$oc_30cm[kssl_points_30cm$SHR7name==x] ~ kssl_points_30cm$clay_30cm[kssl_points_30cm$SHR7name==x] + kssl_points_30cm$annual.P[kssl_points_30cm$SHR7name==x] + kssl_points_30cm$annual.T[kssl_points_30cm$SHR7name==x]))
  }
)

summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='1. Coarse & no restrictions',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='2. Loamy & no restrictions',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='3. Coarse-loamy & restrictive layers',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='4. Loamy & restrictive layers',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='5. Coarse-loamy saline-sodic',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='6. Fine saline-sodic',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='7. Shrink-swell',]))

summary(lm(oc_30cm ~ annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='6. Fine saline-sodic',]))

summary(lm(oc_30cm ~ annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='1. Coarse & no restrictions',]))
summary(lm(oc_30cm ~ annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='2. Loamy & no restrictions',]))
summary(lm(oc_30cm ~ annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='3. Coarse-loamy & restrictive layers',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='4. Loamy & restrictive layers',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='5. Coarse-loamy saline-sodic',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='6. Fine saline-sodic',]))
summary(lm(oc_30cm ~ clay_30cm + annual.T + annual.P, data = kssl_points_30cm[kssl_points_30cm$oc_30cm < 3.5 & kssl_points_30cm$SHR7name=='7. Shrink-swell',]))