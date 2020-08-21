laptop <- TRUE
account_for_compost <- FALSE
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
  FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/Figures/valley_final/CalAg' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kerri data'
  CalAgDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/napa lodi analysis'
} else { #on UCD desktop
  dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/summaries/valley_final' #was valley_trial
  FiguresDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/Figures/valley_final' #was valley_trial
  ksslDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kssl'
  kerriDir <- 'C:/Users/smdevine/Desktop/PostDoc/soil health/kerri data'
}

clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')

soil_data <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', "")) #colClasses = columnClass,
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))

#CDFA C validation for 7-region model
#compare 3 management types in cluster 6 across each depth segment
pts_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
pts_30cm$cc_type <- soil_data$perennial.Covercrop.vs.annual.covercrop[match(pts_30cm$Concatenate, soil_data$Concatenate)]
pts_30cm$SHR7name <- clus_7_names[pts_30cm$cluster_7]
pts_10cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_10cm.csv'), stringsAsFactors = FALSE)
pts_10cm$SHR7name <- clus_7_names[pts_10cm$cluster_7]
pts_50cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_50cm.csv'), stringsAsFactors = FALSE)
pts_50cm$SHR7name <- clus_7_names[pts_50cm$cluster_7]

# columnClass <- c(Concatenate='character')
unique(pts_30cm$SHR7name)

#count up
#test overall effects in AOV
napa_lodi_data <- data.frame(ID=pts_30cm$Concatenate[!is.na(pts_30cm$SHR7name)], tillage=pts_30cm$tillage[!is.na(pts_30cm$SHR7name)], irrigation=pts_30cm$irrigated_vs_dryfarm[!is.na(pts_30cm$SHR7name)], compost=pts_30cm$compost_added[!is.na(pts_30cm$SHR7name)], cc_type=pts_30cm$cc_type[!is.na(pts_30cm$SHR7name)], stringsAsFactors = FALSE)
head(napa_lodi_data)
dim(napa_lodi_data)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
napa_lodi_data$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
napa_lodi_data$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
napa_lodi_data$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
napa_lodi_data$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
napa_lodi_data$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(napa_lodi_data$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
napa_lodi_data$SHR7name <- pts_30cm$SHR7name[match(napa_lodi_data$ID, pts_30cm$Concatenate)]
colnames(napa_lodi_data)
sum(apply(napa_lodi_data[,5:9], 1, function(x) {all(!is.na(x))})) #101 rows have all data


table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='yes' & napa_lodi_data$cc_type !='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='yes' & napa_lodi_data$cc_type=='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='yes' & napa_lodi_data$cc_type!='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='yes' & napa_lodi_data$cc_type!='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='yes' & napa_lodi_data$cc_type=='resident')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='none' & napa_lodi_data$cc_type!='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='none' & napa_lodi_data$cc_type=='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='none' & napa_lodi_data$cc_type!='perennial')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='none')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='none' & napa_lodi_data$cc_type!='perennial')])

table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm')])

summary(aov(soil_C_0_5cm ~ tillage + irrigation + compost + cc_type + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + compost + cc_type + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name + tillage:SHR7name, data = napa_lodi_data))
boxplot(soil_C_0_5cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ cc_type, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ SHR7name, data = napa_lodi_data)


summary(aov(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_5_10cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
boxplot(soil_C_5_10cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ cc_type, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ SHR7name, data = napa_lodi_data)


summary(aov(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_10_30cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
boxplot(soil_C_10_30cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ cc_type, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ SHR7name, data = napa_lodi_data)

summary(aov(soil_C_30_50cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_30_50cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_30_50cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_30_50cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
boxplot(soil_C_30_50cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_30_50cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_30_50cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_30_50cm ~ SHR7name, data = napa_lodi_data)

summary(aov(soil_C_50_100cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_50_100cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_50_100cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_50_100cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
boxplot(soil_C_50_100cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_50_100cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_50_100cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_50_100cm ~ SHR7name, data = napa_lodi_data)

napa_lodi_data[which(napa_lodi_data$compost=='none' & napa_lodi_data$SHR7name=='2. Loamy & no restrictions'),]
napa_lodi_data[which(napa_lodi_data$compost=='none' & napa_lodi_data$SHR7name=='1. Coarse & no restrictions'),]
napa_lodi_data[which(napa_lodi_data$compost=='none' & napa_lodi_data$SHR7name=='3. Coarse-loamy & restrictive layers'),]
napa_lodi_data[which(is.na(napa_lodi_data$compost) & napa_lodi_data$SHR7name=='3. Coarse-loamy & restrictive layers'),]
napa_lodi_data[which(napa_lodi_data$compost=='none' & napa_lodi_data$SHR7name=='4. Loamy & restrictive layers'),]

soil_data_shrink_swell <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="7. Shrink-swell")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="7. Shrink-swell")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="7. Shrink-swell")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="7. Shrink-swell")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="7. Shrink-swell")], stringsAsFactors = FALSE)
head(soil_data_shrink_swell)
dim(soil_data_shrink_swell)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_shrink_swell$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_shrink_swell$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_shrink_swell$soil_C_0_10cm <- apply(soil_data_shrink_swell[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_shrink_swell$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_shrink_swell$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_shrink_swell$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_shrink_swell$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_shrink_swell$kgOrg.m2_10cm <- pts_10cm$kgOrg.m2_10cm[match(soil_data_shrink_swell$ID, pts_10cm$Concatenate)]
soil_data_shrink_swell$kgOrg.m2_30cm <- pts_30cm$kgOrg.m2_30cm[match(soil_data_shrink_swell$ID, pts_30cm$Concatenate)]
soil_data_shrink_swell$kgOrg.m2_50cm <- pts_50cm$kgOrg.m2_50cm[match(soil_data_shrink_swell$ID, pts_50cm$Concatenate)]
colnames(soil_data_shrink_swell)
soil_data_shrink_swell <- soil_data_shrink_swell[order(soil_data_shrink_swell$tillage, soil_data_shrink_swell$irrigation, soil_data_shrink_swell$cc_type, soil_data_shrink_swell$compost), ]
soil_data_shrink_swell$SHR <- 7
# write.csv(soil_data_shrink_swell, file.path(CalAgDir, 'soil_data_shrink_swell.csv'), row.names = FALSE)
shrink_swell <- data.frame(C_means=sapply(soil_data_shrink_swell[,8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_shrink_swell[,8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_shrink_swell[,8:11], function(x) {length(x[!is.na(x)])}), C_q25=sapply(soil_data_shrink_swell[,8:11], quantile, probs=0.25, na.rm=TRUE), C_medians=sapply(soil_data_shrink_swell[,8:11], median, na.rm=TRUE), C_q75=sapply(soil_data_shrink_swell[,8:11], quantile, probs=0.75, na.rm=TRUE))
shrink_swell
# write.csv(shrink_swell, file.path(CalAgDir, 'shrink_swell_soilC.csv'), row.names = FALSE)


soil_data_loamy_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], stringsAsFactors = FALSE)
head(soil_data_loamy_w_no_res)
dim(soil_data_loamy_w_no_res)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_loamy_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_no_res$soil_C_0_10cm <- apply(soil_data_loamy_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_loamy_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_loamy_w_no_res$kgOrg.m2_10cm <- pts_10cm$kgOrg.m2_10cm[match(soil_data_loamy_w_no_res$ID, pts_10cm$Concatenate)]
soil_data_loamy_w_no_res$kgOrg.m2_30cm <- pts_30cm$kgOrg.m2_30cm[match(soil_data_loamy_w_no_res$ID, pts_30cm$Concatenate)]
soil_data_loamy_w_no_res$kgOrg.m2_50cm <- pts_50cm$kgOrg.m2_50cm[match(soil_data_loamy_w_no_res$ID, pts_50cm$Concatenate)]
colnames(soil_data_loamy_w_no_res)
soil_data_loamy_w_no_res <- soil_data_loamy_w_no_res[order(soil_data_loamy_w_no_res$tillage, soil_data_loamy_w_no_res$irrigation, soil_data_loamy_w_no_res$cc_type, soil_data_loamy_w_no_res$compost), ]
soil_data_loamy_w_no_res$SHR <- 2
# write.csv(soil_data_loamy_w_no_res, file.path(CalAgDir, 'soil_data_loamy_w_no_res.csv'), row.names = FALSE)

loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[,8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[,8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[,8:11], function(x) {length(x[!is.na(x)])}), C_q25=sapply(soil_data_loamy_w_no_res[,8:11], quantile, probs=0.25, na.rm=TRUE), C_medians=sapply(soil_data_loamy_w_no_res[,8:11], median, na.rm=TRUE), C_q75=sapply(soil_data_loamy_w_no_res[,8:11], quantile, probs=0.75, na.rm=TRUE))
loamy_w_no_res
# write.csv(loamy_w_no_res, file.path(CalAgDir, 'loamy_w_no_res_soilC.csv'), row.names = FALSE)

till_irr_ann_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])})) #2a in plot
till_irr_ann_loamy_w_no_res

notill_irr_ann_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])})) #2b in plot
notill_irr_ann_loamy_w_no_res

notill_irr_per_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & soil_data_loamy_w_no_res$cc_type == 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])}))
notill_irr_per_loamy_w_no_res

till_dry_ann_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])})) #2c in plot
till_dry_ann_loamy_w_no_res

notill_dry_ann_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & soil_data_loamy_w_no_res$cc_type != 'perennial' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='none'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])}))
notill_dry_ann_loamy_w_no_res

soil_data_coarse_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], stringsAsFactors = FALSE)
soil_data_coarse_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_w_no_res$soil_C_0_10cm <- apply(soil_data_coarse_w_no_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_coarse_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_coarse_w_no_res$kgOrg.m2_10cm <- pts_10cm$kgOrg.m2_10cm[match(soil_data_coarse_w_no_res$ID, pts_10cm$Concatenate)]
soil_data_coarse_w_no_res$kgOrg.m2_30cm <- pts_30cm$kgOrg.m2_30cm[match(soil_data_coarse_w_no_res$ID, pts_30cm$Concatenate)]
soil_data_coarse_w_no_res$kgOrg.m2_50cm <- pts_50cm$kgOrg.m2_50cm[match(soil_data_coarse_w_no_res$ID, pts_50cm$Concatenate)]
soil_data_coarse_w_no_res <- soil_data_coarse_w_no_res[order(soil_data_coarse_w_no_res$tillage, soil_data_coarse_w_no_res$irrigation, soil_data_coarse_w_no_res$cc_type, soil_data_coarse_w_no_res$compost), ]
soil_data_coarse_w_no_res$SHR <- 1
# write.csv(soil_data_coarse_w_no_res, file.path(CalAgDir, 'soil_data_coarse_w_no_res.csv'), row.names = FALSE)

coarse_w_no_res <- data.frame(C_means=sapply(soil_data_coarse_w_no_res[,8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_w_no_res[,8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_w_no_res[,8:11], function(x) {length(x[!is.na(x)])}), C_q25=sapply(soil_data_coarse_w_no_res[,8:11], quantile, probs=0.25, na.rm=TRUE), C_medians=sapply(soil_data_coarse_w_no_res[,8:11], median, na.rm=TRUE), C_q75=sapply(soil_data_coarse_w_no_res[,8:11], quantile, probs=0.75, na.rm=TRUE))
coarse_w_no_res
# write.csv(coarse_w_no_res, file.path(CalAgDir, 'coarse_w_no_res_soilC.csv'), row.names = FALSE)

till_irr_coarse_w_no_res <- data.frame(C_means=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])}))
till_irr_coarse_w_no_res #excludes 3 dryfarm coarse w/no res in Napa

soil_data_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], stringsAsFactors = FALSE)
soil_data_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_res$soil_C_0_5cm[soil_data_loamy_w_res$soil_C_0_5cm==0] <- NA #convert 0 value to NA
soil_data_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_res$soil_C_0_10cm <- apply(soil_data_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_loamy_w_res$soil_C_0_10cm[soil_data_loamy_w_res$ID==853] <- soil_data_loamy_w_res$soil_C_5_10cm[soil_data_loamy_w_res$ID==853]
soil_data_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_loamy_w_res <- soil_data_loamy_w_res[order(soil_data_loamy_w_res$tillage, soil_data_loamy_w_res$irrigation, soil_data_loamy_w_res$cc_type, soil_data_loamy_w_res$compost), ]
soil_data_loamy_w_res$kgOrg.m2_10cm <- pts_10cm$kgOrg.m2_10cm[match(soil_data_loamy_w_res$ID, pts_10cm$Concatenate)]
soil_data_loamy_w_res$kgOrg.m2_30cm <- pts_30cm$kgOrg.m2_30cm[match(soil_data_loamy_w_res$ID, pts_30cm$Concatenate)]
soil_data_loamy_w_res$kgOrg.m2_50cm <- pts_50cm$kgOrg.m2_50cm[match(soil_data_loamy_w_res$ID, pts_50cm$Concatenate)]
soil_data_loamy_w_res$SHR <- 4
# write.csv(soil_data_loamy_w_res, file.path(CalAgDir, 'soil_data_loamy_w_res.csv'), row.names = FALSE)
colnames(soil_data_loamy_w_res)
loamy_w_res <- data.frame(C_means=sapply(soil_data_loamy_w_res[,8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_res[,8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[,8:11], function(x) {length(x[!is.na(x)])}), C_q25=sapply(soil_data_loamy_w_res[,8:11], quantile, probs=0.25, na.rm=TRUE), C_medians=sapply(soil_data_loamy_w_res[,8:11], median, na.rm=TRUE), C_q75=sapply(soil_data_loamy_w_res[,8:11], quantile, probs=0.75, na.rm=TRUE))
loamy_w_res
# write.csv(loamy_w_res, file.path(CalAgDir, 'loamy_w_res_soilC.csv'), row.names = FALSE)

notill_irr_ann_loamy_w_res <- data.frame(C_means=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type!='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])})) #4 in plot
notill_irr_ann_loamy_w_res

notill_irr_per_loamy_w_res <- data.frame(C_means=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & soil_data_loamy_w_res$cc_type=='perennial' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),8:11], function(x) {length(x[!is.na(x)])}))
notill_irr_per_loamy_w_res

soil_data_coarse_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], cc_type=pts_30cm$cc_type[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], stringsAsFactors = FALSE)

soil_data_coarse_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_loamy_w_res$soil_C_0_10cm <- apply(soil_data_coarse_loamy_w_res[,c('soil_C_0_5cm', 'soil_C_5_10cm')], 1, mean)
soil_data_coarse_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm[soil_data_coarse_loamy_w_res$soil_C_50_100cm==0] <- NA
soil_data_coarse_loamy_w_res$kgOrg.m2_10cm <- pts_10cm$kgOrg.m2_10cm[match(soil_data_coarse_loamy_w_res$ID, pts_10cm$Concatenate)]
soil_data_coarse_loamy_w_res$kgOrg.m2_30cm <- pts_30cm$kgOrg.m2_30cm[match(soil_data_coarse_loamy_w_res$ID, pts_30cm$Concatenate)]
soil_data_coarse_loamy_w_res$kgOrg.m2_50cm <- pts_50cm$kgOrg.m2_50cm[match(soil_data_coarse_loamy_w_res$ID, pts_50cm$Concatenate)]
soil_data_coarse_loamy_w_res <- soil_data_coarse_loamy_w_res[order(soil_data_coarse_loamy_w_res$tillage, soil_data_coarse_loamy_w_res$irrigation, soil_data_coarse_loamy_w_res$cc_type, soil_data_coarse_loamy_w_res$compost), ]
soil_data_coarse_loamy_w_res$SHR <- 3
# write.csv(soil_data_coarse_loamy_w_res, file.path(CalAgDir, 'soil_data_coarse_loamy_w_res.csv'), row.names = FALSE)
coarse_loamy_w_res <- data.frame(C_means=sapply(soil_data_coarse_loamy_w_res[,8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_loamy_w_res[,8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_loamy_w_res[,8:11], function(x) {length(x[!is.na(x)])}), C_q25=sapply(soil_data_coarse_loamy_w_res[,8:11], quantile, probs=0.25, na.rm=TRUE), C_medians=sapply(soil_data_coarse_loamy_w_res[,8:11], median, na.rm=TRUE), C_q75=sapply(soil_data_coarse_loamy_w_res[,8:11], quantile, probs=0.75, na.rm=TRUE))
coarse_loamy_w_res
# write.csv(coarse_loamy_w_res, file.path(CalAgDir, 'coarse_loamy_w_res_soilC.csv'), row.names = FALSE)

till_irr_ann_coarse_loamy_w_res <- data.frame(C_means=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated' & soil_data_coarse_loamy_w_res$cc_type!='perennial'),8:11], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'& soil_data_coarse_loamy_w_res$cc_type!='perennial'), 8:11], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}), n=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'& soil_data_coarse_loamy_w_res$cc_type!='perennial'), 8:11], function(x) {length(x[!is.na(x)])}))
till_irr_ann_coarse_loamy_w_res

soil_data_allSHR <- rbind(soil_data_coarse_w_no_res, soil_data_loamy_w_no_res, soil_data_coarse_loamy_w_res, soil_data_loamy_w_res, soil_data_shrink_swell)
soil_data_allSHR$SHRcolor <- ifelse(soil_data_allSHR$SHR==1, 'lightgoldenrod', ifelse(soil_data_allSHR$SHR==2, 'tan4', ifelse(soil_data_allSHR$SHR==3, 'gold', ifelse(soil_data_allSHR$SHR==4, 'firebrick3', ifelse(soil_data_allSHR$SHR==7, 'violetred', 'black')))))
table(soil_data_allSHR$SHR)
soil_data_allSHR$Area <- soil_data$Area[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$Vineyard_Mgmt <- soil_data$Vineyard.Management[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$yrs_in_vineyard <- soil_data$years.in.vineyard[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$slope <- soil_data$Cluster.SLOPE[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$great_group <- soil_data$Cluster.GREATGROUP[match(soil_data_allSHR$ID, soil_data$Concatenate)]
soil_data_allSHR$soil_series <- soil_data$Soil.Series.SSURGO[match(soil_data_allSHR$ID, soil_data$Concatenate)]
write.csv(soil_data_allSHR, file.path(CalAgDir, 'soil_data_allSHR.csv'), row.names = FALSE)

lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$Area[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$soil_series[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$great_group[soil_data_allSHR$SHR==x])})
lapply(c(1,2,3,4,7), function(x) {table(soil_data_allSHR$Vineyard_Mgmt[soil_data_allSHR$SHR==x])})
soil_data_allSHR[soil_data_allSHR$SHR==1,]
soil_data_allSHR[soil_data_allSHR$SHR==2,]


depths_cm <- c(5, 20, 40, 75)
clus_7_names
tiff(file = file.path(FiguresDir, 'Napa_Lodi_soilC_profile_0_100cm_SHR_means.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(loamy_w_no_res$C_means, -depths_cm, type='b', xlim=c(0.05,2.7), pch=21, bg='tan4', col='tan4', cex=1, xlab='', ylab='', yaxt='n')
axis(side = 2, at=c(-50, -30, -10), labels=c(50, 30, 10))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(coarse_w_no_res$C_means, -depths_cm, type = 'b', pch=21, bg='lightgoldenrod', col='lightgoldenrod')
points(coarse_loamy_w_res$C_means, -depths_cm, type = 'b', pch=21, bg='gold', col='gold')
points(loamy_w_res$C_means, -depths_cm, type='b', pch=21, bg='firebrick3', col='firebrick3')
points(shrink_swell$C_means, -depths_cm, type='b', pch=21, bg='violetred', col='violetred')
legend('bottomright', legend=c('1. Coarse & no restrictions', '2. Loamy & no restrictions', '3. Coarse-loamy & restrictive horizons','4. Loamy & restrictive horizons','7. Shrink-swell'), pch=21, lty=1, pt.bg=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'), col=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'))
dev.off()

?rgb
gray_transparent <- rgb(red=242, green=242, blue = 242, max =255, alpha=130)
tiff(file = file.path(FiguresDir, 'Napa_Lodi_soilC_profile_0_100cm_SHR_medianIQR.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(loamy_w_no_res$C_medians, -depths_cm, type='n', xlim=c(0.05,3), pch=21, bg='tan4', col='tan4', cex=1, xlab='', ylab='', yaxt='n')
polygon(x=c(loamy_w_no_res$C_q25, rev(loamy_w_no_res$C_q75)), y=c(-depths_cm, -rev(depths_cm)), col=gray_transparent, border = NA)
polygon(x=c(coarse_w_no_res$C_q25, rev(coarse_w_no_res$C_q75)), y=c(-depths_cm, -rev(depths_cm)), col=gray_transparent, border = NA)
polygon(x=c(coarse_loamy_w_res$C_q25, rev(coarse_loamy_w_res$C_q75)), y=c(-depths_cm, -rev(depths_cm)), col=gray_transparent, border = NA)
polygon(x=c(loamy_w_res$C_q25, rev(loamy_w_res$C_q75)), y=c(-depths_cm, -rev(depths_cm)), col=gray_transparent, border = NA)
polygon(x=c(shrink_swell$C_q25, rev(shrink_swell$C_q75)), y=c(-depths_cm, -rev(depths_cm)), col=gray_transparent, border = NA)
lines(x=loamy_w_no_res$C_q25, y=-depths_cm, col='tan4', lty = 2, lwd=0.5)
lines(x=loamy_w_no_res$C_q75, y=-depths_cm, col='tan4', lty = 2, lwd=0.5)
lines(x=coarse_w_no_res$C_q25, y=-depths_cm, col='lightgoldenrod', lty = 2, lwd=0.5)
lines(x=coarse_w_no_res$C_q75, y=-depths_cm, col='lightgoldenrod', lty = 2, lwd=0.5)
lines(x=coarse_loamy_w_res$C_q25, y=-depths_cm, col='gold', lty =2, lwd=0.5)
lines(x=coarse_loamy_w_res$C_q75, y=-depths_cm, col='gold', lty =2, lwd=0.5)
lines(x=loamy_w_res$C_q25, y=-depths_cm, col='firebrick3', lty =2, lwd=0.5)
lines(x=loamy_w_res$C_q75, y=-depths_cm, col='firebrick3', lty =2, lwd=0.5)
lines(x=shrink_swell$C_q25, y=-depths_cm,  col='violetred', lty =2, lwd=0.5)
lines(x=shrink_swell$C_q75, y=-depths_cm,  col='violetred', lty =2, lwd=0.5)
axis(side = 2, at=c(-50, -30, -10), labels=c(50, 30, 10))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(loamy_w_no_res$C_medians, -depths_cm, type='b', pch=21, bg='tan4', col='tan4', lty=1, lwd=0.6, cex=0.8)
points(coarse_w_no_res$C_medians, -depths_cm, type = 'b', pch=21, bg='lightgoldenrod', col='lightgoldenrod', lty=1, lwd=0.6, cex=0.8)
points(coarse_loamy_w_res$C_medians, -depths_cm, type = 'b', pch=21, bg='gold', col='gold', lty=1, lwd=0.6, cex=0.8)
points(loamy_w_res$C_medians, -depths_cm, type='b', pch=21, bg='firebrick3', col='firebrick3', lty=1, lwd=0.6, cex=0.8)
points(shrink_swell$C_medians, -depths_cm, type='b', pch=21, bg='violetred', col='violetred', lty=1, lwd=0.6, cex=0.8)
legend('bottomright', legend=c('1. Coarse & no restrictions', '2. Loamy & no restrictions', '3. Coarse-loamy & restrictive horizons','4. Loamy & restrictive horizons','7. Shrink-swell'), pch=21, lty=1, pt.bg=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'), col=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'))
dev.off()

tiff(file = file.path(FiguresDir, 'Napa_Lodi_soilC_profile_0_100cm_SHR_alldata.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(loamy_w_no_res$C_means, -depths_cm, type='b', xlim=c(0.05,4.5), pch=21, bg='tan4', col='tan4', cex=1, xlab='', ylab='', yaxt='n')
axis(side = 2, at=c(-50, -30, -10), labels=c(50, 30, 10))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(coarse_w_no_res$C_means, -depths_cm, type = 'b', pch=21, bg='lightgoldenrod', col='lightgoldenrod')
points(coarse_loamy_w_res$C_means, -depths_cm, type = 'b', pch=21, bg='gold', col='gold')
points(loamy_w_res$C_means, -depths_cm, type='b', pch=21, bg='firebrick3', col='firebrick3')
points(shrink_swell$C_means, -depths_cm, type='b', pch=21, bg='violetred', col='violetred')
# segments(x0=loamy_w_res$C_means + qt(0.975,df=loamy_w_res$n)*loamy_w_res$C_se, y0=-depths_cm, x1=loamy_w_res$C_means + qt(0.025,df=loamy_w_res$n)*loamy_w_res$C_se, y1=-depths_cm, lwd = 1)
points(x=soil_data_allSHR$soil_C_0_10cm[!is.na(soil_data_allSHR$soil_C_0_10cm)], y=rep(-depths_cm[1], sum(!is.na(soil_data_allSHR$soil_C_0_10cm))), col=soil_data_allSHR$SHRcolor[!is.na(soil_data_allSHR$soil_C_0_10cm)], pch=1, cex=0.6)
points(x=soil_data_allSHR$soil_C_10_30cm[!is.na(soil_data_allSHR$soil_C_10_30cm)], y=rep(-depths_cm[2], sum(!is.na(soil_data_allSHR$soil_C_10_30cm))), col=soil_data_allSHR$SHRcolor[!is.na(soil_data_allSHR$soil_C_10_30cm)], pch=1, cex=0.6)
points(x=soil_data_allSHR$soil_C_30_50cm[!is.na(soil_data_allSHR$soil_C_30_50cm)], y=rep(-depths_cm[3], sum(!is.na(soil_data_allSHR$soil_C_30_50cm))), col=soil_data_allSHR$SHRcolor[!is.na(soil_data_allSHR$soil_C_30_50cm)], pch=1, cex=0.6)
points(x=soil_data_allSHR$soil_C_50_100cm[!is.na(soil_data_allSHR$soil_C_50_100cm)], y=rep(-depths_cm[4], sum(!is.na(soil_data_allSHR$soil_C_50_100cm))), col=soil_data_allSHR$SHRcolor[!is.na(soil_data_allSHR$soil_C_50_100cm)], pch=1, cex=0.6)
legend('bottomright', legend=c('1. Coarse & no restrictions', '2. Loamy & no restrictions', '3. Coarse-loamy & restrictive horizons','4. Loamy & restrictive horizons','7. Shrink-swell'), pch=21, lty=1, pt.bg=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'), col=c('lightgoldenrod', 'tan4', 'gold', 'firebrick3', 'violetred'))
dev.off()
#option 2 names: c('1. Very coarse w/no restrictions', '2a. Coarse w/no resrictions', '2b. Coarse w/no resrictions: No-till', '2c. Coarse w/no resrictions: Dry farm', '3. Very coarse w/pans', '4. Coarse w/pans: No-till')

#make barplot of 30cm SOC kg m^-2
#see confidence interval formula for unknown mean and standard deviation at http://www.stat.yale.edu/Courses/1997-98/101/confint.htm
plot_da_error <- function(barnum, barcenters, df) {
  segments(x0=barcenters[barnum,], y0=df$means[barnum] + qt(p=0.975, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), x1=barcenters[barnum,], y1=df$means[barnum] + qt(p=0.025, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), lwd = 1.2)
}

#SOC concentration barplots

Cmeans_SHRmgmt <- rbind(till_irr_coarse_w_no_res$C_means, till_irr_ann_loamy_w_no_res$C_means, notill_irr_ann_loamy_w_no_res$C_means, notill_irr_per_loamy_w_no_res$C_means, till_dry_ann_loamy_w_no_res$C_means, notill_dry_ann_loamy_w_no_res$C_means, till_irr_ann_coarse_loamy_w_res$C_means, notill_irr_ann_loamy_w_res$C_means, notill_irr_per_loamy_w_res$C_means, shrink_swell$C_means)
row.names(Cmeans_SHRmgmt) <- c('1', '2a', '2b', '2c', '2d', '2e', '3', '4a', '4b', '7')
colnames(Cmeans_SHRmgmt) <- c('0-10 cm', '10-30 cm', '30-50 cm', '50-100 cm')
C_SE_SHRmgmt <- rbind(till_irr_coarse_w_no_res$C_se, till_irr_ann_loamy_w_no_res$C_se, notill_irr_ann_loamy_w_no_res$C_se, notill_irr_per_loamy_w_no_res$C_se, till_dry_ann_loamy_w_no_res$C_se, notill_dry_ann_loamy_w_no_res$C_se, till_irr_ann_coarse_loamy_w_res$C_se, notill_irr_ann_loamy_w_res$C_se, notill_irr_per_loamy_w_res$C_se, shrink_swell$C_se)
C_n_SHRmgmt <- rbind(till_irr_coarse_w_no_res$n, till_irr_ann_loamy_w_no_res$n, notill_irr_ann_loamy_w_no_res$n, notill_irr_per_loamy_w_no_res$n, till_dry_ann_loamy_w_no_res$n, notill_dry_ann_loamy_w_no_res$n, till_irr_ann_coarse_loamy_w_res$n, notill_irr_ann_loamy_w_res$n, notill_irr_per_loamy_w_res$n, shrink_swell$n)

tiff(file = file.path(FiguresDir, "Napa Lodi", 'Napa_Lodi_soilC_barchart_alldepths_v3.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 4.25, units = 'in', res=800, compression = 'lzw')
par(mar=c(2.4, 3.5, 1.25, 0.5))
Cmeans_bardims <- barplot(Cmeans_SHRmgmt, beside=TRUE, legend.text=c('1. Coarse & no restrictions (n=9)', '2a. Loamy & no restrictions (n=18)', '2b. Loamy & no restrictions: no-till (n=5)', '2c. Loamy & no restrictions: no-till & perennial (n=13)', '2d. Loamy & no restrictions: dryfarm (n=12)', '2e. Loamy & no restrictions: no-till & dryfarm (n=6)', '3. Coarse-loamy & restrictive horizons (n=15)', '4a. Loamy & restrictive horizons: no-till (n=5)', '4b. Loamy & restrictive horizons: no-till & perennial (n=6)', '7. Shrink-swell (n=4)'), col=c('lightgoldenrod', 'tan4', 'tan4', 'tan4',  'tan4', 'tan4', 'gold', 'firebrick3', 'firebrick3', 'violetred'), ylab='', ylim=c(0,4.2), args.legend=list(x='topright', cex=1, pt.cex=1, bty='n', ncol=1, inset=c(0,-0.05), x.intersp=0.4), cex.names=1)
mtext(text = row.names(Cmeans_SHRmgmt), side = 1, line=0.1, at=Cmeans_bardims, cex=0.9)
mtext('Soil organic carbon (%)', side=2, line=2.25)
segments(x0=Cmeans_bardims[1,], y0=Cmeans_SHRmgmt[1,] + qt(0.975,df=C_n_SHRmgmt[1,]) * C_SE_SHRmgmt[1,], x1=Cmeans_bardims[1,], y1=Cmeans_SHRmgmt[1,] + qt(0.025,df=C_n_SHRmgmt[1,]) * C_SE_SHRmgmt[1,], lwd = 1)
segments(x0=Cmeans_bardims[2,], y0=Cmeans_SHRmgmt[2,] + qt(0.975,df=C_n_SHRmgmt[2,]) * C_SE_SHRmgmt[2,], x1=Cmeans_bardims[2,], y1=Cmeans_SHRmgmt[2,] + qt(0.025,df=C_n_SHRmgmt[2,]) * C_SE_SHRmgmt[2,], lwd = 1)
segments(x0=Cmeans_bardims[3,], y0=Cmeans_SHRmgmt[3,] + qt(0.975,df=C_n_SHRmgmt[3,]) * C_SE_SHRmgmt[3,], x1=Cmeans_bardims[3,], y1=Cmeans_SHRmgmt[3,] + qt(0.025,df=C_n_SHRmgmt[3,]) * C_SE_SHRmgmt[3,], lwd = 1)
segments(x0=Cmeans_bardims[4,], y0=Cmeans_SHRmgmt[4,] + qt(0.975,df=C_n_SHRmgmt[4,]) * C_SE_SHRmgmt[4,], x1=Cmeans_bardims[4,], y1=Cmeans_SHRmgmt[4,] + qt(0.025,df=C_n_SHRmgmt[4,]) * C_SE_SHRmgmt[4,], lwd = 1)
segments(x0=Cmeans_bardims[5,], y0=Cmeans_SHRmgmt[5,] + qt(0.975,df=C_n_SHRmgmt[5,]) * C_SE_SHRmgmt[5,], x1=Cmeans_bardims[5,], y1=Cmeans_SHRmgmt[5,] + qt(0.025,df=C_n_SHRmgmt[5,]) * C_SE_SHRmgmt[5,], lwd = 1)
segments(x0=Cmeans_bardims[6,], y0=Cmeans_SHRmgmt[6,] + qt(0.975,df=C_n_SHRmgmt[6,]) * C_SE_SHRmgmt[6,], x1=Cmeans_bardims[6,], y1=Cmeans_SHRmgmt[6,] + qt(0.025,df=C_n_SHRmgmt[6,]) * C_SE_SHRmgmt[6,], lwd = 1)
segments(x0=Cmeans_bardims[7,], y0=Cmeans_SHRmgmt[7,] + qt(0.975,df=C_n_SHRmgmt[7,]) * C_SE_SHRmgmt[7,], x1=Cmeans_bardims[7,], y1=Cmeans_SHRmgmt[7,] + qt(0.025,df=C_n_SHRmgmt[7,]) * C_SE_SHRmgmt[7,], lwd = 1)
segments(x0=Cmeans_bardims[8,], y0=Cmeans_SHRmgmt[8,] + qt(0.975,df=C_n_SHRmgmt[8,]) * C_SE_SHRmgmt[8,], x1=Cmeans_bardims[8,], y1=Cmeans_SHRmgmt[8,] + qt(0.025,df=C_n_SHRmgmt[8,]) * C_SE_SHRmgmt[8,], lwd = 1)
segments(x0=Cmeans_bardims[9,], y0=Cmeans_SHRmgmt[9,] + qt(0.975,df=C_n_SHRmgmt[9,]) * C_SE_SHRmgmt[9,], x1=Cmeans_bardims[9,], y1=Cmeans_SHRmgmt[9,] + qt(0.025,df=C_n_SHRmgmt[9,]) * C_SE_SHRmgmt[9,], lwd = 1)
segments(x0=Cmeans_bardims[10,], y0=Cmeans_SHRmgmt[10,] + qt(0.975,df=C_n_SHRmgmt[10,]) * C_SE_SHRmgmt[10,], x1=Cmeans_bardims[10,], y1=Cmeans_SHRmgmt[10,] + qt(0.025,df=C_n_SHRmgmt[10,]) * C_SE_SHRmgmt[10,], lwd = 1)
dev.off()
#0-30 cm data aggregation


#simplified classes for comparison updated 7/21/20
kgOrgC_30cm_data <- data.frame(means=c(mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type!='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$cc_type=='perennial' & pts_30cm$SHR7name=="4. Loamy & restrictive layers")])))) #, row.names = c('1. Coarse & no restrictions', '2a. Loamy & no restrictions', '2b. Loamy & no restrictions: no-till', '2c. Loamy & no restrictions: no-till & perennial', '3. Coarse-loamy & restrictive horizons','4a. Loamy & restrictive horizons: no-till', '4b. Loamy & restrictive horizons: no-till & perennial'))
row.names(kgOrgC_30cm_data) <- c('1. Coarse & no restrictions', '2a. Loamy & no restrictions', '2b. Loamy & no restrictions: no-till', '2c. Loamy & no restrictions: no-till & perennial', '3. Coarse-loamy & restrictive horizons','4a. Loamy & restrictive horizons: no-till', '4b. Loamy & restrictive horizons: no-till & perennial')
kgOrgC_30cm_data


tiff(file = file.path(FiguresDir, 'CDFA_soilC_kg_0_30cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_30cm <- barplot(height=kgOrgC_30cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,8.5))
mtext(text=expression('0-30 cm SOC (kg'~m^-2*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(2, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(3, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(4, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(5, barcenters_30cm, kgOrgC_30cm_data)
plot_da_error(6, barcenters_30cm, kgOrgC_30cm_data)
dev.off()

#30 cm test
pts_30cm$aov_code <- ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions", 'Coarse w/no restrictions', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions', ifelse(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: no-till', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: dry farm', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers", 'Coarse w/pans', ifelse(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers", 'Loamy w/pans: no-till', NA))))))

tapply(pts_30cm$kgOrg.m2_30cm, pts_30cm$aov_code, mean, na.rm=TRUE)
kgOrgC_30cm_data
summary(aov(kgOrg.m2_30cm ~ aov_code, data=pts_30cm))
TukeyHSD(aov(kgOrg.m2_30cm ~ aov_code, data=pts_30cm))

