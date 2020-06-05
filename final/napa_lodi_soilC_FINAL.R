laptop <- TRUE
account_for_compost <- TRUE
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

clus_7_names <- c('6. Fine saline-sodic', '3. Coarse-loamy & restrictive layers', '4. Loamy & restrictive layers', '1. Coarse & no restrictions', '2. Loamy & no restrictions', '7. Shrink-swell', '5. Coarse-loamy saline-sodic')

#CDFA C validation for 7-region model
#compare 3 management types in cluster 6 across each depth segment
pts_30cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_30cm.csv'), stringsAsFactors = FALSE)
pts_30cm$SHR7name <- clus_7_names[pts_30cm$cluster_7]
pts_10cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_10cm.csv'), stringsAsFactors = FALSE)
pts_10cm$SHR7name <- clus_7_names[pts_10cm$cluster_7]
pts_50cm <- read.csv(file.path(kerriDir, 'FINAL', 'CDFA_samples_cluster_50cm.csv'), stringsAsFactors = FALSE)
pts_50cm$SHR7name <- clus_7_names[pts_50cm$cluster_7]

# columnClass <- c(Concatenate='character')
soil_data <- read.csv(file.path(kerriDir, 'CDFA Soil Survey All Data_copy.csv'), stringsAsFactors = FALSE, na.strings = c('#N/A', 'pOOR DATA', 'POOR DATA', 'NO Sample', 'Missing Data', '?', "")) #colClasses = columnClass,
dim(soil_data)
soil_data$GPS.N <- as.numeric(gsub("N", "", soil_data$GPS.N))
soil_data$GPS.W <- -as.numeric(gsub("W", "", soil_data$GPS.W))

unique(pts_30cm$SHR7name)
soil_data_loamy_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="2. Loamy & no restrictions")], stringsAsFactors = FALSE)
head(soil_data_loamy_w_no_res)
dim(soil_data_loamy_w_no_res)
length(soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])])
# "0_5"    "5_10"   "10_30"  "30_50"  "50_100"
soil_data_loamy_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
colnames(soil_data_loamy_w_no_res)
lapply(soil_data_loamy_w_no_res[,5:9], summary)
lapply(soil_data_loamy_w_no_res[,5:9], function(x) tapply(x, soil_data_loamy_w_no_res$tillage, mean, na.rm=TRUE))
lapply(soil_data_loamy_w_no_res[,5:9], function(x) tapply(x, soil_data_loamy_w_no_res$irrigation, mean, na.rm=TRUE))
lapply(soil_data_loamy_w_no_res[,5:9], function(x) tapply(x, soil_data_loamy_w_no_res$compost, mean, na.rm=TRUE))

#test overall effects in AOV
napa_lodi_data <- data.frame(ID=pts_30cm$Concatenate[!is.na(pts_30cm$SHR7name)], tillage=pts_30cm$tillage[!is.na(pts_30cm$SHR7name)], irrigation=pts_30cm$irrigated_vs_dryfarm[!is.na(pts_30cm$SHR7name)], compost=pts_30cm$compost_added[!is.na(pts_30cm$SHR7name)], stringsAsFactors = FALSE)
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
sum(apply(napa_lodi_data[,5:9], 1, function(x) {all(!is.na(x))})) #100 rows have all data

table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='yes')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='yes')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='yes')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='yes')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='none')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated' & napa_lodi_data$compost=='none')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='none')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm' & napa_lodi_data$compost=='none')])

table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='irrigated')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='irrigated')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='till' & napa_lodi_data$irrigation=='dryfarm')])
table(napa_lodi_data$SHR7name[which(napa_lodi_data$tillage=='no till' & napa_lodi_data$irrigation=='dryfarm')])

summary(aov(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_0_5cm ~ tillage + irrigation + compost + SHR7name + tillage:SHR7name, data = napa_lodi_data))
boxplot(soil_C_0_5cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_0_5cm ~ SHR7name, data = napa_lodi_data)

summary(aov(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_5_10cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_5_10cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
boxplot(soil_C_5_10cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ compost, data = napa_lodi_data)
boxplot(soil_C_5_10cm ~ SHR7name, data = napa_lodi_data)


summary(aov(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
TukeyHSD(aov(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_10_30cm ~ tillage + irrigation + compost + SHR7name, data = napa_lodi_data))
summary(lm(soil_C_10_30cm ~ tillage + irrigation + SHR7name, data = napa_lodi_data))
boxplot(soil_C_10_30cm ~ tillage, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ irrigation, data = napa_lodi_data)
boxplot(soil_C_10_30cm ~ compost, data = napa_lodi_data)
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

#make a barplot showing compost, irrigation, and tillage by SHR
#compost application has 6 NAs
practicesbySHR <- data.frame(SHR7name=unique(napa_lodi_data$SHR7name)[order(unique(napa_lodi_data$SHR7name))], compost_yes=as.numeric(tapply(napa_lodi_data$compost, napa_lodi_data$SHR7name, function(x){sum(x=='yes', na.rm = TRUE)})), compost_none=as.numeric(tapply(napa_lodi_data$compost, napa_lodi_data$SHR7name, function(x){sum(x=='none', na.rm = TRUE)})), compost_NA=as.numeric(tapply(napa_lodi_data$compost, napa_lodi_data$SHR7name, function(x){sum(is.na(x))})), irrigated=as.numeric(tapply(napa_lodi_data$irrigation, napa_lodi_data$SHR7name, function(x){sum(x=='irrigated')})), dryfarm=as.numeric(tapply(napa_lodi_data$irrigation, napa_lodi_data$SHR7name, function(x){sum(x=='dryfarm')})), tillage=as.numeric(tapply(napa_lodi_data$tillage, napa_lodi_data$SHR7name, function(x){sum(x=='till')})), no_till=as.numeric(tapply(napa_lodi_data$tillage, napa_lodi_data$SHR7name, function(x){sum(x=='no till')})))
compostbySHR <- data.frame(SHR7name=practicesbySHR$SHR7name, compost_yes=practicesbySHR$compost_yes/apply(practicesbySHR[,2:4], 1, sum), compost_NA=practicesbySHR$compost_NA/apply(practicesbySHR[,2:4], 1, sum), compost_no=practicesbySHR$compost_none/apply(practicesbySHR[,2:4], 1, sum))
tillagebySHR <- data.frame(SHR7name=practicesbySHR$SHR7name, tillage=practicesbySHR$tillage/apply(practicesbySHR[,7:8], 1, sum), no_till=practicesbySHR$no_till/apply(practicesbySHR[,7:8], 1, sum))
irrigationbySHR <- data.frame(SHR7name=practicesbySHR$SHR7name, irrigated=practicesbySHR$irrigated/apply(practicesbySHR[,5:6], 1, sum), dryfarm=practicesbySHR$dryfarm/apply(practicesbySHR[,5:6], 1, sum))

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'compost_by_SHR.tif'), pointsize = 11, family = 'Times New Roman', width = 3, height = 2.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(0.75, 4, 0.75, 0.5))
barplot(t(as.matrix(compostbySHR[1:4,2:4])*100), beside=FALSE, axisnames=FALSE, ylab='Vineyards applying compost (%)')
text('yes', x=0.7, y=40)
text('none', x=0.7, y=92)
text('NA', x=3.1, y=67)
text('a', x=4.3, y=97)
dev.off()

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'tillage_by_SHR.tif'), pointsize = 11, family = 'Times New Roman', width = 3, height = 2.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(0.75, 4, 0.75, 0.5))
barplot(t(as.matrix(tillagebySHR[1:4,2:3])*100), beside=FALSE, axisnames=FALSE, ylab='Vineyard tillage (%)')
text('yes', x=0.7, y=40)
text('no-till', x=0.7, y=92)
text('b', x=4.3, y=97)
dev.off()

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'irrigation_by_SHR.tif'), pointsize = 11, family = 'Times New Roman', width = 3, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 4, 0.75, 0.5))
barplot(t(as.matrix(irrigationbySHR[1:4,2:3])*100), beside=FALSE, ylab='Vinyeard irrigation (%)', xlab='', axisnames=FALSE)
mtext(1:4, side=1, line=0.4, at=c(0.7, 1.9, 3.1, 4.3))
mtext('Soil health region', side=1, line=1.5)
text('yes', x=0.7, y=40)
text('dry-\nfarm', x=0.7, y=90)
text('c', x=4.3, y=97)
dev.off()

till_irr_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}))

notill_irr_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}))#2b in plot

till_dry_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='till' & soil_data_loamy_w_no_res$irrigation=='dryfarm' & if(account_for_compost){soil_data_loamy_w_no_res$compost=='yes'}else{TRUE}),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))})) #2c in plot

if(!account_for_compost){notill_dry_loamy_w_no_res <- data.frame(C_means=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm'),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_no_res[which(soil_data_loamy_w_no_res$tillage=='no till' & soil_data_loamy_w_no_res$irrigation=='dryfarm'),5:9], function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}))} #2d in plot

soil_data_coarse_w_no_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="1. Coarse & no restrictions")], stringsAsFactors = FALSE)
soil_data_coarse_w_no_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_w_no_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_w_no_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_w_no_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_w_no_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_w_no_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]

till_irr_coarse_w_no_res <- data.frame(C_means=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='irrigated' & if(account_for_compost){soil_data_coarse_w_no_res$compost=='yes'}else{TRUE}),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}))

if(!account_for_compost){till_dry_coarse_w_no_res <- data.frame(C_means=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='dryfarm'),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_w_no_res[which(soil_data_coarse_w_no_res$tillage=='till' & soil_data_coarse_w_no_res$irrigation=='dryfarm'),5:9], function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}))} #1b in plot

soil_data_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="4. Loamy & restrictive layers")], stringsAsFactors = FALSE)
soil_data_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_loamy_w_res$soil_C_0_5cm[soil_data_loamy_w_res$soil_C_0_5cm==0] <- NA #convert 0 value to NA

notill_irr_loamy_w_res <- data.frame(C_means=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_loamy_w_res[which(soil_data_loamy_w_res$tillage=='no till' & soil_data_loamy_w_res$irrigation=='irrigated' & if(account_for_compost){soil_data_loamy_w_res$compost=='yes'}else{TRUE}),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))})) #4 in plot

soil_data_coarse_loamy_w_res <- data.frame(ID=pts_30cm$Concatenate[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], tillage=pts_30cm$tillage[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], irrigation=pts_30cm$irrigated_vs_dryfarm[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], compost=pts_30cm$compost_added[which(pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], stringsAsFactors = FALSE)
soil_data_coarse_loamy_w_res$soil_C_0_5cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='0_5'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='0_5'])]
soil_data_coarse_loamy_w_res$soil_C_5_10cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='5_10'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='5_10'])]
soil_data_coarse_loamy_w_res$soil_C_10_30cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='10_30'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='10_30'])]
soil_data_coarse_loamy_w_res$soil_C_30_50cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='30_50'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='30_50'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm <-  soil_data$C....g.g..Oven.dry.converted[soil_data$Depth..cm.=='50_100'][match(soil_data_coarse_loamy_w_res$ID, soil_data$Concatenate[soil_data$Depth..cm.=='50_100'])]
soil_data_coarse_loamy_w_res$soil_C_50_100cm[soil_data_coarse_loamy_w_res$soil_C_50_100cm==0] <- NA

till_irr_coarse_loamy_w_res <- data.frame(C_means=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}))

notill_irr_coarse_loamy_w_res <- data.frame(C_means=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='no till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'),5:9], mean, na.rm=TRUE), C_se=sapply(soil_data_coarse_loamy_w_res[which(soil_data_coarse_loamy_w_res$tillage=='no till' & soil_data_coarse_loamy_w_res$irrigation=='irrigated'),5:9], function(x) {sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))}))


depths_cm <- c(2.5, 7.5, 20, 40, 75)
clus_7_names
tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'CDFA_soilC_profile_v2.tif'), pointsize = 11, family = 'Times New Roman', width = 9, height = 3.5, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 3.5, 0.5, 0.5))
plot(till_irr_loamy_w_no_res$C_means, -depths_cm, type='b', xlim=c(0.1,3.1), pch=24, bg='tan4', col='tan4', cex=1, xlab='', ylab='', yaxt='n')
axis(side = 2, at=seq(from=-100, to=0, by=20), labels=-seq(from=-100, to=0, by=20))
mtext(text='Soil organic carbon (%)', side = 1, line=2.25)
mtext(text='Depth (cm)', side = 2, line=2.25)
points(till_dry_loamy_w_no_res$C_means, -depths_cm, type='b', pch=24, bg='tan4', col='tan4', lty=2)
points(notill_irr_loamy_w_no_res$C_means, -depths_cm, type='b', pch=23, bg='tan4', col='tan4')
points(notill_dry_loamy_w_no_res$C_means, -depths_cm, type='b', pch=23, bg='tan4', col='tan4', lty=2)
points(till_irr_coarse_w_no_res$C_means, -depths_cm, type = 'b', pch=24, bg='tan2', col='tan2')
points(till_dry_coarse_w_no_res$C_means, -depths_cm, type = 'b', pch=24, bg='tan2', col='tan2', lty=2)
points(till_irr_coarse_loamy_w_res$C_means, -depths_cm, type = 'b', pch=24, bg='gold', col='gold')
points(notill_irr_coarse_loamy_w_res$C_means, -depths_cm, type = 'b', pch=23, bg='gold', col='gold')
points(notill_irr_loamy_w_res$C_means, -depths_cm, type='b', pch=23, bg='firebrick3', col='firebrick3')
legend('bottomright', legend=c('1. Coarse & no restrictions', '1b. Coarse & no restrictions: Dryfarm', '2a. Loamy & no restrictions', '2b. Loamy & no restrictions: No-till', '2c. Loamy & no restrictions: Dryfarm', '2d. Loamy & no restrictions: No-till & Dryfarm', '3. Coarse-loamy & restrictive layers', '3b. Coarse-loamy & restrictive layers: No-till', '4. Loamy & restrictive layers: No-till'), pch=c(24,24,24,23,24,23,24,23,23), lty=c(1,2,1,1,2,2,1,1,1), pt.bg=c('tan2', 'tan2', 'tan4', 'tan4', 'tan4', 'tan4', 'gold', 'gold', 'firebrick3'), col=c('tan2', 'tan2', 'tan4', 'tan4', 'tan4', 'tan4', 'gold', 'gold', 'firebrick3'))
dev.off()
#option 2 names: c('1. Very coarse w/no restrictions', '2a. Coarse w/no resrictions', '2b. Coarse w/no resrictions: No-till', '2c. Coarse w/no resrictions: Dry farm', '3. Very coarse w/pans', '4. Coarse w/pans: No-till')

#make barplot of 30cm SOC kg m^-2
#see confidence interval formula for unknown mean and standard deviation at http://www.stat.yale.edu/Courses/1997-98/101/confint.htm
plot_da_error <- function(barnum, barcenters, df) {
  segments(x0=barcenters[barnum,], y0=df$means[barnum] + qt(p=0.975, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), x1=barcenters[barnum,], y1=df$means[barnum] + qt(p=0.025, df=df$n[barnum]-1) * df$sd[barnum] / sqrt(df$n[barnum]), lwd = 1.2)
}
kgOrgC_30cm_data <- data.frame(means=c(mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_30cm$kgOrg.m2_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'CDFA_soilC_kg_0_30cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
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
kgOrgC_50cm_data <- data.frame(means=c(mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'CDFA_soilC_kg_0_50cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
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
kgOrgC_10cm_data <- data.frame(means=c(mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_10cm$kgOrg.m2_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))

tiff(file = file.path(FiguresDir, 'FINAL', 'Napa-Lodi', 'CDFA_soilC_kg_0_10cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
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
pts_50cm$aov_code <- ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="1. Coarse & no restrictions", 'Coarse w/no restrictions', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: no-till', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$compost_added=='yes' & pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: dry farm', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers", 'Coarse w/pans', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$compost_added=='yes'& pts_50cm$SHR7name=="4. Loamy & restrictive layers", 'Loamy w/pans: no-till', NA))))))

tapply(pts_50cm$kgOrg.m2_50cm, pts_50cm$aov_code, mean, na.rm=TRUE)
kgOrgC_50cm_data
summary(aov(kgOrg.m2_50cm ~ aov_code, data=pts_50cm))
TukeyHSD(aov(kgOrg.m2_50cm ~ aov_code, data=pts_50cm))

#30 cm test
pts_30cm$aov_code <- ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions", 'Coarse w/no restrictions', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions', ifelse(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: no-till', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: dry farm', ifelse(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers", 'Coarse w/pans', ifelse(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers", 'Loamy w/pans: no-till', NA))))))

tapply(pts_30cm$kgOrg.m2_30cm, pts_30cm$aov_code, mean, na.rm=TRUE)
kgOrgC_30cm_data
summary(aov(kgOrg.m2_30cm ~ aov_code, data=pts_30cm))
TukeyHSD(aov(kgOrg.m2_30cm ~ aov_code, data=pts_30cm))

#10 cm test
pts_10cm$aov_code <- ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions", 'Coarse w/no restrictions', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions', ifelse(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: no-till', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: dry farm', ifelse(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers", 'Coarse w/pans', ifelse(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers", 'Loamy w/pans: no-till', NA))))))

tapply(pts_10cm$kgOrg.m2_10cm, pts_10cm$aov_code, mean, na.rm=TRUE)
kgOrgC_10cm_data
summary(aov(kgOrg.m2_10cm ~ aov_code, data=pts_10cm))
TukeyHSD(aov(kgOrg.m2_10cm ~ aov_code, data=pts_10cm))

#bulk density plots
#0-10 cm
BD_10cm_data <- data.frame(means=c(mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='dryfarm' & pts_10cm$compost_added=='yes' & pts_10cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_10cm$bd_10cm[which(pts_10cm$tillage=='no till' & pts_10cm$irrigated_vs_dryfarm=='irrigated' & pts_10cm$compost_added=='yes'& pts_10cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))
summary(aov(bd_10cm ~ aov_code, data=pts_10cm))
TukeyHSD(aov(bd_10cm ~ aov_code, data=pts_10cm))
tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_BD_0_10cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_10cm <- barplot(height=BD_10cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,3.9))
mtext(text=expression('0-10 cm bulk density (g'~cm^-3*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_10cm, BD_10cm_data)
plot_da_error(2, barcenters_10cm, BD_10cm_data)
plot_da_error(3, barcenters_10cm, BD_10cm_data)
plot_da_error(4, barcenters_10cm, BD_10cm_data)
plot_da_error(5, barcenters_10cm, BD_10cm_data)
plot_da_error(6, barcenters_10cm, BD_10cm_data)
#text(x=barcenters_10cm, y=0.2, labels = c('A', 'ABC', 'BC', 'C', 'AB', 'BC'), adj=0.5)
dev.off()

#bulk density 0-30 cm
BD_30cm_data <- data.frame(means=c(mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='dryfarm' & pts_30cm$compost_added=='yes' & pts_30cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_30cm$bd_30cm[which(pts_30cm$tillage=='no till' & pts_30cm$irrigated_vs_dryfarm=='irrigated' & pts_30cm$compost_added=='yes'& pts_30cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))
summary(aov(bd_30cm ~ aov_code, data=pts_30cm))
TukeyHSD(aov(bd_30cm ~ aov_code, data=pts_30cm))
tiff(file = file.path(FiguresDir, 'v2', 'CDFA validation', '7 region', 'CDFA_BD_0_30cm.tif'), pointsize = 11, family = 'Times New Roman', width = 4.5, height = 3, units = 'in', res=800, compression = 'lzw')
par(mar=c(3, 4, 0.5, 0.5))
barcenters_30cm <- barplot(height=BD_30cm_data$means, col=c('tan2', 'tan4', 'tan4', 'tan4', 'gold', 'firebrick3'), ylab='', names.arg=c('1', '2a', '2b', '2c', '3', '4'), ylim=c(0,3.9))
mtext(text=expression('0-30 cm bulk density (g'~cm^-3*')'), side = 2, line = 2.5)
plot_da_error(1, barcenters_30cm, BD_30cm_data)
plot_da_error(2, barcenters_30cm, BD_30cm_data)
plot_da_error(3, barcenters_30cm, BD_30cm_data)
plot_da_error(4, barcenters_30cm, BD_30cm_data)
plot_da_error(5, barcenters_30cm, BD_30cm_data)
plot_da_error(6, barcenters_30cm, BD_30cm_data)
#text(x=barcenters_30cm, y=0.2, labels = c('A', 'ABC', 'BC', 'C', 'AB', 'BC'), adj=0.5)
dev.off()

#test including 'no compost'
#analysis of variance
pts_50cm$aov_code_test <- ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="1. Coarse & no restrictions", 'Coarse w/no restrictions', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: no-till', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$SHR7name=="2. Loamy & no restrictions", 'Loamy w/no restrictions: dry farm', ifelse(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers", 'Coarse w/pans', ifelse(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="4. Loamy & restrictive layers", 'Loamy w/pans: no-till', NA))))))

tapply(pts_50cm$kgOrg.m2_50cm, pts_50cm$aov_code_test, mean, na.rm=TRUE)
tapply(pts_50cm$kgOrg.m2_50cm, pts_50cm$aov_code, mean, na.rm=TRUE)
kgOrgC_50cm_data
summary(aov(kgOrg.m2_50cm ~ aov_code_test, data=pts_50cm))
TukeyHSD(aov(kgOrg.m2_50cm ~ aov_code_test, data=pts_50cm))

kgOrgC_50cm_data_test <- data.frame(means=c(mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), mean(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), sd=c(sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="1. Coarse & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$SHR7name=="2. Loamy & no restrictions")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")], na.rm=TRUE), sd(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="4. Loamy & restrictive layers")], na.rm=TRUE)), n=c(sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="1. Coarse & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='dryfarm' & pts_50cm$SHR7name=="2. Loamy & no restrictions")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="3. Coarse-loamy & restrictive layers")])), sum(!is.na(pts_50cm$kgOrg.m2_50cm[which(pts_50cm$tillage=='no till' & pts_50cm$irrigated_vs_dryfarm=='irrigated' & pts_50cm$SHR7name=="4. Loamy & restrictive layers")]))), row.names = c('1. Coarse w/no restrictions', '2a. Loamy w/no resrictions', '2b. Loamy w/no restrictions: No-till', '2c. Loamy w/no restrictions: Dry farm', '3. Coarse w/pans', '4. Loamy w/pans: No-till'))
