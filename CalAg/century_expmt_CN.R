dataDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/citations/century experiment/RR_2014_datasets'
list.files(dataDir)


ltras_cn <- read.csv(file.path(dataDir, 'soil_total_CN_all.csv'), stringsAsFactors = FALSE)
colnames(ltras_cn)
unique(ltras_cn$date)
unique(ltras_cn$soil_location_description)
sum(is.na(ltras_cn$soil_location_description))
unique(ltras_cn$plot)[order(unique(ltras_cn$plot))]
ltras_cn$depth_class <- paste0(ltras_cn$depth_top_cm, '_', ltras_cn$depth_bottom_cm)

ltras_cn_1993 <- ltras_cn[ltras_cn$soil_location_description=='whole_plot_compilation' & ltras_cn$date=='9/1/1993',]
unique(ltras_cn_1993$plot)[order(unique(ltras_cn_1993$plot))] #has e, w, and & compilation per plot
length(unique(ltras_cn_1993$plot)) #63 plots
unique(ltras_cn_1993$depth_class)
table(ltras_cn_1993$system_name[ltras_cn_1993$depth_class=='0_15'])

ltras_cn_1997 <- ltras_cn[ltras_cn$date=='9/1/1997', ] #ltras_cn$soil_location_description=='whole_plot_compilation' & 
unique(ltras_cn_1997$depth_class)
unique(ltras_cn_1997$plot)[order(unique(ltras_cn_1997$plot))] #only has compilations
length(unique(ltras_cn_1997$plot)) #62 plots

ltras_cn_2003 <- ltras_cn[ ltras_cn$date=='9/1/2003', ] #ltras_cn$soil_location_description=='whole_plot_compilation' &
unique(ltras_cn_2003$depth_class)
unique(ltras_cn_2003$plot)[order(unique(ltras_cn_2003$plot))] #very inconsistent dataset
length(unique(ltras_cn_2003$plot)) #108 e & w subplots
ltras_cn_2003$plot_overall <- substr(ltras_cn_2003$plot, start = 1, stop = 3)

ltras_cn_2012 <- ltras_cn[ltras_cn$date=='10/8/2012',] #ltras_cn$soil_location_description=='whole_plot_compilation' & 
unique(ltras_cn_2012$plot)[order(unique(ltras_cn_2012$plot))]
length(unique(ltras_cn_2012$plot)) #54 plots
ltras_cn_2012$plot[ltras_cn_2012$depth_class=='0_15']
ltras_cn_2012$plot[ltras_cn_2012$depth_class=='15_30']
ltras_cn_2012$plot[ltras_cn_2012$depth_class=='30_60']
ltras_cn_2012$plot[ltras_cn_2012$depth_class=='60_100']
ltras_cn_2012$plot[ltras_cn_2012$depth_class=='100_200']
table(ltras_cn_2012$depth_class)

soilC_0_15_1993 <- data.frame(plot=names(tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='0_15'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='0_15'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='0_15'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='0_15'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_1993)
soilC_0_15_1993$system_name <- ltras_cn_1993$system_name[match(soilC_0_15_1993$plot, ltras_cn_1993$plot)]
unique(soilC_0_15_1993$plot)[order(unique(soilC_0_15_1993$plot))]
table(soilC_0_15_1993$system_name)

soilC_0_15_1997 <- data.frame(plot=names(tapply(ltras_cn_1997$total_carbon_.[ltras_cn_1997$depth_class=='0_15'], ltras_cn_1997$plot[ltras_cn_1997$depth_class=='0_15'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1997$total_carbon_.[ltras_cn_1997$depth_class=='0_15'], ltras_cn_1997$plot[ltras_cn_1997$depth_class=='0_15'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_1997)
soilC_0_15_1997$system_name <- ltras_cn_1997$system_name[match(soilC_0_15_1997$plot, ltras_cn_1997$plot)]
unique(soilC_0_15_1997$plot)[order(unique(soilC_0_15_1997$plot))]
table(soilC_0_15_1997$system_name) #same as 1993 except rainfed wheat has one less

soilC_0_15_2003 <- data.frame(plot=names(tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='0_15'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='0_15'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='0_15'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='0_15'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_2003)
soilC_0_15_2003$system_name <- ltras_cn_2003$system_name[match(soilC_0_15_2003$plot, ltras_cn_2003$plot_overall)]
unique(soilC_0_15_2003$plot)[order(unique(soilC_0_15_2003$plot))]
table(soilC_0_15_2003$system_name)

soilC_0_15_2012 <- data.frame(plot=names(tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='0_15'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='0_15'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='0_15'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='0_15'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_0_15_2012)
soilC_0_15_2012$system_name <- ltras_cn_2012$system_name[match(soilC_0_15_2012$plot, ltras_cn_2012$plot)]
unique(soilC_0_15_2012$plot)[order(unique(soilC_0_15_2012$plot))]
table(soilC_0_15_2012$system_name)

soilC_0_15 <- soilC_0_15_1993
colnames(soilC_0_15)[2] <- 'SOC_1993'
soilC_0_15 <- soilC_0_15[,c(1,3,2)]
soilC_0_15$SOC_1997 <- soilC_0_15_1997$soilC[match(soilC_0_15$plot, soilC_0_15_1997$plot)]
soilC_0_15$SOC_2003 <- soilC_0_15_2003$soilC[match(soilC_0_15$plot, soilC_0_15_2003$plot)]
soilC_0_15$SOC_2012 <- soilC_0_15_2012$soilC[match(soilC_0_15$plot, soilC_0_15_2012$plot)]
dim(soilC_0_15)
soilC_0_15
soilC_0_15[soilC_0_15$system_name=='Legume/Corn/Tomato',]
soilC_0_15[soilC_0_15$system_name=='Organic Corn/Tomato',]
soilC_0_15[soilC_0_15$system_name=='Conventional Corn/Tomato',]

#15-30 cm SOC by date
soilC_15_30_1993 <- data.frame(plot=names(tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='15_30'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='15_30'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='15_30'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='15_30'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_1993)
soilC_15_30_1993$system_name <- ltras_cn_1993$system_name[match(soilC_15_30_1993$plot, ltras_cn_1993$plot)]
unique(soilC_15_30_1993$plot)[order(unique(soilC_15_30_1993$plot))]
table(soilC_15_30_1993$system_name)

soilC_15_30_1997 <- data.frame(plot=names(tapply(ltras_cn_1997$total_carbon_.[ltras_cn_1997$depth_class=='15_30'], ltras_cn_1997$plot[ltras_cn_1997$depth_class=='15_30'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1997$total_carbon_.[ltras_cn_1997$depth_class=='15_30'], ltras_cn_1997$plot[ltras_cn_1997$depth_class=='15_30'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_1997)
soilC_15_30_1997$system_name <- ltras_cn_1997$system_name[match(soilC_15_30_1997$plot, ltras_cn_1997$plot)]
unique(soilC_15_30_1997$plot)[order(unique(soilC_15_30_1997$plot))]
table(soilC_15_30_1997$system_name) #same as 1993 except rainfed wheat has one less

soilC_15_30_2003 <- data.frame(plot=names(tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='15_30'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='15_30'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='15_30'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='15_30'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_2003)
soilC_15_30_2003$system_name <- ltras_cn_2003$system_name[match(soilC_15_30_2003$plot, ltras_cn_2003$plot_overall)]
unique(soilC_15_30_2003$plot)[order(unique(soilC_15_30_2003$plot))]
table(soilC_15_30_2003$system_name)

soilC_15_30_2012 <- data.frame(plot=names(tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='15_30'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='15_30'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='15_30'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='15_30'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_15_30_2012)
soilC_15_30_2012$system_name <- ltras_cn_2012$system_name[match(soilC_15_30_2012$plot, ltras_cn_2012$plot)]
unique(soilC_15_30_2012$plot)[order(unique(soilC_15_30_2012$plot))]
table(soilC_15_30_2012$system_name)

soilC_15_30 <- soilC_15_30_1993
colnames(soilC_15_30)[2] <- 'SOC_1993'
soilC_15_30 <- soilC_15_30[,c(1,3,2)]
soilC_15_30$SOC_1997 <- soilC_15_30_1997$soilC[match(soilC_15_30$plot, soilC_15_30_1997$plot)]
soilC_15_30$SOC_2003 <- soilC_15_30_2003$soilC[match(soilC_15_30$plot, soilC_15_30_2003$plot)]
soilC_15_30$SOC_2012 <- soilC_15_30_2012$soilC[match(soilC_15_30$plot, soilC_15_30_2012$plot)]
dim(soilC_15_30)
soilC_15_30
soilC_15_30[soilC_15_30$system_name=='Legume/Corn/Tomato',]
soilC_15_30[soilC_15_30$system_name=='Organic Corn/Tomato',]
soilC_15_30[soilC_15_30$system_name=='Conventional Corn/Tomato',]

#30-60 cm SOC by date
soilC_30_60_1993 <- data.frame(plot=names(tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='30_60'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='30_60'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='30_60'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='30_60'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_1993)
soilC_30_60_1993$system_name <- ltras_cn_1993$system_name[match(soilC_30_60_1993$plot, ltras_cn_1993$plot)]
unique(soilC_30_60_1993$plot)[order(unique(soilC_30_60_1993$plot))]
table(soilC_30_60_1993$system_name)

#no data in 1997

soilC_30_60_2003 <- data.frame(plot=names(tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='30_60'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='30_60'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='30_60'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='30_60'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_2003)
soilC_30_60_2003$system_name <- ltras_cn_2003$system_name[match(soilC_30_60_2003$plot, ltras_cn_2003$plot_overall)]
unique(soilC_30_60_2003$plot)[order(unique(soilC_30_60_2003$plot))]
table(soilC_30_60_2003$system_name)

soilC_30_60_2012 <- data.frame(plot=names(tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='30_60'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='30_60'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='30_60'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='30_60'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_30_60_2012)
soilC_30_60_2012$system_name <- ltras_cn_2012$system_name[match(soilC_30_60_2012$plot, ltras_cn_2012$plot)]
unique(soilC_30_60_2012$plot)[order(unique(soilC_30_60_2012$plot))]
table(soilC_30_60_2012$system_name)

soilC_30_60 <- soilC_30_60_2003
colnames(soilC_30_60)[2] <- 'SOC_2003'
soilC_30_60$SOC_1993 <- soilC_30_60_1993$soilC[match(soilC_30_60$plot, soilC_30_60_1993$plot)]

colnames(soilC_30_60)
soilC_30_60 <- soilC_30_60[,c(1,3,4,2)]
soilC_30_60$SOC_2012 <- soilC_30_60_2012$soilC[match(soilC_30_60$plot, soilC_30_60_2012$plot)]
dim(soilC_30_60)
soilC_30_60
soilC_30_60[soilC_30_60$system_name=='Legume/Corn/Tomato',]
soilC_30_60[soilC_30_60$system_name=='Organic Corn/Tomato',]
soilC_30_60[soilC_30_60$system_name=='Conventional Corn/Tomato',]

#60-100 cm SOC by date
soilC_60_100_1993 <- data.frame(plot=names(tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='60_100'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='60_100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_1993$total_carbon_.[ltras_cn_1993$depth_class=='60_100'], ltras_cn_1993$plot[ltras_cn_1993$depth_class=='60_100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_1993)
soilC_60_100_1993$system_name <- ltras_cn_1993$system_name[match(soilC_60_100_1993$plot, ltras_cn_1993$plot)]
unique(soilC_60_100_1993$plot)[order(unique(soilC_60_100_1993$plot))]
table(soilC_60_100_1993$system_name)

#no data in 1997

soilC_60_100_2003 <- data.frame(plot=names(tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='60_100'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='60_100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2003$total_carbon_.[ltras_cn_2003$depth_class=='60_100'], ltras_cn_2003$plot_overall[ltras_cn_2003$depth_class=='60_100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_2003)
soilC_60_100_2003$system_name <- ltras_cn_2003$system_name[match(soilC_60_100_2003$plot, ltras_cn_2003$plot_overall)]
unique(soilC_60_100_2003$plot)[order(unique(soilC_60_100_2003$plot))]
table(soilC_60_100_2003$system_name)

soilC_60_100_2012 <- data.frame(plot=names(tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='60_100'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='60_100'], mean, na.rm=TRUE)), soilC=tapply(ltras_cn_2012$total_carbon_.[ltras_cn_2012$depth_class=='60_100'], ltras_cn_2012$plot[ltras_cn_2012$depth_class=='60_100'], mean, na.rm=TRUE), stringsAsFactors = FALSE)
dim(soilC_60_100_2012)
soilC_60_100_2012$system_name <- ltras_cn_2012$system_name[match(soilC_60_100_2012$plot, ltras_cn_2012$plot)]
unique(soilC_60_100_2012$plot)[order(unique(soilC_60_100_2012$plot))]
table(soilC_60_100_2012$system_name)

soilC_60_100 <- soilC_60_100_2003
colnames(soilC_60_100)[2] <- 'SOC_2003'
soilC_60_100$SOC_1993 <- soilC_60_100_1993$soilC[match(soilC_60_100$plot, soilC_60_100_1993$plot)]

colnames(soilC_60_100)
soilC_60_100 <- soilC_60_100[,c(1,3,4,2)]
soilC_60_100$SOC_2012 <- soilC_60_100_2012$soilC[match(soilC_60_100$plot, soilC_60_100_2012$plot)]
dim(soilC_60_100)
soilC_60_100
soilC_60_100[soilC_60_100$system_name=='Legume/Corn/Tomato',]
soilC_60_100[soilC_60_100$system_name=='Organic Corn/Tomato',]
soilC_60_100[soilC_60_100$system_name=='Conventional Corn/Tomato',]

#look at 1993 to 2012 SOC changes in Legume/Corn/Tomato (aka Conv-WCC) reported by Tautges et al.
soilC_0_15_1993[soilC_0_15_1993$system_name=='Legume/Corn/Tomato', ]
ltras_cn_1997[ltras_cn_1997$system_name=='Legume/Corn/Tomato' & ltras_cn_1997$depth_class=='0_15',]
ltras_cn_2003[ltras_cn_2003$system_name=='Legume/Corn/Tomato' & ltras_cn_2003$depth_class=='0_15',]
soilC_0_15_2012[soilC_0_15_2012$system_name=='Legume/Corn/Tomato',]
mean(soilC_0_15_1993$soilC[soilC_0_15_1993$system_name=='Legume/Corn/Tomato']) - mean(soilC_0_15_2012$soilC[soilC_0_15_2012$system_name=='Legume/Corn/Tomato'])

ltras_cn_1993[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='15_30',]
ltras_cn_1997[ltras_cn_1997$system_name=='Legume/Corn/Tomato' & ltras_cn_1997$depth_class=='15_30',]
ltras_cn_2003[ltras_cn_2003$system_name=='Legume/Corn/Tomato' & ltras_cn_2003$depth_class=='15_30',]
ltras_cn_2012[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='15_30',]
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='15_30']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='15_30'])

ltras_cn_1993[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='30_60',]
ltras_cn_1997[ltras_cn_1997$system_name=='Legume/Corn/Tomato' & ltras_cn_1997$depth_class=='30_60',]
ltras_cn_2003[ltras_cn_2003$system_name=='Legume/Corn/Tomato' & ltras_cn_2003$depth_class=='30_60',]
ltras_cn_2012[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='30_60',]
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='30_60']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='30_60'])

ltras_cn_1993[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='60_100',]
ltras_cn_1997[ltras_cn_1997$system_name=='Legume/Corn/Tomato' & ltras_cn_1997$depth_class=='60_100',]
ltras_cn_2003[ltras_cn_2003$system_name=='Legume/Corn/Tomato' & ltras_cn_2003$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$plot=='2-4',] #60-100 cm data is big drop off from 30-60 cm and 250-300 cm is very high
ltras_cn_1993[ltras_cn_1993$plot=='2-4',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='4-4',] #250-300 cm very high
ltras_cn_1993[ltras_cn_1993$plot=='4-4',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='6-9',] #250-300 cm VERY high
ltras_cn_1993[ltras_cn_1993$plot=='6-9',] #looks ok
10*(mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='60_100']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='60_100'])) #declined by 1.09 g kg^-1 but plot 2-4 change is VERY suspicious

ltras_cn_1993[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='100_200',]
ltras_cn_1997[ltras_cn_1997$system_name=='Legume/Corn/Tomato' & ltras_cn_1997$depth_class=='100_200',]
ltras_cn_2003[ltras_cn_2003$system_name=='Legume/Corn/Tomato' & ltras_cn_2003$depth_class=='100_200',]
ltras_cn_2012[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='100_200',]
10*(mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Legume/Corn/Tomato' & ltras_cn_1993$depth_class=='100_200']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Legume/Corn/Tomato' & ltras_cn_2012$depth_class=='100_200'])) #0.23 g kg^-1 decline less than that reported by Tautges et al.


#same for Corn/Tomato (aka Conv)
ltras_cn_1993[ltras_cn_1993$system_name=='Conventional Corn/Tomato' & ltras_cn_1993$depth_class=='60_100',]
ltras_cn_1997[ltras_cn_1997$system_name=='Conventional Corn/Tomato' & ltras_cn_1997$depth_class=='60_100',]
ltras_cn_2003[ltras_cn_2003$system_name=='Conventional Corn/Tomato' & ltras_cn_2003$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$system_name=='Conventional Corn/Tomato' & ltras_cn_2012$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$plot=='1-4',] #250-300 cm data VERY high
ltras_cn_1993[ltras_cn_1993$plot=='1-4',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='5-5',] #60-100 cm greater than 30-60 cm
ltras_cn_1993[ltras_cn_1993$plot=='5-5',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='7-8',] #looks ok
ltras_cn_1993[ltras_cn_1993$plot=='7-8',] #looks ok
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Conventional Corn/Tomato' & ltras_cn_1993$depth_class=='60_100']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Conventional Corn/Tomato' & ltras_cn_2012$depth_class=='60_100']) #increased by 0.038 g kg^-1

ltras_cn_1993[ltras_cn_1993$system_name=='Conventional Corn/Tomato' & ltras_cn_1993$depth_class=='100_200',]
ltras_cn_1997[ltras_cn_1997$system_name=='Conventional Corn/Tomato' & ltras_cn_1997$depth_class=='100_200',]
ltras_cn_2003[ltras_cn_2003$system_name=='Conventional Corn/Tomato' & ltras_cn_2003$depth_class=='100_200',]
ltras_cn_2012[ltras_cn_2012$system_name=='Conventional Corn/Tomato' & ltras_cn_2012$depth_class=='100_200',]
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Conventional Corn/Tomato' & ltras_cn_1993$depth_class=='100_200']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Conventional Corn/Tomato' & ltras_cn_2012$depth_class=='100_200']) #0.031 g kg^-1

#same for Organic Corn/Tomato
ltras_cn_1993[ltras_cn_1993$system_name=='Organic Corn/Tomato' & ltras_cn_1993$depth_class=='60_100',] #only 2 plots
ltras_cn_1997[ltras_cn_1997$system_name=='Organic Corn/Tomato' & ltras_cn_1997$depth_class=='60_100',]
ltras_cn_2003[ltras_cn_2003$system_name=='Organic Corn/Tomato' & ltras_cn_2003$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$system_name=='Organic Corn/Tomato' & ltras_cn_2012$depth_class=='60_100',]
ltras_cn_2012[ltras_cn_2012$plot=='1-4',] #250-300 cm data VERY high
ltras_cn_1993[ltras_cn_1993$plot=='1-4',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='5-5',] #60-100 cm greater than 30-60 cm
ltras_cn_1993[ltras_cn_1993$plot=='5-5',] #looks ok
ltras_cn_2012[ltras_cn_2012$plot=='7-8',] #looks ok
ltras_cn_1993[ltras_cn_1993$plot=='7-8',] #looks ok
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Organic Corn/Tomato' & ltras_cn_1993$depth_class=='60_100']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Organic Corn/Tomato' & ltras_cn_2012$depth_class=='60_100']) #increased by 0.038 g kg^-1

ltras_cn_1993[ltras_cn_1993$system_name=='Organic Corn/Tomato' & ltras_cn_1993$depth_class=='100_200',]
ltras_cn_1997[ltras_cn_1997$system_name=='Organic Corn/Tomato' & ltras_cn_1997$depth_class=='100_200',]
ltras_cn_2003[ltras_cn_2003$system_name=='Organic Corn/Tomato' & ltras_cn_2003$depth_class=='100_200',]
ltras_cn_2012[ltras_cn_2012$system_name=='Organic Corn/Tomato' & ltras_cn_2012$depth_class=='100_200',]
mean(ltras_cn_1993$total_carbon_.[ltras_cn_1993$system_name=='Organic Corn/Tomato' & ltras_cn_1993$depth_class=='100_200']) - mean(ltras_cn_2012$total_carbon_.[ltras_cn_2012$system_name=='Organic Corn/Tomato' & ltras_cn_2012$depth_class=='100_200']) #0.031 g kg^-1