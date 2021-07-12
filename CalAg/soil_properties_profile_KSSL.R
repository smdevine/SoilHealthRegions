SOC_to_SOM <- 1.72
library(aqp)
library(soilDB)
library(lattice)
library(extrafont)
library(extrafontdb)
loadfonts(device = 'win')
ksslDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/kssl'
FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/soil health/publication/California Agriculture/Figures'
# Newville_data <- fetchKSSL(series='Newville')
# Newville_data
# site(Newville_data)
clus_7_names <- c('6. Fine salt-affected', '3. Low OM with restrictive horizons', '4. High OM with restrictive horizons', '1. Coarse with no restrictions', '2. Loamy with no restrictions', '7. Shrink-swell', '5. Coarse-loamy salt-affected')
clus_7_colors <- c('deepskyblue', 'olivedrab3', 'firebrick3', 'lightgoldenrod', 'tan4', 'violetred', 'lightblue1')
kssl_horizons_SHR <- read.csv(file.path(ksslDir, 'kssl_horizons_SHRonly.csv'), stringsAsFactors = FALSE)
colnames(kssl_horizons_SHR)
summary(kssl_horizons_SHR$oc) #497 NAs
summary(kssl_horizons_SHR$oc_est) #this was was produced in kssl_validation_FINAL.R
kssl_horizons_SHR$SOM <- kssl_horizons_SHR$oc_est*SOC_to_SOM
summary(kssl_horizons_SHR$SOM)
kssl_horizons_SHR$SHR7name <- clus_7_names[kssl_horizons_SHR$SHR7code]
kssl_horizons_SHR$SHR7name <- as.factor(kssl_horizons_SHR$SHR7name)
depths(kssl_horizons_SHR) <- pedon_key ~ hzn_top + hzn_bot
site(kssl_horizons_SHR) <- ~ SHR7name
kssl_horizons_SHR.slab <- slab(kssl_horizons_SHR, SHR7name ~ clay + ph_h2o + ec_12pre + SOM, slab.structure = 1)
dim(kssl_horizons_SHR.slab)
str(kssl_horizons_SHR.slab)
levels(kssl_horizons_SHR.slab$variable)
levels(kssl_horizons_SHR.slab$variable) <- c('Clay (%)', 'Soil pH', 'Soil EC', 'Soil organic matter (%)')
tps <- list(superpose.line=list(col=clus_7_colors[order(clus_7_names)], lwd=2))

tiff(file = file.path(FiguresDir, 'soil_properties_KSSL_profiles_revised.tif'), family = 'Times New Roman', width = 9, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
xyplot(top ~ p.q50 | variable, groups=SHR7name, data=kssl_horizons_SHR.slab, ylab='Depth (cm)', xlab='KSSL median bounded by 25th and 75th percentiles', lower=kssl_horizons_SHR.slab$p.q25, upper=kssl_horizons_SHR.slab$p.q75, ylim=c(155,-5), xlim = list(c(2,55), c(5.5, 8.9), c(0,12), c(0,3.5)), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(4,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE, lwd=2))
dev.off()
#key=list(text=list(levels(dframe$Z)), space='top', points=list(pch=1:nlevels(dframe$Z), col=col), lines=list(col=col), columns=nlevels(dframe$Z))

#make eps version
setEPS()
postscript(file = file.path(FiguresDir, 'EPS versions', 'Figure_2_FINAL.eps'), pointsize = 12, family = 'Times New Roman', width = 9, height = 5)
xyplot(top ~ p.q50 | variable, groups=SHR7name, data=kssl_horizons_SHR.slab, ylab='Depth (cm)', xlab='KSSL median bounded by 25th and 75th percentiles', lower=kssl_horizons_SHR.slab$p.q25, upper=kssl_horizons_SHR.slab$p.q75, ylim=c(155,-5), xlim = list(c(2,55), c(5.5, 8.9), c(0,12), c(0,3.5)), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(4,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE, lwd=2))
dev.off()

#make pdf version
cairo_pdf(file = file.path(FiguresDir, 'PDF versions', 'Figure_2_FINAL.pdf'), pointsize = 12, family = 'Times New Roman', width = 9, height = 5)
xyplot(top ~ p.q50 | variable, groups=SHR7name, data=kssl_horizons_SHR.slab, ylab='Depth (cm)', xlab='KSSL median bounded by 25th and 75th percentiles', lower=kssl_horizons_SHR.slab$p.q25, upper=kssl_horizons_SHR.slab$p.q75, ylim=c(155,-5), xlim = list(c(2,55), c(5.5, 8.9), c(0,12), c(0,3.5)), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(4,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE, lwd=2))
dev.off()

kssl_horizons_SHR.slab_v2 <- slab(kssl_horizons_SHR, SHR7name ~ clay + ph_h2o + SOM, slab.structure = 1)
levels(kssl_horizons_SHR.slab_v2$variable) <- c('Clay (%)', 'Soil pH', 'Soil organic matter (%)')
tiff(file = file.path(FiguresDir, 'soil_properties_KSSL_profiles_revised_v2.tif'), family = 'Times New Roman', width = 9, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
xyplot(top ~ p.q50 | variable, groups=SHR7name, data=kssl_horizons_SHR.slab_v2, ylab='Depth (cm)', xlab='KSSL data median bounded by 25th and 75th percentiles', lower=kssl_horizons_SHR.slab_v2$p.q25, upper=kssl_horizons_SHR.slab_v2$p.q75, ylim=c(155,-5), xlim = list(c(2,55), c(5.5, 8.9), c(0,4)), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(3,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE))
dev.off()

mean.and.sd <- function(values) {
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  upper <- m + s
  lower <- m - s
  res <- c(mean=m, lower=lower, upper=upper)
  return(res)
}
kssl_horizons_SHR.slab_v3 <- slab(kssl_horizons_SHR, SHR7name ~ clay + ph_h2o + oc_est, slab.structure = 1, slab.fun = mean.and.sd)
levels(kssl_horizons_SHR.slab_v3$variable) <- c('Clay (%)', 'Soil pH', 'Soil organic carbon (%)')
tiff(file = file.path(FiguresDir, 'CalAg', 'KSSL', 'soil_properties_KSSL_profiles_v3.tif'), family = 'Times New Roman', width = 9, height = 5, pointsize = 12, units = 'in', res=800, compression='lzw')
xyplot(top ~ mean | variable, groups=SHR7name, data=kssl_horizons_SHR.slab_v3, ylab='Depth (cm)', xlab='KSSL validation data mean \u00B1 1 standard deviation', lower=kssl_horizons_SHR.slab_v3$lower, upper=kssl_horizons_SHR.slab_v3$upper, ylim=c(155,-5), xlim = list(c(2,55), c(5.5, 8.9), c(0,2)), panel=panel.depth_function, alpha=0.4, sync.colors=TRUE, prepanel=prepanel.depth_function, par.strip.text=list(cex=0.8), strip=strip.custom(bg=grey(0.85)), layout=c(3,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=3, lines=TRUE, points=FALSE))
dev.off()
