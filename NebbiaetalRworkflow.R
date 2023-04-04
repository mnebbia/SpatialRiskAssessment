## Code for the paper "Spatial risk assessment and the protection of cultural heritage in southern Tajikistan" 
# Authors: Marco Nebbia, Federica Cilio, Bobomullo Bobomulloyev
# Date: 08/09/2020

## Libraries
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(spatstat)
#library(sparr)
library(viridis)

## Set up
#setwd("~/....../NebbiaetalSuppMaterial") # save and exctract the folder on local drive and set it as working directory

## Read in 
region <- readOGR("shp/tj_region.shp", layer="tj_region")
elev <- raster("raster/tj_elev.grd")
rivers <- readOGR("shp/rivers_tj.shp", layer="rivers_tj")
high.riv <- readOGR("shp/high_riv.shp", layer="high_riv") # rivers with strahler order >= 3
tj.sites <- readOGR("shp/sites.shp", layer="sites")
towns <- readOGR("shp/tj_towns.shp", layer="tj_towns")
sites.eros <- subset(tj.sites, rivErosion > 0) # sites affected by riverine erosion
sites.hum <- subset(tj.sites, urban > 0 | construct > 0 | agricult > 0) # sites affected by human activities
logAll <- readOGR("shp/logAllCov.shp", layer="logAllCov") # all sites + random non-site sample
mod <- raster("raster/mod.grd") # modified catchment areas

## Read in covariates
twi <- raster("raster/twi.grd") # topographic wetness index calculated with SAGA-GIS
towndist <- raster("raster/tj_towndist.tif") # distance from modern settlements
slope <- raster("raster/tj_slope.tif")
hirivdist <- raster("raster/tj_rivdist.tif") #distance from high strahler order (>=3) rivers

## Convert to spastat formats
tjsitesppp <- ppp(x=coordinates(tj.sites)[,1],y=coordinates(tj.sites)[,2], window=as.owin(region))
erosppp <- ppp(x=coordinates(sites.eros)[,1],y=coordinates(sites.eros)[,2], window=as.owin(region))
humppp <- ppp(x=coordinates(sites.hum)[,1],y=coordinates(sites.hum)[,2], window=as.owin(region))
twi.im <- as.im(twi)
slope.im <- as.im(slope)
towndist.im <- as.im(towndist)
hirivdist.im <- as.im(hirivdist)
riverpsp <- as.psp(rivers, window=as.owin(region))
rivdist <- distmap(riverpsp) # distance from ALL rivers
rivdist.im <- as.im(rivdist)

## First order intensities all sites
twi.rho <- rhohat(tjsitesppp, twi.im, confidence=0.95)
towndist.rho <- rhohat(tjsitesppp, towndist.im, confidence=0.95)
slope.rho <- rhohat(tjsitesppp, slope.im, confidence=0.95)
hirivdist.rho <- rhohat(tjsitesppp, hirivdist.im, confidence=0.95)
rivdist.rho <- rhohat(tjsitesppp, rivdist.im, confidence=0.95)

## First order logistic model GLM
#all sites
glmAll <- glm(class ~ rivdist + slope + towndist + twi, data=logAll, family=binomial)
summary(glmAll)
predAll <- 3.881e+00 + (hirivdist.im*-3.331e-05) + (slope.im*-1.372e-01) + (towndist.im*7.985e-05) + (twi.im*-2.218e-01)
#erosion sites - all other site non affected by erosion as non-site locations
#predEreg <- -3.073e+00 + (hirivdist.im*2.110e-04) + (towndist.im*8.490e-05)
glmE <- glm(class_E ~ rivdist + slope + towndist + twi, data=tj.sites, family=binomial)
summary(glmE)
predE <- -8.996e+00+(twi.im*2.963e-01)
#human activities sites - all other site non affected by human activities as non-site locations
glmH <- glm(class_H ~ rivdist + slope + towndist + twi, data=tj.sites, family=binomial) 
summary(glmH)
predH <- -1.379e+00 +(hirivdist.im*-1.010e-04)+(slope.im*-1.044e-01)+(towndist.im*4.294e-05)

## Relative risk surfaces 
#bandwidth selection
#bw.sites <- bw.ppl(tjsitesppp)
#bw.sites # sigma=912.9484
#bw.eros <- bw.ppl(erosppp)
#bw.eros # sigma=5121.654
#bw.hum <- bw.ppl(humppp)
#bw.hum # sigma=1357.358
bw <- 5121.654 # coarser sigma of the three ppl adjustments
cellres <- 100
sitesdens <- density(tjsitesppp, sigma=bw, eps=cellres, edge=F)
erosdens <- density(erosppp, sigma=bw, eps=cellres, edge=F)
humdens <- density(humppp, sigma=bw, eps=cellres, edge=F)
erosRR <- erosdens/sitesdens *100
humRR <- humdens/sitesdens *100
ramperos <- plot(erosRR, main="", box=FALSE, do.plot=FALSE) # save colour ramp
dev.off()
ramphum <- plot(humRR, main="", box=FALSE, do.plot=FALSE) # save colour ramp
dev.off()

## Relative risk with tolerance p-value contours
hum_rr <- risk(humppp, tjsitesppp, h0=bw, adapt=TRUE, tolerate=TRUE, pilot.symmetry="pooled", davies.baddeley=0.05)
eros_rr <- risk(erosppp, tjsitesppp, h0=bw, adapt=TRUE, tolerate=TRUE, pilot.symmetry="pooled", davies.baddeley=0.05)

## Plot site types
pdf(paper = "a4", file="outputs/figA1.pdf")
par(mar=c(13,3,3,3))
plot(tj.sites$MonType, las=3)
dev.off()

## Plot covariates
pdf(height=9, width = 12, paper="a4r", file ="outputs/fig6.pdf")
par(mfrow=c(2,2), mar=c(1,1,1,8))
plot(twi, box=F, axes=F, main="a")
plot(towndist, box=F, axes=F, main="c")
plot(slope, box=F, axes=F, main="b")
plot(hirivdist, box=F, axes=F, main="d")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
dev.off()

## Plot first order intensities
pdf(height=11, width = 9, file="outputs/fig7.pdf")
par(mfrow=c(3,2))
plot(twi.rho, legend=F, main="a.Topographic wetness index", xlab="index", ylab="")
plot(slope.rho, legend=F, main="b.Slope", xlab="degrees", ylab="")
plot(towndist.rho, legend=F, main="c.Distance from towns", xlab="metres", ylab="")
plot(hirivdist.rho, legend=F, main="d.Distance from high order rivers", xlab="metres", ylab="")
plot(rivdist.rho, legend=F, main="e.Distance from rivers", xlab="metres", ylab="")
dev.off()

## Plot first order logistic model (all sites)
pdf(height=7, width=9, paper = "a4r", file="outputs/fig8.pdf")
plot(predAll, box=F, main="")
points(tjsitesppp, pch=19, cex=0.1, col="black")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
dev.off()

## Plot first order logistic model (sites affected by erosion)
pdf(height=11, width=9, file="outputs/figA4.pdf")
par(mfrow=c(2,1), mar=c(1,1,1,1))
plot(predE, box=F, main="")
points(erosppp, pch=19, cex=0.3, col="black")
plot(predH, box=F, main="")
points(humppp, pch=19, cex=0.3, col="black")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
dev.off()

## Plot kde and rel risk erosion
pdf(height=11, width=9, file="outputs/figA5_1.pdf")
par(mfrow=c(2,1), mar=c(1,1,1,1))
plot(erosdens, main="", box=F)
points(tj.sites, pch=19, cex=0.1, col="grey")
points(sites.eros, pch=1, cex=1, col="green")
plot(erosRR, main="", box=F, col=viridis)
points(tj.sites, pch=19, cex=0.1, col="grey")
points(sites.eros, pch=1, cex=1, col="green")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
dev.off()

## Plot erosRR vs mod catchment areas
erosRR[as.matrix(erosRR) < 20] <- NA # truncating values below 20%
pdf(height=7, width=9, paper = "a4r", file="outputs/fig9.pdf")
plot(mod, box=F, axes=F, legend=F, col=gray.colors(n=20))
plot(region, lwd=0.5, add=T)
lines(high.riv, col="cyan", lwd=0.5)
plot(erosRR, main="", box=F, add=T)
points(tj.sites, pch=19, cex=0.2, col="green")
points(sites.eros, pch=1, cex=1, col="red")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
plot(ramperos, vertical=T, ylim=c(4100000, 4120000), xlim=c(990000, 995000), add=T, las=1)
text(992500, 4122000, labels="%")
legend(1017000, 4125000, legend = c("All sites", "Sites affected by erosion"), col = c("green", "red"), pch = c(19, 1), bty = "n")
dev.off()

## Plot kde and rel risk human activities
pdf(height=11, width=9, file="outputs/figA5_2.pdf")
par(mfrow=c(2,1), mar=c(1,1,1,1))
plot(humdens, main="", box=F)
points(tj.sites, pch=19, cex=0.1, col="grey")
points(sites.hum, pch=1, cex=1, col="cyan")
plot(humRR, main="", box=F, col=viridis)
points(tj.sites, pch=19, cex=0.1, col="grey")
points(sites.hum, pch=1, cex=1, col="cyan")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
dev.off()

## Plot humRR
humRR[as.matrix(humRR) < 20] <- NA #truncating values below 20%
pdf(height=7, width=9, paper = "a4r", file="outputs/fig10.pdf")
plot(elev, box=F, axes=F, col=gray.colors(7))
plot(region, lwd=0.5, add=T)
plot(humRR, main="", box=F, add=T)
points(tj.sites, pch=19, cex=0.2, col="green")
points(towns, pch=19, cex=0.3, col="black")
scalebar(50000, xy=c(1020000, 4100000), type="bar", divs=2, below="km", label=c(0,25,50))
plot(ramphum, vertical=T, ylim=c(4100000, 4120000), xlim=c(990000, 995000), add=T, las=1)
text(992500, 4122000, labels="%")
legend(1017000, 4125000, legend = c("Modern towns", "All sites"), col = c("black", "green"), pch = c(19, 19), bty = "n")
dev.off()

## Plot RelRisk with tolerance contour
pdf(height = 11, width = 7, paper = "a4", file = "outputs/fig11.pdf")
par(mfrow=c(2,1))
plot(eros_rr, xlab="Easting", ylab="Northing", main="Sites affected by riverine erosion")
points(erosppp, pch=19, cex=0.2, col="black")
plot(hum_rr, xlab="Easting", ylab="Northing", main="Sites affected by human activties")
points(humppp, pch=19, cex=0.2, col="black")
points(towns, pch=19, cex=0.3, col="white")
dev.off()
