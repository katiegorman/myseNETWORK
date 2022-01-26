

################################################################################################################################
################################################################################################################################
#                                       Indiana Bat Colony Foraging

# Code by Alex Silvis, 14 September 2012, using R version 2.15
# revised 27 September 2012
# revised 21 February 2013
# revised 11 April 2013
# revised 25 April 2013
# revised 27 September 2013

################################################################################################################################
################################################################################################################################

# set the working directory where all of the data files are located

setwd("C:/Users/silvis/My Documents/R working directory/andrew ibats")

# load data from previous run
load("ibat colony foraging.RData", envir=parent.frame())
library(adehabitatHR)
library(ks)
library(raster)

# read in datashet in csv format
# this should be the sequential list of locations with animal/group identifier
tnfd <- data.frame(read.table(file='ibatforagingtelemetry.csv', sep=',', header=TRUE))
# get the names of the variables in the data and check the first few lines of data
t(names(tnfd))
head(tnfd)
tail(tnfd)
dim(tnfd)

length(unique(tnfd$name[tnfd$year=='2009']))
length(unique(tnfd$name[tnfd$year=='2010']))

length(tnfd$name[tnfd$year=='2009'])
length(tnfd$name[tnfd$year=='2010'])

# getting a summary of the data for the animals/groups
summary(as.factor(tnfd$year))
# note that for kernel estimation here there must be > 4 locations for every animal/group

################################################################################################################################
################################################################################################################################

# some data processing to make sure locations data are in the correct/useful format

# combine the location data into a single dataframe; it is best to use UTM values
locs<-cbind(tnfd$x,tnfd$y)
plot(locs)

####
# center the locations for data archiving purposes while preserving privacy
mx<-mean(tnfd$x)
my<-mean(tnfd$y)
locs2<-cbind(tnfd$x-mx,tnfd$y-my)
plot(locs2)
head(locs2)
####

# turn the output into a data frame
locs<-data.frame(locs)
names(locs)<-c('UTM.X','UTM.Y')
# plot the points to make sure they look correct
plot(locs)

# make a spatial points dataframe - required for home range calculations
library(sp)
id<-as.factor(tnfd$year)
idsp <- data.frame(id)
coordinates(idsp) <- locs
class(idsp)
library(sp)
loc<-SpatialPoints(locs, CRS(as.character("+proj=utm +zone=17 ellps=WGS84", data=locs)))
plot(loc)
summary(idsp)

################################################################################################################################
################################################################################################################################

# create home ranges and utilization distributions for tracked animals

# load the R library that contains the functions used below
library(adehabitatHR)

# calculate home range area using the kernel method; manually select fixed or adaptive
ud <- kernelUD(idsp, h = "href", grid = 200, same4all = T, kern = c("bivnorm"), extent = 1)

# the UD of the animals/groups
#windows()
image(ud)
# calculation of the home range area at requested percentage values
ka<-kernel.area(ud, percent = seq(20, 95, by = 5),unin = c("m"),unout = c("ha"), standardize = FALSE)
ka
plot(ka)
# get home range area at a single percentage value- will be used to transfer home ranges to a GIS later
ver <- getverticeshr(ud, 50)
# plot the home ranges
plot(ver)

# export the home range shapefile created in the previous step so that it can be used in a GIS
library(maptools)
#writePolyShape(ver, "foraging 50% range reference method")
# the shapefiles will be written to the working directory set; i.e., where the datafile is

# make a pretty plot of the UD's
volud<-getvolumeUD(ud, standardize = F)
#mage(volud, col=heat.colors(1000))

# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 11, height = 6)
par(mfrow=c(1,2))
par(mar=c(1,1,2,1))

image(volud[[1]], col=heat.colors(1000),xlim=c(314500, 328418.1), ylim=c(4387085, 4399461))
# must add singly
xyzA <- as.image.SpatialGridDataFrame(volud[[1]])
contour(xyzA, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
box()
title(main = "2009")

image(volud[[2]], col=heat.colors(1000),xlim=c(314500, 328418.1), ylim=c(4387085, 4399461))
# must add singly
xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
contour(xyzB, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
box()
title(main = "2010")

# scale bar and label will be placed interactively
SpatialPolygonsRescale(layout.scale.bar(), offset = c(locator(1)),
                       scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)

text(locator(1), "5 km")

############################

# look for spatial similarity of the roosting areas between years at selected percent home range areas
# note that this is calculated from the points, not the calculated UD!

wro<-kerneloverlaphr(ud, meth="BA", conditional=TRUE, percent = 95) # entire home range
wro

cro<-kerneloverlaphr(ud, meth="BA", conditional=TRUE, percent = 50) # core home range
cro

################################################################################################################################
################################################################################################################################

# export to GIS

Gud<-estUDm2spixdf(volud)
image(Gud[1])
image(Gud[2])

library(GDAL)
#writeGDAL(Gud, 'foraging ud')

#hri<-c(25,50,75,95)
#for(i in 1:length(hri)){
#  hris[i]<-getverticeshr(ud,hri[i])
#}

library(maptools)
ver <- getverticeshr(ud, 25)
#writePolyShape(ver, "foraging 25% homerange")

################################################################################################################################
################################################################################################################################

# make some quality graphics for publication

#Arial <- Type1Font(family = "Arial", metrics = c("ArialMT.afm", "arial-BoldMT.afm", "Arial-ItalicMT.afm", 
#                                                 "Arial-BoldItalicMT.afm"))
#names(postscriptFonts())
#postscriptFonts(Arial=Arial)

#http://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.07/bin/gswin64c.exe")

library(extrafont)
font_import()
loadfonts(device="postscript")
fonts()

# for Windows - in each session
# adjust the path to match the installation of Ghostscript

##################

postscript(file="foraging uds.eps",family='sans',horizontal=FALSE, onefile=TRUE,
           paper="special",width=11,height=6)

par(mfrow=c(1,2))
par(mar=c(1,1,2,1))

image(volud[[1]], col=heat.colors(1000),xlim=c(314500, 328418.1), ylim=c(4387085, 4399461))
# must add singly
xyzA <- as.image.SpatialGridDataFrame(volud[[1]])
contour(xyzA, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
box()
title(main = "2009")

image(volud[[2]], col=heat.colors(1000),xlim=c(314500, 328418.1), ylim=c(4387085, 4399461))
# must add singly
xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
contour(xyzB, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
box()
title(main = "2010")

# scale bar and label will be placed interactively
SpatialPolygonsRescale(layout.scale.bar(), offset = c(locator(1)),
                       scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)

dev.off()
#embedFonts(file="test1.eps",format="pswrite", outfile="test2.eps")


################################################################################################################################
################################################################################################################################

# calculate foraging area UD using a different bandwidth method

t(names(tnfd))

y09 <- subset(tnfd, year == '2009', select=c(1:11))
dim(y09)
y09<-droplevels(y09)

y10 <- subset(tnfd, year == '2010', select=c(1:11))
dim(y10)
y10<-droplevels(y10)

locs09<-data.frame(y09$x, y09$y)
locs10<-data.frame(y10$x, y10$y)

plot(locs09)
plot(locs10)

library(ks)

#####################

## bivariate bandwidth selection 
##using the plug-in bandwidt selector of Wand, M.P. & Jones, M.C. (1994) Multivariate plugin bandwidth selection. Computational Statistics, 9, 97-116.

#####################

# 2009 foraging area

H.scv09 <- Hpi.diag(locs09, nstage=2, binned=FALSE)
fhat09 <- kde(locs09, H=H.scv09, compute.cont=TRUE)
names(fhat09)
plot(fhat09, drawpoints=TRUE, drawlabels=TRUE, col='black', lwd=3, cex=0.5)
plot(fhat09, display="persp", border=NA, shade=0.3)
plot(fhat09, display="image", col=rev(heat.colors(1000)))

x11()
plot(fhat09, display="filled.contour2", cont=seq(5,95,by=5),drawpoints=T,drawlabels=T, xlab="Easting", ylab="Northing")
plot(fhat09, display='slice', cont=c(95,75,50,25), approx.cont=T, drawlabels=T, drawpoints=F, add=T, pch=19, ptcol='black')

contourSizes(fhat09)
# area given in square meters
contourSizes(fhat09, cont=seq(5,95))
contourSizes(fhat09, cont=95)


#####################

# 2010 foraging area

H.scv10 <- Hpi.diag(locs10, nstage=2, binned=FALSE)
fhat10 <- kde(locs10, H=H.scv10, compute.cont=TRUE)
names(fhat10)
plot(fhat10, drawpoints=TRUE, drawlabels=TRUE, col='black', lwd=3, cex=0.5)
plot(fhat10, display="persp", border=NA, shade=0.3)
plot(fhat10, display="image", col=rev(heat.colors(1000)))

x11()
plot(fhat10, display="filled.contour2", cont=seq(5,95,by=5),drawpoints=T,drawlabels=T, xlab="Easting", ylab="Northing")
plot(fhat10, display='slice', cont=c(95,75,50,25), approx.cont=T, drawlabels=T, drawpoints=F, add=T, pch=19, ptcol='black')

contourSizes(fhat10)
# area given in square meters
contourSizes(fhat10, cont=seq(5,95))
contourSizes(fhat10, cont=95)

###########################################################################

library(raster)

ud09<-raster(fhat09)
ud09
ud10<-raster(fhat10)
ud10
ud09 <- as(ud09, 'SpatialPixelsDataFrame')
ud10 <- as(ud10, 'SpatialPixelsDataFrame')

image(ud09)
image(ud10)

### conversion to adehabitat format

fullgrid(ud09) <- FALSE
hli <- list(h = 1, meth="specified")
ud09 <- new("estUD", ud09)
slot(ud09, "h") <- hli
slot(ud09, "vol") <- FALSE

fullgrid(ud10) <- FALSE
hli2 <- list(h = 1, meth="specified")
ud10 <- new("estUD", ud10)
slot(ud10, "h") <- hli2
slot(ud10, "vol") <- FALSE

# combine the home ranges into a single estUDm object for further analysis
liud <- list(animal1=ud09, animal2=ud10)
class(liud) <- "estUDm"
liud
image(liud)

kal<-kernel.area(liud, percent = seq(20, 95, by = 5),unin = c("m"),unout = c("ha"), standardize = FALSE)
kal
plot(kal)
# get home range area at a single percentage value- will be used to transfer home ranges to a GIS later
verl <- getverticeshr(liud, 95)
# plot the home ranges
plot(verl)

# export the home range shapefile created in the previous step so that it can be used in a GIS
library(maptools)
#writePolyShape(ver, "foraging 95% range plug- method")
# the shapefiles will be written to the working directory set; i.e., where the datafile is

# pretty plot of the UD's
volud2<-getvolumeUD(liud, standardize = F)

# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 9, height = 6)
par(mfrow=c(1,2))
par(mar=c(3,1,2,1))

image(volud2[[1]], col=heat.colors(1000),xlim=c(315000, 325700), ylim=c(4385000, 4405500))
# must add singly
#xyzA <- as.image.SpatialGridDataFrame(volud[[1]])
#contour(xyzA, add=TRUE)
plot(fhat09, display='slice', cont=c(95,75,50,25), approx.cont=T, drawlabels=T, drawpoints=F, add=T, pch=19, ptcol='black')
box()
title(main = "2009")
#points(tnfd$x[tnfd$year=='2009'],tnfd$y[tnfd$year=='2009'])

image(volud2[[2]], col=heat.colors(1000),xlim=c(315000, 325700), ylim=c(4385000, 4405500))
# must add singly
#xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
#contour(xyzB, add=TRUE)
plot(fhat10, display='slice', cont=c(95,75,50,25), approx.cont=T, drawlabels=T, drawpoints=F, add=T, pch=19, ptcol='black')
box()
title(main = "2010")
#points(tnfd$x[tnfd$year=='2010'],tnfd$y[tnfd$year=='2010'])

SpatialPolygonsRescale(layout.scale.bar(), offset=c(locator(1)),
                       scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)
text(locator(1), "5 km")


#####################

# check for spatial overlap

wro2<-kerneloverlaphr(liud, method = c("BA"), percent = 95, conditional = T)
wro2
cro2<-kerneloverlaphr(liud, method = c("BA"), percent = 50, conditional = T)
cro2


################################################################################################################################
################################################################################################################################


# save workspace
save.image(file = "ibat colony foraging.RData")

# delete saved workspace
#unlink("ibat colony foraging.RData")


