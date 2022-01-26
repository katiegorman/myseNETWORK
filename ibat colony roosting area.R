

################################################################################################################################
################################################################################################################################
#                                       Indiana Bat COlony Roosting

# Code by Alex Silvis, 14 September 2012, using R version 2.15
# revised 27 September 2012
# revised 21 February 2013
# revised 11 April 2013

################################################################################################################################
################################################################################################################################

# set the working directory where all of the data files are located

setwd("//minnow.cc.vt.edu/cnre1/silvis/My Documents/R working directory/andrew ibats")

# data for home range estimation should be a list of locations with animal/group identifier
# for location dates, day of year is generally much easier to use than date in dd/mm/yyyy or similar format
# so data should include a column for day of year, which can be generated easily in Excel using
#       =TEXT((A2-DATEVALUE("1/1/" &TEXT(A2,"yy"))+1),"000")
#       where A2 is the cell number containing the date value

# load data from previous run
load("ibat colony roosting.RData", envir=parent.frame())
library(adehabitatLT)
library(adehabitatHR)
library(sp)
library(maptools)
library(raster)
library(GDAL)

# read in datashet in csv format
# this should be the sequential list of locations with animal/group identifier
#rts <- data.frame(read.table(file='ibatnetworkdata.csv', sep=',', header=TRUE))
# get the names of the variables in the data and check the first few lines of data
t(names(rts))
head(rts)
tail(rts)
dim(rts)

# remove bats not in the network
#rts <- subset(rts, Bat != 15450, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)
#rts <- subset(rts, Bat != 15531, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)
#rts <- subset(rts, Bat != 18858, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)

length(unique(rts$Bat[rts$Year=='2009']))
length(unique(rts$Bat[rts$Year=='2010']))

length(rts$Bat[rts$Year=='2009'])
length(rts$Bat[rts$Year=='2010'])

# getting a summary of the data for the animals/groups
summary(as.factor(rts$Year))
# note that for kernel estimation here there must be > 4 locations for every animal/group

################################################################################################################################
################################################################################################################################

# some data processing to make sure locations data are in the correct/useful format

# combine the location data into a single dataframe; it is best to use UTM values
locs<-cbind(rts$X,rts$Y)
plot(locs)

# turn the output into a data frame
locs<-data.frame(locs)
names(locs)<-c('UTM.X','UTM.Y')
# plot the points to make sure they look correct
plot(locs)

# make a spatial points dataframe - required for home range calculations
library(sp)
id<-as.factor(rts$Year)
idsp <- data.frame(id)
coordinates(idsp) <- locs
class(idsp)
library(sp)
loc<-SpatialPoints(locs, CRS(as.character("+proj=utm +zone=17 ellps=WGS84", data=locs)))
plot(loc)
summary(idsp)
plot(idsp)

################################################################################################################################
################################################################################################################################

# look at the tracking histories

# some summary statistics

str(rts$Bat)
str(rts$Tree)

summary(rts$daynum)
summary(as.factor(rts$ID))

# table containing the tracking histories
#table(rts$daynum, rts$Bat)

plot.min.max(as.factor(rts$Bat[rts$Year=='2009']), rts$daynum[rts$Year=='2009'], dots=F)
plot.min.max(as.factor(rts$Bat[rts$Year=='2010']), rts$daynum[rts$Year=='2010'], dots=F)

################################################################################################################################
################################################################################################################################

# look at the trajectories of tracked animals

# load the library with the trajectory functions
library(adehabitatLT)

# this will plot the trajectories of all animals/groups
# if all do not appear, then the window frame is too small
trax<-as.ltraj(locs, id = rts$Bat, typeII = FALSE)
trax
# generic quick plotting
plot(trax)

# look at the trajectory data for a single animal
trax[[1]]
# view summary for distance data for a single animal
summary(trax[[1]]$dist)

# get the mean distance moved by all animals
# only use if there are no missing values in locations
tr<-ld(trax)
mean(tr$dist, na.rm=T)
sd(tr$dist, na.rm=T)

# plot the trajectory of a specific animal
plot(trax, perani=F, final=F, burst=rts$Bat[rts$Bat=='13452'])

# create an interactive trajectory that shows the individual movements
# we open a new window specifically for the interactive trajectory
# use the commands in the results window to interact with the trajectory
# (n/p on the keyboard control next/previous relocation respectively)
# note that any animal can be displayed, but this starts with the first in
# alphabetical or numerical order
#windows()
#trajdyn(trax, burst = attr(trax[[1]], "burst"), hscale = 1, vscale = 1,
        #recycle = TRUE, display = c("guess", "windows", "tk")) 

# convert trajectories into a spatial line object that shows all trajectories
# load library with conversion function
library(sp)
# save trajectory as a spatial lines object
pl<-ltraj2sldf(trax) 
# check the raw dataa for the spatial lines object
pl
# see the summary data for the trajectories
summary(pl)
# plot all the trajectories simultaneously in a single window
plot(pl)

# write the spatial lines object to a GIS friendly format
# load library with needed function and export
library(maptools)
# here, pl is the name of the spatial lines object, and trajectories is the name of the output file
#writeSpatialShape(pl, "trajectories", factor2char = TRUE, max_nchar=254)

################################################################################################################################
################################################################################################################################

# create home ranges and utilization distributions for tracked animals

# load the R library that contains the functions used below
library(adehabitatHR)

# calculate home range area using the kernel method; manually select fixed or adaptive
ud <- kernelUD(idsp, h = "href", grid = 200,same4all = T, hlim = c(0.1, 1.5),
               kern = c("bivnorm"), extent = 1)
# the UD of the animals/groups
#windows()
image(ud)
# calculation of the home range area at requested percentage values
ka<-kernel.area(ud, percent = seq(20, 95, by = 5),unin = c("m"),unout = c("ha"), standardize = FALSE)
ka
plot(ka)
# get home range area at a single percentage value- will be used to transfer home ranges to a GIS later
ver <- getverticeshr(ud, 95)
# plot the home ranges
plot(ver)

# export the home range shapefile created in the previous step so that it can be used in a GIS
library(maptools)
#writePolyShape(ver, "roost 95% homerange")
# the shapefiles will be written to the working directory set; i.e., where the datafile is

# make a pretty plot of the UD's
volud<-getvolumeUD(ud, standardize = F)
#image(volud, col=heat.colors(1000))


# find centroids for each year
m1<-(mean(rts$X[rts$Year=='2009']))
m2<-(mean(rts$Y[rts$Year=='2009']))

m3<-(mean(rts$X[rts$Year=='2010']))
m4<-(mean(rts$Y[rts$Year=='2010']))

# distance between centroids
a<-(m1-m3)^2
b<-(m2-m4)^2
sqrt(abs(a-b))


# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 11, height = 6)
par(mfrow=c(1,2))
par(mar=c(1,1,2,1))

image(volud[[1]], col=heat.colors(1000),xlim=c(314751.4, 323257.1), ylim=c(4391025, 4397669))
# must add singly
xyzA <- as.image.SpatialGridDataFrame(volud[[1]])
contour(xyzA, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc[rts$Year=="2009"],col='black',pch=19)
box()
title(main = "2009")
#points(m1,m2, pch=20, col="blue", cex=2)

image(volud[[2]], col=heat.colors(1000),xlim=c(314751.4, 323257.1), ylim=c(4391025, 4397669))
# must add singly
xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
contour(xyzB, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc[rts$Year=="2010"],col='black',pch=19)
box()
title(main = "2010")
#points(m3,m4, pch=20, col="blue", cex=2)

# scale bar and label will be placed interactively
SpatialPolygonsRescale(layout.scale.bar(), offset = c(locator(1)),
                       scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)
text(locator(1), "5 km")

############################

# look for spatial similarity of the roosting areas between years at selected percent home range areas
wro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 95) # entire home range
wro
cro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 50) # "core" home range
cro

############################

# export home range estimates to GIS format

Gud<-estUDm2spixdf(volud)
image(Gud[1])
image(Gud[2])

library(rgdal)
writeGDAL(Gud, 'roosting ud')


#hri<-c(25,50,75,95)
#for(i in 1:length(hri)){
#  hris[i]<-getverticeshr(ud,hri[i])
#}

library(maptools)
ver <- getverticeshr(ud, 25)
#writePolyShape(ver, "roost 25% homerange")


################################################################################################################################
################################################################################################################################

# save workspace
save.image(file = "ibat colony roosting.RData")

# delete saved workspace
#unlink("ibat colony roosting.RData")

