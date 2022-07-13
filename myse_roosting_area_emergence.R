

## THIS IS THE ROOSTING AREA CODE BUT RUN THROUGH WITH MYSTERY BATS PUT IN
## TO REFLECT THE NUMBER OF BATS THAT WERE IN THE ROOSTS ON THE NIGHTS WE
## TRACKED -- 


library(adehabitatLT)
library(adehabitatHR)
library(sp)
library(maptools)
library(raster)
library(GDAL)
library(rgdal)

# read in dataset in csv format
# this should be the sequential list of locations with animal/group identifier

rts <- read.csv("myse_emergence.csv", header=TRUE, sep = ",", 
                quote = "\"", dec = ".", numerals = "no.loss", 
                fill = TRUE, comment.char = "")
# get the names of the variables in the data and check the first few lines of data
## Silvis has it set up Bat, Tree, X, Y, Date, daynum, ReproByBat, Age, Sex, Year
t(names(rts))
head(rts)
tail(rts)
dim(rts)



length(unique(rts$Bat)) ## 32 bats, but this means nothing because we don't 
  ## know the identities of our mystery bats

length(rts$Bat) ## 132 obs, still kinda means nothing


# getting a summary of the data for the animals/groups
# summary(as.factor(rts$Bat)) ## number of locations per bat

# note that for kernel estimation here there must be > 4 locations for every animal/group

################################################################################################################################
################################################################################################################################

# some data processing to make sure locations data are in the correct/useful format

# combine the location data into a single dataframe; it is best to use UTM values
## kg added print command to make sure we didn't lose decimal places
locs<-cbind(print(rts$X, digits = 14),print(rts$Y, digits = 14))
plot(locs)

# turn the output into a data frame
locs<-data.frame(locs)
names(locs)<-c('UTM.X','UTM.Y')
# plot the points to make sure they look correct
plot(locs)

# make a spatial points dataframe - required for home range calculations
library(sp)
# id<-as.factor(rts$Year)
# idsp <- data.frame(id)
# coordinates(idsp) <- locs
# class(idsp)

coordinates(locs) <- locs
class(locs)

loc<-SpatialPoints(locs, CRS(as.character("+proj=utm +zone=18 ellps=WGS84", data=locs)))
plot(loc)

# summary(idsp)
# plot(idsp)
summary(locs)
plot(locs)

## save work so far
save.image("myse_roosting_area_emergence.RDS")
################################################################################################################################
################################################################################################################################

# look at the tracking histories

# some summary statistics

str(rts$Bat)
str(rts$Tree)

# summary(rts$daynum)
# summary(as.factor(rts$ID)) ## kg doesn't know what this col could be, it isn't in
## original silvis image file

# table containing the tracking histories
# table(rts$daynum, rts$Bat)

## the following 2 lines don't work for kg, despite running network analysis functions code
# plot.min.max(as.factor(rts$Bat[rts$Year=='2018']), rts$daynum[rts$Year=='2019'], dots=F)
# plot.min.max(as.factor(rts$Bat[rts$Year=='2018']), rts$daynum[rts$Year=='2019'], dots=F)

################################################################################################################################
################################################################################################################################

# look at the trajectories of tracked animals

# load the library with the trajectory functions
library(adehabitatLT)

# this will plot the trajectories of all animals/groups
# if all do not appear, then the window frame is too small
# trax<-as.ltraj(locs, id = rts$Bat, typeII = FALSE)
# trax
# # generic quick plotting
# plot(trax)
# 
# # look at the trajectory data for a single animal
# trax[[1]]
# # view summary for distance data for a single animal
# summary(trax[[16]]$dist)
# 
# # get the mean distance moved by all animals
# # only use if there are no missing values in locations
# tr<-ld(trax)
# mean(tr$dist, na.rm=T)
# sd(tr$dist, na.rm=T)
# 
# # plot the trajectory of a specific animal
# plot(trax, perani=F, final=F, burst=rts$Bat[rts$Bat=='72407'])
# 
# # create an interactive trajectory that shows the individual movements
# # we open a new window specifically for the interactive trajectory
# # use the commands in the results window to interact with the trajectory
# # (n/p on the keyboard control next/previous relocation respectively)
# # note that any animal can be displayed, but this starts with the first in
# # alphabetical or numerical order
# windows()
# trajdyn(trax, burst = attr(trax[[3]], "burst"), hscale = 1, vscale = 1,
#         recycle = TRUE, display = c("guess", "windows", "tk")) 
# 
# # convert trajectories into a spatial line object that shows all trajectories
# # load library with conversion function
# # library(sp)
# # save trajectory as a spatial lines object
# pl<-ltraj2sldf(trax) 
# # check the raw data for the spatial lines object
# pl
# # see the summary data for the trajectories
# summary(pl)
# # plot all the trajectories simultaneously in a single window
# plot(pl)

# write the spatial lines object to a GIS friendly format
# load library with needed function and export
library(maptools)
library(rgdal)
# here, pl is the name of the spatial lines object, and trajectories is the name of the output file

# writeSpatialShape(pl, "trajectories", factor2char = TRUE, max_nchar=254)

sfdir <- file.path("G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/GIS") 
# writeOGR(pl, sfdir, "trajectories", driver = "ESRI Shapefile")
# 
# save.image("myse_roosting_area_2yr.RDS")

################################################################################################################################
################################################################################################################################

# create home ranges and utilization distributions for tracked animals

# load the R library that contains the functions used below
library(adehabitatHR)

# calculate home range area using the kernel method; manually select fixed or adaptive
ud <- kernelUD(loc, h = "href", grid = 200, same4all = T, hlim = c(0.1, 1.5),
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

## kg updated saving hr shapefiles, and saved to designated GIS folder for this proj
writeOGR(ver, sfdir, "roost 95% homerange 2yr emergence", driver = "ESRI Shapefile")

# make a pretty plot of the UD's
volud<-getvolumeUD(ud, standardize = F)
image(volud, col=heat.colors(1000))


# find centroids for each year
m1<-(mean(rts$X))
m2<-(mean(rts$Y))


# # distance between centroids
# a<-(m1-m3)^2
# b<-(m2-m4)^2
# sqrt(abs(a-b)) ## 300.1671


# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 6, height = 6)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))

image(volud, col=heat.colors(1000),xlim=c(683715.5, 682582.5), ylim=c(4514639, 4516529))
# must add singly
xyzA <- as.image.SpatialGridDataFrame(volud)
contour(xyzA, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc[rts],col='black',pch=19)
box()
# plot(WIFL_sf, add = TRUE)
title(main = "2018–2019 Utilization Distribution")
# points(m1,m2, pch=20, col="blue", cex=2)

# image(volud[[2]], col=heat.colors(1000),xlim=c(682715.5, 684082.5), ylim=c(4514239, 4516529))
# # must add singly
# xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
# contour(xyzB, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
# points(loc[rts$Year=="2019"],col='black',pch=19)
# box()
# title(main = "2019")
# points(m3,m4, pch=20, col="blue", cex=2)

# scale bar and label will be placed interactively
## scale bar code makes r crash for kg
# SpatialPolygonsRescale(layout.scale.bar(), offset = c(locator(1)),
#                        scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)
# text(locator(1), "5 km")

############################

# look for spatial similarity of the roosting areas between years at selected percent home range areas
# wro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 95) # entire home range
# wro
# cro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 50) # "core" home range
# cro

############################

# export home range estimates to GIS format
class(volud)
# image(volud)

writeOGR(volud, sfdir, "volud_emergence", driver = "ESRI Shapefile")

# Gud<-estUDm2spixdf(volud)
# # image(Gud[1])
# # image(Gud[2])
# 
# library(rgdal)
# # writeGDAL(Gud, 'roosting ud')
# 
# writeOGR(volud, sfdir, "roosting UD both years", driver = "ESRI Shapefile")


#hri<-c(25,50,75,95)
#for(i in 1:length(hri)){
#  hris[i]<-getverticeshr(ud,hri[i])
#}
# 
# library(maptools)
# ver <- getverticeshr(ud, 25)
# #writePolyShape(ver, "roost 25% homerange")
# writeOGR(ver, sfdir, "roost 25% homerange", driver = "ESRI Shapefile")



################################################################################################################################
################################################################################################################################

# save workspace
# save.image(file = "ibat colony roosting.RData")
save.image("myse_roosting_area_emergence.RDS")

# load("myse_roosting_area2.RDS")
## try to map it here
library(rgdal)
library(sf)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(raster)
library(mapview)


## load in WIFL outline shapefile
WIFL_sf <- st_read(dsn = "G:/My Drive/FIRE, FIRE ISLAND!/GIS/wfe.shp", 
                   crs = "+proj=utm +zone=18 ellps=WGS84")

volud_sf <- st_read(dsn = "G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/GIS/volud_emergence.shp", 
                    crs = "+proj=utm +zone=18 ellps=WGS84")

# 
# 
# class(volud_sf)
# 
# volud_sp <- as(volud_sf, "Spatial")

volud2 <- volud
volud2$n <- ifelse(volud2$n > 99, NA, volud2$n)

## 
kareas <- getverticeshr(ud, 95)
kareas50 <- getverticeshr(ud, 50)
kareas25 <- getverticeshr(ud, 25)
kdareas <- fortify(kareas)

## mapview
mapview(WIFL_sf, lwd = 2, alpha.regions = 0.01, legend = FALSE, map.types = "OpenStreetMap") + 
  mapview(volud2, col.regions = rev(brewer.pal(9, "YlOrRd")), na.color = "transparent") +
  mapview(loc, col.regions = "black", cex = 5, alpha = 0.01) + 
  mapview(kareas, lw = 2, alpha.regions = 0.01, legend = FALSE) +
  mapview(kareas50, lwd = 2, alpha.regions = 0.01, legend = FALSE) + 
  mapview(kareas25, lwd = 2, alpha.regions = 0.01, legend = FALSE)



