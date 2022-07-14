



#                                       Indiana Bat COlony Roosting

# Code by Alex Silvis, 14 September 2012, using R version 2.15
# revised 27 September 2012
# revised 21 February 2013
# revised 11 April 2013

## Code revised by Katie Gorman, 24 January 2022, using R version 4.0.3 for MYSE
## revised 13 July 2022 by KG to explore other mapping packages (tmap), 
  ## since mapview is stupid and won't let me add graticules




## working directory is already set because working within a github project

# data for home range estimation should be a list of locations with animal/group identifier
# for location dates, day of year is generally much easier to use than date in dd/mm/yyyy or similar format
# so data should include a column for day of year, which can be generated easily in Excel using
#       =TEXT((A2-DATEVALUE("1/1/" &TEXT(A2,"yy"))+1),"000")
#       where A2 is the cell number containing the date value

# load data from previous run
load("myse_roosting_area_20220713.RDS")

#### libraries ####

library(adehabitatLT) ## trajectories
library(adehabitatHR) ## create homerange and UD
library(sp) ## make spatialpointsDF
library(maptools) ## trajectories to GIS friendly
library(raster)
library(rgdal) ## trajectories to GIS friendly
library(ggplot2)
library(dplyr) ## piping



# read in dataset in csv format
# this should be the sequential list of locations with animal/group identifier

rts <- read.csv("myse_roosts.csv", header=TRUE, sep = ",", 
                           quote = "\"", dec = ".", numerals = "no.loss", 
                           fill = TRUE, comment.char = "")


# get the names of the variables in the data and check the first few lines of data
## Silvis has it set up Bat, Tree, X, Y, Date, daynum, ReproByBat, Age, Sex, Year
t(names(rts))
head(rts)
tail(rts)
dim(rts)
unique(rts$Tree) ## 30 unique roosts

# remove bats not in the network -- 
## katie has none, I think they are all in the network even if not connected
## however, unconnected female and juv male were removed from network analysis
## so on this run I'm taking them out

rts <- subset(rts, Bat != 72414, select=c(1:10)) ## gets rid of 2 roosts used by A-F
dim(rts)
unique(rts$Tree) ## 28 roosts
rts<-droplevels(rts)
rts <- subset(rts, Bat != 72418, select=c(1:10)) ## gets rid of 3 more roosts used by J-M
dim(rts)
rts<-droplevels(rts)


length(unique(rts$Bat[rts$Year=='2018'])) ## 14
length(unique(rts$Bat[rts$Year=='2019'])) ## 2

length(rts$Bat[rts$Year=='2018']) ## 59
length(rts$Bat[rts$Year=='2019']) ## 12

# getting a summary of the data for the animals/groups
summary(as.factor(rts$Year))
# note that for kernel estimation here there must be > 4 locations for every animal/group


## ------------------------------------------------------------------------------------- ##
#### some data processing to make sure locations data are in the correct/useful format ####
## ------------------------------------------------------------------------------------- ##


# combine the location data into a single dataframe; it is best to use UTM values
## kg added print command to make sure we didn't lose decimal places
locs<-cbind(print(rts$X, digits = 14), print(rts$Y, digits = 14))
plot(locs)

# turn the output into a data frame and rename the cols
locs<-data.frame(locs)
names(locs)<-c('UTM.X','UTM.Y')

# plot the points to make sure they look correct
plot(locs) 
## yup

# make a spatial points dataframe - required for home range calculations
id<-as.factor(rts$Year)
idsp <- data.frame(id)
coordinates(idsp) <- locs
class(idsp)

loc<-SpatialPoints(locs, CRS(as.character("+proj=utm +zone=18 ellps=WGS84", data=locs)))
plot(loc)
summary(idsp)
plot(idsp)


## also do spatial points df for factor female so both years are in one
id.f <- as.factor(rts$Sex)
idsp.f <- data.frame(id.f)
coordinates(idsp.f) <- locs
class(idsp.f)


## save work so far
save.image("myse_roosting_area_20220713.RDS")


## ---------------------------------- ##
#### look at the tracking histories ####
## ---------------------------------- ##


# some summary statistics
str(rts$Bat)
str(rts$Tree)

## which days were recorded?
summary(rts$daynum)
# summary(as.factor(rts$ID)) ## kg doesn't know what this col could be, it isn't in
                            ## original silvis image file

# table containing the tracking histories
table(rts$daynum, rts$Bat)
t(names(rts))

## the following 2 lines don't work for kg, despite running network analysis functions code
# plot.min.max(as.factor(rts$Bat[rts$Year=='2018']), rts$daynum[rts$Year=='2019'], dots=F)
# plot.min.max(as.factor(rts$Bat[rts$Year=='2018']), rts$daynum[rts$Year=='2019'], dots=F)

## here is KG attempt
rts.track <- rts
rts.track$Date <- as.Date(rts.track$Date, format = "%m/%d/%Y")
head(rts.track)

rts.track <- rts.track %>%
  group_by(Bat, ReproByBat, Year) %>%
  summarize(minDate = min(Date),
            maxDate = max(Date))

## the plot below is fine for now
## troubleshooting: 72411 has two repro stages because 2 years of tracking,
## and bats 72412/72403 only had one data point (capture roosts?) so they are blank
ggplot(data = rts.track) +
  geom_segment(aes(xend = strftime(maxDate,format="%d-%b"), x = strftime(minDate,format="%d-%b"),
                 y = Bat, yend = Bat,
                 color = ReproByBat)) +
  theme_classic() +
  labs(x="Date") 

save.image("myse_roosting_area_20220713.RDS")


## ----------------------------------------------- ##
#### look at the trajectories of tracked animals ####
## ----------------------------------------------- ##


# this will plot the trajectories of all animals/groups
# if all do not appear, then the window frame is too small
trax<-as.ltraj(locs, id = rts$Bat, typeII = FALSE)
trax
# generic quick plotting
plot(trax)

# look at the trajectory data for a single animal
trax[[11]]

# view summary for distance data for a single animal
summary(trax[[11]]$dist)

# get the mean distance moved by all animals
# only use if there are no missing values in locations
tr<-ld(trax)
mean(tr$dist, na.rm=T)
sd(tr$dist, na.rm=T)

# plot the trajectory of a specific animal
plot(trax, perani=F, final=F, burst=rts$Bat[rts$Bat=='72411'])

# create an interactive trajectory that shows the individual movements
# we open a new window specifically for the interactive trajectory
# use the commands in the results window to interact with the trajectory
# (n/p on the keyboard control next/previous relocation respectively)
# note that any animal can be displayed, but this is just showing the 
## most interesting one
windows()
trajdyn(trax, burst = attr(trax[[11]], "burst"), hscale = 1, vscale = 1,
        recycle = TRUE, display = c("guess", "windows", "tk")) 

# convert trajectories into a spatial line object that shows all trajectories
# load library with conversion function

# save trajectory as a spatial lines object
pl<-ltraj2sldf(trax) 
# check the raw data for the spatial lines object
pl
# see the summary data for the trajectories
summary(pl)
# plot all the trajectories simultaneously in a single window
plot(pl)

# write the spatial lines object to a GIS friendly format
# here, pl is the name of the spatial lines object, and trajectories is the name of the output file
## the 2 in front of the filename indicates this is the second run though

sfdir <- file.path("G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/GIS"); dir.create(sfdir) 
writeOGR(pl, sfdir, "2 trajectories", driver = "ESRI Shapefile")

save.image("myse_roosting_area_20220713.RDS")


## ------------------------------------------------------------------------ ##
#### create home ranges and utilization distributions for tracked animals ####
## ------------------------------------------------------------------------ ##


# calculate home range area using the kernel method; manually select fixed or adaptive
ud <- kernelUD(idsp, h = "href", grid = 200, same4all = T, hlim = c(0.1, 1.5),
               kern = c("bivnorm"), extent = 1)

# the UD of the animals/groups
#windows()
image(ud)
print(ud)
ud

# calculation of the home range area at requested percentage values
ka<-kernel.area(ud, percent = seq(20, 95, by = 5),unin = c("m"),unout = c("ha"), standardize = FALSE)
ka
plot(ka)

# get home range area at a single percentage value- will be used to transfer home ranges to a GIS later
ver <- getverticeshr(ud, 95)

# plot the home ranges
plot(ver)

# export the home range shapefile created in the previous step so that it can be used in a GIS
## saved to designated GIS folder for this proj
## 2 in front of file name indicates this is the second time running through the code
writeOGR(ver, sfdir, "2 roost 95% homerange", driver = "ESRI Shapefile")

# make a pretty plot of the UD's
volud<-getvolumeUD(ud, standardize = F)
image(volud, col=heat.colors(1000))


# find centroids for each year
m1<-(mean(rts$X[rts$Year=='2018']))
m2<-(mean(rts$Y[rts$Year=='2018']))

m3<-(mean(rts$X[rts$Year=='2019']))
m4<-(mean(rts$Y[rts$Year=='2019']))

# distance between centroids
a<-(m1-m3)^2
b<-(m2-m4)^2
sqrt(abs(a-b)) ## 172.8974


# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 11, height = 6)
par(mfrow=c(1,2))
par(mar=c(1,1,2,1))

image(volud[[1]], col=heat.colors(1000),xlim=c(682715.5, 683582.5), ylim=c(4514439, 4516529))
# must add singly
xyzA <- as.image.SpatialGridDataFrame(volud[[1]])
contour(xyzA, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc[rts$Year=="2018"],col='black',pch=19)
box()
title(main = "2018")
# points(m1,m2, pch=20, col="blue", cex=2)

image(volud[[2]], col=heat.colors(1000),xlim=c(682715.5, 684082.5), ylim=c(4514239, 4516529))
# must add singly
xyzB <- as.image.SpatialGridDataFrame(volud[[2]])
contour(xyzB, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc[rts$Year=="2019"],col='black',pch=19)
box()
title(main = "2019")
# points(m3,m4, pch=20, col="blue", cex=2)

# scale bar and label will be placed interactively
## scale bar code makes r crash for kg
# SpatialPolygonsRescale(layout.scale.bar(), offset = c(locator(1)),
#                         scale = 5000, fill=c("transparent","black"), plot.grid=FALSE)
#  text(locator(1), "5 km")

# look for spatial similarity of the roosting areas between years at selected percent home range areas
wro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 95) # entire home range
wro

#           2018      2019
# 2018 0.9499957 0.3554692
# 2019 0.3554692 0.9499837

cro<-kerneloverlap(idsp,grid=500, meth="BA", conditional=TRUE, percent = 50) # "core" home range
cro

#           2018       2019
# 2018 0.49997627 0.01619581
# 2019 0.01619581 0.49984414


# export home range estimates to GIS format

Gud<-estUDm2spixdf(volud)
image(Gud[1])
image(Gud[2])

writeOGR(Gud, sfdir, "2 roosting UD", driver = "ESRI Shapefile")


# hri<-c(25,50,75,95)
# for(i in 1:length(hri)){
#  hri[i]<-getverticeshr(ud,hri[i])
# }


ver <- getverticeshr(ud, 25)
writeOGR(ver, sfdir, "2 roost 25% homerange", driver = "ESRI Shapefile")

save.image("myse_roosting_area_20220713.RDS")


## ------------------------------------------------ ##
#### do it all again but pool both years together ####
## ------------------------------------------------ ##


## remember that spdf.f you made earlier? time to use it

## a refresher:
# id.f <- as.factor(rts$Sex)
# idsp.f <- data.frame(id.f)
# coordinates(idsp.f) <- locs
# class(idsp.f)


summary(idsp.f)
plot(idsp.f)

summary(locs)
plot(locs)

## ---------------------------------------------------------------------------- ##
#### create home ranges and utilization distributions for both years together ####
## ---------------------------------------------------------------------------- ##



# calculate home range area using the kernel method; manually select fixed or adaptive
ud.f <- kernelUD(idsp.f, h = "href", grid = 200, same4all = T, hlim = c(0.1, 1.5),
               kern = c("bivnorm"), extent = 1)

# the UD of the animals/groups
#windows()
image(ud.f)
print(ud.f)
ud.f

# calculation of the home range area at requested percentage values
ka.f <- kernel.area(ud.f, percent = seq(20, 95, by = 5),unin = c("m"),unout = c("ha"), standardize = FALSE)
ka.f
plot(ka.f)

# get home range area at a single percentage value- will be used to transfer home ranges to a GIS later
ver.f <- getverticeshr(ud.f, 95)

# plot the home ranges
plot(ver.f)

# export the home range shapefile created in the previous step so that it can be used in a GIS
## saved to designated GIS folder for this proj
## 2 in front of file name indicates this is the second time running through the code
writeOGR(ver.f, sfdir, "2 roost 95% homerange", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# make a pretty plot of the UD's
volud.f <- getvolumeUD(ud.f, standardize = F)
image(volud.f, col=heat.colors(1000))


# find centroid 
m1.f<-(mean(rts$X))
m2.f<-(mean(rts$Y))


# make a better plot with the contour lines corresponding to the percent utilization distribution
x11(width = 5, height = 6)
# par(mfrow=c(1,2))
# par(mar=c(1,1,2,1))

image(volud.f[[1]], col=heat.colors(1000),xlim=c(682715.5, 683582.5), ylim=c(4514439, 4516529))
# must add singly
xyzA.f <- as.image.SpatialGridDataFrame(volud.f[[1]])
contour(xyzA.f, add=TRUE, levels=(c(25,50,75,95)), nlevels=2)
points(loc,col='black',pch=19)
box()
title(main = "UD two years")
points(m1.f,m2.f, pch=20, col="blue", cex=2)


# export home range estimates to GIS format
class(volud.f)
image(volud.f)

writeOGR(volud.f, sfdir, "2 volud all bats", driver = "ESRI Shapefile")

Gud.f <- estUDm2spixdf(volud.f)
image(Gud.f[1])


ver.f.25 <- getverticeshr(ud.f, 25)
writeOGR(ver.f.25, sfdir, "2 roost 25% homerange all bats", driver = "ESRI Shapefile")

save.image("myse_roosting_area_20220713.RDS")

## ----------- ##
#### MAP?!?! ####
## ----------- ##

library(rgdal)
library(sf)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(raster)
library(mapview)
library(graticule)


## load in WIFL outline shapefile
WIFL_sf <- st_read(dsn = "G:/My Drive/FIRE, FIRE ISLAND!/GIS/wfe.shp", 
                   crs = "+proj=utm +zone=18 ellps=WGS84")

volud_sf <- st_read(dsn = "G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/GIS/volud.shp", 
                    crs = "+proj=utm +zone=18 ellps=WGS84")



class(volud_sf)

volud_sp <- as(volud_sf, "Spatial")

volud2 <- volud
volud2$n <- ifelse(volud2$n > 95, NA, volud2$n)

## 
kareas <- getverticeshr(ud, 95)
kareas50 <- getverticeshr(ud, 50)
kareas25 <- getverticeshr(ud, 25)
kdareas <- fortify(kareas)

## define projection
prj <- "+proj=utm +zone=18 ellps=WGS84"

## specify where we want meridians and parallels
lats <- seq(4513600, 4516600, length = 6)
lons <- seq(652000, 685000, length = 6)

## push them out a little bit on each side
x1 <- range(lons) + c(-1000, 1000)
y1 <- range(lats) + c(-1000, 1000)

## build the lines with precise locations and ranges
grat <- graticule(lons, lats, proj = prj, xlim = x1, ylim = y1)

## build the labels, here they sit exactly on the western and northern extent
## of our line ranges
labs <- graticule_labels(lons, lats, xline = min(xl), yline = max(yl), proj = prj)


## mapview
mapview(WIFL_sf, alpha.regions = 0.01) + 
  mapview(volud2, col.regions = rev(brewer.pal(9, "YlOrRd")), na.color = "transparent") +
  mapview(loc, col.regions = "black", cex = 5, alpha = 0.01) + mapview(kareas, alpha.regions = 0.01) +
  mapview(kareas50, alpha.regions = 0.01) + mapview(kareas25, alpha.regions = 0.01)




# 
# ## try with leaflet
# library(leaflet)
# 



ggplot(kdareas) + 
  geom_polygon(aes(x=long, y=lat, group = group, fill = id, color = id)) +


  geom_sf(data = WIFL_sf, fill = "NA", color = "grey") +
  
  

  # geom_sf(data = lines_HARV, aes(color = TYPE), size = 1) +
  # geom_sf(data = point_HARV) +
  # ggtitle("NEON Harvard Forest Field Site") + 
  coord_sf()
