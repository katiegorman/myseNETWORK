

################################################################################################################################
################################################################################################################################
#                                       Indiana Bat Individual Foraging Areas
#                                          - tracking histories
#                                          - trajectories
#                                          - home ranges

# Code by Alex Silvis, 12 August 2013

################################################################################################################################
################################################################################################################################

# set the working directory where all of the data files are located

setwd("//minnow.cc.vt.edu/cnre1/silvis/My Documents/R working directory/andrew ibats")


# load data from previous run
load("ibat individual foraging.RData", envir=parent.frame())
library(adehabitatLT)
library(adehabitatHR)
library(sp)
library(maptools)

# read in datashet in csv format
# this should be the sequential list of locations with animal/group identifier
#fds <- data.frame(read.table(file='ibatforagingsubset.csv', sep=',', header=TRUE))
# this data is a subset of the total foraging data; only individuals with > 39 locations included

# get the names of the variables in the data and check the first few lines of data
t(names(fds))
head(fds)
tail(fds)
dim(fds)

# getting a summary of the data for the animals/groups
(s1<-summary(as.factor(fds$name)))
mean(s1)
stderr(s1)

# how many bats with ~40 foraging points in 2009?
length(unique(fds$name[fds$year=='2009']))

# how many bats with ~40 foraging points in 2010?
length(unique(fds$name[fds$year=='2010']))

################################################################################################################################
################################################################################################################################

# some data processing to make sure locations data are in the correct/useful format

# combine the location data into a single dataframe; it is best to use UTM values
locs<-cbind(fds$x,fds$y)
plot(locs)

# turn the output into a data frame
locs<-data.frame(locs)
names(locs)<-c('UTM.X','UTM.Y')
# plot the points to make sure they look correct
plot(locs)

# make a spatial points dataframe - required for home range calculations
library(sp)
id<-as.factor(fds$name)
idsp <- data.frame(id)
coordinates(idsp) <- locs
class(idsp)
library(sp)
loc<-SpatialPoints(locs, CRS(as.character("+proj=utm +zone=17 ellps=WGS84", data=locs)))
plot(loc)
summary(idsp)


################################################################################################################################
################################################################################################################################

# look at the trajectories of tracked animals

# load the library with the trajectory functions
library(adehabitatLT)

# this will plot the trajectories of all animals/groups
# if all do not appear, then the window frame is too small

### conversion of the date to the format POSIX so the location date and time can be used in the trajectory
head(fds$DATETIME4)
da <- as.POSIXct(strptime(as.character(fds$DATETIME4),"%m%d%y %I%M %p"))
head(da)

trax<-as.ltraj(locs, id = fds$name, typeII = TRUE, date=da)
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
plot(trax, perani=F, final=F, burst=fds$name[fds$name=='8270'])

# create an interactive trajectory that shows the individual movements
# we open a new window specifically for the interactive trajectory
# use the commands in the results window to interact with the trajectory
# (n/p on the keyboard control next/previous relocation respectively)
# note that any animal can be displayed, but this stafds with the first in
# alphabetical or numerical order
#windows()
#trajdyn(trax, burst = attr(trax[[1]], "burst"), hscale = 1, vscale = 1,
#        recycle = TRUE, display = c("guess", "windows", "tk")) 

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
windows();plot(pl)

# write the spatial lines object to a GIS friendly format
# load library with needed function and export
library(maptools)
# here, pl is the name of the spatial lines object, and trajectories is the name of the output file
#writeSpatialShape(pl, "individtrajectories", factor2char = TRUE, max_nchar=254)

################################################################################################################################
################################################################################################################################

# create home ranges and utilization distributions for tracked animals

# load the R library that contains the functions used below
library(adehabitatHR)

# calculate the home range using biased random bridges

# estimate the diffusion component using the plug-in method
vv <- BRB.D(trax, Tmax = 60*60, Lmin = 50)
vv

# note that the values are given here as m^2/s
# values in m^2 per min are:
vv[[1]][,2]*60

# note that an alternative estimation of the Diffusion coefficient
# could be found using maximum likelihood
vv2 <- BRB.likD(trax, Tmax = 60*60, Lmin = 50)
vv2
#vv[[1]][,2]*60

# hmin is the error of the location data that will be used to smooth

## Estimation of the UD 
BRBud <- BRB(trax, D = vv, Tmax = 60*60, Lmin = 50, hmin = 88, filtershort = FALSE,
             b = FALSE, grid = 250, same4all=F, extent=1)

# show the UD.
image(BRBud, col=heat.colors(1000))

# BRBver <- getverticeshr(BRBud, 95) # fails, not sure why
# so calling the UDs one by one (which works for some reason)
# getting the 95% UD and plotting them all together
# first, get the areas and spatial polygons
ver<-
  lapply(1:length(BRBud), 
         function(i){
    cat(i, ' ')
    #cat(i/length(BRBud)," % complete   ")
    getverticeshr(BRBud[[i]],95)
  })
#ver
# now make the plot
plot(ver[[1]], xlim=c(315500, 325770), ylim=c(4385000, 4405500))
for(i in 2:length(ver)){
  plot(ver[[i]],add=T)
}

# get BRB foraging home range area
FHR<-rep(NA,length(ver))
for(i in 1:length(ver)){
  FHR[i]<-ver[[i]]$area
}
FHR
mean(FHR)
sd(FHR)


# get the probability contours from the BRB objects
BRBvolud<-getvolumeUD(BRBud, standardize = F)

# plot the UDs with probability contours
x11(width = 6, height = 6); par(mfrow=n2mfrow(length(BRBvolud)), mar=c(0,0,2,0))
lapply(1:length(BRBvolud), function(i) {
  image(BRBvolud[[i]],col=heat.colors(1000), main=names(BRBud)[i], xlim=c(315000, 325770), ylim=c(4385000, 4405500))
  xyz <- as.image.SpatialGridDataFrame(BRBvolud[[i]])
  contour(xyz, add=TRUE)
  box()
})

############################

# see an individual foraging range with the trajectory overlain

x11()
image(BRBvolud[[11]], col=heat.colors(1000),xlim=c(312621.2, 325773.3), ylim=c(4387950, 4397939))
contour(as.image.SpatialGridDataFrame(BRBvolud[[11]]),add=TRUE)
box()
par(new=T)
lines(trax[[11]])



############################

# look for spatial similarity of the roosting areas between years at selected percent home range areas
# note that this is calculated from the points, not the calculated UD!

wro<-kerneloverlaphr(BRBud, meth="UDOI", conditional=TRUE, percent = 95) # entire home range
wro

wro2<-as.matrix(wro)
wro2
m<-unique(fds$name)

rownames(wro2)=c(m);wro2 #assign row names
colnames(wro2)=c(m);wro2 #assign column names

#write.csv(wro2, file="ibat foraging overlap.csv")

las<-array(dim=dim(wro2))
for(i in 1:nrow(wro2)){
  las[i,]<-ifelse(wro2[i,]>0.99,1,0)
}

las
rownames(las)=c(m);las #assign row names
colnames(las)=c(m);las #assign column names

#write.csv(las, file="ibat foraging overlap.csv")

################################################################################################################################
################################################################################################################################

# save workspace
save.image(file = "ibat individual foraging.RData")

# delete saved workspace
#unlink("ibat individual foraging.RData")





