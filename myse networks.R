



#                                       Indiana Bat Network Analysis
#                                          - bipartite networks
#                                          - unipartite networks

# Code by Alex Silvis, 14 September 2012, using R version 2.15
# revised 27 September 2012
# revised 21 February 2013
# revised 11 April 2013
# revised 30 July 2013
# revised 12 August 2013
# revised 6 September 2013

## Code revised by Katie Gorman, 26 January 2022, using R version 4.0.3 for MYSE
## revised by KG, 28 June 2022, to run without juvenile





# set the working directory where all of the data files are located
## working directory set because working with a GIT project

#setwd("C:/Users/silvis/My Documents/R working directory/andrew ibats")
# setwd("//minnow.cc.vt.edu/cnre1/silvis/My Documents/R working directory/andrew ibats")

# load data from previous run
# load("ibat networks.RData", envir=parent.frame())
library(igraph)
library(tnet)
library(ggplot2)

save.image("myse networks.RDS")
save.image("myse networks adults.RDS")



load("myse networks.RDS")


# read in datashet in csv format
# this should be the sequential list of locations with animal/group identifier
#rts <- data.frame(read.table(file='ibatnetworkdata.csv', sep=',', header=TRUE))
# get the names of the variables in the data and check the first few lines of data

rts <- read.csv("myse_roosts.csv", header=TRUE, sep = ",", 
                quote = "\"", dec = ".", numerals = "no.loss", 
                fill = TRUE, comment.char = "")

t(names(rts))
head(rts)
tail(rts)

# remove bats not in the network
# determined on a previous run
# so that some calculations that follow are correct
#rts <- subset(rts, Bat != 15450, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)
#rts <- subset(rts, Bat != 15531, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)
#rts <- subset(rts, Bat != 18858, select=c(1:10))
#dim(rts)
#rts<-droplevels(rts)

#attach(rts)

# split data by year
#y09 <- subset(rts, Year == '2009', select=c(1:10))
#dim(y09)
#y09<-droplevels(y09)

#y10 <- subset(rts, Year == '2010', select=c(1:10))
#dim(y10)
#y10<-droplevels(y10)



## ------------------------------ ##
#### get some summary data first####
## ------------------------------ ##



# get tree use counts
# and merge with locations

# (tc09<-data.frame(table(y09$Tree)))
# names(tc09)<-c('Tree','Uses')
# 
# (tc10<-data.frame(table(y10$Tree)))
# names(tc10)<-c('Tree','Uses')

## kg does both years together
(tc09<-data.frame(table(rts$Tree)))
names(tc09)<-c('Tree','Uses')

(my09<-data.frame(rts$Tree,rts$X,rts$Y,rts$Year))
names(my09)<-c('Tree','X','Y','Year')
(my09<-unique(my09))

## how many times was each tree used, each year?
(rul09<-merge(tc09,my09, by.y='Tree'))
# 
# (my10<-data.frame(y10$Tree,y10$X,y10$Y,y10$Year))
# names(my10)<-c('Tree','X','Y','Year')
# (my10<-unique(my10))
# 
# (rul10<-merge(tc10,my10, by.y='Tree'))

(ttu<-rul09)

loguses<-log(ttu$Uses)

(ttu<-cbind(ttu,loguses))

# write merged files to .csv format for use elsewhere (like a GIS)

write.csv(ttu, file="tree uses.csv")

# roost uses in individual years
# mean(ttu$Uses[ttu$Year=='2009'])
# sd(ttu$Uses[ttu$Year=='2009'])
# 
# mean(ttu$Uses[ttu$Year=='2010'])
# sd(ttu$Uses[ttu$Year=='2010'])


mean(ttu$Uses) ## 2.967742
sd(ttu$Uses) ## 4.11083

## ------------------- ##
#### remove juvenile ####
## ------------------- ##

## just save rts df as rts.all but then create a new df for rts wihout juvie
rts.all <- rts

rts <- rts[rts$Age == "adult",] ## gets rid of 4 entries for juvenile male
unique(rts$Bat)

## rename RRR034 to 034034 because need all bats to be numeric for later
rts$Bat <- ifelse(rts$Bat == "RRR034", "34034", rts$Bat)


(tc09<-data.frame(table(rts$Tree)))
names(tc09)<-c('Tree','Uses')

(my09<-data.frame(rts$Tree,rts$X,rts$Y,rts$Year))
names(my09)<-c('Tree','X','Y','Year')
(my09<-unique(my09))

## how many times was each tree used, each year?
(rul09<-merge(tc09,my09, by.y='Tree'))
# 
# (my10<-data.frame(y10$Tree,y10$X,y10$Y,y10$Year))
# names(my10)<-c('Tree','X','Y','Year')
# (my10<-unique(my10))
# 
# (rul10<-merge(tc10,my10, by.y='Tree'))

(ttu<-rul09)

loguses<-log(ttu$Uses)

(ttu<-cbind(ttu,loguses))

# write merged files to .csv format for use elsewhere (like a GIS)

write.csv(ttu, file="tree uses adults.csv")

# roost uses in individual years
# mean(ttu$Uses[ttu$Year=='2009'])
# sd(ttu$Uses[ttu$Year=='2009'])
# 
# mean(ttu$Uses[ttu$Year=='2010'])
# sd(ttu$Uses[ttu$Year=='2010'])


mean(ttu$Uses) ## 3.142857
sd(ttu$Uses) ## 4.292235


## --------------------------- ##
#### get a bat x roost table ####
## --------------------------- ##

n09<-table(rts$Bat,rts$Tree)
n09<-as.matrix(n09)
n09

#n10<-table(y10$Bat,y10$Tree)
#n10<-as.matrix(n10)
#n10



t(names(rts))
y09 <- rts


# look at tracking periods
# plot.min.max(as.factor(y09$Bat), y09$daynum, dots=F)
# plot.min.max(as.factor(y10$Bat), y10$daynum, dots=F)
# 
# plot.min.max.cond(as.factor(y09$Bat), y09$daynum, y09$ReproByBat, exact=T)
# plot.min.max.cond(as.factor(y10$Bat), y10$daynum, y10$ReproByBat, exact=T)
# 
# plot.min.max.cond(as.factor(y09$Bat), y09$daynum, y09$Age, exact=T)
# plot.min.max.cond(as.factor(y10$Bat), y10$daynum, y10$Age, exact=T)

# mean number relocations
# 2009
(rs9<-rowSums(n09)) # get the number of locations by bat
mean(rs9) ## 5
sd(rs9) ## 3.273268

# 2010
#(rs10<-rowSums(n10))
# mean(rs10)
# sd(rs10)



## ------------------------- ##
#### look at tracking data ####
## ------------------------- ##


# 2009

#rowSums(t1) # this gives the same result as tabulate on the raw data
rc09<-array(dim=dim(n09))
for(i in 1:nrow(rc09)){
  rc09[i,]<-ifelse(n09[i,]>0,1,0)
}
# number of roosts used by a bat
(n.roosts09<-rowSums(rc09))
mean(n.roosts09) ## 3.13333
sd(n.roosts09) ## 1.684665

# number of relocations of a bat
(n.relocations09<-data.frame(table(y09$Bat)))
mean(n.relocations09$Freq) ## 5
sd(n.relocations09$Freq) ## 3.273268

cbind(n.relocations09,n.roosts09)

# roost switching frequency
rdr09<-n.relocations09$Freq/n.roosts09
(RSF09<-mean(rdr09)) ## 1.571111
sd(rdr09) ## 0.4631871


# 2010

#rowSums(t1) # this gives the same result as tabulate on the raw data
# rc10<-array(dim=dim(n10))
# for(i in 1:nrow(rc10)){
#   rc10[i,]<-ifelse(n10[i,]>0,1,0)
# }
# # number of roosts used by a bat
# (n.roosts10<-rowSums(rc10))
# mean(n.roosts10)
# sd(n.roosts10)
# 
# # number of relocations of a bat
# (n.relocations10<-data.frame(table(y10$Bat)))
# mean(n.relocations10$Freq)
# sd(n.relocations10$Freq)
# 
# cbind(n.relocations10,n.roosts10)
# 
# # roost switching frequency
# rdr10<-n.relocations10$Freq/n.roosts10
# (RSF10<-mean(rdr10))
# sd(rdr10)




# is the number of roosts used a function of the number of locations obtained?
rd09<-data.frame(n.roosts09, rs9, names(rs9))
names(rd09)<-c('NumRoosts','NumLocations','Bat')

brc1<-data.frame(y09$Bat, y09$ReproByBat)
brc1<-unique(brc1)
names(brc1)<-c('Bat','Repro')
brc1
summary(brc1[2])
(mrd09<-merge(brc1,rd09))

# rd10<-data.frame(n.roosts10, rs10, names(rs10), rep(2010, length(n.roosts10)))
# names(rd10)<-c('NumRoosts','NumLocations','Bat','Year')
# 
# brc2<-data.frame(y10$Bat, y10$ReproByBat)
# brc2<-unique(brc2)
# names(brc2)<-c('Bat','Repro')
# brc2
# summary(brc2[2])
# (mrd10<-merge(brc2,rd10))

# arsd<-rbind(mrd09,mrd10)

arsd <- mrd09

rum<-glm(NumRoosts~NumLocations+Repro, family=poisson(link='log'), data=arsd)
summary(rum)
Dsquared(rum)

#lm9<-glm(n.roosts09~rs9, family=poisson(link='log'))
#summary(lm9)
#Dsquared(lm9)
#confint(lm9)

# cross-validation
#library(boot)
#(cv.10.err9 <- cv.glm(data.frame(n.roosts09,rs9), lm9, K = 10)$delta) # 10-fold cross-validation




## ------------------------------------------- ##
#### prep data and create bipartite networks ####
## ------------------------------------------- ##


library(igraph)
t(names(rts))

n09<-table(y09$Bat,y09$Tree)
n09<-as.matrix(n09)
n09

# n10<-table(y10$Bat,y10$Tree)
# n10<-as.matrix(n10)
# n10

# 2009 network

## ----------------- ##
#### graph network ####
## ----------------- ##

n09<-graph.incidence(n09, directed=F, multiple=T)
n09
V(n09) # vertices
E(n09) # edges
#V(n09)$label <- V(n09)$name
V(n09)$label[1:length(V(n09))] <- NA
length(V(n09)$type)
vtypen<-get.vertex.attribute(n09, 'type')
V(n09)$color <- c(rep('#f1a340',length(vtypen[vtypen=='FALSE'])), rep('#998ec3',length(vtypen[vtypen=='TRUE'])))
E(n09)$color <- 'black'
plot(n09, layout=layout.reingold.tilford)
ew9<-count.multiple(n09)
E(n09)$width <- ew9
V(n09)$size <- 11.5
n09<-simplify(n09, edge.attr.comb=min)

tkplot(n09)
coord09<-tkplot.getcoords(1)

x11()
plot(n09, layout=coord09)
legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
       col=c('#f1a340','#998ec3'), pch=16, bty='n')

#V(n09)$type2 <- c(rep('bat',length(vtypen[vtypen=='FALSE']))), rep('roost', length(vtypen[vtypen=='TRUE'])))
#V(bbn)$size <- degree(bbn)+1
#tkplot(n09)

# 2010 network
# 
# n10<-graph.incidence(n10, directed=F, multiple=T)
# n10
# V(n10) # vertices
# E(n10) # edges
# vtypet<-get.vertex.attribute(n10, 'type')
# #V(n10)$label <- V(n10)$name
# V(n10)$label[1:length(V(n10))] <- NA
# vtypet<-get.vertex.attribute(n10, 'type')
# V(n10)$color <- c(rep('slateblue',length(vtypet[vtypet=='FALSE'])), rep('yellowgreen',length(vtypet[vtypet=='TRUE'])))
# E(n10)$color <- 'black'
# plot(n10, layout=layout.reingold.tilford)
# ew10<-count.multiple(n10)
# E(n10)$width <- ew10
# V(n10)$size <- 11.5
# 
# n10<-simplify(n10, edge.attr.comb=min)
# tkplot(n10)
# coord10<-tkplot.getcoords(2)
# 
# x11()
# plot(n10, layout=coord10)
# 
# #x11()
# plot(n10, layout=coord10)
# legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
#        col=c('slateblue','yellowgreen'), pch=16, bty='n')

#V(n10)$size <- degree(n10)+1
#tkplot(n10)




## ------------ ##
#### analysis ####
## ------------ ##

# 2009 analysis


# separate graphs from a single file so that only connected components remain in each group
dn09<-decompose.graph(n09, mode = 'weak',
                      max.comps = NA, min.vertices = 0)
dn09

# extract a single graph
n09<-(dn09[[1]])
n09
plot(n09)

n09.2<-simplify(n09)
plot(n09.2)
#V(n09)$size <- degree(n09)
#tkplot(n09)

# are there multiple edges between individual nodes
has.multiple(n09.2) # are there?
is.multiple(n09) # which ones?
count.multiple(n09) # how many edges per node?

################################################

# analyze the bipartite graph of the network

# look at the degree distribution
degd<-degree(n09)
hist(degd)

# but what about top and bottom node differences?

hist(degd[vtypen=='FALSE']) # bats
hist(degd[vtypen=='TRUE']) # roosts

mean(degd[vtypen=='FALSE']) # bats; 3.670666
sd(degd[vtypen=='FALSE']) # bats; 1.75119

mean(degd[vtypen=='TRUE'], na.rm = TRUE) # roosts; 1.833333
sd(degd[vtypen=='TRUE'], na.rm = TRUE) # roosts; 2.407717




# Number of top and bottom nodes
top<-length(V(n09)[type==FALSE]); top ## 14 bats
bottom<-length(V(n09)[type==TRUE]); bottom ## 25 roosts

# so, FALSE/top == bats
# and TRUE/bottom == roosts

# Number of edges
m<-ecount(n09); m ## 45

# Mean degree for top and bottom nodes
ktop<-m/top; ktop # bats; 3.214286
kbottom<-m/bottom; kbottom # roosts; 1.8

# Density for bipartite network
# must remove the multiple connections between bats and roosts
bipartite.density(simplify(n09)) ## 0.1285714

# Largest connected component for top and bottom nodes (used for mean distance for nodes below):
gclust<-clusters(n09, mode='weak')
lcc<-induced.subgraph(n09, V(n09)[gclust$membership==1])
lcctop<-length(V(lcc)[type==FALSE])
lccbottom<-length(V(lcc)[type==TRUE])
# Mean distance for top and bottom nodes
distop<-mean(shortest.paths(lcc, v=V(lcc)[type==FALSE], to=V(lcc)[type==FALSE], mode = 'all'))
disbottom<-mean(shortest.paths(lcc, v=V(lcc)[type==TRUE], to=V(lcc)[type==TRUE], mode = 'all'))
distop # bats; 2.346939
disbottom # roosts; 3.8976

# clustering as suggested in Tore Opsahl. Triadic closure in two-mode networks: Redefining 
# the global and local clustering coefficients. arXiv:1006.0887
library(tnet)
# the node set of primary interest must be the first list
# select the appropriate one as needed
nel<-data.frame(get.edgelist(n09))
ncl<-unique(nel)
clustering_tm(ncl) ## didn't work for kg

#### centrality and centralization as given in ####
# Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks 19, 243--269. 
# get centrality values for nodes in the bipartite graph
bipartite.betweenness.centrality(n09)
bipartite.closeness.centrality(n09)
bipartite.degree.centrality(n09)

# get centralization values for the bipartite graph
# so, FALSE/top == bats
# and TRUE/bottom == roosts
single.mode.betweenness.centralization(n09) ## bats = 0.1880565; roosts = 0.5763666
single.mode.closeness.centralization(n09) ## bats = 0.2650857; roosts =  0.750348
#single.mode.degree.centralization(nd) # nonsensical b/c normalization becomes meaningless
single.mode.degree.centralization(simplify(n09)) # simplify removes multi-edges; roosts = 0.7371795; bats = 0.125

## ------------------ ##
#### random network ####
## ------------------ ##

# do the values we observed differ from what we would observe in a random network?

(bgrpr09<-bi.graph.random.perm(n09,perms=10000,full=F))
bgrpr09
tkplot(n09)

################################################################################################################################
################################################################################################################################

# convert into unipartite graphs

unid<-bipartite.projection(n09, multiplicity = F, probe1 = NULL)
unid

#tkplot(unid[[1]]) # bats
#tkplot(unid[[2]]) # roosts

# select the network
ndt<-unid[[1]]
ndt
ndt<-simplify(ndt)
E(ndt)$color <- 'black'
plot(ndt)

#centralization.betweenness(ndt, normalized=T, loops=F)
centralization.degree(ndt, normalized=T, loops=F)
mean(degree(simplify(ndt)))

centralization.betweenness(ndt, normalized = TRUE)

graph.density(ndt,loops=F)
transitivity(ndt)
#V(ndt)$label

# select the network
ndt2<-unid[[2]]
ndt2
ndt2<-simplify(ndt2)
E(ndt2)$color <- 'black'
plot(ndt2)

#centralization.betweenness(ndt, normalized=T, loops=F)
centralization.degree(ndt2, normalized=T, loops=F)
mean(degree(simplify(ndt2)))

centralization.betweenness(ndt2, normalized = TRUE)

graph.density(ndt2,loops=F)
transitivity(ndt2)


# calculate shortest path length
# average
(aspl09<-average.path.length(ndt))
# all shortest paths
(spl09<-shortest.paths(ndt, mode='all', weight=NULL))
(spl09<-data.frame(spl09))

# determine whether dyads are more or less than the average distance
abapl09<-spl09-aspl09;abapl09

# convert distance to binary scale
las1<-array(dim=dim(abapl09))
for(i in 1:nrow(abapl09)){
  las1[i,]<-ifelse(abapl09[i,]>0,1,0)
}

las1
colnames(las1)=colnames(abapl09)
rownames(las1)=colnames(abapl09)
las1

#write.csv(las1, file="ibat network dist09.csv")


# define a few properties for plotting
# weight is turned off by the 'multiplicity = F' argument in bipartite.projection above
#E(ndt)$weight # the multiplicity of the edge
#E(ndt)$width <- E(ndt)$weight # make edge width equal to the edge weight
#E(ndt)$color <- 'wheat3'

plot(ndt)
#tkplot(ndt)

# measuring graph modularity - looking for subgroups within the graph

# modularity based on eigenvectors
lec<-leading.eigenvector.community(ndt, steps = -1, start = NULL, options = igraph.arpack.default, 
                                   callback = NULL, extra = NULL, env = parent.frame())
lec

V(ndt)$membership<-lec$membership

#colors<-ifelse(lec$membership=='1','black',ifelse(lec$membership=='2','red','blue'))
#V(ndt)$color <- colors

#plot(ndt)
#tkplot(ndt)


################################################

# do the values we observed differ from what we would observe in a random network?

grpr09<-graph.random.perm(ndt,perms=10000,full=F)
grpr09
#tkplot(ndt)

################################################################################################################################
################################################################################################################################

# test homophily

# show that vertex list of bats doesn't match the order of unique bats
bats<-unique(y09$Bat)
bats
V(n09)

t(names(y09))

# get the bat-reproductive condition combination data
al<-plot.min.max.cond(as.factor(y09$Bat), y09$daynum, y09$Age, exact=F)
al<-data.frame(al)
names(al) <- c('bat', 'minday', 'maxday', 'repro')
al
# need to order the data by bat; this will match the vertex id's
order(al$bat,al$repro)
A=data.frame(al$bat,al$repro)
A
B=A[order(al$bat,al$repro),]
B
# get vector of reproductive conditions
repro<-B$al.repro
repro

# set reproductive condition = 0 for trees (necessary to complete vertex attribute list)
(ts<-rep(0,length(vtypen[vtypen=='TRUE'])))
# combine bat and tree repro cond
rcl<-c(repro,ts)
# set repro cond as a vertex attribute
V(n09)$rcl <- rcl
V(n09)$rcl
# project bipartite graph
i95<-bipartite.projection(n09, multiplicity = F, probe1 = NULL)
i95
print(i95)
# select the bat-bat graph
abnet<-i95[[1]]
abnet
V(abnet)$rcl
# set node color based on reproductive condition
#colors<-ifelse(rcl=='1','black',ifelse(rcl=='2','red',ifelse(rcl=='4','blue','green')))
colors<-ifelse(V(abnet)$rcl=='1','slateblue4','lightblue')
V(abnet)$color <- colors
# red = non-reproductive, blue = pregnant, black = lactating, green = post-lactation
# matches the plot.min.max.cond plot

# set edge color
E(abnet)$color <- 'black'

# delete vertices with degree=0 (no connections) that will not be useful for assortativity
abnet <- delete.vertices(abnet, V(abnet)[ degree(abnet)==0 ])
abnet
V(abnet)$rcl

# check the graphs
x11(width = 8, height = 6.5)
plot(abnet, layout=layout.kamada.kawai)
legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Adult','Juvenile'), cex=1.2, 
       col=c('slateblue4','lightblue'), pch=19, bty='n')


# look at assortativity / look for differences in homiphily by age class
(ass09<-assortativity.nominal (abnet, V(abnet)$rcl, directed = F) )
(homt09<-homophily.perm.test(abnet, V(abnet)$rcl, 10000))


################################################################################################################################
################################################################################################################################
################################################################################################################################


# 2010 analysis

# separate graphs from a single file so that only connected components remain in each group
dn10<-decompose.graph(n10, mode = 'weak',
                      max.comps = NA, min.vertices = 0)
dn10

# extract a single graph
n10<-(dn10[[1]])
n10
plot(n10)
#V(n10)$size <- degree(n10)
#tkplot(n10)

# are there multiple edges between individual nodes
has.multiple(n10) # are there?
is.multiple(n10) # which ones?
count.multiple(n10) # how many edges per node?

################################################

# analyze the bipartite graph of the network

# look at the degree distribution
dege<-degree(n10)
hist(dege)

# but what about top and bottom node differences?

hist(dege[vtypet=='FALSE']) # bats
hist(dege[vtypet=='TRUE']) # roosts

mean(dege[vtypet=='FALSE']) # bats
sd(dege[vtypet=='FALSE']) # bats

mean(dege[vtypet=='TRUE']) # roosts
sd(dege[vtypet=='TRUE']) # roosts

# Number of top and bottom nodes
top<-length(V(n10)[type==FALSE]); top
bottom<-length(V(n10)[type==TRUE]); bottom

# so, FALSE/top == bats
# and TRUE/bottom == roosts

# Number of edges
m<-ecount(n10); m

# Mean degree for top and bottom nodes
ktop<-m/top; ktop # bats
kbottom<-m/bottom; kbottom # roosts

# Density for bipartite network
# must remove the multiple connections between bats and roosts
bipartite.density(simplify(n10))

# Largest connected component for top and bottom nodes (used for mean distance for nodes below):
gclust<-clusters(n10, mode='weak')
lcc<-induced.subgraph(n10, V(n10)[gclust$membership==1])
lcctop<-length(V(lcc)[type==FALSE])
lccbottom<-length(V(lcc)[type==TRUE])
# Mean distance for top and bottom nodes
distop<-mean(shortest.paths(lcc, v=V(lcc)[type==FALSE], to=V(lcc)[type==FALSE], mode = 'all'))
disbottom<-mean(shortest.paths(lcc, v=V(lcc)[type==TRUE], to=V(lcc)[type==TRUE], mode = 'all'))
distop # bats
disbottom # roosts

# clustering as suggested in Tore Opsahl. Triadic closure in two-mode networks: Redefining 
# the global and local clustering coefficients. arXiv:1006.0887
library(tnet)
# the node set of primary interest must be the first list
# select the appropriate one as needed
tel<-data.frame(get.edgelist(n10))
tcl<-unique(tel)
clustering_tm(tcl)

# centrality and centralization as given in 
# Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks 19, 243--269. 
# get centrality values for nodes in the bipartite graph
bipartite.betweenness.centrality(n10)
bipartite.closeness.centrality(n10)
bipartite.degree.centrality(n10)

# get centralization values for the bipartite graph
# so, FALSE/top == bats
# and TRUE/bottom == roosts
single.mode.betweenness.centralization(n10)
single.mode.closeness.centralization(n10)
#single.mode.degree.centralization(nd) # nonsensical b/c normalization becomes meaningless
single.mode.degree.centralization(simplify(n10)) # simplify removes multi-edges

################################################################################################################################
################################################################################################################################

# do the values we observed differ from what we would observe in a random network?

bgrpr10<-bi.graph.random.perm(n10,perms=10000,full=F)
bgrpr10
tkplot(n10)

################################################################################################################################
################################################################################################################################

# convert into unipartite graphs

unidt<-bipartite.projection(n10, multiplicity = F, probe1 = NULL)
unidt

#tkplot(unid[[1]]) # bats
#tkplot(unid[[2]]) # roosts

# select the network
tdt<-unidt[[1]]
tdt
tdt<-simplify(tdt)
E(tdt)$color <- 'black'
plot(tdt)

#centralization.betweenness(tdt, normalized=T)
centralization.degree(tdt, normalized=T, loops=F)
mean(degree(simplify(tdt)))
graph.density(tdt, loops=F)
transitivity(tdt)
#V(tdt)$label

average.path.length(tdt)
shortest.paths(tdt, mode='all', weight=NULL)

aspl10<-average.path.length(tdt);aspl10
spl10<-shortest.paths(tdt, mode='all', weight=NULL);spl10
spl10<-data.frame(spl10)

abapl10<-spl10-aspl10;abapl10

las2<-array(dim=dim(abapl10))
for(i in 1:nrow(abapl10)){
  las2[i,]<-ifelse(abapl10[i,]>0,1,0)
}

las2

colnames(las2)=colnames(abapl10)
rownames(las2)=colnames(abapl10)
las2

#write.csv(las2, file="ibat network dist10.csv")

# define a few properties for plotting
# weight is turned off by the 'multiplicity = F' argument in bipartite.projection above
#E(tdt)$weight # the multiplicity of the edge
#E(tdt)$width <- E(tdt)$weight # make edge width equal to the edge weight
#E(ndt)$color <- 'wheat3'

plot(tdt)
#tkplot(tdt)

# measuring graph modularity - looking for subgroups within the graph

# modularity based on eigenvectors
lec<-leading.eigenvector.community(tdt, steps = -1, start = NULL, options = igraph.arpack.default, 
                                   callback = NULL, extra = NULL, env = parent.frame())
lec

V(tdt)$membership<-lec$membership

#colors<-ifelse(lec$membership=='1','black',ifelse(lec$membership=='2','red','blue'))
#V(tdt)$color <- colors

#plot(tdt)
#tkplot(tdt)


################################################

# do the values we observed differ from what we would observe in a random network?

grpr10<-graph.random.perm(tdt,perms=10000,full=F)
grpr10
tkplot(tdt)

################################################################################################################################
################################################################################################################################

# show that vertex list of bats doesn't match the order of unique bats
bats2<-unique(y10$Bat)
bats2
V(n10)

t(names(y10))

# get the bat-reproductive condition combination data
al2<-plot.min.max.cond(as.factor(y10$Bat), y10$daynum, y10$Age, exact=F)
al2<-data.frame(al2)
names(al2) <- c('bat', 'minday', 'maxday', 'repro')
al2
# need to order the data by bat; this will match the vertex id's
order(al2$bat,al2$repro)
A2=data.frame(al2$bat,al2$repro)
A2
B2=A2[order(al2$bat,al2$repro),]
B2
# get vector of reproductive conditions
repro2<-B2$al2.repro
repro2

# set reproductive condition = 0 for trees (necessary to complete vertex attribute list)
(ts2<-rep(0,length(vtypen[vtypet=='TRUE'])))
# combine bat and tree repro cond
rcl2<-c(repro2,ts2)
# set repro cond as a vertex attribute
V(n10)$rcl <- rcl2
V(n10)$rcl
# project bipartite graph
i96<-bipartite.projection(n10, multiplicity = F, probe1 = NULL)
i96
print(i96)
# select the bat-bat graph
abnet2<-i96[[1]]
abnet2
V(abnet2)$rcl
# set node color based on reproductive condition
#colors2<-ifelse(rcl=='1','black',ifelse(rcl=='2','red',ifelse(rcl=='4','blue','green')))
colors2<-ifelse(V(abnet2)$rcl=='1','slateblue4','lightblue')
V(abnet2)$color <- colors2
# red = non-reproductive, blue = pregnant, black = lactating, green = post-lactation
# matches the plot.min.max.cond plot

# set edge color
E(abnet2)$color <- 'black'

# delete vertices with degree=0 (no connections) that will not be useful for assortativity
abnet2 <- delete.vertices(abnet2, V(abnet2)[ degree(abnet2)==0 ])
abnet2
V(abnet2)$rcl

# check the graphs
# check the graphs
x11(width = 8, height = 6.5)
plot(abnet2, layout=layout.kamada.kawai)
legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Adult','Juvenile'), cex=1.2, 
       col=c('slateblue4','lightblue'), pch=19, bty='n')


# look at assortativity / look for differences in homiphily by age class
assortativity.nominal (abnet2, V(abnet2)$rcl, directed = F) 

hom10<-homophily.perm.test(abnet2, V(abnet2)$rcl, 10000)
hom10

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

# plot the 2009 and 2010 bipartite networks together

x11(width = 9, height = 6)
par(mfrow=c(1,2))
par(mar=c(1,1,2,1))

plot(n09)
title('2009'); box()
plot(n10)
title('2010'); box()


################################################################################################################################
################################################################################################################################

# random node removal 2009

# to perform a single random removal experiment

# removal.sim(input network, number of simulations, number nodes to remove)

ndtt<-unid[[2]]

removal.sim(ndtt, 3, 3)

# conduct iterative random removal simulations and plot the results
# batch.rm.sim(netowrk, number simulations, vector of the number of nodes to remove)
y=rep(1:23)
b1<-batch.rm.sim(ndtt,500,y)

b1$Num.Clusters.Per.Sim
mn<-b1$Mean.No.Clusters
se<-b1$SE.No.Clusters
prp<-b1$Proportion.Removed
# prp has too many digits, so we round to 2
prp<-round(prp, digits = 2)


upr<-mn+se
lwr<-mn-se

x11();plot(prp,mn, type='b', xlab="Proportion of Roosts Removed", ylab="Number of Networks", font.lab=2, pch=19, cex=0.75)
arrows(prp, mn, prp, lwr, angle = 90, length = 0.1)
arrows(prp, mn, prp, upr, angle = 90, length = 0.1)


################################################################################################################################
################################################################################################################################

# random node removal 2010

# to perform a single random removal experiment

# removal.sim(input network, number of simulations, number nodes to remove)

tdtt<-unidt[[2]]

removal.sim(tdtt, 3, 3)

# conduct iterative random removal simulations and plot the results
# batch.rm.sim(netowrk, number simulations, vector of the number of nodes to remove)
y2=rep(1:11)
b2<-batch.rm.sim(tdtt,500,y2)

b2$Num.Clusters.Per.Sim
mn2<-b2$Mean.No.Clusters
se2<-b2$SE.No.Clusters
prp2<-b2$Proportion.Removed
# prp has too many digits, so we round to 2
prp2<-round(prp2, digits = 2)


upr2<-mn2+se2
lwr2<-mn2-se2

x11();plot(prp2, mn2, type='b', xlab="Proportion Roosts Removed", ylab="Number of Networks", font.lab=2, pch=19, cex=0.75)
arrows(prp2, mn2, prp2, lwr2, angle = 90, length = 0.1)
arrows(prp2, mn2, prp2, upr2, angle = 90, length = 0.1)

################################################################################################################################
################################################################################################################################

# figure with both removal simulations on the same graph

x11();plot(prp,mn, type='b', xlab="Proportion of Roosts Removed", ylab="Number of Networks", font.lab=2, pch=19, cex=0.75)
arrows(prp, mn, prp, lwr, angle = 90, length = 0.1)
arrows(prp, mn, prp, upr, angle = 90, length = 0.1)

points(prp2,mn2,pch=1)
lines(prp2,mn2)
arrows(prp2, mn2, prp2, lwr2, angle = 90, length = 0.1)
arrows(prp2, mn2, prp2, upr2, angle = 90, length = 0.1)

legend("topleft",legend=c("2009","2010"),pch=c(19,1))

# now make a better looking version

library(ggplot2)

(prps<-c(prp,prp2))
(mnvs<-c(mn,mn2))
(Year<-c(rep('2009',length(mn)),rep('2010',length(mn2))))
(uper<-c(upr,upr2))
(lwer<-c(lwr,lwr2))
(simpl<-data.frame(prps,mnvs,Year,uper,lwer))

cbPalette <- c("#56B4E9", "#E69F00")


x11(width=6, height=5)
ggplot(simpl, aes(x=prps, y=mnvs, colour=Year)) + 
  geom_errorbar(aes(ymin=lwer, ymax=uper), width=.01) +  geom_line() +  geom_point() + 
  labs(fill='Year') + theme_bw(base_size=12, base_family="") + xlab('Proportion of Roosts Removed') +
  ylab('Number of Networks') + theme(legend.justification=c(1,0), legend.position=c(0.2,0.75)) +
  scale_color_manual(values=cbPalette) 

################################################################################################################################
################################################################################################################################

# targeted node removal
# removing the primary roost (tree 6)

plot(ndtt)
rpn09<-delete.vertices(ndtt,1)
plot(rpn09)

plot(tdtt)
rpn10<-delete.vertices(tdtt,1)
plot(rpn10)

################################################################################################################################
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
#font_import()
#loadfonts(device="postscript")
fonts()
library(grDevices)
library(grid)


##################

# the removal simulation

#cbPalette <- c("#56B4E9", "#E69F00")

pdf(file="removalsimulation.pdf",family='sans', onefile=TRUE,useDingbats=F,bg='white',
           paper="special",width=4.86,height=3.5, pointsize=12)

#x11(width=4.86, height= 3.5)

ggplot(simpl, aes(x=prps, y=mnvs, colour=Year)) + 
  geom_errorbar(aes(ymin=lwer, ymax=uper), width=.01) +  geom_line() +  geom_point() + 
  labs(fill='Year') + theme_bw(base_size=10, base_family="") + xlab('Proportion of Roosts Removed') +
  ylab('Number of Networks') + 
  theme(legend.key=element_rect(colour='white'), legend.justification=c(1,0), legend.position=c(0.22,0.62),
        legend.text = element_text(colour="black", size = 8), 
        legend.background = element_rect(colour = "black", size=0.5), 
        legend.title = element_text(colour="black", size = 8)) +
  scale_color_manual(values=cbPalette) 

dev.off()
#library(grDevices)

##################

# bipartite networks

#postscript(file="bipartitenetworks.eps",family='sans',horizontal=FALSE, onefile=TRUE,
#           paper="special",width=6.83,height=4.4)


pdf(file="bipartitenetworks.pdf",family='sans', onefile=TRUE,useDingbats=F,bg='white',
    paper="special",width=6.83,height=4.4,pointsize=12)

#x11(width=6.83, height= 4.4)

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

plot(n09, layout=coord09)
box()
text(-1.005411, 1.366885, expression(bold('A')))
#title('2009')

#legend("bottomleft", inset=c(0.02,0.02), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
#       col=c('slateblue','yellowgreen'), pch=16, bty='n')
#legend("bottomleft", inset=c(0.02,0.02), "(x,y)", legend=c('',''), cex=1.2, 
#       col=c('black','black'), pch=1, bty='n')

plot(n10, layout=coord10)
box()
text(-1.005411, 1.366885, expression(bold('B')))
#title('2010')
legend("bottomright", inset=c(0.02,0.0000000005), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
       col=c('slateblue','yellowgreen'), pch=16, bty='n')
#legend("bottomright", inset=c(0.02,0.02), "(x,y)", legend=c('',''), cex=1.2, 
#       col=c('black','black'), pch=1, bty='n')


dev.off()
#embedFonts(file="bipartitenetworks.eps",format="pswrite", outfile="bipartitenetworksEF.eps")

##################

#postscript(file="unipartitenetworks.eps",family='sans',horizontal=FALSE, onefile=TRUE,
#           paper="special",width=6.83,height=4.4)

pdf(file="unipartitenetworks.pdf",family='sans', onefile=TRUE,useDingbats=F,bg='white',
    paper="special",width=6.83,height=4.4,pointsize=12)

#x11(width=6.83, height= 4.4)

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

plot(abnet)
text(-1.005411, 1.366885, expression(bold('A')))
#title('2009')
box()

#legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Adult','Juvenile'), cex=1.2, 
#       col=c('slateblue4','lightblue'), pch=19, bty='n')
#legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('',''), cex=1.2, 
#       col=c('black','black'), pch=1, bty='n')

plot(abnet2)
text(-1.005411, 1.366885, expression(bold('B')))
#title('2010')
box()

legend("bottomright", inset=c(0.02,0.0000000005), "(x,y)", legend=c('Adult','Juvenile'), cex=1.2, 
       col=c('slateblue4','lightblue'), pch=19, bty='n')
#legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('',''), cex=1.2, 
#       col=c('black','black'), pch=1, bty='n')

dev.off()
#embedFonts(file="unipartitenetworks.eps",format="pswrite", outfile="unipartitenetworksEF.eps")

################################################################################################################################
################################################################################################################################

# save workspace
save.image(file = "ibat networks.RData")

# delete saved workspace
#unlink("ibat networks.RData")






