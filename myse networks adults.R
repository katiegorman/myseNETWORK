



#                                       Indiana Bat Network Analysis
#                                          - bipartite networks
#                                          - unipartite networks

## FIRST RUN *NETWORK ANALYSIS FUNCTIONS.R*

# Code by Alex Silvis, 14 September 2012, using R version 2.15
# revised 27 September 2012
# revised 21 February 2013
# revised 11 April 2013
# revised 30 July 2013
# revised 12 August 2013
# revised 6 September 2013

## Code revised by Katie Gorman, 26 January 2022, using R version 4.0.3 for MYSE
## revised by KG, 28 June 2022, to run without juvenile
## also to make all nodes numeric so functions run properly for whatever dumb reason
## revised by KG, 29 June 2022, to remove UNCONNECTED adult female


## leaving?
save.image("myse networks adults.RDS")

## coming back?
load("myse networks adults.RDS")

## working directory set because working with a GIT project

library(igraph)
library(tnet)
library(ggplot2)

## some pathnames
path.root <- "G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks"  # network analysis folder
path.git <- paste(path.root, "/myseNETWORK/myseNETWORK", sep = "") ## actual working directory for git proj
path.figs <- paste(path.root, "/NEW figures June 2022", sep = "") ## new figs


# read in datashet in csv format
# this should be the sequential list of locations with animal/group identifier

# get the names of the variables in the data and check the first few lines of data

rts <- read.csv("myse_roosts_numbers.csv", header=TRUE, sep = ",", 
                quote = "\"", dec = ".", numerals = "no.loss", 
                fill = TRUE, comment.char = "")

t(names(rts))
## cols are: Bat|Tree|X|Y|Date|daynum|ReproByBat|Age|Sex|Year
head(rts)
tail(rts)

bat.days <- data.frame(table(rts$Bat))
bat.days

names(bat.days)<-c('Bat','Days tracked')

bat.days.year <- data.frame(rts$Bat, rts$Year)
names(bat.days.year)<-c('Bat','Year')
bat.days.year <- unique(bat.days.year)

tot.bat.days <- merge(bat.days, bat.days.year, by.y = "Bat")

## fix the bats that were tracked for 2 years
## 72406 used 6 in 2018 and 2 in 2019
## 72411 used 4 in 2018 and 10 in 2019

tot.bat.days$`Days tracked`[7] <- 6
tot.bat.days$`Days tracked`[8] <- 2
tot.bat.days$`Days tracked`[13] <- 4
tot.bat.days$`Days tracked`[14] <- 10

## check
tot.bat.days

mean(tot.bat.days$`Days tracked`) ## 4.388889
sd(tot.bat.days$`Days tracked`) ## 2.304443

## take away unconnected bats
tot.bat.days.connect <- tot.bat.days[(1:16),]
mean(tot.bat.days.connect$`Days tracked`) ## 4.4375
sd(tot.bat.days.connect$`Days tracked`) ##2.448639

## ------------------------------- ##
#### get some summary data first ####
## ------------------------------- ##

# get tree use counts
# and merge with locations

## kg does both years together
tc09<-data.frame(table(rts$Tree))
names(tc09)<-c('Tree','Uses')

my09<-data.frame(rts$Tree,rts$X,rts$Y,rts$Year)
names(my09)<-c('Tree','X','Y','Year')
my09<-unique(my09)

## how many times was each tree used, each year?
rul09<-merge(tc09,my09, by.y='Tree')

## that didn't work for kg
## manually go in and change 1806 to be used 12x in 2018, 1x in 2019
rul09$Uses[6] <- 12
rul09$Uses[7] <- 1

## check
rul09


ttu<-rul09

loguses<-log(ttu$Uses)

ttu<-cbind(ttu,loguses)

# write merged files to .csv format for use elsewhere (like a GIS)

write.csv(ttu, file="tree uses all bats.csv")



mean(ttu$Uses) ## 2.548387
sd(ttu$Uses) ## 3.585514

## --------------------------------------- ##
#### remove juvenile and unconnected bat ####
## --------------------------------------- ##

## just save rts df as rts.all but then create a new df for rts wihout juvie
rts.all <- rts

rts <- rts[rts$Age == "adult",] ## gets rid of 4 entries for juvenile male
unique(rts$Bat)

## also get rid of 72414 because she wasn't connected
rts <- rts[!rts$Bat == "72414",] ## gets rid of 4 entries for that bat

tc09<-data.frame(table(rts$Tree))
names(tc09)<-c('Tree','Uses')

my09<-data.frame(rts$Tree,rts$X,rts$Y,rts$Year)
names(my09)<-c('Tree','X','Y','Year')
my09<-unique(my09)

## how many times was each tree used, each year?
rul09<-merge(tc09,my09, by.y='Tree')


## again: manually go in and change 1806 to be used 12x in 2018, 1x in 2019
rul09$Uses[6] <- 12
rul09$Uses[7] <- 1

## check
rul09

ttu<-rul09

loguses<-log(ttu$Uses)

ttu<-cbind(ttu,loguses)

# write merged files to .csv format for use elsewhere (like a GIS)

write.csv(ttu, file="tree uses connected adults.csv")


mean(ttu$Uses) ## 3.142857
sd(ttu$Uses) ## 4.292235
min(ttu$Uses) ## 1
max(ttu$Uses) ## 17

## --------------------------- ##
#### get a bat x roost table ####
## --------------------------- ##

n09<-table(rts$Bat,rts$Tree)
n09<-as.matrix(n09)
n09

## write matrix to csv for use in ucinet
write.csv(n09, 
          file = "G:/My Drive/FIRE, FIRE ISLAND!/Analysis/Networks/UCINET2/myse_connect_matrix.csv", 
          row.names = TRUE)


t(names(rts))
y09 <- rts


rs9<-rowSums(n09) # get the number of locations by bat
mean(rs9) ## 5.071429
sd(rs9) ## 3.384678
min(rs9) ## 1-14



## ------------------------- ##
#### look at tracking data ####
## ------------------------- ##


## rowSums(t1) # this gives the same result as tabulate on the raw data
rc09<-array(dim=dim(n09))
for(i in 1:nrow(rc09)){
  rc09[i,]<-ifelse(n09[i,]>0,1,0)
}
# number of roosts used by a bat
n.roosts09<-rowSums(rc09)
mean(n.roosts09) ## 3.214286
sd(n.roosts09) ## 1.717716

# number of relocations of a bat
n.relocations09<-data.frame(table(y09$Bat))
mean(n.relocations09$Freq) ## 5.071429
sd(n.relocations09$Freq) ## 3.384678

cbind(n.relocations09,n.roosts09)

# roost switching frequency
rdr09<-n.relocations09$Freq/n.roosts09
RSF09<-mean(rdr09) 
RSF09 ## 1.540476
sd(rdr09) ## 0.4646345


# is the number of roosts used a function of the number of locations obtained?
rd09<-data.frame(n.roosts09, rs9, names(rs9))
names(rd09)<-c('NumRoosts','NumLocations','Bat')

brc1<-data.frame(y09$Bat, y09$ReproByBat)
brc1<-unique(brc1)
names(brc1)<-c('Bat','Repro')
brc1
summary(brc1[2])
mrd09<-merge(brc1,rd09)


arsd <- mrd09

rum<-glm(NumRoosts~NumLocations+Repro, family=poisson(link='log'), data=arsd)
summary(rum)
Dsquared(rum)
confint(rum)


## YES!!
## p = 0.007 for NumLocations but not repro


## ------------------------------------------- ##
#### prep data and create bipartite networks ####
## ------------------------------------------- ##


library(igraph)
t(names(rts))

n09<-table(y09$Bat,y09$Tree)
n09<-as.matrix(n09)
n09



## ----------------- ##
#### graph network ####
## ----------------- ##

n09<-graph.incidence(n09, directed=F, multiple=T)
n09
V(n09) # vertices
E(n09) # edges
V(n09)$label <- V(n09)$name
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
coord09<-tkplot.getcoords(15)

x11()
plot(n09, layout=coord09)
legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
       col=c('#f1a340','#998ec3'), pch=16, bty='n')

setwd(path.figs)

pdf( "bivariate_network_weights.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white"          # Background color
)

plot(n09, layout=coord09)
legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Bat','Roost'), cex=1.2, 
       col=c('#f1a340','#998ec3'), pch=16, bty='n')

dev.off()

setwd(path.git)

## ------------ ##
#### analysis ####
## ------------ ##


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


# are there multiple edges between individual nodes
has.multiple(n09) # are there?
is.multiple(n09) # which ones?
count.multiple(n09) # how many edges per node?


## --------------------------------------------- ##
#### analyze the bipartite graph of the network ###
## --------------------------------------------- ##

# look at the degree distribution
degd<-degree(n09)
hist(degd)

# but what about top and bottom node differences?

hist(degd[vtypen=='FALSE']) # bats
hist(degd[vtypen=='TRUE']) # roosts

mean(degd[vtypen=='FALSE']) # bats; 3.214286
sd(degd[vtypen=='FALSE']) # bats; 1.717716

mean(degd[vtypen=='TRUE'], na.rm = TRUE) # roosts; 1.8
sd(degd[vtypen=='TRUE'], na.rm = TRUE) # roosts; 2.362908




# Number of top and bottom nodes
top<-length(V(n09)[type==FALSE]); top ## 14 bats
bottom<-length(V(n09)[type==TRUE]); bottom ## 25 roosts

# so, FALSE/top == bats
# and TRUE/bottom == roosts

# Number of edges
m<-ecount(n09); m ## 45

# Mean degree for top and bottom nodes
## these values are the same as degree distribution above
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

## here you are getting two cols according to all 45 connections
nel<-data.frame(get.edgelist(n09))

## then just taking unique ones, so no back and forth bats
ncl<-unique(nel)

## doesn't look like there were any
clustering_tm(ncl) ## 0.1168831

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

set.seed(1234)
bgrpr09<-bi.graph.random.perm(n09,perms=10000,full=F)
bgrpr09
n09


save.image("myse networks adults.RDS")

## ---------------------------------- ##
#### convert into unipartite graphs ####
## ---------------------------------- ##

unid<-bipartite.projection(n09, multiplicity = F, probe1 = NULL)
unid

# tkplot(unid[[1]]) # bats
# tkplot(unid[[2]]) # roosts

# select the network
ndt<-unid[[1]]
ndt
ndt<-simplify(ndt)
E(ndt)$color <- 'black'
plot(ndt)

tkplot(ndt)
coordndt<-tkplot.getcoords(16) ## change this number according to window

x11()
plot(ndt, layout=coordndt)
legend("topleft", inset=c(0.02,0.02), "(x)", legend=c('Bats'), cex=1.2, 
       col=c('#f1a340'), pch=16, bty='n')



setwd(path.figs)

pdf( "univariate_network_bats_legend.pdf",         # File name
     width = 8, height = 7, # Width and height in inches
     bg = "white"          # Background color
)


plot(ndt, layout=coordndt)
legend("topleft", inset=c(0.02,0.02), "(x)", legend=c('Bats'), cex=1.2, 
       col=c('#f1a340'), pch=16, bty='n')

dev.off()

setwd(path.git)


#centralization.betweenness(ndt, normalized=T, loops=F)
centralization.degree(ndt, normalized=T, loops=F) ## centralization = 0.3076923
mean(degree(simplify(ndt))) ## 9.571429

centralization.betweenness(ndt, normalized = TRUE) ## centralization = 0.8678501

graph.density(ndt,loops=F) ## 0.7362637
transitivity(ndt) ## 0.8846154
#V(ndt)$label

# select the network
ndt2<-unid[[2]]
ndt2
ndt2<-simplify(ndt2)
E(ndt2)$color <- 'black'
plot(ndt2)


tkplot(ndt2)
coordndt2<-tkplot.getcoords(17) ## change this number according to window

x11()
plot(ndt2, layout=coordndt2)
legend("topleft", inset=c(0.02,0.02), "(x)", legend=c('Roosts'), cex=1.2, 
       col=c('#998ec3'), pch=16, bty='n')



setwd(path.figs)

## do it once with legened and once without
## change file name accordingly
pdf( "univariate_network_roosts_legend.pdf",         # File name
     width = 8, height = 7, # Width and height in inches
     bg = "white"          # Background color
)


plot(ndt2, layout=coordndt2)
legend("topleft", inset=c(0.02,0.02), "(x)", legend=c('Roosts'), cex=1.2,
       col=c('#998ec3'), pch=16, bty='n')

dev.off()

setwd(path.git)




#centralization.betweenness(ndt, normalized=T, loops=F)
centralization.degree(ndt2, normalized=T, loops=F) ## 0.548913
mean(degree(simplify(ndt2))) ## 4.88

centralization.betweenness(ndt2, normalized = TRUE) ## 0.4722977

graph.density(ndt2,loops=F) ## 0.2033333
transitivity(ndt2) ## 0.4454756


# calculate shortest path length
# average
aspl09<-average.path.length(ndt) 
aspl09 ## 1.26 3736

# all shortest paths
spl09<-shortest.paths(ndt, mode='all', weight=NULL)
spl09<-data.frame(spl09)
spl09

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

## write.csv if you want
## i think 0s are bats that are closer than avg


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




# save workspace
save.image(file = "myse networks adults.RDS")

## delete saved workspace
unlink("myse networks adults.RDS")






