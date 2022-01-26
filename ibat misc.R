

setwd("C:/Users/silvis/My Documents/R working directory/andrew ibats")

# load saved workspace from previous run
load("ibat misc.RData", envir=parent.frame())

library(ggplot2)
library(igraph)

################################################################################################################################
################################################################################################################################

# check for correlation between exit count and bat roost telemetry data

#ED<-read.csv("emergence3.csv", header=TRUE)

t(names(ED))
dim(ED)

#ED <- subset(ED, treecode != 't24', select=c(1:10))
#ED <- subset(ED, treecode != 't39', select=c(1:10))
#ED <- subset(ED, treecode != 't41', select=c(1:10))
#ED <- subset(ED, treecode != 't42', select=c(1:10))
#ED <- subset(ED, treecode != 't37', select=c(1:10))
#ED <- subset(ED, treecode != 't51', select=c(1:10))
#ED <- subset(ED, treecode != 't52', select=c(1:10))
#ED <- subset(ED, treecode != 't53', select=c(1:10))
#dim(ED)

ED

n09<-table(y09$Bat,y09$Tree)
n09

#rowSums(t1) # this gives the same result as tabulate on the raw data
rc09<-array(dim=dim(n09))
for(i in 1:nrow(rc09)){
  rc09[i,]<-ifelse(n09[i,]>0,1,0)
}
rc09

# number of bats that used a roost
(n.bats09<-colSums(rc09))
mean(n.bats09)
sd(n.bats09)

cbind(unique(y09$Tree),n.bats09)

n10<-table(y10$Bat,y10$Tree)
n10

#rowSums(t1) # this gives the same result as tabulate on the raw data
rc10<-array(dim=dim(n10))
for(i in 1:nrow(rc10)){
  rc10[i,]<-ifelse(n10[i,]>0,1,0)
}
rc10

# number of bats that used a roost
(n.bats10<-colSums(rc10))
mean(n.bats10)
sd(n.bats10)

cbind(unique(y10$Tree),n.bats10)



rbu09<-data.frame(cbind(unique(y09$Tree),n.bats09))
names(rbu09)<-c('treenum','batsusingtree09')
rbu10<-data.frame(cbind(unique(y10$Tree),n.bats10))
names(rbu10)<-c('treenum','batsusingtree10')

dim(ED)
(new<-merge(ED, rbu09, by='treenum', all.x=T, all.y=F))
dim(test)

(new2<-merge(test, rbu10, by='treenum', all.x=T, all.y=F))
dim(new2)

new2
ED


cbind(c(ED$treenum,ED$treenum),c(ED$roostdays09,ED$roostdays10),c(ED$maxemerge09,ED$Maxemerge10))

################################################################################################################################
################################################################################################################################

# plot emergence vs several observed items

#plot(ED$maxemerge,ED$totroostdays)
#plot(ED$maxemerge,ED$batsusingtree)

summary(ED$maxemerg09)
summary(as.factor(ED$maxemerg09))
summary(ED$roostdays09)
summary(as.factor(ED$roostdays09))

summary(ED$maxemerge10)
summary(as.factor(ED$maxemerge10))
summary(ED$roostdays10)
summary(as.factor(ED$roostdays10))

# correlation between max emergence and bats using tree is unfair, b/c the bats in the tree obviously
# are related to the number in the tree, so disregard these correlations

#cor.test(ED$totroostdays,ED$maxemerge, method='pear')
#cor.test(ED$maxemerge,ED$batsusingtree, method='pear')

# does the number of days a roost was used reflect the number of bats in the tree? a fair question

cor.test(ED$roostdays09,ED$maxemerg09,method='pear')
cor.test(ED$roostdays10,ED$maxemerge10, method='pear')


################################################################################################################################
################################################################################################################################

# create  roost tree data graph
# but only for roosts with max emergence counts
ED<-subset(ED, maxemerge != '0', select=c(1:10))
ED

t(names(ED))
attach(ED)
ED2<-data.frame(treecode,totroostdays,maxemerge,batsusingtree)
detach(ED)
head(ED2)
names(ED2)<-c('Roost','Total Roost Days','Maximum Emergence','Number Tagged Bats')
formula(ED2)
(ED2<-stack(ED2))
head(ED2)
ED2<-cbind(rep(ED$treenum-5,3),ED2)
head(ED2)
names(ED2)<-c('Roost','Count','CountType')
head(ED2)

library(ggplot2)

t(names(ED2))

# create custom color palette
#colors()
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cc<-colorRampPalette(c('darkseagreen4','goldenrod','seashell4')) (3)
cbPalette <- c('darkslategray', '#E69F00', '#56B4E9')

# check to see if we like the colors
pie(rep(1,3), col=cbPalette) 

x11(width=7.5,height=7)
ggplot(data=ED2, aes(x=factor(Roost), y=Count, fill=CountType)) + geom_bar(stat="identity",position='dodge') + 
  xlab('Roost number (in order of observed use)') + scale_fill_manual(values=cbPalette) + 
  labs(fill='Count Type') + theme_bw(base_size=12, base_family="") + 
  theme(legend.justification=c(1,0), legend.position=c(1,0.72)) 


################################################################################################################################
################################################################################################################################

t(names(new2))

# create  roost tree data graph
# but only for roosts with max emergence counts
new3<-subset(new2, maxemerg09 != '0', select=c(1:15))
new3

t(names(new3))
attach(new3)
new3<-data.frame(treecode,roostdays09,maxemerg09,batsusingtree09.x)
detach(new3)
head(new3)
names(new3)<-c('RoostC','Total Roost Days','Maximum Emergence','Number Tagged Bats')
formula(new3)
head(new3)
(new4<-stack(new3))
head(new4)
new4<-cbind(rep(new3$Roost,3),new4)
head(new4)
new4<-cbind(new4, rep(c(6,7,8,11,12,13,15,19,25,26,40),3), rep('2009',33))
names(new4)<-c('RoostC','Count','CountType','Roost','Year')
head(new4)

library(ggplot2)

t(names(new3))

# create custom color palette
#colors()
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cc<-colorRampPalette(c('darkseagreen4','goldenrod','seashell4')) (3)
cbPalette <- c('darkslategray', '#E69F00', '#56B4E9')

# check to see if we like the colors
pie(rep(1,3), col=cbPalette) 

x11(width=7.5,height=7)
ggplot(data=new4, aes(x=factor(Roost), y=Count, fill=CountType)) + geom_bar(stat="identity",position='dodge') + 
  xlab('Roost number (in order of observed use)') + scale_fill_manual(values=cbPalette) + 
  labs(fill='Count Type') + theme_bw(base_size=12, base_family="") + 
  theme(legend.justification=c(1,0), legend.position=c(1,0.72)) 

################################################################################################################################
################################################################################################################################


t(names(new2))

# create  roost tree data graph
# but only for roosts with max emergence counts
new5<-subset(new2, maxemerge10 != '0', select=c(1:15))
new5

t(names(new5))
attach(new5)
new5<-data.frame(treecode,roostdays10,maxemerge10,batsusingtree10.x)
detach(new5)
head(new5)
names(new5)<-c('RoostC','Total Roost Days','Maximum Emergence','Number Tagged Bats')
formula(new5)
head(new5)
(new6<-stack(new5))
head(new6)
new6<-cbind(rep(new5$Roost,3),new6)
head(new6)
new6<-cbind(new6, rep(c(6,7,8,26,49,50,55),3), rep('2010',21))
names(new6)<-c('RoostC','Count','CountType','Roost','Year')
head(new6)

library(ggplot2)

t(names(new6))

# create custom color palette
#colors()
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#cc<-colorRampPalette(c('darkseagreen4','goldenrod','seashell4')) (3)
cbPalette <- c('darkslategray', '#E69F00', '#56B4E9')

# check to see if we like the colors
pie(rep(1,3), col=cbPalette) 

x11(width=7.5,height=7)
ggplot(data=new6, aes(x=factor(Roost), y=Count, fill=CountType)) + geom_bar(stat="identity",position='dodge') + 
  xlab('Roost number (in order of observed use)') + scale_fill_manual(values=cbPalette) + 
  labs(fill='Count Type') + theme_bw(base_size=12, base_family="") + 
  theme(legend.justification=c(1,0), legend.position=c(1,0.72)) 


################################################################################################################################
################################################################################################################################

head(new4)
head(new6)

new7<-rbind(new4,new6)

x11(width=4.86,height=4.86)
ggplot(data=new7, aes(x=factor(Roost), y=Count, fill=CountType)) + geom_bar(stat="identity",position='dodge') + 
  xlab('Roost number (in order of observed use)') + scale_fill_manual(values=cbPalette) + 
  labs(fill='Count Type') + theme_bw(base_size=12, base_family="") + 
  facet_grid(Year ~ ., scales = "fixed") +
  theme(legend.justification=c(1,0), legend.position=c(1,0.75)) 


################################################################################################################################
################################################################################################################################

################################################################################################################################
################################################################################################################################

# create publication quality image

library(grid)
library(grDevices)
library(extrafont)
#font_import()
loadfonts(device="postscript")
fonts()

# for Windows - in each session
# adjust the path to match the installation of Ghostscript
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.07/bin/gswin64c.exe")

pdf(file="roostcounts.pdf",family='sans', onefile=TRUE,useDingbats=F,bg='white',
    paper="special",width=4.86,height=4.86, pointsize=12)

#x11(width=4.86, height=4.86)
#par(mar=c(1,1,0,1))

ggplot(data=new7, aes(x=factor(Roost), y=Count, fill=CountType)) + geom_bar(stat="identity",position='dodge') + 
  xlab('Roost number (in order of observed use)') + scale_fill_manual(values=cbPalette) + 
  labs(fill='Count Type') + theme_bw(base_size=10, base_family="") + 
  facet_grid(Year ~ ., scales = "fixed") +
  theme(legend.justification=c(1,0), 
        legend.text = element_text(colour="black", size = 8), legend.background = element_rect(), 
        legend.title = element_text(colour="white", size = 0), legend.position="top",
        legend.key.size= unit(.25, 'cm'))

dev.off()



################################################################################################################################
################################################################################################################################

# check to see if bats that had foraging areas that overlapped more than expected by chance were also of 
# closer than average social proximity in the bat-bat social network
# social proximity and UDOI calculated in other scripts; data reduced in MS Access (sorry!) to appropriate pairs

#(pl09<-read.csv("09pathlength.csv", header=TRUE,row.names=1))
#(pl10<-read.csv("10pathlength.csv", header=TRUE, row.names=1))
#(UDOI09<-read.csv("09UDOI.csv", header=TRUE, row.names=1))
#(UDOI10<-read.csv("10UDOI.csv", header=TRUE, row.names=1))

dim(pl09)
dim(UDOI09)
dim(pl10)
dim(UDOI10)

# in imported path length data, 0 represents social proximity closer than average
# in imported UDOI data, 0 represents less foraging area overlap than expected
# for the following analysis, we want the path length values and UDOI values to be of the same interpretation
# so we subtract 1 from the path length dataframes, then multiply by -1 to accomplish this

(pl09<-(pl09-1)*-1)
(pl10<-(pl10-1)*-1)

# diagonals show an individual bats overlap/proximity to itself, so we need to remove the diagonal
# similarly, because we have a square matrix, every combination is represented twice, so we
# need to drop one side of the matrix

pl09
upper.tri(pl09,diag=T)
pl09[upper.tri(pl09,diag=T)]=NA
pl09

pl10[upper.tri(pl10,diag=T)]=NA
pl10

UDOI09[upper.tri(UDOI09,diag=T)]=NA
UDOI09

UDOI10[upper.tri(UDOI10,diag=T)]=NA
UDOI10

# do some simple addition of the corresponding path length and UDOI data 
# 2 represents more foraging area overlap than expected AND closer than average social proximity
# 1 represents EITHER more foraging area overlap OR closer than average social proximity
# 0 represents NEITHER more foraging area overlap OR closer than average social proximity

(d09<-as.matrix(pl09+UDOI09))
(d10<-as.matrix(pl10+UDOI10))

# look at the counts of individual instances
table(d09)
table(d10)

# get the proportion of the instances
table(d09)/45
table(d10)/171

# check to see how many bats had more overlap than expected in foraging area
table(as.matrix(UDOI09))
table(as.matrix(UDOI10))

# the following test isn't relevant
#chisq.test(table(d09))
#chisq.test(table(d10))

################################################################################################################################
################################################################################################################################

# plot a network representation of the foraging association network

library(igraph)

# 2009
fn09<-simplify(graph.adjacency(UDOI09, mode='undirected'))
V(fn09)$label[1:length(V(fn09))] <- NA
E(fn09)$color <- 'black'
V(fn09)$size <- 11.5
V(fn09)$color <- c(rep('slateblue1',length(V(fn09))))
plot(fn09)
tkplot(fn09)
coord09<-tkplot.getcoords(1)
plot(fn09, layout=coord09)

# remove unconnected nodes
#fn09 <- delete.vertices(fn09, V(fn09)[ degree(fn09)==0 ])
#plot(fn09)

# 2010
fn10<-simplify(graph.adjacency(UDOI10, mode='undirected'))
V(fn10)$label[1:length(V(fn10))] <- NA
E(fn10)$color <- 'black'
V(fn10)$size <- 11.5
V(fn10)$color <- c(rep('slateblue1',length(V(fn10))))
plot(fn10)
tkplot(fn10)
coord10<-tkplot.getcoords(2)
plot(fn10, layout=coord10)

# remove unconnected nodes
#fn10 <- delete.vertices(fn10, V(fn10)[ degree(fn10)==0 ])
#plot(fn10)

# to prove that graph.adjacency works for this
# degree should match column sums minus 1 (minus 1 removes the bat
# connection with itself)
degree(fn09)
colSums(UDOI09)-1

degree(fn10)
colSums(UDOI10)-1

############################
############################

# network analysis of the foraging networks
# don't do unless the unconnected nodes are removed

# 2009
mean(degree(simplify(fn09)))
centralization.degree(fn09, normalized=T, loops=F)
graph.density(fn09,loops=F)
transitivity(fn09)

graph.random.perm(fn09,500,full=F)

# 2010
mean(degree(simplify(fn10)))
centralization.degree(fn10, normalized=T, loops=F)
graph.density(fn10,loops=F)
transitivity(fn10)

graph.random.perm(fn10,500,full=F)


################################################################################################################################
################################################################################################################################

# create publication quality image

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.07/bin/gswin64c.exe")

library(grid)
library(grDevices)
library(extrafont)
#font_import()
loadfonts(device="postscript")
fonts()

#postscript(file="foragingnetworks.eps",family='sans',horizontal=FALSE, onefile=TRUE,
#           paper="special",width=6.83,height=4.4)

pdf(file="foragingnetworks.pdf",family='sans', onefile=TRUE,useDingbats=F,bg='white',
    paper="special",width=6.83,height=4.4, pointsize=12)

#x11(width=6.83, height= 4.4)

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

plot(fn09, layout=coord09)
text(-1.005411, 1.366885, expression(bold('A')))
#title('2009')
box()

#legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('Adult','Juvenile'), cex=1.2, 
#       col=c('slateblue4','lightblue'), pch=19, bty='n')
#legend("topleft", inset=c(0.02,0.02), "(x,y)", legend=c('',''), cex=1.2, 
#       col=c('black','black'), pch=1, bty='n')

plot(fn10, layout=coord10)
text(-1.005411, 1.366885, expression(bold('B')))
#title('2010')
box()

dev.off()

################################################################################################################################
################################################################################################################################


# save workspace
save.image(file = "ibat misc.RData")

# delete saved workspace
#unlink("ibat misc.RData")


















