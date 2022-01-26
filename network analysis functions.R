
##############################################################################################################
# calculate the proportion of deviance explained in a glm
# this is the R squared value, aka "D squared"
# see Guisan et al. 2000. Predictive habitat distribution models in ecology. Ecological Modelling 135: 147-186 
Dsquared<-function(fit){
  D2<-((fit$null.deviance-fit$deviance)/fit$null.deviance)
  MAD2<-(1-((fit$df.null)/(fit$df.residual))*(1-D2))
  return((c("Max.Adj.D^2"=MAD2, "D^2"=D2, "Null deviance"=fit$null.deviance,"Null df"=fit$df.null, 
            "Residual deviance"=fit$deviance, "Residual df"=fit$df.residual)))
}
##############################################################################################################
##############################################################################################################




##############################################################################################################
##############################################################################################################

# calculate standard error
# regular formulation, required below

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

##############################################################################################################
##############################################################################################################

library(igraph)

bipartite.betweenness.centrality<-function(g){
  if (!is.null(V(g)$type)){
    # determine maximal raw scores for both vertex subsets
    if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
      mrs_TRUE <- 2*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
    }
    else{
      mrs_TRUE <- 0.5*length(V(g)[type==FALSE])*(length(V(g)[type==FALSE])-1)+0.5*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-2)+(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
    }
    if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
      mrs_FALSE <- 0.5*length(V(g)[type==TRUE])*(length(V(g)[type==TRUE])-1)+0.5*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-2)+(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
      
    }
    else{
      mrs_FALSE <- 2*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
    }
    
    # get raw betweenness centrality scores from igraph
    betweenness_rs <- betweenness(g,directed=FALSE)
    # "bipartite" normalization of scores
    for (i in V(g)){
      if (V(g)[i]$type==TRUE){
        V(g)[i]$betweenness.centrality <- betweenness_rs[i]/mrs_TRUE
      }
      else{
        V(g)[i]$betweenness.centrality <- betweenness_rs[i]/mrs_FALSE
      }
    }
    # return value as list
    return(list("Bipartite.Betweenness.Centrality"=V(g)$betweenness.centrality))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

#bipartite.betweenness.centrality(netwd)

##############################################################################################################
##############################################################################################################


bipartite.closeness.centrality<-function(g){
  if (!is.null(V(g)$type)){
    # determine maximal raw scores for both vertex subsets
    mrs_TRUE <- length(V(g)[type==FALSE]) + 2*length(V(g)[type==TRUE]) - 2
    mrs_FALSE <- length(V(g)[type==TRUE]) + 2*length(V(g)[type==FALSE]) - 2
    # get sum of all geodesic paths for each vertex
    rowsums_shortest_paths <- rowSums(shortest.paths(g))
    # "bipartite" normalization of scores
    for (i in V(g)){
      if (V(g)[i]$type==TRUE){
        V(g)[i]$closeness.centrality <- mrs_TRUE/rowsums_shortest_paths[i]
      }
      else{
        V(g)[i]$closeness.centrality <- mrs_FALSE/rowsums_shortest_paths[i]
      }
    }
    # return value as list
    return(list("Bipartite.Closeness.Centrality"=V(g)$closeness.centrality))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

#bipartite.closeness.centrality(netwd)

##############################################################################################################
##############################################################################################################


bipartite.degree.centrality<-function(g,loops=FALSE){
  # if boolean vertex attribute <type> is present, calculate bipartite degree centrality, otherwise monopartite degree centrality
  if (!is.null(V(g)$type)){
    for (i in V(g)){
      V(g)[i]$degree.centrality <- degree(g,v=i)/length(V(g)[type==!V(g)[i]$type])
    }
    # return value vector as list item
    return(list("Bipartite.Degree.Centrality"=V(g)$degree.centrality))
  }
  else {
    for (i in V(g)){
      if (!loops){
        V(g)[i]$degree.centrality <- degree(g,v=i,loops=FALSE)/(length(V(g))-1)
      }
      else{
        V(g)[i]$degree.centrality <- degree(g,v=i,loops=TRUE)/(length(V(g)))
      }
    }
    # return value vector as list item
    return(list("Monopartite.Degree.Centrality"=V(g)$degree.centrality))
  }
}

#bipartite.degree.centrality(netwd)

##############################################################################################################
##############################################################################################################


bipartite.density<-function(g){
  if (!is.null(V(g)$type)){
    # return value as list item
    if (!is.directed(g)){
      return(list("Density"=length(E(g))/(length(which(V(g)$type))*length(which(!V(g)$type)))))
    }
    else{
      return(list("Density"=length(E(g))/(2*length(which(V(g)$type))*length(which(!V(g)$type)))))
    }
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

#bipartite.density(netwd)

####

# # an equivalent way to calculate density in a bipartite network
# 
# # Number of top and bottom nodes
# top<-length(V(g)[type==FALSE])
# bottom<-length(V(g)[type==TRUE])
# # Number of edges
# m<-ecount(g)
# # Mean degree for top and bottom nodes
# ktop<-m/top
# kbottom<-m/bottom
# # Density for bipartite network
# bidens<-m/(top*bottom)


##############################################################################################################
##############################################################################################################


single.mode.betweenness.centralization<-function(g){
  if (!is.null(V(g)$type)){
    # determine denominators for both vertex subsets
    if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
      denom_TRUE <- 2*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
    }
    else{
      denom_TRUE <- (length(V(g)[type==TRUE])-1) * ( 0.5*length(V(g)[type==FALSE])*(length(V(g)[type==FALSE])-1) + 0.5*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-2) + (length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1) )
    }
    if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
      denom_FALSE <- 2*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
    }
    else{
      denom_FALSE <- (length(V(g)[type==FALSE])-1) * ( 0.5*length(V(g)[type==TRUE])*(length(V(g)[type==TRUE])-1) + 0.5*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-2) + (length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1) )
    }
    
    # determine raw betweenness centrality scores from igraph
    V(g)$betweenness.raw <- betweenness(g,directed=FALSE)
    # get maximum scores
    max_c_TRUE <- max(V(g)[type==TRUE]$betweenness.raw)
    max_c_FALSE <- max(V(g)[type==FALSE]$betweenness.raw)
    
    # determine centralization for TRUE vertex subset
    g$single.mode.betweenness.centralization.true <-sum(max_c_TRUE - V(g)[type==TRUE]$betweenness.raw)/denom_TRUE
    # determine centralization for FALSE vertex subset
    g$single.mode.betweenness.centralization.false <-sum(max_c_FALSE - V(g)[type==FALSE]$betweenness.raw)/denom_FALSE
    # return both values as list
    return(list("Single.Mode.Betweenness.Centralization.TRUE"=g$single.mode.betweenness.centralization.true,"Single.Mode.Betweenness.Centralization.FALSE"=g$single.mode.betweenness.centralization.false))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

#single.mode.betweenness.centrality(netwd)

##############################################################################################################
##############################################################################################################


single.mode.closeness.centralization<-function(g){
  if (!is.null(V(g)$type)){
    # determine bipartite closeness centrality scores
    V(g)$Bipartite.Closeness.Centrality <- bipartite.closeness.centrality(g)[[1]]
    
    # get maximum scores
    max_c_TRUE <- max(V(g)[type==TRUE]$Bipartite.Closeness.Centrality)
    max_c_FALSE <- max(V(g)[type==FALSE]$Bipartite.Closeness.Centrality)
    
    # determine denominators for both vertex subsets
    if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
      denom_TRUE <- ((length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-2))/(2*length(V(g)[type==TRUE])-3) + ((length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-length(V(g)[type==FALSE])))/(length(V(g)[type==TRUE])+length(V(g)[type==FALSE])-2)
    }
    else{
      denom_TRUE <- ((length(V(g)[type==TRUE])-2)*(length(V(g)[type==TRUE])-1))/(2*length(V(g)[type==TRUE])-3)
    }
    if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
      denom_FALSE <- ((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-2))/(2*length(V(g)[type==FALSE])-3) + ((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-length(V(g)[type==TRUE])))/(length(V(g)[type==FALSE])+length(V(g)[type==TRUE])-2)
    }
    else{
      denom_FALSE <- ((length(V(g)[type==FALSE])-2)*(length(V(g)[type==FALSE])-1))/(2*length(V(g)[type==FALSE])-3)
    }
    
    # determine centralization for TRUE vertex subset
    g$single.mode.closeness.centralization.true <-sum(max_c_TRUE - V(g)[type==TRUE]$Bipartite.Closeness.Centrality)/denom_TRUE
    # determine centralization for FALSE vertex subset
    g$single.mode.closeness.centralization.false <-sum(max_c_FALSE - V(g)[type==FALSE]$Bipartite.Closeness.Centrality)/denom_FALSE
    # return both values as list
    return(list("Single.Mode.Closeness.Centralization.TRUE"=g$single.mode.closeness.centralization.true,"Single.Mode.Closeness.Centralization.FALSE"=g$single.mode.closeness.centralization.false))
  }
  else {
    # boolean vertex attribute 'type' is required
    cat("vertex attribute <type> is missing")
  }
}

#single.mode.closeness.centrality(netwd)

##############################################################################################################
##############################################################################################################


single.mode.degree.centralization <-
  function(g){
    if (!is.null(V(g)$type)){
      V(g)$degree <- degree(g)
      # determine maximum degrees for each vertex set
      max_d_TRUE <- max(V(g)[type==TRUE]$degree)
      max_d_FALSE <- max(V(g)[type==FALSE]$degree)
      # determine centralization for TRUE vertex subset
      g$single.mode.degree.centralization.true <-sum(max_d_TRUE - V(g)[type==TRUE]$degree)/((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1))
      # determine centralization for FALSE vertex subset
      g$single.mode.degree.centralization.false <-sum(max_d_FALSE - V(g)[type==FALSE]$degree)/((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1))
      # return both values as list
      return(list("Single.Mode.Degree.Centralization.TRUE"=g$single.mode.degree.centralization.true,"Single.Mode.Degree.Centralization.FALSE"=g$single.mode.degree.centralization.false))
    }
    else {
      # boolean vertex attribute 'type' is required
      cat("vertex attribute <type> is missing")
    }
  }

#single.mode.degree.centralization(netwd)


##############################################################################################################
##############################################################################################################


# compare (specific) observed values to values from random bipartite ER graphs

bi.graph.random.perm<-function(graph,perms,full){
  require(igraph)
  require(tnet)
  require(tcltk)
  graph<-simplify(graph)
  perms<-perms
  # get probability of connection for randomization
  x<-get.edgelist(graph)
  x<-data.frame(unique(x))
  n1<-length(unique(x$X1))
  n2<-length(unique(x$X2))
  ptie<-nrow(x)/(n1*n2)
  
  # create empty slots for permutations
  degreetop<-rep(NA, length(perms))
  degreebottom<-rep(NA, length(perms))
  density<-rep(NA, length(perms))
  centralization<-rep(NA, length(perms))
  clustering<-rep(NA, length(perms))
  
  # calculate observed values
  # degree
  o.ktop<-(ecount(graph))/(length(V(graph)[type==FALSE]))
  o.kbottom<-(ecount(graph))/(length(V(graph)[type==TRUE]))
  # density
  ds<-bipartite.density(simplify(graph))
  o.density<-ds$Density
  # degree centralization
  t1<-single.mode.degree.centralization(simplify(graph))
  o.dcentralization<-t1$Single.Mode.Degree.Centralization.TRUE
  # clustering
  grel<-get.edgelist(graph)
  grel<-data.frame(unique(grel))
  o.clustering<-clustering_tm(grel)
  
  # create list to hold randomized networks
  df.list <- vector("list", perms)
  
  # create progress bar
  pb <- tkProgressBar(title = "Progress", min = 0,
                      max = perms, width = 300)
  
  set.seed(1142) # for reproducibility
  for(i in 1:perms){
    
    # get progress bar updates
    setTkProgressBar(pb, i, label=paste( round(i/perms*100,0),"% done"))
    
    # randomize graph
    df.list[[i]]<-rg_tm(ni=n1,np=n2,ties=ptie)
    require(igraph)
    tb<-table(df.list[[i]])
    pg<-graph.incidence(tb, directed=F)
    pg<-simplify(pg)
    
    # calculate values
    # degree
    degreetop[i]<-(ecount(pg))/(length(V(pg)[type==FALSE]))
    degreebottom[i]<-(ecount(pg))/(length(V(pg)[type==TRUE]))
    
    # density
    density[i]<-bipartite.density(pg)
    
    # centralization
    tmp<-single.mode.degree.centralization(pg)
    centralization[i]=tmp$Single.Mode.Degree.Centralization.TRUE
    
    # clustering
    require(tnet)
    clustering[i]<-clustering_tm(df.list[[i]])
    
  }
  close(pb)
  
  # calculate Monte Carlo two-tailes p values
  kt<-mean(as.numeric(degreetop))
  kt2<-as.numeric(degreetop)-as.numeric(kt)
  o.ktop2<-as.numeric(o.ktop)-as.numeric(kt)
  pktop<-(1+sum(abs(kt2) >= abs(o.ktop2)))/(perms+1)
  
  kb<-mean(as.numeric(degreebottom))
  kb2<-as.numeric(degreebottom)-as.numeric(kb)
  o.kbottom2<-as.numeric(o.kbottom)-as.numeric(kb)
  pkbottom<-(1+sum(abs(kb2) >= abs(o.kbottom2)))/(perms+1)
  
  d<-mean(as.numeric(density))
  density2<-as.numeric(density)-as.numeric(d)
  o.density2<-as.numeric(o.density)-as.numeric(d)
  pdensity<-(1+sum(abs(density2) >= abs(o.density2)))/(perms+1)
  
  ct<-mean(as.numeric(centralization))
  cent2<-as.numeric(centralization)-as.numeric(ct)
  o.cent2<-o.dcentralization-ct
  pcent<-(1+sum(abs(cent2) >= abs(o.cent2)))/(perms+1)
  
  cl<-mean(as.numeric(clustering),na.rm=T)
  clust2<-as.numeric(clustering)-as.numeric(cl)
  o.clust2<-o.clustering-cl
  pclust<-(1+sum(abs(clust2) >= abs(o.clust2)))/(perms+1)
  
  # get standard error for simulated values
  std.ktop<-stderr(as.numeric(degreetop))
  std.kbottom<-stderr(as.numeric(degreebottom))
  std.density<-stderr(as.numeric(density))
  std.centralization<-stderr(as.numeric(centralization))
  std.clustering<-stderr(as.numeric(clustering))
  
  # plot the simulation value density plot with the observed values
  # corresponds to p value but also shows direction of difference
  x11(width = 14, height = 3)
  par(mfrow=c(1,5))
  (plot(density(as.numeric(degreetop)),main="Mean Top Degree", xlim=c(0,15)))
  (abline(v=o.ktop, lty=2))
  (plot(density(as.numeric(degreebottom)),main="Mean Bottom Degree", xlim=c(0,15)))
  (abline(v=o.kbottom, lty=2))
  (plot(density(as.numeric(density)),main="Density", xlim=c(0,1)))
  (abline(v=o.density, lty=2))
  (plot(density(as.numeric(centralization)),main="Degree Centralization",xlim=c(0,1)))
  (abline(v=o.dcentralization, lty=2))
  (plot(density(as.numeric(clustering),na.rm=T),main="Clustering",xlim=c(0,1)))
  (abline(v=o.clustering, lty=2)) 
  
  if(full == T){  
    return((c("Obs Degree Top"=o.ktop, "Mean SIM ktop"=kt, "SIM ktop SE"=std.ktop, "ktop p.val"=pktop,
              "obs Degree Bottom"=o.kbottom, "Mean SIM kbottom"=kb, "SIM kbottom SE"=std.kbottom, "kbottom p.val"=pkbottom,
              "Obs Density"=o.density, "Mean SIM Density"=d,"SIM Density SE"=std.density, "Density p.val"=pdensity,
              "Obs D.Centralization"=o.dcentralization, "Mean SIM Centralization"=ct, "SIM D.Centralization SE"=std.centralization, "D.Centralization p.val"=pcent,
              "Obs Clustering"=o.clustering, "Mean SIM Clustering"=cl, "SIM Clustering SE"=std.clustering,  "Clustering p.val"=pclust
    )))  
  }
  
  
  if(full == F){  
    return((c("Obs Degree Top"=o.ktop, "ktop p.val"=pktop,
              "Obs Degree Bottom"=o.kbottom, "kbottom p.val"=pkbottom,
              "Obs Density"=o.density, "Density p.val"=pdensity,
              "Obs D.Centralization"=o.dcentralization, "D.Centralization p.val"=pcent,
              "Obs Clustering"=o.clustering, "Clustering p.val"=pclust
    )))
  }
  
}
#test
#bi.graph.random.perm(,perms=50)

##############################################################################################################
##############################################################################################################


# compare (specific) observed values to values from random unipartite ER graphs

graph.random.perm<-function(graph,perms,full){
  require(igraph)
  require(tnet)
  require(tcltk)
  graph<-simplify(graph)
  x<-get.edgelist(graph)
  x<-data.frame(unique(x))
  perms<-perms
  
  # calculate probabilities for randomizations
  # set to graph density
  ptie<-graph.density(graph, loops=F)
  
  # calculate observed values
  o.density<-graph.density(graph)
  t1<-centralization.degree(graph, normalized=T, loops=F)
  o.dcentralization<-t1$centralization
  o.clustering<-transitivity(graph)
  o.path<-average.path.length(graph)
  o.degree<-mean(degree(graph))
  
  t2<-leading.eigenvector.community(graph, steps = -1, start = NULL, options = igraph.arpack.default, 
                                    callback = NULL, extra = NULL, env = parent.frame())
  o.modularity<-t2$modularity
  
  # create empty slots for simulation values
  density<-rep(NA, length(perms))
  centralization<-rep(NA, length(perms))
  clustering<-rep(NA, length(perms))
  path<-rep(NA, length(perms))
  modularity<-rep(NA, length(perms))
  degree<-rep(NA, length(perms))
  
  # create list
  df.list <- vector("list", perms)
  
  # create progress bar
  pb <- tkProgressBar(title = "Progress", min = 0,
                      max = perms, width = 300)
  
  set.seed(1142) # makes the results reproducible
  for(i in 1:perms){
    
    # get data for the progress bar
    setTkProgressBar(pb, i, label=paste( round(i/perms*100,0),"% done"))
    
    # randomize network
    df.list[[i]]<-rg_w(nodes = length(V(graph)), arcs = ptie, directed = FALSE)
    df.list
    #pg<-graph.edgelist(as.matrix(df.list[[i]]), directed=F)
    pg<-tnet_igraph(df.list[[i]],type=NULL, directed=FALSE)
    pg<-simplify(pg)
    
    # calculate randomized values
    # density
    density[i]<-graph.density(pg)
    
    # degree centralization
    tmp<-centralization.degree(pg, normalized=T, loops=F)
    centralization[i]=tmp$centralization
    
    # clustering
    clustering[i]<-transitivity(pg)
    
    # average path length
    path[i]<-average.path.length(pg)
    
    # mean degree
    degree[i]<-mean(degree(pg))
    
    # leading eigenvector modularity
    tmp2<-leading.eigenvector.community(pg, steps = -1, start = NULL, options = igraph.arpack.default, 
                                        callback = NULL, extra = NULL, env = parent.frame())
    modularity[i]<-tmp2$modularity
    
  }
  close(pb)
  
  # calculate Monte Carlo two-tailed P values
  g<-mean(as.numeric(degree))
  degree2<-as.numeric(degree)-as.numeric(g)
  o.degree2<-as.numeric(o.degree)-as.numeric(g)
  pdegree<-(1+sum(abs(degree2) >= abs(o.degree2)))/(perms+1)
  
  d<-mean(as.numeric(density))
  density2<-as.numeric(density)-as.numeric(d)
  o.density2<-as.numeric(o.density)-as.numeric(d)
  pdensity<-(1+sum(abs(density2) >= abs(o.density2)))/(perms+1)
  
  ct<-mean(as.numeric(centralization))
  cent2<-as.numeric(centralization)-as.numeric(ct)
  o.cent2<-o.dcentralization-ct
  pcent<-(1+sum(abs(cent2) >= abs(o.cent2)))/(perms+1)
  
  cl<-mean(as.numeric(clustering))
  clust2<-as.numeric(clustering)-as.numeric(cl)
  o.clust2<-o.clustering-cl
  pclust<-(1+sum(abs(clust2) >= abs(o.clust2)))/(perms+1)
  
  pl<-mean(as.numeric(path))
  path2<-as.numeric(path)-as.numeric(pl)
  o.path2<-as.numeric(o.path)-as.numeric(pl)
  ppath<-(1+sum(abs(path2) >= abs(o.path2)))/(perms+1)
  
  ml<-mean(as.numeric(modularity))
  modularity2<-as.numeric(modularity)-as.numeric(ml)
  o.modularity2<-as.numeric(o.modularity)-as.numeric(ml)
  pmodularity<-(1+sum(abs(modularity2) >= abs(o.modularity2)))/(perms+1)
  
  # standard errors of the simulation values
  std.degree<-stderr(as.numeric(degree))
  std.density<-stderr(as.numeric(density))
  std.centralization<-stderr(as.numeric(centralization))
  std.clustering<-stderr(as.numeric(clustering))
  std.path<-stderr(as.numeric(path))
  std.modularity<-stderr(as.numeric(modularity))
  
  # plot the simulation value density plot with the observed values
  # corresponds to p value but also shows direction of difference
  x11(width = 15, height = 3)
  par(mfrow=c(1,6))
  (plot(density(as.numeric(degree)),main="Mean Degree",xlim=c(0,10)))
  (abline(v=o.degree, lty=2))
  (plot(density(as.numeric(path)),main="Average Path Length"))
  (abline(v=o.path, lty=2))
  (plot(density(as.numeric(density)),main="Density",xlim=c(0,1)))
  (abline(v=o.density, lty=2))
  (plot(density(as.numeric(centralization)),main="Degree Centralization",xlim=c(0,1)))
  (abline(v=o.dcentralization, lty=2))
  (plot(density(as.numeric(clustering)),main="Clustering",xlim=c(0,1)))
  (abline(v=o.clustering, lty=2))
  (plot(density(as.numeric(modularity)),main="Modularity",xlim=c(0,1)))
  (abline(v=o.modularity, lty=2))
  
  
  # return results
  if(full == T){  
    return((c("Obs Mean Degree"=o.degree, "Mean SIM Degree"=g,"SIM Degree SE"=std.degree,"Degree p.val"=pdegree,
              "Obs Avg Path"=o.path, "Mean SIM Avg Path"=pl, "SIM Avg Path SE"=std.path, "Path p.val"=ppath,
              "Obs Density"=o.density, "Mean SIM Density"=d,"SIM Density SE"=std.density, "Density p.val"=pdensity,
              "Obs D.Centralization"=o.dcentralization, "Mean SIM Centralization"=ct, "SIM D.Centralization SE"=std.centralization, "D.Centralization p.val"=pcent,
              "Obs Clustering"=o.clustering, "Mean SIM Clustering"=cl, "SIM Clustering SE"=std.clustering,  "Clustering p.val"=pclust,
              "Obs L.Eig Modularity"=o.modularity, "Mean L.Eig Modularity"=ml, "SIM Modularity SE"=std.modularity, "L.Eig Modularity p.val"=pmodularity
    )))
  }
  
  if(full == F){  
    return((c("Obs Degree"=o.degree, "Degree p.val"=pdegree,
              "Obs Avg Path"=o.path, "Path p.val"=ppath,
              "Obs Density"=o.density, "Density p.val"=pdensity,
              "Obs D.Centralization"=o.dcentralization, "D.Centralization p.val"=pcent,
              "Obs Clustering"=o.clustering, "Clustering p.val"=pclust,
              "Obs L.Eig Modularity"=o.modularity, "L.Eig Modularity p.val"=pmodularity
    )))
  }
  
}
# test
#graph.random.perm(tdt,500,full=T)

##############################################################################################################
##############################################################################################################


# a function to assess whether observed homophily for our network is significant
# by permuting the node label, NOT the network

homophily.perm.test<-function(graph,cond,perms){
  require(igraph)
  require(tcltk)
  
  graph = simplify(graph)
  
  o.hom <- assortativity.nominal (graph, cond, directed = F)
  
  hom<-rep(NA, perms)
  tmp<-rep(NA,length(cond))
  
  # create progress bar
  pb <- tkProgressBar(title = "Progress", min = 0,
                      max = perms, width = 300)
  
  set.seed(1142) # makes the results reproducible
  for(i in 1:perms){
    # get progress bar updates
    setTkProgressBar(pb, i, label=paste( round(i/perms*100,0),"% done"))
    
    tmp <- sample(cond, length(cond), replace=F)
    hom[i] <- assortativity.nominal (graph, tmp, directed = F)
  }
  close(pb)
  
  # calculate Monte Carlo two-tailed p values
  mh<-mean(as.numeric(hom))
  mh2<-as.numeric(hom)-as.numeric(mh)
  o.hom2<-as.numeric(o.hom)-as.numeric(mh)
  phom<-(1+sum(abs(mh2) >= abs(o.hom2)))/(perms+1)
  
  # get standard error for simulated values
  std.hom<-stderr(as.numeric(hom))
  
  # plot the simulation value density plot with the observed values
  # corresponds to p value but also shows direction of difference
  x11()
  (plot(density(as.numeric(hom)),main="Mean Homophily", xlim=c(-1,1)))
  (abline(v=o.hom, lty=2))
  
  return(c("Obs Homophily"=o.hom, "Mean SIM Homophily"=mh, "SIM Homophily SE"=std.hom, "Homophily p.val"=phom))  
  
}

#homophily.perm.test(abnet, V(abnet)$rcl, 10000)

##############################################################################################################
##############################################################################################################


# a function to randomly remove nodes in a unipartite and return values from the resulting networks
# set.seed is used to generate reproducible results
# a function to randomly remove nodes and return values from the resulting networks
# set.seed is used to generate reproducible results

removal.sim<-function(network, simulations, remove){
  net<-network
  sim<-simulations
  rem<-remove
  nodes<-array(dim=c(sim,rem))
  number.nodes<-rep(1:length(V(net)$type)) 
  components<-rep(NA,length(sim))
  max.cluster.size<-rep(NA,length(sim))
  set.seed(1142)
  for (i in 1:sim){
    nodes[i,]<-sample(number.nodes, rem)
    tmp<-delete.vertices(net, nodes[i,])
    cl<-clusters(tmp)
    components[i]<-cl$no
    max.cluster.size[i]<-max(cl$csize)
  }
  tmp.out<-list("Components"=components,"Largest.Component"=max.cluster.size,"Removed.Nodes"=nodes)
  return(tmp.out)
}

# removal.sim(input network, number of simulations, number nodes to remove)

# a test
#removal.sim(ndt, 3, 3)

##############################################################################################################
##############################################################################################################

# to perform iterative repeated random removals across an increasing proportion of nodes

# batch mode simulations
# useful when wanting to do iterative calculations
# here, removing an increasing proportion of roosts

batch.rm.sim<-function(network, simulations, remove){
  net<-network
  sim<-simulations
  rem<-remove
  prp<-rep(NA,length(rem))
  mn<-rep(NA, length(rem))
  csse<-rep(NA, length(rem))
  noc<-vector('list', length(rem))
  for (i in 1:length(rem)){
    noc[[i]]<-rep(NA,sim)
    set.seed(1142)
    tmp<-removal.sim(net, sim, rem[i])
    noc[[i]]<-tmp$Components
    prp[i]<-rem[i]/length(V(net)$type)
    mn[i]<-mean(noc[[i]])
    csse[i]<-stderr(noc[[i]])
  }
  out<-list("Num.Clusters.Per.Sim"=noc,"Proportion.Removed"=prp,"Mean.No.Clusters"=mn,"SE.No.Clusters"=csse)
  return(out)
}

# batch.rm.sim(network, number simulations, vector of the number of nodes to remove)
# test
#y=rep(1:5)
#batch.rm.sim(ndt,10,y)

##############################################################################################################
##############################################################################################################

# add shapes for vertex plotting

shapes <- setdiff(vertex.shapes(), "")
g <- graph.ring(length(shapes))
set.seed(42)
#plot(g, vertex.shape=shapes, vertex.label=shapes, vertex.label.dist=1,
#     vertex.size=15, vertex.size2=15,
#     vertex.pie=lapply(shapes, function(x) if (x=="pie") 2:6 else 0),
#     vertex.pie.color=list(heat.colors(5)))

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                 plot=mytriangle)
#plot(g, vertex.shape="triangle", vertex.color=rainbow(vcount(g)),
#     vertex.size=seq(10,20,length=vcount(g)))

# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add.vertex.shape("star", clip=igraph.shape.noclip,
                 plot=mystar, parameters=list(vertex.norays=5))
#plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#     vertex.size=seq(10,20,length=vcount(g)))
#plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#     vertex.size=seq(10,20,length=vcount(g)),
#     vertex.norays=rep(4:8, length=vcount(g)))


##############################################################################################################
##############################################################################################################

# # Largest connected component for top and bottom nodes (used for mean distance for nodes below):
# gclust<-clusters(, mode='weak')
# lcc<-induced.subgraph(, V()[gclust$membership==1])
# lcctop<-length(V(lcc)[type==FALSE])
# lccbottom<-length(V(lcc)[type==TRUE])
# # Mean distance for top and bottom nodes
# distop<-mean(shortest.paths(lcc, v=V(lcc)[type==FALSE], to=V(lcc)[type==FALSE], mode = 'all'))
# disbottom<-mean(shortest.paths(lcc, v=V(lcc)[type==TRUE], to=V(lcc)[type==TRUE], mode = 'all'))
# #distop # bats
# #disbottom # roosts





