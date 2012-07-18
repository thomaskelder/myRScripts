##########################################
## Functions dealing with igraph graphs ##
##                                      ##
## Thomas Kelder, 2012                  ##
##########################################
library(igraph)

###############################################################
## Export igraph graph to GML format compatible to Cytoscape ##
###############################################################
saveGML = function(g, fileName, title = "untitled") {
  attrToString = function(x) {
    m = mode(x)
    if(m == "numeric") {
      xc = sprintf("%.12f", x)
      
      xc[is.na(x)] = "NaN"
      xc[x == "Infinity"]= "Infinity"
      xc[x == "-Infinity"]= "-Infinity"
      x = xc
    } else {
      x = paste("\"", x , "\"", sep="")
    }
    x
  }
  
  vAttrNames = list.vertex.attributes(g)
  vAttrNames = vAttrNames[vAttrNames != "id"]
  vAttr = lapply(vAttrNames, function(x) attrToString(get.vertex.attribute(g, x)))
  names(vAttr) = vAttrNames
  eAttrNames = list.edge.attributes(g)
  eAttrNames = eAttrNames[eAttrNames != "id"]
  eAttr = lapply(eAttrNames, function(x) attrToString(get.edge.attribute(g, x)))
  names(eAttr) = eAttrNames
  
  f = file(fileName, "w")
  cat("graph\n[", file=f)
  cat(" directed ", as.integer(is.directed(g)), "\n", file=f)
  for(i in seq_len(vcount(g))) {
    cat(" node\n [\n", file=f)
    cat("    id", i, "\n", file=f)
    for(n in names(vAttr)) {
      cat("   ", gsub("[\\._]", "", n), vAttr[[n]][i], "\n", file=f)
    }
    cat(" ]\n", file=f)
  }
  
  el = get.edgelist(g, names=FALSE)
  for (i in seq_len(nrow(el))) { 
    cat(" edge\n  [\n", file=f) 
    cat("  source", el[i,1], "\n", file=f) 
    cat("  target", el[i,2], "\n", file=f) 
    for(n in names(eAttr)) {
      cat("   ", gsub("[\\._]", "", n), eAttr[[n]][i], "\n", file=f)
    }
    cat(" ]\n", file=f) 
  }
  
  cat("]\n", file=f)
  cat("Title \"", title, '"', file=f, sep="")
  close(f)
}

################################################################
## Combine undirected interaction networks into a             ##
## single graph. Attributes will be merged when possible.     ##
################################################################
mergeGraphs = function(graphs, setSourceAttr = T, firstAsBase = F) {
  interactions = graph.empty()
  if(firstAsBase) {
    interactions = graphs[[1]]
    graphs = graphs[2:length(graphs)]
    vas = list.vertex.attributes(interactions)
    if(!("name" %in% vas) && "identifier" %in% vas) V(interactions)$name = V(interactions)$identifier
  }
  
  for(graphName in names(graphs)) {
    message("processing ", graphName)
    g = graphs[[graphName]]
    vas = list.vertex.attributes(g)
    if(!("name" %in% vas) && "identifier" %in% vas) {
      print("Setting identifier attribute as name")
      V(g)$name = V(g)$identifier
    }
    
    #Add nodes
    newv = V(g)[!(name %in% V(interactions)$name)]
    newnames = V(g)[newv]$name
    attr = list()
    attr$name = newnames
    
    interactions = add.vertices(
      interactions, length(newv), attr = attr
      )
    
    ## Merge node attributes
    newAttr = setdiff(list.vertex.attributes(g), list.vertex.attributes(interactions))
    for(an in newAttr) {
      interactions = set.vertex.attribute(interactions, an, V(interactions), "")
      interactions = set.vertex.attribute(interactions, an, V(interactions)[V(g)$name], get.vertex.attribute(g, an))
    }
    mergeAttr = intersect(list.vertex.attributes(g), list.vertex.attributes(interactions))
    for(an in mergeAttr) {
      if(length(newnames) > 0)
        interactions = set.vertex.attribute(interactions, an, V(interactions)[newnames], get.vertex.attribute(g, an, V(g)[newnames]))
      
      ovl = setdiff(V(g)$name, newnames)
      if(length(ovl) > 0) {
        ovl.attr = cbind(
          get.vertex.attribute(interactions, an, V(interactions)[ovl]),
          get.vertex.attribute(g, an, V(g)[ovl])
          )
        ovl.attr = apply(ovl.attr, 1, function(x) {
          xs = gsub("([\\.\\(\\)]{1})", "\\\\\\1", x[2])
          if(is.na(x[1])) x[2]
          else if(
            x[1] != x[2] && 
              length(grep(paste(xs, ", ", sep=""), x[1])) == 0 &&
              length(grep(paste(xs, "$", sep=""), x[1])) == 0) {
            paste(x, collapse = ", ")
          } else {
            x[1]
          }
        })
        interactions = set.vertex.attribute(interactions, an, V(interactions)[ovl], ovl.attr)
      }
    }
    
    #Add edges
    newNodeIndex = as.numeric(V(interactions)[V(g)$name])
    edges = numeric(ecount(g)*2)
    el = get.edgelist(g, names = F)
    for(n in 1:nrow(el)) {
      if(n %% 1000 == 0) message(n, " out of ", ecount(g))
      edges[c(n*2-1,n*2)] = c(newNodeIndex[el[n,1]], newNodeIndex[el[n,2]])
    }
    attr = list()
    for(n in list.edge.attributes(g)) {
      attr[[n]] = as.character(get.edge.attribute(g, n))
    }
    if(setSourceAttr) attr$SourceFile = rep(graphName, ecount(g))
    attr$sourceEdge = as.numeric(E(g))
    interactions = add.edges(interactions, edges, attr = attr)
  }
  interactions
}

##################################################
## Transfer data matrix to node/edge attributes ##
##################################################
dataToGraph = function(graph, data, cols, matchCol, edges = T, nodes = T, edge.function = max, edge.na = 1, combine.function = function(rows, col) mean(rows[,col], na.rm=T)) {
  for(col in cols) {
    message("Processing ", col)
    
    if(nodes) {
      nodeWeights = sapply(V(graph)$name, function(x) {
        rows = data[data[,matchCol] == x,]
        if(!is.null(nrow(rows)) && nrow(rows) > 0) w = combine.function(rows, col)
        else w = rows[col]
        if(is.nan(w)) w = NA
        w
      })
      graph = set.vertex.attribute(graph, col, V(graph), as.numeric(nodeWeights))
    }
    
    if(edges) {
      nodeWeights = get.vertex.attribute(graph, col)
      
      edgeWeights = sapply(E(graph), function(e) {
        edge = get.edge(graph, e)
        v1 = nodeWeights[edge[1] + 1]
        v2 = nodeWeights[edge[2] + 1]
        if(is.na(v1)) v1 = edge.na
        if(is.na(v2)) v2 = edge.na
        edge.function(v1, v2)
      })
      graph = set.edge.attribute(graph, col, E(graph), edgeWeights)
    }
  }
  graph
}

#########################################
## Remove unconnected nodes from graph ##
#########################################
removeLonelyNodes = function(g) {
  induced.subgraph(g, which(igraph::degree(g) != 0) - 1)
}

#########################################
## Find communities per component      ##
#########################################
communitiesPerComponent = function(g, components, comFun = walktrap.community, ...) {
  lapply(components, function(comp) {
    cl = comFun(comp, ...)
    if(max(cl$membership) == (vcount(comp) -1)) { ##Invalid clustering (each node separate)
      cl$membership = rep(0, length(cl$membership))
    }
    cl
  })
}

assignMembership = function(g, clusterings, components, attr = "community", min.edges.per.node = 1, componentsAsClusterSize = vcount(components[[1]])) {
  if(length(components) == 0) return(g)
  
  ## Assign cluster membership to nodes
  last = 0
  for(i in 1:length(components)) {
    component = components[[i]]
    
    members = clusterings[[i]]$membership + last
    last = max(members) + 1
    names(members) = V(component)$name
    index = V(g)[name %in% names(members)]
    
    g = set.vertex.attribute(g, attr, index, members[V(g)[index]$name])
  }
  g = set.vertex.attribute(g, attr, V(g)[is.na(get.vertex.attribute(g, attr))], -1)
  
  ## Assign cluster membership to edges within cluster
  members.edge = apply(get.edgelist(g, names=F), 1, function(e) {
    m = get.vertex.attribute(g, attr, e)
    if(m[1] == m[2]) m[1]
    else -1
  })
  g = set.edge.attribute(g, attr, value = members.edge)
  
  g = filterCommunitiesByEdgeCount(g, attr, min.edges.per.node)
  
  # Set cluster membership to component number if component is too small to be clustered
  for(i in 1:length(components)) {
    component = components[[i]]
    if(vcount(component) < componentsAsClusterSize) {
      nodes = V(g)[name %in% V(component)$name]
      edges = E(g)[adj(nodes)]
      g = set.edge.attribute(g, attr, edges, paste("comp", i, sep=""))
      g = set.edge.attribute(g, paste("epn_", attr, sep=""), edges, paste("comp", i, sep=""))
      g = set.vertex.attribute(g, attr, nodes, paste("comp", i, sep=""))
    }
  }
  
  g
}

filterCommunitiesByEdgeCount = function(g, attr, min.edges.per.node) {
  members.edge = get.edge.attribute(g, attr)
  ## Add edge attribute for edge/node ratio filter
  members.edge.count = table(get.edge.attribute(g, attr))
  members.node.count = table(get.vertex.attribute(g, attr))
  edges.per.node = members.edge.count / members.node.count[names(members.edge.count)]
  inv = members.edge %in% names(which(edges.per.node <= min.edges.per.node))
  #g = set.edge.attribute(g, attr, which(inv) -1, -2)
  g = set.edge.attribute(g, paste("epn_", attr, sep=""), value = members.edge)
  g = set.edge.attribute(g, paste("epn_", attr, sep=""), which(inv), -2)
  g
}

removeUndirectedMultipleEdges = function(g, dirAttr = "Directed", undirValue = "false") {
  rmEdges = c()
  for(e in E(g)) {
    if(!(e %in% rmEdges)) {
      en = get.edge(g, e)
      er = E(g)[en[2] %->% en[1]]
      if(length(er) > 0 && get.edge.attribute(g, dirAttr, er) == undirValue) {
        rmEdges = c(rmEdges, er)
      }
    }
  }
  gn = delete.edges(g, rmEdges)
}

################################################
## Convenience function to perform community  ##
## finding and annotate communities with GO   ##
################################################
detectAndAnnotateCommunities = function(g, clustAttr = "cluster", organism = "mouse", componentsAsClusterSize = NULL, ...) {
  components = decompose.graph(g, "weak")
  components = components[order(sapply(components, vcount), decreasing=T)]
  components = components[sapply(components, vcount) > 2]
  communities = communitiesPerComponent(g, components, ...)
  if(is.null(componentsAsClusterSize)) componentsAsClusterSize = vcount(components[[1]])
  g = assignMembership(
    g, 
    communities, 
    components, attr = clustAttr,
    min.edges.per.node = 0,
    componentsAsClusterSize = componentsAsClusterSize
  )
  
  library(WGCNA)
  vClusters = get.vertex.attribute(g, clustAttr)
  clusters = unique(vClusters)
  go = GOenrichmentAnalysis(vClusters, V(g)$entrez, organism=organism)
  go.tab = go$bestPTerms[[4]]$enrichment
  
  g = set.vertex.attribute(g, paste0("GO_", clustAttr), value = rep("", vcount(g)))
  g = set.vertex.attribute(g, paste0("GOCount_", clustAttr), value = rep(-1, vcount(g)))
  g = set.vertex.attribute(g, paste0("GOSig_", clustAttr), value = rep("", vcount(g)))
  
  for(mod in clusters) {
    value = rep("", 3)
    modGO = which(go.tab[,'module'] == mod)
    if(length(modGO) == 0) {
      value = c("", "", -1)
    } else {
      sel = go.tab[modGO, ]
      sig = sel[sel[,"enrichmentP"] < 1E-3,]
      
      value = c(
        paste(sig[,"termName"], collapse = "; "), 
        paste(sig[,"nModGenesInTerm"], collapse = "; "), 
        as.numeric(nrow(sig) > 0)
      )
    }
    g = set.vertex.attribute(g, paste0("GO_", clustAttr), index = V(g)[vClusters == mod], value = value[1])
    g = set.vertex.attribute(g, paste0("GOCount_", clustAttr), index = V(g)[vClusters == mod], value = value[2])
    g = set.vertex.attribute(g, paste0("GOSig_", clustAttr), index = V(g)[vClusters == mod], value = value[3])
  }
  
  g
}
