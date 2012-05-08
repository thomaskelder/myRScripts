##########################################
## Functions dealing with igraph graphs ##
##                                      ##
## Thomas Kelder, 2012                  ##
##########################################

###############################################################
## Export igraph graph to GML format compatible to Cytoscape ##
###############################################################
saveGML = function(g, fileName, title = "untitled") {
	attrToString = function(x) {
		m = mode(x)
		if(m == "numeric") {
			x = sprintf("%.12f", x)
			
			x[is.na(x)] = "NaN"
			x[x == "Infinity"]= "Infinity"
			x[x == "-Infinity"]= "-Infinity"	
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
		cat("    id", i-1, "\n", file=f)
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
      edges[c(n*2-1,n*2)] = c(newNodeIndex[el[n,1]+1], newNodeIndex[el[n,2]+1])
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
dataToGraph = function(graph, data, cols, matchCol, edges = T, nodes = T, edge.function = max, edge.na = 1) {
  for(col in cols) {
    message("Processing ", col)
    
    if(nodes) {
      nodeWeights = sapply(V(graph)$name, function(x) {
        w = data[data[,matchCol] == x, col]
        w = mean(w, na.rm=T)
        if(is.nan(w)) w = NA
        w
      })
      graph = set.vertex.attribute(graph, col, V(graph), nodeWeights)
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
  subgraph(g, which(igraph::degree(g) != 0) - 1)
}

#########################################
## Find communities per component      ##
## using WalkTrap                      ##
#########################################
communitiesPerComponent = function(g, components) {
  lapply(components, function(comp) {
    cl = walktrap.community(comp)
    if(max(cl$membership) == (vcount(comp) -1)) { ##Invalid clustering (each node separate)
      cl$membership = rep(0, length(cl$membership))
    }
    cl
  })
}