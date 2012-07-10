pruneBranches = function(g, ignore = c()) {
	## Iteratively remove all nodes with only one connection
	## until there are none left
	checkDegree = function(g) {
		d = igraph::degree(g)
		dns = V(g)[d == 1]$name
		setdiff(dns, ignore)
	}
	rm = checkDegree(g)
	while(length(rm) > 0) {
		g = subgraph(g, V(g)[!(name %in% rm)])
		rm = checkDegree(g)
	}
	g
}

doCutoff = function(g, er, cutoff, removeNodes = T) {
	gc = delete.edges(g, E(g)[er < cutoff])
	if(removeNodes) removeLonelyNodes(gc)
	else gc
}

findEdgeRelevanceCutoff = function(g.kw, toInclude, optional, minOptional = 1, attr = "kwalks", incr = 0.001) {
  # Determine minimum edge relevance so that
  # both the drugs and module nodes are in same
  # component
  er = get.edge.attribute(g.kw, attr)
  
  g.kw.cut = NULL
  if(length(incr) == 1) incr = seq(max(er, na.rm=T), 0, -incr)
  for(cutoff in incr) {
    g.kw.cut = doCutoff(g.kw, er, cutoff, F)
    cl = clusters(g.kw.cut, mode = "weak")$membership
    cl.incl = cl[as.numeric(V(g.kw.cut)[toInclude]) + 1]
    cl.opt = cl[as.numeric(V(g.kw.cut)[optional]) + 1]
    if(length(unique(cl.incl)) == 1) {
    	if(sum(cl.opt == unique(cl.incl)) >= minOptional) break
    }
  }
  g.kw.cut = removeLonelyNodes(g.kw.cut)
  list(graph = g.kw.cut, cutoff = cutoff)
}

## Convert graph to a format readable for kwalks
kwalkMatrix = function(g, file, ...) {
  adj = get.adjacency(g, ...)
  smx = sapply(1:nrow(adj), function(n) {
    nbs = which(adj[n,] > 0)
    if(length(nbs > 0)) {
      e = paste(nbs, sprintf("%.32f", adj[n,nbs]), sep=":")
      paste(n, paste(e, collapse = " "))
    } else {
      n
    }
  })
  smx = c(vcount(g), smx)
}

runKwalk = function(g, smx, a, b, lmax = 50, tmpPath = "/tmp/", scriptPath = getwd(), attr = "kwalks", iterations = 1) {
  # Write the graph file
  fin = paste(tmpPath, "kwalks.graph.txt", sep="")
  write.table(smx, fin, quote=F, row.names = F, col.names=F)
  
  fout = paste(tmpPath, "kwalks.out", sep="")

  # Run kwalks
  for(i in seq_len(iterations)) {
	  rel = paste(paste(a, collapse=":"), paste(b, collapse=":"), sep="#")
	  cmd = paste(scriptPath, "/lkwalk", sep="")
	  cmd = paste(cmd, "-g", fin, "-l", lmax, "-k", rel, "-o", fout)
	  message(cmd)
	  system(cmd)
	  
	  fin = paste(fout, ".dif", sep="")
  }
  message("Kwalks finished!")
  # Parse results
  g.w = g
  
  fout.dif = paste(fout, ".dif", sep="")
  dif = unlist(read.delim(fout.dif, sep="\t", stringsAsFactor = F, header=F))
  dif = dif[2:length(dif)]
  
  edges = c()
  values = c()
  
  weights = matrix(0, nrow = vcount(g.w), ncol = vcount(g.w))
  
  for(i in 1:length(dif)) {
    if(i %% 100 == 0) {
      message(i, " out of ", length(dif))
    }
    x = strsplit(dif[i], " ")[[1]]
    n1 = as.numeric(x[1])
    for(y in x[2:length(x)]) {
      y = as.numeric(strsplit(y, ":")[[1]])
      n2 = y[1]
      w = y[2]
      weights[n1, n2] = weights[n2, n1] = w
    }
  }
  values = apply(get.edgelist(g.w, names = F), 1, function(x) {
    weights[x[1] + 1, x[2] + 1]
  })
  g.w = set.edge.attribute(g.w, attr, E(g.w), values)
  nw = as.character(unlist(read.delim(paste(fout, ".N", sep=""), stringsAsFactor = F, header = F)))
  nw = sapply(nw, function(x) as.numeric(strsplit(x, ":")[[1]]))
  g.w = set.vertex.attribute(g.w, attr, V(g.w)[nw[1,] - 1], nw[2,])
  g.w
}
