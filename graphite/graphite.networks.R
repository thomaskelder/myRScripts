#########################################################
## Methods to deal with networks from graphite package ##
#########################################################
pathwayToSpecies = function(pathway, homology, outProt) {
  library("org.Hs.eg.db")
  
  pathwayEntrez = convertIdentifiers(pathway, 'entrez')

  ## Find homologs for human ids in pathway
  entrezHuman = nodes(pathwayEntrez)
  human2hom = lapply(entrezHuman, mapGeneByHomology, homology = homology, outProt = outProt, inProt = org.Hs.egENSEMBLPROT)
  names(human2hom) = entrezHuman
  human2hom[is.null(human2hom)] = c()
  
  ## Convert edges to homologs
  edgesEntrez = edges(pathwayEntrez)
  edgesHom = edgesEntrez[F, , drop=F]
  for(i in seq_len(nrow(edgesEntrez))) {
    hrow = edgesEntrez[i,]
    h1 = as.character(hrow[1])
    h2 = as.character(hrow[2])
    for(x1 in human2hom[[h1]]) {
      for(x2 in human2hom[[h2]]) {
        xrow = hrow
        xrow[1] = x1
        xrow[2] = x2
        edgesHom = rbind(edgesHom, xrow)
      }
    }
  }
  
  pathwayHom = pathwayEntrez
  pathwayHom@nodes = as.character(unlist(human2hom))
  pathwayHom@edges = edgesHom
  pathwayHom
}

mapGeneByHomology = function(g, homology, outProt, inProt) {
  hp = unlist(mget(g, inProt, ifnotfound=NA))
  hp = hp[!is.na(hp)]
  if(length(hp) == 0) return(c())
  mp = unlist(mget(hp, homology, ifnotfound=NA))
  mp = mp[!is.na(mp)]
  if(length(mp) == 0) return(c())
  me = unlist(mget(mp, outProt, ifnotfound=NA))
  me = me[!is.na(me)]
  unique(me)
}

pathwaysToIGraph = function(pathways) {
  library(devtools)
  source_url("https://raw.github.com/thomaskelder/myRScripts/master/network-analysis/igraph.functions.R")
  
  graph = graph.empty(directed = T)
  for(pn in names(pathways)) {
    p = pathways[[pn]]
    d = graphite::edges(p)
    if(nrow(d) == 0) next
    
    ## Double edges for undirected
    d = rbind(d, d[d$direction == 'undirected', c(2, 1, 3:ncol(d))])
    
    pgraph = graph.data.frame(d)
    m = list(graph)
    m[[pn]] = pgraph
    graph = mergeGraphs(m, firstAsBase=T)
  }
  
  graph
}
