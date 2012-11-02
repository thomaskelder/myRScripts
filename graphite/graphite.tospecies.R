################################################
## Convert graphite pathways to other species ##
################################################

library(foreach)
library(doMC)
registerDoMC(cores = 25)

source("graphite.networks.R")

## Assemble pathway list in other R process, because
## RSQLite has problems with next (parallel) step if
## already loaded in parent process.
if(!file.exists("graphitePathways.RData")) {
  message("Preparing graphite pathways list")
  library(graphite)
  appendPathways = function(l, p, name) {
    names(p) = paste0(names(p), " (", name, ")")
    append(l, p)
  }
  graphitePathways = biocarta
  names(graphitePathways) = paste(names(biocarta), "(biocarta)")
  graphitePathways = appendPathways(graphitePathways, kegg, "kegg")
  graphitePathways = appendPathways(graphitePathways, reactome, "reactome")
  graphitePathways = appendPathways(graphitePathways, nci, "nci")
  save(graphitePathways, file = "graphitePathways.RData")
} else {
  message("Loading graphite pathways list")
  load("graphitePathways.RData")
}

## Mouse
graphiteMouse = foreach(p=graphitePathways , .verbose = T) %dopar% {
  library(graphite)
  library(hom.Hs.inp.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  pathwayToSpecies(p, homology = hom.Hs.inpMUSMU, outProt = org.Mm.egENSEMBLPROT2EG)
}
names(graphiteMouse) = names(graphitePathways)
save(graphiteMouse, file = "graphitePathwaysMouse.RData");

graphiteMouseGraph = pathwaysToIGraph(graphiteMouse)
save(graphiteMouseGraph, file = "graphiteGraphMouse.RData")
saveGML(graphiteMouseGraph, fileName = "graphiteGraphMouse.gml", "Graphite-mouse")

## Rat
graphiteRat = foreach(p=graphitePathways , .verbose = T) %dopar% {
  library(graphite)
  library(hom.Hs.inp.db)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
  pathwayToSpecies(p, homology = hom.Hs.inpRATNO, outProt = org.Rn.egENSEMBLPROT2EG)
}
names(graphiteRat) = names(graphitePathways)
save(graphiteRat, file = "graphitePathwaysRat.RData");

graphiteRatGraph = pathwaysToIGraph(graphiteRat)
save(graphiteRatGraph, file = "graphiteGraphRat.RData")
saveGML(graphiteRatGraph, fileName = "graphiteGraphRat.gml", "Graphite-rat")