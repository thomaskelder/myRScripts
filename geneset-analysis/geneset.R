#####################################
## Functions for gene set analysis ##
#####################################
library(GSEABase)

defaultSetFilter = function(set) {
  bcCategory(collectionType(set)) %in% c("c2", "c3", "c5") &
  bcSubCategory(collectionType(set)) != "CGP"
}

loadBroadSets = function(broadXml = NULL, cache = NULL, filterFun = defaultSetFilter) {
  if(!is.null(cache)) {
    message("Using cache ", cache)
    load(cache)
    sets
  } else {
    setsAll = getBroadSets(broadXml)
    selectSets = sapply(setsAll, defaultSetFilter)
    setsAll[selectSets]
  }
}

mapSetsToSpecies = function(sets, inProt, outProt, homology, setsEntrez = sets) {
  i = 1
  setsMapped = lapply(names(sets), function(setName) {
    message(setName, "\t\t", i, " out of ", length(sets))
    i <<- i + 1
    
    s = sets[[setName]]
    se = setsEntrez[[setName]]
    hp = unlist(mget(geneIds(se), inProt, ifnotfound=NA))
    hp = hp[!is.na(hp)]
    mp = unlist(mget(hp, homology, ifnotfound=NA))
    mp = mp[!is.na(mp)]
    me = unlist(mget(mp, outProt, ifnotfound=NA))
    me = me[!is.na(me)]
    
    geneIds(s) = unique(me)
    s
  })
  names(setsMapped) = names(sets)
  setsMapped
}

saveSetCache = function(sets, file) {
  save(sets, file = file)
}

performEnrichment = function(tscores, sets, ids = rownames(tscores), ...) {
  i = 0
  gsea = t(sapply(sets, function(s) {
    message(i, " out of ", length(sets))
    i <<- i + 1
    
    set = intersect(geneIds(s), rownames(tscores))
    p = sapply(colnames(tscores), function(g) {
      ts = tscores[,g]
      geneSetTest(ids %in% set, ts, ...)
    })
    c(description = description(s), size = length(set), p)
  }))
  rownames(gsea) = sapply(sets, setName)
  
  gsea.list = list()
  gsea.list$description = gsea[, 'description']
  
  setLabels = rownames(gsea)
  setLabels = gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", setLabels, perl=TRUE)
  setLabels = gsub("_", " ", setLabels)
  gsea.list$name = setLabels

  gsea.num  = gsea[, 2:ncol(gsea)]
  mode(gsea.num) = 'numeric'
  
  gsea.list$size = gsea.num[, 'size']
  gsea.list$pvalue = gsea.num[, 2:ncol(gsea.num)]
  
  gsea.adj = sapply(colnames(tscores), function(col) p.adjust(gsea.num[,col], method="BH"))
  gsea.list$adj.pvalue = gsea.adj
  as.data.frame(gsea.list, stringsAsFactors=F)
}

topEnrichmentGenes = function(gsea, sets, data, tCol, ids = rownames(data), pCutoff = 0.01, nrTop = nrow(data)) {
  pcol = grep("^adj.pvalue.", colnames(gsea))
  gsea = gsea[gsea[,pcol] < pCutoff,]
  top = gsea[order(gsea[, pcol]), ]
  
  datSorted = data[order(abs(data[, tCol]), decreasing=T),]
  idsSorted = ids[order(abs(data[, tCol]), decreasing=T)]
  datSummary = rbind(colnames(data))
  for(sn in rownames(top)) {
    set = intersect(geneIds(sets[[sn]]), ids)
    setData = as.matrix(datSorted[idsSorted %in% set, , drop = F])
    setData = setData[1:min(nrow(setData), nrTop),]
    header = rep("", ncol(datSummary))
    header[1] = sn
    header[2] = length(set)
    datSummary = rbind(datSummary, header)
    datSummary = rbind(datSummary, setData)
  }
  colnames(datSummary) = datSummary[1,]
  rownames(datSummary) = NULL
  datSummary = datSummary[2:nrow(datSummary),]
  datSummary
}