#####################################
## Functions for gene set analysis ##
#####################################
library(GSEABase)
library(GO.db)

defaultSetFilter = function(set) {
  bcCategory(collectionType(set)) %in% c("c2", "c3", "c5") &
  bcSubCategory(collectionType(set)) != "CGP"
}

noGOnoCGPSetFilter = function(set) {
     bcCategory(collectionType(set)) %in% c("c2", "c3", "c5") &
         !(bcSubCategory(collectionType(set)) %in% c("CGP", "CC", "MF", "BP"))
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

performEnrichment = function(tscores, sets, ids = rownames(tscores), signed = F, ...) {
  i = 0
  gsea = t(sapply(sets, function(s) {
    message(i, " out of ", length(sets))
    i <<- i + 1
    
    set = intersect(geneIds(s), rownames(tscores))
    p = sapply(colnames(tscores), function(g) {
      ts = tscores[,g]
      pval = geneSetTest(ids %in% set, ts, ...)
      if(signed) {
        m = mean(ts[set], na.rm=T)
        if(!is.na(m) && m != 0) pval = sign(m) * pval
      }
      pval
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
  
  gsea.adj = sapply(colnames(tscores), function(col) {
    padj = p.adjust(abs(gsea.num[,col]), method="BH")
    padj * sign(gsea.num[,col])
  })
  gsea.list$adj.pvalue = gsea.adj
  as.data.frame(gsea.list, stringsAsFactors=F)
}

buildEnrichmentTable = function(gsea.signed, gsea.mixed) {
  #gsea.signed = performEnrichment(tscores.bygene, sets, alternative = "either", signed = T)
  #gsea.mixed = performEnrichment(tscores.bygene, sets, alternative = "mixed", signed = F)

  gsea = gsea.signed
  colnames(gsea)[grep("pvalue", colnames(gsea))] = paste("abs signed", grep("pvalue", colnames(gsea), value = T))
  signed.fortable = gsea.signed[, grep("pvalue", colnames(gsea.signed))]
  colnames(signed.fortable) = paste("signed ", colnames(signed.fortable))
  gsea = cbind(gsea, signed.fortable, gsea.mixed[, grep("pvalue", colnames(gsea.mixed))])
  gsea[, grep("abs", colnames(gsea))] = abs(gsea[, grep("abs", colnames(gsea))])

  gsea
}

topEnrichmentGenes = function(gsea, sets, data, tCol, ids = rownames(data), pCutoff = 0.01, nrTop = nrow(data), pCol = "^adj.pvalue.") {
  print(rownames(gsea))
  pcol = grep(pCol, colnames(gsea))
  print(pcol)
  gsea = gsea[gsea[,pcol] < pCutoff,]
  top = gsea[order(gsea[, pcol]), ]
  
  if(nrow(top) == 0) {
    message("No significant sets!")
    return(c())
  }
  
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

enrichmentHeatmap = function(gsea, signed = T, rowNameTruncate = 25, pCutoff = 0.001, minSig = 1, minmax = 6, setNames = rownames(gsea), colNames = colnames(gsea), parseGroups = NULL, cluster_cols = T, colDist = function(x) as.dist(1 - cor(x)), colClustMeth = "complete", ...) {
  gsea.signed = gsea
  if(signed) {
    gsea = abs(gsea.signed)
  }
  gsea.log10 = -log10(gsea)
  for(cn in colnames(gsea.log10)) {
    gsea.log10[, cn] = sign(gsea.signed[, cn]) * gsea.log10[, cn]
  }
  library(pheatmap)
  shortSetNames = sapply(setNames, function(sn) {
    if(nchar(sn) > rowNameTruncate) sn = paste0(substr(sn, 1, rowNameTruncate), '...')
    sn
  })
  
  gsea.log10 = as.matrix(gsea.log10)
  rownames(gsea.log10) = shortSetNames
  
  hdata = gsea.log10[rowSums(abs(gsea.log10) > -log10(pCutoff)) >= minSig,]
  colnames(hdata) = colNames
  
  if(!is.null(parseGroups)) {
    groups = parseGroups(colnames(hdata))
    print(groups)
    matlist = lapply(unique(groups), function(g) {
      cols = colnames(hdata)[groups == g]
      if(cluster_cols) colorder = hclust(colDist(hdata[,cols]), method=colClustMeth)$order
      else colorder = 1:length(cols)
      cbind(hdata[, cols[colorder]], matrix(0, nrow = nrow(hdata), ncol = 1))
    })
    mat = matlist[[1]]
    for(m in 2:length(matlist)) mat = cbind(mat, matlist[[m]])
    mat = as.matrix(mat)
    colnames(mat)[grep("^matrix\\(", colnames(mat))] = ""
    mat = mat[, 1:(ncol(mat)-1)]
    hdata = mat
    cluster_cols = F
  }
  
  nc = 256
  colors = colorRampPalette(c('blue', 'white', 'red'))(nc+2)
  breaks = seq(-minmax, minmax, (2*minmax)/(nc-1))
  breaks = c(min(hdata, na.rm=T), breaks, max(hdata, na.rm=T))
  pheatmap(
    hdata, cluster_cols = cluster_cols, color = colors, scale = "none", breaks = breaks, ...
  )
}

getGOSets = function(goAnn, minSize = c(BP = 15, MF = 15, CC = 15), maxSize = c(BP = 500, MF = 500, CC = 500), ontologies = c("BP", "MF", "CC"), evidence = c("IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "TAS")) {
  annList = as.list(goAnn)
  print(length(annList))
  sel = sapply(names(annList), function(x) {
    xt = GOTERM[[x]]
    xo = Ontology(xt)
    xo %in% ontologies &&
      length(annList[[x]]) >= minSize[xo] &&
      length(annList[[x]]) <= maxSize[xo]
  })
  annList = annList[sel]
  print(length(annList))
  terms = names(annList)
  gosets = lapply(terms, function(t) {
    geneIds = annList[[t]][names(annList[[t]]) %in% evidence]
    geneIds = unique(as.character(geneIds))
    geneSet = GeneSet(geneIds, geneIdType=EntrezIdentifier(), 
                      setName=Term(GOTERM[[t]]))
  })
  names(gosets) = sapply(gosets, setName)
  gosets
}

getKeggDiseaseIds = function() {
	keggDisease = c(
		"Transcriptional misregulation in cancer",
		"Chemical carcinogenesis",
		"Viral carcinogenesis",
		"Colorectal cancer",
		"Pancreatic cancer",
		"Glioma",
		"Thyroid cancer",
		"Acute myeloid leukemia",
		"Chronic myeloid leukemia",
		"Basal cell carcinoma",
		"Melanoma",
		"Renal cell carcinoma",
		"Bladder cancer",
		"Prostate cancer",
		"Endometrial cancer",
		"Small cell lung cancer",
		"Non-small cell lung cancer",
		"Asthma",
		"Systemic lupus erythematosus",
		"Rheumatoid arthritis",
		"Autoimmune thyroid disease",
		"Allograft rejection",
		"Graft-versus-host disease",
		"Primary immunodeficiency",	
		"Alzheimer's disease",
		"Parkinson's disease",
		"Amyotrophic lateral sclerosis (ALS)",
		"Huntington's disease",
		"Prion diseases",	
		"Cocaine addiction",
		"Amphetamine addiction",
		"Morphine addiction",
		"Nicotine addiction",
		"Alcoholism",	
		"Hypertrophic cardiomyopathy (HCM)",
		"Arrhythmogenic right ventricular cardiomyopathy (ARVC)",
		"Dilated cardiomyopathy (DCM)",
		"Viral myocarditis",	
		"Type I diabetes mellitus",
		"Type II diabetes mellitus",
		"Maturity onset diabetes of the young",
		"Vibrio cholerae infection",
		"Vibrio cholerae pathogenic cycle",
		"Epithelial cell signaling in Helicobacter pylori infection",
		"Pathogenic Escherichia coli infection",
		"Salmonella infection",
		"Shigellosis",
		"Pertussis",
		"Legionellosis",
		"Staphylococcus aureus infection",
		"Tuberculosis",
		"Bacterial invasion of epithelial cells",
		"HTLV-I infection",
		"Measles",
		"Influenza A",
		"Hepatitis C",
		"Herpes simplex infection",
		"Epstein-Barr virus infection",
		"Amoebiasis",
		"Malaria",
		"Toxoplasmosis",
		"Leishmaniasis",
		"Chagas disease (American trypanosomiasis)",
		"African trypanosomiasis",
                "Leishmania infection"
	)
	
	keggDiseaseID = toupper(keggDisease)
	keggDiseaseID = gsub("[\\(\\)',-]{1}", "", keggDiseaseID)
	keggDiseaseID = paste0("KEGG_", keggDiseaseID)
	keggDiseaseID = gsub(" ", "_", keggDiseaseID)
	
	keggDiseaseID
}
