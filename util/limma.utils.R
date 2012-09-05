##########################################
## Utility functions for limma analysis ##
##########################################
createStatsTable = function(fit, contrasts, annot = NULL, file = NULL, adjust.method = "BH") {
  deg = NULL
  for(g in contrasts) {
    d = fit[,g]
    d = cbind(
      d$coefficients,
      d$t,
      d$p.value,
      apply(d$p.value, 2, p.adjust, method = adjust.method)
    )
    colnames(d) = paste(c("logFC", "T", "p.value", "adj.p.value"), colnames(d), sep='_')
    if(is.null(deg)) deg = d
    else deg = cbind(deg, d)
  }
  deg = as.data.frame(deg)
  if(!is.null(annot)) {
    annot = as.data.frame(annot, stringsAsFactors=F)
    deg = cbind(annot, deg)
  }
  if(!is.null(file)) write.table(deg, file, sep="\t", quote = F, row.names = F)
  deg
}

createCleanDataTable = function(normData, fit, contrasts, groups, annot = NULL, file = NULL, adjust.method = "BH") {
  ## Normalized intensities
  cleanTable = normData
  
  ## Average per group
  for(g in unique(groups)) {
    avg = rowMeans(normData[, groups == g], na.rm = T)
    cleanTable = cbind(cleanTable, avg)
    colnames(cleanTable)[ncol(cleanTable)] = paste0("mean_", g)
  }

  ## Limma stats
  cleanTable = cbind(cleanTable, createStatsTable(fit, contrasts, annot = NULL, file = NULL, adjust.method = "BH"))
  
  ## Annotation
  if(!is.null(annot)) {
    annot = as.data.frame(annot, stringsAsFactors=F)
    cleanTable = cbind(annot, cleanTable)
  }
  
  if(!is.null(file)) write.table(cleanTable, file, sep="\t", quote = F, row.names = F)
  cleanTable
}

sigCounts = function(fit, ...) {
  sigTest = decideTests(fit, ...)
  numSig = summary(sigTest)
  numSig = rbind(numSig, numSig[1,] + numSig[3,])
  rownames(numSig) = c("Down", "Not significant", "Up", "Significant")
  numSig = numSig[c("Up", "Down", "Significant", "Not significant"),]
  numSig
}

sigCountsRange = function(fit, p.range = c(0.05, 0.01, 0.005, 0.001), adj.range = c("none", "BH"), ...) {
  sigCounts = NULL
  for(adj in adj.range) {
    for(p in p.range) {
      counts = sigCounts(fit, adjust.method = adj, p.value = p, ...)["Significant", , drop=F]
      adjLabel = ""
      if(adj != "none") adjLabel = paste0(adj, " adj. ")
      rownames(counts) = paste0(adjLabel, "p < ", p)
      if(is.null(sigCounts)) sigCounts = counts
      else sigCounts = rbind(sigCounts, counts)
    }
  }
  sigCounts
}