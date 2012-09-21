######################################
## Utility functions for T-profiler ##
######################################
require(pheatmap)
	
tprofilerHeatmap = function(
  gctFile, minSigSamples = 3, cutoffT = 4, minmax = 5, parseGroupsForSigSamples = NULL,
  colors = c('blue', 'white', 'red'), modifyColnames = function(x) x,
  cellwidth = 8, cellheight = 8, selectCols = function(x) x,
  clustering_distance_rows="correlation", clustering_distance_cols="correlation",
  ...
) {
  
  tprof = read.delim(gctFile, as.is=T, skip = 3)
  rownames(tprof) = tprof[,1]
  tprof = tprof[,3:ncol(tprof)]
  
  tprof = tprof[, selectCols(colnames(tprof))]

  if(!is.null(parseGroupsForSigSamples)) {
    groupsForSigSamples = parseGroupsForSigSamples(colnames(tprof))
    groupSig = sapply(unique(groupsForSigSamples), function(g) {
      rowSums(abs(tprof[, groupsForSigSamples == g]) > cutoffT) >= minSigSamples
    })
    tprof.filt = tprof[rowSums(groupSig) > 0,]
  } else {
    tprof.filt = tprof[rowSums(abs(tprof) > cutoffT) >= minSigSamples,]
  }
  
  hdata = tprof.filt
  
  colnames(hdata) = modifyColnames(colnames(hdata))
  nc = 256
  colors = colorRampPalette(colors)(nc+2)
  breaks = seq(-minmax, minmax, (2*minmax)/(nc-1))
  breaks = c(min(hdata, na.rm=T)-1, breaks, max(hdata, na.rm=T)+1)
  pheatmap(
    hdata, color = colors, scale = "none", breaks = breaks,
    cellwidth = cellwidth, cellheight = cellheight, 
    clustering_distance_rows=clustering_distance_rows, clustering_distance_cols=clustering_distance_cols,
    ...
  )
}
