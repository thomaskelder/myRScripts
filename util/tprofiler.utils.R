######################################
## Utility functions for T-profiler ##
######################################
require(pheatmap)
require(devtools)
source_url("https://raw.github.com/thomaskelder/myRScripts/master/util/heatmap.R")

loadTprofilerResults = function(file, skip = 3) {
  tprof = read.delim(file, as.is=T, skip = skip)
  rownames(tprof) = tprof[,1]
  tprof = tprof[,3:ncol(tprof)]
  tprof
}

tprofilerHeatmap = function(
  gctFile = NULL, gct = NULL, minSigSamples = 3, cutoffT = 4, minmax = 5,
  colors = c('blue', 'white', 'red'), modifyColnames = function(x) x,
  cellwidth = 8, cellheight = 8, selectCols = function(x) x,
  clustering_distance_rows="correlation", clustering_distance_cols="correlation", excludeSigSampleGroups = c(),
  parseGroups = NULL, parseGroupsForSigSamples = parseGroups, clusterPerGroup = F, cluster_cols = T, colDist = function(x) as.dist(1 - cor(x)), colClustMeth = "complete", ...
) {
  
  if(!is.null(gctFile)) {
    tprof = loadTprofilerResults(gctFile)
  } else {
    tprof = gct
  }
  
  tprof = tprof[, selectCols(colnames(tprof))]
  if(!is.null(parseGroupsForSigSamples)) {
    groupsForSigSamples = parseGroupsForSigSamples(colnames(tprof))
    groups = unique(groupsForSigSamples)
    groups = setdiff(groups, excludeSigSampleGroups)
    groupSig = sapply(groups, function(g) {
      rowSums(abs(tprof[, groupsForSigSamples == g]) > cutoffT) >= minSigSamples
    })
    tprof.filt = tprof[rowSums(groupSig) > 0,]
  } else {
    tprof.filt = tprof[rowSums(abs(tprof) > cutoffT) >= minSigSamples,]
  }
  
  groupedHeatmap(
    tprof.filt, minmax = minmax, colors = colors, modifyColnames = modifyColnames,
    cellwidth = cellwidth, cellheight = cellheight,
    clustering_distance_rows = clustering_distance_rows, clustering_distance_cols = clustering_distance_cols,
    parseGroups = parseGroups, clusterPerGroup = clusterPerGroup, cluster_cols = cluster_cols, 
    colDist = colDist, colClustMeth = colClustMeth, ...
  )
}
