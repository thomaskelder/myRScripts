######################################
## Utility functions for T-profiler ##
######################################
require(pheatmap)
	
tprofilerHeatmap = function(
  gctFile, minSigSamples = 3, cutoffT = 4, minmax = 5,
  colors = c('blue', 'white', 'red'), modifyColnames = function(x) x,
  cellwidth = 8, cellheight = 8, selectCols = function(x) x,
  clustering_distance_rows="correlation", clustering_distance_cols="correlation",
  parseGroups = NULL, parseGroupsForSigSamples = parseGroups, clusterPerGroup = F, cluster_cols = T, colDist = function(x) as.dist(1 - cor(x)), colClustMeth = "complete", ...
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

  if(clusterPerGroup & !is.null(parseGroups)) {
   groups = parseGroups(colnames(hdata))
  
	matlist = lapply(unique(groups), function(g) {
		cols = colnames(hdata)[groups == g]
		colorder = hclust(colDist(hdata[,cols]), method=colClustMeth)
		cbind(hdata[, cols[colorder$order]], matrix(0, nrow = nrow(hdata), ncol = 1))
	})
	mat = matlist[[1]]
	for(m in 2:length(matlist)) mat = cbind(mat, matlist[[m]])
	mat = as.matrix(mat)
   colnames(mat)[grep("^matrix\\(", colnames(mat))] = ""
   mat = mat[, 1:(ncol(mat)-1)]
	hdata = mat
	cluster_cols = F
  }
  
  colnames(hdata) = modifyColnames(colnames(hdata))
  nc = 256
  colors = colorRampPalette(colors)(nc+2)
  breaks = seq(-minmax, minmax, (2*minmax)/(nc-1))
  breaks = c(min(hdata, na.rm=T)-1, breaks, max(hdata, na.rm=T)+1)
  pheatmap(
    hdata, color = colors, scale = "none", breaks = breaks,
    cellwidth = cellwidth, cellheight = cellheight, 
    clustering_distance_rows=clustering_distance_rows, clustering_distance_cols=clustering_distance_cols,
    cluster_cols = cluster_cols, ...
  )
}
