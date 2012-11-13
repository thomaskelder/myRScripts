#################################################################
## Functions for co-expression analysis with the WGCNA package ##
#################################################################

plotTraitMod = function(datTraits, trait, MEs, outPath, group, mod, r = "", title = "", xlab = mod, ylab = trait) {
	library(ggplot2)
	x = datTraits[,trait]
	y = MEs[,mod]
	dta = data.frame(
		mouse = rownames(datTraits), trait = x, eg = y, group = group
		)
	
	p = ggplot(dta, aes(eg, trait))
	p = p + geom_point(aes(colour = group), size = 4)
	#p = p + geom_text(aes(label = mouse), hjust=0, vjust=0, size=3)
	p = p + opts(title = title)
	p = p + xlab(xlab) + ylab(ylab)
	
	pdf(paste(outPath, "cor.", r, ".", trait, ".", mod, ".pdf", sep=""))
	print(p)
	dev.off()
}

plotTraitCorrelations = function(moduleTraitCor, moduleTraitPvalue, MEs, datTraits, file = NULL, plotPvalue = T, windowSize = c(10,6), colors = greenWhiteRed(50)) {
  sizeGrWindow(windowSize[1], windowSize[2])
  textMatrix = signif(moduleTraitCor, 2)
  if(plotPvalue) {
  	textMatrix = paste(textMatrix, "\n(",
  	  signif(moduleTraitPvalue, 1), ")", sep = "")
  }
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
    xLabels = colnames(datTraits),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colors,
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1,1),
    main = paste("Module-trait relationships"))
  
  if(!is.null(file)) dev.copy2pdf(file = file)
}

plotSoftThreshold = function(sft, rs, file = NULL, powers = NULL) {
  ## Evaluate soft threshold
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red"
  )
  abline(h=rs,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red"
  )
  
  if(!is.null(file)) dev.copy2pdf(file = file)
}


plotHardThreshold = function(hrd, rs, cutoffs, file = NULL) {
  ## Evaluate soft threshold
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(hrd$fitIndices[,1], -sign(hrd$fitIndices[,4])*hrd$fitIndices[,5],
       xlab="Hard Threshold (TOM)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(hrd$fitIndices[,1], -sign(hrd$fitIndices[,4])*hrd$fitIndices[,5],
       labels=cutoffs,cex=cex1,col="red"
       )
  abline(h=rs,col="red")
  plot(hrd$fitIndices[,1], hrd$fitIndices[,6],
       xlab="Hard Threshold (TOM)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(hrd$fitIndices[,1], hrd$fitIndices[,6], labels=cutoffs, cex=cex1,col="red"
       )
  
  if(!is.null(file)) dev.copy2pdf(file = file)
}

plotNetAsDendro = function(net, file = NULL) {
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
      "Module colors",
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05)
  
  if(!is.null(file)) dev.copy2pdf(file = file)
}