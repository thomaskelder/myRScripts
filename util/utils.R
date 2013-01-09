toExcel = function(path, xls) {
  require(devtools)
  options(java.parameters = "-Xmx2048m")
  source_url("https://raw.github.com/thomaskelder/myRScripts/master/util/excel.R")
  files = listFiles(path, ext = "txt")
  file.remove(xls)
  textFilesToExcel(files, xls)
}

performORA = function(set, all) {
	set = as.logical(set)
	all = as.logical(all)

	a = sum(set)
	b = sum(all) - a
	c = sum(!set)
	d = sum(!all) - c
	fisher.test(matrix(c(a, c, b, d),2,2), alternative='greater')$p.value
}

## Recursively lits files with a given extension ##
## from multiple paths.                          ##
listFiles = function(paths, ext) {
  dirs = file.info(paths)$isdir
  files = dir(paths[dirs], pattern=paste("\\.", ext, "$", sep=""), 
              recursive=T, full.names=T, ignore.case=T)
  files = append(files, paths[!dirs])
  files
}

## Map a list of entrez ids to another species by homology
## Use one of the hom.* bioconductor packages to provide homology information and the org.*
## packages to provide mappings from gene to ensembl protein
#inProt = org.Hs.egENSEMBLPROT, outProt = org.Mm.egENSEMBLPROT2EG, homology = hom.Hs.inpMUSMU
mapEntrezToSpecies = function(ids, inProt, outProt, homology) {
  message("mapping to proteins")
  ip = unlist(mget(ids, inProt, ifnotfound=NA))
  ip = ip[!is.na(ip)]
  if(length(ip) == 0) return(c())
  message("mapping to homolog proteins")
  op = unlist(mget(ip, homology, ifnotfound=NA))
  op = op[!is.na(op)]
  if(length(op) == 0) return(c())
  message("mapping to homolog genes")
  oe = unlist(mget(op, outProt, ifnotfound=NA))
  oe = oe[!is.na(oe)]
  unique(oe)
}


getLumiAnnot = function(nu.ids, pkg = "lumiMouseAll.db") {
  require(pkg, character.only = T)
  nu.annot = nuID2IlluminaID(nu.ids, lib.mapping = "lumiMouseIDMapping", idType = "All")
  
  cbind(
    probeId = nu.annot[,"Probe_Id"],
    entrez = as.character(mget(nu.ids, lumiMouseAllENTREZID)),
    symbol = as.character(mget(nu.ids, lumiMouseAllSYMBOL)),
    geneName = as.character(mget(nu.ids, lumiMouseAllGENENAME))
  )
}

cor.test.p = function(x, y = NULL, ...) {
  if(is.null(y)) y = x
  ps = matrix(1, ncol(x), ncol(y), dimnames=list(colnames(x), colnames(y)))
  for(i in 1:ncol(x)) {
    for(j in 1:ncol(y)) {
      ps[i,j] = cor.test(x[,i], y[,j], ...)$p.value
    }
  }
  ps
}
