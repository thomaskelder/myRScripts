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
  message("mapping to homolog proteins")
  op = unlist(mget(ip, homology, ifnotfound=NA))
  op = op[!is.na(op)]
  message("mapping to homolog genes")
  oe = unlist(mget(op, outProt, ifnotfound=NA))
  oe = oe[!is.na(oe)]
  unique(oe)
}