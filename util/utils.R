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