#options(java.parameters = "-Xmx2048m")
library(XLConnect)

textFilesToExcel = function(files, out, freeze = T) {
  wb = loadWorkbook(out, create = TRUE)
  for(f in files) {
    message(f)
    d = read.delim(f, sep="\t", as.is=T)
    n = basename(f)
    n = gsub("\\.txt$", "", n)
    if(nchar(n) > 31) n = substr(n, 1, 31)
    createSheet(wb, name=n)
    writeWorksheet(wb, data=d, sheet=n, header = T)
    
    if(freeze) createFreezePane(wb, n, "B", 2)
  }
  saveWorkbook(wb)
}