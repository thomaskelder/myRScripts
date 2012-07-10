
## Create colored pie chart based on data table using the Google image chart API
## The data argument should be a numerical matrix. A pie chart url will be generated for each row.
## Each column will become a slice in the pie chart and each slice will be colored according to the value in the matrix,
## based on the given colorRamp function.
##
## http://chart.apis.google.com/chart?chs=300x300&cht=p&chco=00FF00|0048FD|FF0000&chd=t:33.33,33.33,33.33
pieChartUrls = function(data, pieColors = colorRamp(c("blue", "white", "red")), lim = 2, width = 100, height = 100, rotation = 0) {
  chartUrlBase = "http://chart.apis.google.com/chart?cht=p&chf=bg,s,FFFFFF00"
  chartUrlBase = paste(chartUrlBase, '&chs=', width, 'x', height, sep='')
  chartUrlBase = paste(chartUrlBase, '&chp=', rotation, sep='')
  charts = apply(data, 1, function(x) {
    parts = paste(round(rep(100/length(x), length(x)), 3), collapse=',')
    url = paste(chartUrlBase, '&chd=t:', parts, sep='')
    xco = x
    xco[which(x > lim)] = lim
    xco[which(x < -lim)] = -lim
    cols = ""
    if(!is.na(sum(xco))) {
      xco = rgb(pieColors((xco + lim) / (2*lim)), max=255)
      xco = gsub('#', '', xco)
      chco = paste(xco, collapse='|')
      paste(url, '&chco=', chco, sep='')
    } else {
      ''
    }
  })
}