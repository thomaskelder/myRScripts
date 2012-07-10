setwd("~/code/scripts/10OAD")

library(igraph0)

source("systemVariables.R")
source("kwalks.functions.R")

outPath = paste(oadPath, "output/kwalks/", sep="")

## Create test graph to reproduce figure 3 from Dupont et al.
graph = graph.formula(
  n1:n2:n3:n4:n5 -- n1:n2:n3:n4:n5, n3:n11:n6 -- n3:n11:n6, n6:n7:n8:n9:n10 -- n6:n7:n8:n9:n10  
)
graph = simplify(graph)
plot.igraph(graph, vertex.label=V(graph.kw)$name, layout=layout.fruchterman.reingold)

smx = kwalkMatrix(graph)
graph.kw = runKwalk(graph, smx, V(graph)["n1"] + 1, V(graph)["n9"] + 1)

plot.igraph(
  graph.kw, vertex.label=V(graph.kw)$name, layout=layout.fruchterman.reingold, edge.width=E(graph.kw)$kwalks * 20,
  vertex.size = V(graph.kw)$kwalks * 10
)

graph.kw.it = runKwalk(graph, smx, V(graph)["n1"] + 1, V(graph)["n9"] + 1, iterations = 5)

plot.igraph(
	graph.kw.it, vertex.label=V(graph.kw)$name, layout=layout.fruchterman.reingold, edge.width=E(graph.kw.it)$kwalks * 20,
	vertex.size = V(graph.kw.it)$kwalks * 10
)

plot(E(graph.kw)$kwalks, E(graph.kw.it)$kwalks)

graph.kw = runKwalk(graph, smx, V(graph)[c("n1", "n2")] + 1, V(graph)[c("n9", "n10")] + 1)

plot.igraph(
  graph.kw, vertex.label=V(graph.kw)$name, layout=layout.fruchterman.reingold, edge.width=E(graph.kw)$kwalks * 20,
  vertex.size = V(graph.kw)$kwalks * 10
)