library(igraph)

id <- 6687479
g <- read.graph(sprintf('4-%d.gml',id), format='gml')

## filter out singletons
#singleton.idx <- which(degree(g) > 0)
#g <- induced.subgraph(g, vids=singleton.idx)

vcol <- rep('black', length(V(g)))
vcol[ which(V(g)$klass == 'doc') ] <- 'red'

pdf(sprintf("%d.pdf", id), 20,20)
par(mar=rep(0,4))
plot(g, layout=layout.fruchterman.reingold, vertex.size=2, vertex.label=NA, vertex.color=vcol)
dev.off()
