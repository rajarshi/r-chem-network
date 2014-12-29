library(rcdk)
library(igraph)
library(fingerprint)
mols <- load.molecules('data/mipe-std.smi')
fps <- lapply(mols, get.fingerprint)
smat <- fp.sim.matrix(fps)
smat[ smat < 0.75 ] <- 0
g <- graph.adjacency(smat, mode='undirected', weight=TRUE, diag=FALSE)

par(mar=rep(0,4))
pdf('img/sim-network.pdf', 6,6)
plot(g, layout=layout.fruchterman.reingold,
     margin=0,
     vertex.label=NA, vertex.size=2)
dev.off()
