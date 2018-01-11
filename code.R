library(rcdk)
library(igraph)
library(fingerprint)

library(org.Hs.eg.db)
library(DBI)
library(RMySQL)

m<-dbDriver("MySQL");
con<-dbConnect(m,user='XXX',password='XXX',host='XXX',dbname='XXX');

mipe <- read.csv('data/mipe.csv', header=TRUE, as.is=TRUE)
mipe <- subset(mipe, GENE_SYMBOL != '')

## Convert syms to uniprot ID's and then pull in Panther class
## but only keep the first Panther class at node_level = 1
eids <- select(org.Hs.eg.db, keys=mipe$GENE, columns=c("SYMBOL", "ENTREZID", "UNIPROT"), keytype="SYMBOL")
sql <- "select a.pclass_id, node_level, class_name from panther_uniprot_map a, panther_class b where accession in (%s) and a.pclass_id = b.pclass_id"
tmp <- by(eids, eids$SYMBOL, function(x) {
  accs <- join(sprintf("'%s'", x$UNIPROT), ',')
  q <- sprintf(sql, accs)
  r <- dbGetQuery(con, q)
  if (nrow(r) == 0) ret <- data.frame(gene=x$SYMBOL[1], pclass_id=NA, node_level=NA, class_name="UNKNOWN")
  else {
    ret <- cbind(gene=x$SYMBOL[1], r)
    ret <- subset(ret, node_level == 1)[1,]
  }
  return(ret)
})
tmp <- do.call(rbind, tmp)
tmp.mipe <- merge(mipe, tmp, by.x='GENE_SYMBOL', by.y='gene')
mipe <- tmp.mipe

mols <- parse.smiles(mipe$SMILES_STD)
fps <- lapply(mols, get.fingerprint, type='extended')
fpsim <- fp.sim.matrix(fps)
colnames(fpsim) <- mipe$SAMPLE_NAME
rownames(fpsim) <- mipe$SAMPLE_NAME
fpsim[ fpsim < 0.6 ] <- 0
summary(as.numeric(fpsim))

## Similarity based network, colored by target
g <- graph.adjacency(fpsim, mode='undirected', weight=TRUE, diag=FALSE)
g <- set.vertex.attribute(g, 'gene', V(g), mipe$GENE)
g <- set.vertex.attribute(g, 'panther', V(g), mipe$class_name)

col21 <- list("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
names(col21) <- unique(V(g)$panther)
col21$UNKNOWN <- 'grey'
vcols <- unlist(sapply(V(g)$panther, function(x) col21[x]))

## Find v's that are not singletons
## genes.in.connected <- V(g)$gene[ which( degree(g) > 0) ]
## tmp <- subset(data.frame(table(genes.in.connected)), Freq >= 10)

pdf('img/sim-network.pdf', 8,7)
layout(matrix(c(1,2,1,2), 2,2, byrow=TRUE), widths=c(3,1), heights=c(1,1))
par(mar=c(0,0,0,0))
plot(g, layout=layout.fruchterman.reingold,
     margin=0,
     vertex.label=NA, vertex.size=2, vertex.color=vcols)
par(mar=c(0,0,0,0.5))
plot.new()
legend('center', legend=names(col21),  fill=unlist(col21),
       border='white', pt.cex=4, bty='n', y.intersp=1.5,
       title='PANTHER Class')
dev.off()

non.singles <- which(degree(g) > 0)
g2 <- induced.subgraph(g, vids = non.singles)
vcols <- unlist(sapply(V(g2)$panther, function(x) col21[x]))
pdf('img/sim-network-connected.pdf', 9,7)
layout(matrix(c(1,2,1,2), 2,2, byrow=TRUE), widths=c(3,1), heights=c(1,1))
par(mar=c(0,0,0,0))
l <- layout.fruchterman.reingold(g2, niter=1000, area=80*vcount(g2)^2)
plot(g2, layout=l,
     margin=0,
     vertex.label=NA, vertex.size=2.5, vertex.color=vcols)
par(mar=c(0,0,0,0.5))
plot.new()
legend('center', legend=names(col21),  fill=unlist(col21),
       border='white', cex=1.2, pt.cex=8, bty='n', y.intersp=1.5,
       title=expression(bold('PANTHER Class')))
dev.off()
