library(igraph)
library(RColorBrewer)
library(Hmisc)

id <- 12477366
##id <- 21568322
g <- read.graph(sprintf('4-%d.xml',id), format='graphml')

## degree of doc nodes - how many fragments did we identify?
doc.idx <- which(V(g)$klass == 'doc')
degree(g)[doc.idx]
frag.idx <- which(V(g)$klass == 'fragment')
degree(g)[frag.idx]


## Year based doc node colors
year <- cut(V(g)$year, 5)
year <- data.frame(year = levels(year))
year.cols <- rev(brewer.pal(5, 'Blues'))[1:nrow(year)]
year$col <- year.cols

tmp <- data.frame(id=1:length(V(g)), year = cut(V(g)$year,5), klass=V(g)$klass)
tmp <- merge(tmp, year, by='year', all.x=TRUE)
tmp$col[ which(tmp$klass == 'fragment') ] <- 'black'
tmp$col[ which(tmp$klass == 'doc' & is.na(tmp$col)) ] <- 'grey'
tmp <- tmp[order(tmp$id),]
vcol <- tmp$col

## MeSH heading based doc node colors
mesh <- data.frame(table(V(g)$mesh))
mesh <- subset(mesh, Freq > 5 & Var1 != '' &  Var1 != 'None')
mesh <-mesh[order(-mesh$Freq),]
mesh.cols <- brewer.pal(nrow(mesh), 'Paired')
mesh$col <- mesh.cols

tmp <- data.frame(id=1:length(V(g)), mesh = V(g)$mesh, klass=V(g)$klass)
tmp <- merge(tmp, mesh, by.x='mesh', by.y='Var1', all.x=TRUE)
tmp$col[ which(tmp$klass == 'fragment') ] <- 'black'
tmp$col[ which(tmp$klass == 'doc' & is.na(tmp$col)) ] <- 'grey'
tmp <- tmp[order(tmp$id),]
vcol <- tmp$col

## Target classification headings
tc <- data.frame(table(V(g)$l2))
tc <- tc[order(tc$Freq),]
tc <- subset(tc, Freq > 3)
tc.cols <- brewer.pal(nrow(tc), 'Paired')
tc$col <- tc.cols

tmp <- data.frame(id=1:vcount(g), tc = V(g)$l2, klass=V(g)$klass)
tmp <- merge(tmp, tc, by.x='tc', by.y='Var1', all.x=TRUE)
tmp$col[ which(tmp$klass == 'fragment') ] <- 'black'
tmp$col[ which(tmp$klass == 'doc' & is.na(tmp$col)) ] <- 'grey'
tmp <- tmp[order(tmp$id),]
vcol <- tmp$col

vframecol <- rep('black', length(V(g)))
vframecol[1] <- 'red'

vsize <- rep(2.5, length(V(g)))
vsize[which(V(g)$klass == 'fragment')] <- 1
vsize[1] <- 8

vshape <- rep('circle', length(V(g)))
vshape[ which(V(g)$klass == 'fragment') ] <- 'square'

#vcol <- rep('black', length(V(g)))
#vcol[ which(V(g)$klass == 'doc') ] <- 'red'

pdf(sprintf("4-%d.pdf", id), 10,10)

layout(matrix(c(1,2,1,2), 2,2, byrow=TRUE), widths=c(3,1), heights=c(1,1))
par(mar=rep(0,4))
l = layout.fruchterman.reingold(g, niter=1000, area=80*vcount(g)^2)
plot(g, layout=l, margin=0,
     vertex.size=vsize, vertex.label=NA, vertex.color=vcol, vertex.shape = vshape,
     vertex.frame.color = vframecol,
     edge.color='light grey')
## Legend
par(mar=c(20,10,15,2))
lut <-  tc
plot(rep(0, nrow(lut)), 1:nrow(lut), xlim=c(0,1),
     ylim=c(1,nrow(lut)), xaxt='n', yaxt='n',  type='n', axes=FALSE, bty='n', xlab='', ylab='')
start.y <- 0
for (i in 1:nrow(lut)) {
  rect(0,start.y, 1,start.y+1, col=lut$col[i], border='light grey')
  start.y <- start.y+1
}
axis(side=2, at=seq(0.5, nrow(lut)-0.5, by=1), labels=capitalize(as.character(lut$Var1)), tick=FALSE, las=2)
dev.off()


