###########################
# F2 Nematodes
# ISA
# Quentin Schorpp
# 19.08.2015
###########################


# Gower dissimilarity matrix

env.x <- env.fin[,c(1,2,6,7,9,13)]
summary(env.x)

env.x$field.ID <- as.factor(env.x$field.ID)

env.gow <- daisy(env.x, "gower")


env.kmeans <- kmeans(env.gow, centers = 5, nstart = 100)
env.kmeans$cluster
iva <- indval(fam,env.kmeans$cluster)

# Table of the significant indicator families
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(fam>0,2,sum)[iva$pval<=0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv); freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]


library(cluster)

fam.ch <- decostand(fam,"norm")

k <- 5
fam.fuz <- fanny(fam.bc, k=k, memb.exp=1.5, maxit=1000)
summary(fam.fuz)

fam.fuz.g <- fam.fuz$clustering

fam.fuz$membership

plot(silhouette(fam.fuz), main="Silhouette plot - Fuzzy clustering", cex.names=0.8, col=fam.fuz$silinfo$widths+1)

dc.pcoa <- cmdscale(fam.bc)
dc.scores <- scores(dc.pcoa, choices=c(1,2))

plot(scores(dc.pcoa), asp=1, type="n", main="Ordination of fuzzy CLusters")
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")

for (i in 1:k) {
  gg <- dc.scores[fam.fuz.g==i,]
  hpts <- chull(gg)
  hpts <- c(hpts, hpts[1])
  lines(gg[hpts,], col=i+1)
}

stars(fam.fuz$membership, location=scores(dc.pcoa),
      draw.segments=TRUE, add=TRUE, scale=FALSE, len=0.1, 
      col.segments=2:(k+1))
legend(locator(1), paste("Cluster", 1:k, sep=" "), pch=15, pt.cex=2, col=2:(k+1), bty="n")


library(labdsv)
indval(fam)
