###########################
# F2 Nematodes
# Fuzzy Clustering
# Quentin Schorpp
# 19.08.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
#env.fin$c <-  env1$c

data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

# Faunal Profile on subsample abundances
fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fuzzy Clusters ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(cluster)

fam.ch <- decostand(fam.org,"norm")
fam.bc <- vegdist(fam.org,"bray")

k <- 5
fam.fuz <- fanny(fam.bc, k=k, memb.exp=1.5, maxit=1000)
summary(fam.fuz)

fam.fuz.g <- fam.fuz$clustering

fam.fuz$membership

windows(record=TRUE)
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
