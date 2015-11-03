###########################
# F2 Nematodes
# PCA Analysis
# Quentin Schorpp
# 14.09.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam.pca <- rda(fam, scale=T)
summary(fam.pca, display=NULL)
ev <- fam.pca$CA$eig
evplot(ev)
screeplot(fam.pca, bstick=T, type="lines")
cleanplot.pca(fam.pca)

par(mfrow=c(1,1))
# Odinationsdiagramm Variablen als Vektoren, Scaling = 2 nach Variablen
pl <- biplot(fam.pca, type = "points")
# ausgew?hlte Variablen labeln
identify(pl,"sp", plot=TRUE, atpen=TRUE, labels=names(fam), cex=1.0)
# Verschiedene Symbole nach Gruppierungsvariable
points(pl, "sites", pch=25, bg=env1$age_class, cex=0.7)
legend("bottomright", legend=c("intensiv", "mittel", "extensiv"), text.col="black", col=c("black", "red", "green", pch=25, pt.cex=1.9))


# Variablen mit signifikantem Beitrag am Ordinationsergebnis, !Auf Samples skalieren!

plot1 <- ordiplot(fam.pca, scaling=1)
ordiequilibriumcircle(fam.pca,plot1)
identify(plot1,"sp", labels=names(fam), cex=1.0)
