###########################
# F2 Nematodes
# NMDS Analysis
# Quentin Schorpp
# 05.08.2015
###########################



rankindex(scale(env1), fam, c("euc","man","bray","jac"))
vare.mds <- metaMDS(fam, trymax=100)
vare.mds1 <- metaMDS(fam,trymax=100, previous.best = vare.mds)

par(mfrow=c(1,1))
plot.varimax <- plot(vare.mds1, display = "sites")
points(plot.varimax, "sites", pch=25, bg=env1$age_class, col="black", cex=1.9)

