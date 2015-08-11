###########################
# F2 Nematodes
# NMDS Analysis
# Quentin Schorpp
# 05.08.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam.pa <- decostand(fam, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)

# remove species that occur at less than 5 sites
fam.fin <- fam[, !fam.sum<5]

str(env1)
# extract only numeric environmental variables
env1.num <- env.fin[sapply(env.fin,is.numeric)]
env1.num <- env1.num[,-c(1,2,5,8,9)] 
env1.num$age_class <- env1$age_class

# Rank correlations between dissimilarity indices and gradient separation
rankindex(env1.num, veg=fam.fin, indices=c("gow","euc","man","bray","jac"), method="spearman")
vare.mds <- metaMDS(fam.fin, trymax=100)
vare.mds1 <- metaMDS(fam.fin,trymax=100, previous.best = vare.mds)
vare.mds1
# Stress: 0.23

par(mfrow=c(1,1))
plot.varimax <- plot(vare.mds1, display = "sites")
plot.varimax <- plot(vare.mds1, type="t")
points(plot.varimax, "sites", pch=25, bg=env1$age_class, col="black", cex=1.9)

stressplot(vare.mds1)

gof  <- goodness(vare.mds1)
plot(vare.mds1, type="t", main="goodness of fit")
points(vare.mds, display="sites", cex=gof*100)


vare.mds <- metaMDS(fam, distance="eu", trymax=100)
vare.mds2 <- metaMDS(fam,trymax=100, distance="eu",previous.best = vare.mds)

pro <- procrustes(vare.mds1, vare.mds2)

plot(pro, cex=1.5)
plot(pro, kind=2)

fit <- envfit(vare.mds1, env.fin, strata="samcam", permu=999)
fit

# function ordisurf fits surfaces of environmental variables to ordinations based on generalized additive models in function gam of package mgcv


plot(vare.mds1, display="sites")
points(vare.mds1, display="sites", pch=25, bg=env1$age_class, col="black", cex=1.9)
plot(fit,p.max=0.05)

env.fin$field.ID <- env1$field.ID
dist.matrix <- dist(fam, method = "manhattan")
adonis(fam ~ age_class, env.fin, permutations = 999,
       method = "bray", strata = NULL, contr.unordered = "contr.sum")

# Repeated measurements PERMANOVA

# Friedman test

# Scheirer-Ray-Hare Test

