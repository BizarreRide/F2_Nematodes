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

# Select abundant species (optional)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fam.pa <- decostand(fam.usc, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)

# remove species that occur at less than 5 sites
fam.fin <- fam.usc[, !fam.sum<5]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Compute NMDS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmds.usc1  <- metaMDS(fam.usc, distance="bray", k=2, trymax=100)
nmds.usc <- metaMDS(fam.usc, distance="bray", k=2,trymax=100, previous.best = nmds.usc1)
nmds.usc # Stress = 0.24
# Stress: 0.24
#**A good rule of thumb: 
# stress > 0.05 provides an excellent representation in reduced dimensions, 
# > 0.1 is great, 
# > 0.2 is good/ok, and stress 
# > 0.3 provides a poor representation.**
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Validation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stressplot(nmds.usc)
 # Large scatter around the line suggests that original dissimilarities are not well preserved 
 # in the reduced number of dimensions. Looks pretty good in this case.
gof  <- goodness(nmds.usc)
plot(nmds.usc, type="t", main="goodness of fit")
points(nmds.usc, display="sites", cex=gof*100)
 # Poorly fitted sites have larger bubbles
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# check dissimilarity measure (useful?)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# extract only numeric environmental variables
env1.num <- env.fin[sapply(env.fin,is.numeric)]
env1.num <- env1.num[,-c(1,2,5,8,9)] 
env1.num$age_class <- env1$age_class
# Rank correlations between dissimilarity indices and gradient separation
rankindex(environmental, veg=fam.usc, indices=c("gow","euc","man","bray","jac"), method="spearman")

# Use procrustes() to compare different dissimilarity measures
pro <- procrustes(vare.mds1, vare.mds2)
plot(pro, cex=1.5) # the plot shows which sites move far
plot(pro, kind=2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# NMDS Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(nmds.usc)
points(nmds.usc, display="sites", pch=25, bg=env1$age_class, col="black", cex=1.2)

ordiplot(nmds.usc,type="n")
ordihull(nmds.usc,groups=env1$age_class,draw="polygon",col="grey90",label=F)
orditorp(nmds.usc,display="species",col="red",air=0.01)
orditorp(nmds.usc,display="sites",cex=1.25,air=0.01)

#ordisurf(nmds.usc,env1$pH,main="",col="forestgreen")
# function ordisurf fits surfaces of environmental variables to ordinations based on generalized additive models in function gam of package mgcv

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Results:
# NMDS with upscaled family abundance data, reveals very similar communities in the Maize plots, 
# caused by the occurence of Pratylenchidae
# SIlphie fields seem to be generally more similar


# Add colours from a clustering result to an NMDS plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Ward clustering of Bray curtis dissimilarity matrix
# and extraction of five (two) groups

# bray curtis dissimilarity matrix of upscaled family data (fam.usc.bray = fubc)

fubc <- vegdist(fam.usc, "bray")

fubc.ward <-  hclust(fubc, "ward.D")
fubcw.groups <- cutree(fubc.ward, k=3)
grp.lev <- levels(factor(fubcw.groups))

# Combination with NMDS result
sit.sc <- scores(nmds.usc)
p <- ordiplot(sit.sc, type="n", main="NMDS/Bray + clusters Ward/Bray")

for(i in 1:length(grp.lev)) {
  points(sit.sc[fubcw.groups==i,], pch=(14+i), cex=2, col=i+1)
}

text(sit.sc, row.names(fam.usc), pos=4, cex=0.7)

# Add the dendrogram
ordicluster(p, fubc.ward, col="dark grey")
legend(locator(1), paste("Group", c(1:length(grp.lev))), 
       pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)),pt.cex=2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Vector overlay and significance test
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fit <- envfit(vare.mds1, env.fin, strata="samcam", permu=999)
fit

plot(nmds.usc, display="sites")
points(nmds.usc, display="sites", pch=25, bg=env1$age_class, col="black", cex=1.9)
plot(fit,p.max=0.05)

# More tests
env.fin$field.ID <- env1$field.ID
dist.matrix <- dist(fam, method = "manhattan")
adonis(fam ~ age_class, env.fin, permutations = 999,
       method = "bray", strata = NULL, contr.unordered = "contr.sum")

# Repeated measurements PERMANOVA

# Friedman test

# Scheirer-Ray-Hare Test
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

