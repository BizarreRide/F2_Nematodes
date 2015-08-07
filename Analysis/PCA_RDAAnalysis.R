###########################
# F2 Nematodes
# PCA and RDA Analysis
# Quentin Schorpp
# 06.08.2015
###########################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
source("analysis/cleanplot.pca.R")
source("analysis/evplot.R")
source("analysis/RequiredPackages.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# multivariate normality ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam.log <- log1p(fam)
fam.z <- data.frame(scale(fam.log))


par(mfrow=c(5,5))
for (i in 1:25) {
  qqnorm(fam.z[,i], xlab=colnames(fam.z[i]))
  qqline(fam.z[,i])
  
  
}

par(mfrow=c(3,4))
for (i in 1:6) {
  hist(fam.z[,i], xlab=colnames(fam.z[i]))
}
summary(fam.z)

par(mfrow=c(1,1))

fam.mnorm <- t(fam.z)
mshapiro.test(fam.mnorm)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PCA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mngmnt.pca <- rda(mngmnt[,-1], scale=T)
summary(mngmnt.pca, display=NULL)
ev <- mngmnt.pca$CA$eig
evplot(ev)
cleanplot.pca(mngmnt.pca)
mngmnt$intensity <- scores(mngmnt.pca, display="sites")[,1] # fertilisation would be a uncorrelated addition to intensity


soil.pca <- rda(soil[,-c(1:3)], scale=T)
cleanplot.pca(soil.pca)
summary(soil.pca, display=NULL)
ev <- soil.pca$CA$eig
evplot(ev) # ommit: cn, silt, sand

climate.pca <- rda(climate1, scale=T)
cleanplot.pca(climate.pca) # ata2 + hum2

climate30.pca <- rda(climate30, scale=T)
cleanplot.pca(climate30.pca) # ata1 + prec1

env.rda <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation"))
env.rda$intensity <- mngmnt$intensity

env.fin.pca <- rda(env.rda[,-c(1:4)], scale=T)
cleanplot.pca(env.fin.pca)

env.rda2 <- subset(env.rda, select=c("field.ID", "age_class", "pH", "mc","ata2", "intensity"))




fam.repmes <- rda(fam.hel ~ age_class + pH + mc + ata2 + n + Condition(field.ID),data=env.rda)

stepping <- ordiR2step(rda(fam.hel ~ 1 + Condition(field.ID),data=env.rda), scope=formula(fam.repmes),direction="forward",pstep=1000,trace=F)

anova(stepping)

fam.repmes <- rda(fam.hel ~ age_class +  mc +  intensity + Condition(field.ID),data=env.rda2)

anova.cca(fam.repmes, step=100)
anova.cca(fam.repmes, step=100, by="axis")

# Partial RDA triplots (with fitted site scores)

# Scaling 1
par(mfrow=c(1,1))
plot(fam.repmes, scaling=1, display=c("sp", "lc","cn"))
plot(fam.repmes, display=c("sp", "lc","cn"))


cleanplot.pca(fam.rda)



# partial RDA ####

fam.hel <- decostand(fam,"hel")
str(fam.hel)
env.rda <- subset(env1, select=c("field.ID", "age_class","crop", "samcam","cn","mc","clay","ata1","pH","prec2"))


fam.rda <- rda(fam.hel)
cleanplot.pca(fam.rda)

fam.repmes <- rda(fam.hel, env1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# RDA after Ralf Schäfer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Here we are interested in families that have a strong relationship with the environmental gradient
# consequently we reove species that occur at 4 or less of the 15(30) study sites.
# Compare the discussion on the removal of rare species in the literature:
# Cao, Y.; Larsen, D. P.; Thorne, R. S., Rare species in multivariate analyysis for bioassessment:
# some considerations
# Marchant, R., Do rare species have any place in multivariate analysis for bioassessment?
# Rare species may especially cause problems in Correspondence Analysis and Canonical Correspondence Analysis
# In this Analysis, we remove them primarily to simplify the data set

# presence-absence transformation to calculate species number per site
fam.pa <- decostand(fam, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)

# remove species that occur at less than 5 sites
fam.fin <- fam[, !fam.sum<5]

sort(apply(fam.fin, 2, max))
# strong differences in order of magnitude of species abundances (?)
sort(apply(fam.fin,2,sd))

# now we look at the environmental variables
env.rda <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation"))
env.rda$intensity <- mngmnt$intensity
# collinearity may hamper interpretation  of the RDA with respect to the relevance of individual
# variables
# we therfore check for collinearity
# note that we could also conduct a PCA and work with orthogonal  variables,
# though this does not necessarily simplify the interpretation.

# function for an overview of environmental variables

panel.cor <- function(x,y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr =c(0,1,0,1))
  r <- abs(cor(x,y, use="complete.obs"))
  txt <- format(c(r,0,123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x,y)
  # borrowed from prinCoefmat
  Signif <- symnum(test$p.value, corr=FALSE, na=FALSE,
                   cutpoints=c(0,0.05,0.1,1),
                   symbols=c("*", "."," "))
  text(0.5, 0.5, txt, cex=cex *r)
  text(.8, .8, Signif, cex=cex, col=2)
}

# since there are too many variables for plotting in one go, we use subsets
pairs(env.rda[,5:10], lower.panel=panel.smooth, upper.panel=panel.cor)
pairs(env.rda[,11:15], lower.panel=panel.smooth, upper.panel=panel.cor)

#c and n should probably not be used together


# we can also check the variance inflation factor
library(faraway)
sort(vif(env.rda[,-c(1:4)]))
# Borcard et al 2011: 175 argue that in the case of RDA, VIFs > 10 should be avoided
# Theres a video of Ralf schäfer, that tells about how to deal with multicollinearity in multiple regression


# checking gradient length with detreded correspondence analysis
decorana(fam.fin)
# gives the average standard deviation
# a completely unimodal gradient has approximately 4 SD
# linear analyses are presumably ok for SD up to 2-3

# check wether Hellinger transformation decreases gradient length
# details on this transformation can be found in:
# Legendre. P; Gallagher, E.D., Ecologically meaningful transformation for ordination of
# species data

fam.hel <- decostand(fam.fin, "hellinger")
decorana(fam.hel)

# Hellinger transforation decreases the axis length, however, it is not necessary here

## RDA ####

# If we had used Hellinger TRansformation we would not need to scale
fam.rda <- rda(fam.fin[1:15]~., data=env.rda[,-c(1,3,4)], scale=TRUE)
summary(fam.rda)
# we obtain 15 RDA (constrained) and 14 PCA (unconstrained) axis,
# Further PCA axes are not displayed, as their eigenvalues are too small and would
# exceed the number of observations

# 67% of our vriance can be constrained and be explained by explanatory variables
# 33% can not be explained

# PCA is done for the unexplained variance only! Refers to the Residuals of the data

# The first two RDA axes explain 30% of the variance, this is less than 50% (indeed 45%) of the explained variability,
# So are these two axes sufficient?

# species scores represent the position of species in the axis in bi- or triplots
# site scores give the coordinates of sites in the spacce of the species
# site constraints give the coordinates of sites in the space of the environmental variables
# (explanatory)
# biplot scores give the coordinates for the environmental variable arrows


# Extraction of canonical coefficients from RDA object
coef(fam.rda)


# Global test of RDA result
set.seed(111)
anova.cca(fam.rda, step=1000)
# The Model is signficant





