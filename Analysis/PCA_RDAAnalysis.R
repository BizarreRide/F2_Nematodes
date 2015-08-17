###########################
# F2 Nematodes
# PCA and RDA Analysis
# Quentin Schorpp
# 06.08.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test for multivariate normality ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Log x+1 transformation of species data
fam.log <- log1p(fam)
# centering species data
fam.z <- data.frame(scale(fam.log))

# Multivariate Shapiro Test of transposed species data
fam.mnorm <- t(fam.z)
mvnormtest::mshapiro.test(fam.mnorm)


# QQ Plots of all variables
par(mfrow=c(5,5))
for (i in 1:25) {
  qqnorm(fam.z[,i], xlab=colnames(fam.z[i]))
  qqline(fam.z[,i]) 
}

# Histrograms of all variables
par(mfrow=c(3,4))
for (i in 1:6) {
  hist(fam.z[,i], xlab=colnames(fam.z[i]))
}
summary(fam.z)
par(mfrow=c(1,1))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# RDA after Ralf Schäfer ####
# https://www.youtube.com/watch?v=gY_iktfpSpQ

## 1. Data Preparation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
env.fin <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation"))
env.fin$intensity <- mngmnt$intensity
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
pairs(env.fin[,5:10], lower.panel=panel.smooth, upper.panel=panel.cor)
pairs(env.fin[,11:15], lower.panel=panel.smooth, upper.panel=panel.cor)

#c and n should probably not be used together


# we can also check the variance inflation factor
library(faraway)
sort(vif(env.fin[,-c(1:4)]))
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

# test for Heteroscedasticity within groups
fam.hel.d1 <- dist(fam.hel)
fam.MHV <- betadisper(fam.hel.d1, env.fin$age_class)
permutest(fam.MHV)

## 2. RDA ####
### for each sampling campaign ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# SAMPLING CAMPAIGN 2 - AUTMUN 2012 ---------------------

# If we had used Hellinger TRansformation we would not need to scale
fam.rda1 <- rda(fam.fin[env.fin$samcam==2,] ~ ., data=env.fin[env.fin$samcam==2,-c(1,3,4,7,10,11)], scale=TRUE)
summary(fam.rda1, display=NULL)

# Reduce environmental variables to significant ones 
stepping <- ordiR2step(rda(fam.fin[env.fin$samcam==2,] ~ 1 ,data=env.fin[env.fin$samcam==2,-c(1,3,4,7,10,11)]), scope=formula(fam.rda1),direction="both",pstep=1000,trace=F)

anova(stepping)

fam.rda1 <- rda(fam.fin[env.fin$samcam==2,]~ age_class + n + pH, data=env.fin[env.fin$samcam==2,], scale=TRUE)
summary(fam.rda1, display=NULL)

RsquareAdj(fam.rda1)


# we obtain 6 RDA (constrained) and 8 PCA (unconstrained) axis,
# Further PCA axes are not displayed, as their eigenvalues are too small and would
# exceed the number of observations

# 54% of our vriance can be constrained and be explained by explanatory variables
# 46% can not be explained
# adjustes Rsquared: 19% can be explained


# PCA is done for the unexplained variance only! Refers to the Residuals of the data

# The first two RDA axes explain 34% of the variance, this is less  63% of the explained variability,
# So are these two axes sufficient?

# species scores represent the position of species in the axis in bi- or triplots
# site scores give the coordinates of sites in the space of the species
# site constraints give the coordinates of sites in the space of the environmental variables
# (explanatory)
# biplot scores give the coordinates for the environmental variable arrows

# Extraction of canonical coefficients from RDA object
coef(fam.rda1)


# Global test of RDA result
set.seed(111)
anova.cca(fam.rda1, step=1000)
# The Model is signficant

# How many axes explain a significant amount of variance in the data?
set.seed(111)
anova.cca(fam.rda1, by="axis", step=1000)
# two axes are significant

# SAMPLING CAMPAIGN 4 - AUTMUN 2013  ---------------------

# If we had used Hellinger TRansformation we would not need to scale
fam.rda2 <- rda(fam.fin[env.fin$samcam==4,]~., data=env.fin[env.fin$samcam==4,-c(1,3,4,7,10,11)], scale=TRUE)
summary(fam.rda2, display=NULL)

# Reduce environmental variables to significant ones 
stepping <- ordiR2step(rda(fam.fin[env.fin$samcam==4,] ~ 1 ,data=env.fin[env.fin$samcam==4,-c(1,3,4,7,10,11)]), scope=formula(fam.rda2),direction="forward",pstep=1000,trace=F)

anova(stepping)
# NO Constrained Component!


# Try Hellinger Distances:
fam.rda2 <- rda(fam.hel[env.fin$samcam==4,]~., data=env.fin[env.fin$samcam==4,-c(1,3,4,7,10,11)], scale=FALSE)
summary(fam.rda2, display=NULL)

# Reduce environmental variables to significant ones 
stepping <- ordiR2step(rda(fam.hel[env.fin$samcam==4,] ~ 1 ,data=env.fin[env.fin$samcam==4,-c(1,3,4,7,10,11)]), scope=formula(fam.rda2),direction="both",pstep=1000,trace=F)

anova(stepping)
# NO Constrained Component!

# AAArrrrrghhH!!!!!!

# Triplots ####

# For Sampling campaign 2 only, s.o.

par(mfrow=c(1,2))
plot(fam.rda1, scaling=1, main="Triplot RDA scaling 1 - wa scores")
fam.sc <- scores(fam.rda1, choices=1:2, scaling=1, display="sp")
arrows(0,0, fam.sc[,1], fam.sc[,2], length =0, lty=1, col="red")
# lines for species can be omitted

plot(fam.rda1, scaling=2, main="Triplot RDA scaling 1 - wa scores")
fam.sc2 <- scores(fam.rda1, choices=1:2, scaling=2, display="sp")
arrows(0,0, fam.sc2[,1], fam.sc2[,2], length =0, lty=1, col="red")



# Scaling 1
# How to interpret this plots?
# Zuur et al 2007, page 211 give the following rules (content modified)
# 1. Sites can be projected perpendicularly on the species of the explanatory lines
# and approximately the fitted value of the site along the variable
# 2. Distances between sites represent a two dimensional approximation of their fitted Euclidean
# distances
# 3. Angles between lines of species and explanatory variables reflect their correlation
# 4. Angles between lines for species or explanatory variables do not represent correlations
# an additional technical explanation is given in Legendre & Legendre 2012: 639-641


# WA scores vs. LC scores

# In the previous plots the position of sites are based in the weighted averages of species positions
# There is an open debate wether the position of the sites should be better expressed as linear combinations of the variables
# see:
# vignette("decision-vegan", package="vegan")

# Scaling 1
x11()
plot(fam.rda1, scaling=1, display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
arrows(0,0,fam.sc[,1], fam.sc[,2],length=0, lty=1, col="red")

# Scaling 2
x11()
plot(fam.rda1, scaling=2, display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")
fam.sc2 <- scores(fam.rda1, choices=1:2, scaling=2, display="sp")
arrows(0,0,fam.sc2[,1], fam.sc2[,2],length=0, lty=1, col="red")
# The rules from above still aplly for interpretation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## 3. partial tb-RDA #### 
# to eliminate the influence of repeated measurements 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fam.repmes <- rda(fam.hel ~ . + Condition(field.ID),data=env.fin[,-c(3,7)])

stepping <- ordiR2step(rda(fam.hel ~ 1 + Condition(field.ID),data=env.fin), scope=formula(fam.repmes),direction="both",pstep=1000,trace=F)

anova(stepping)

fam.repmes <- rda(fam.hel ~ age_class + n + mc +  Condition(field.ID),  data=env.fin[,-c(3,7)])
summary(fam.repmes, display=NULL)
summary(fam.repmes)

RsquareAdj(fam.repmes)

# We obtain 6 RDa axes and 15 PCA axis
# Further PCA axes are not displayed, as their eigenvalues are too small and would
# exceed the number of observations

# 10% of variance is explained by condition 
# 42% of our vriance can be constrained and be explained by explanatory variables
# 48% can not be explained
# adjustes Rsquared: 30% can be explained


# PCA is done for the unexplained variance only! Refers to the Residuals of the data

# The first two RDA axes explain 36% of the variance, this is less  86% of the explained variability,
# So are these two axes sufficient?

# species scores represent the position of species in the axis in bi- or triplots
# site scores give the coordinates of sites in the space of the species
# site constraints give the coordinates of sites in the space of the environmental variables
# (explanatory)
# biplot scores give the coordinates for the environmental variable arrows

# Extraction of canonical coefficients from RDA object
coef(fam.repmes)


anova.cca(fam.repmes, step=100)
anova.cca(fam.repmes, step=100, by="axis")
# Three Axis are significant


# Partial RDA triplots (with fitted site scores)

plot(fam.repmes, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(fam.repmes, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

plot(fam.repmes, scaling=1, main="Triplot RDA - scaling 1 - wa scores")
plot(fam.repmes, scaling=2, main="Triplot RDA - scaling 2 - wa scores")

plot(fam.repmes, type="n")
ordiellipse(fam.repmes, env.fin$age_class)
points(fam.repmes, col=env.fin$age_class)

plot(fam.repmes, type="n", scaling=2)
ordispider(fam.repmes, env.fin$age_class)
ordihull(fam.repmes, env.fin$age_class)
points(fam.repmes, col=env.fin$age_class)


# With ggplot2

fam.rm.sc1 <- scores(fam.repmes, display=c("sp", "lc","cn", "bp"), scaling=1)
fam.rm.sc1
fam.rm.sc2 <- scores(fam.repmes, display=c("sp", "lc","cn", "bp"), scaling=2)

# age_class centroids
ac.cid1 <- data.frame(fam.rm.sc1$centroids)
ac.cid1$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")
ac.cid2 <- data.frame(fam.rm.sc2$centroids)
ac.cid2$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")


# species
spec.cid1 <- data.frame(fam.rm.sc1$species)
spec.cid1$txt <- rownames(spec.cid1)
spec.cid2 <- data.frame(fam.rm.sc2$species)
spec.cid2$txt <- rownames(spec.cid2) 
# you can also put label = rownames() etc in the ggplot code. 
# But sometimes you may want it in a column for later

# arrows for continous environmental
env.ar1 <- data.frame(fam.rm.sc1$biplot)
env.ar1$txt <- rownames(env.ar1)
env.ar2 <- data.frame(fam.rm.sc2$biplot)
env.ar2$txt <- rownames(env.ar2)

mult <- attributes(scores(fam.repmes))$const

#sites
site_scores1 <- data.frame(scores(fam.repmes, scaling=1)$sites)
site_scores1$age_class <- env.fin$age_class

site_scores2 <- data.frame(scores(fam.repmes, scaling=2)$sites)
site_scores2$age_class <- env.fin$age_class

# site constraints 
site_constraints1 <- data.frame(fam.rm.sc1$constraints)
site_constraints1$age_class <- env.fin$age_class
site_constraints2 <- data.frame(fam.rm.sc2$constraints)
site_constraints2$age_class <- env.fin$age_class


#sc.mean=aggregate(site_scores[,1:2],list(group=site_scores$age_class),mean)

df_ell1.1 <- data.frame()
for(g in levels(site_scores1$age_class)){
  df_ell1.1 <- rbind(df_ell1.1, cbind(as.data.frame(with(site_scores1[site_scores1$age_class==g,],
                                                   veganCovEllipse(cov.wt(cbind(RDA1,RDA2),wt=rep(1/length(RDA1),length(RDA1)))$cov,center=c(mean(RDA1),mean(RDA2)))))
                                ,age_class=g))
}

df_ell1.2 <- data.frame()
for(g in levels(site_constraints1$age_class)){
  df_ell1.2 <- rbind(df_ell1.2, cbind(as.data.frame(with(site_constraints1[site_constraints1$age_class==g,],
                                                   veganCovEllipse(cov.wt(cbind(RDA1,RDA2),wt=rep(1/length(RDA1),length(RDA1)))$cov,center=c(mean(RDA1),mean(RDA2)))))
                                ,age_class=g))
}


df_ell2.1 <- data.frame()
for(g in levels(site_scores2$age_class)){
  df_ell2.1 <- rbind(df_ell2.1, cbind(as.data.frame(with(site_scores2[site_scores2$age_class==g,],
                                                   veganCovEllipse(cov.wt(cbind(RDA1,RDA2),wt=rep(1/length(RDA1),length(RDA1)))$cov,center=c(mean(RDA1),mean(RDA2)))))
                                ,age_class=g))
}

df_ell2.2 <- data.frame()
for(g in levels(site_constraints2$age_class)){
  df_ell2.2 <- rbind(df_ell2.2, cbind(as.data.frame(with(site_constraints2[site_constraints2$age_class==g,],
                                                   veganCovEllipse(cov.wt(cbind(RDA1,RDA2),wt=rep(1/length(RDA1),length(RDA1)))$cov,center=c(mean(RDA1),mean(RDA2)))))
                                ,age_class=g))
}


# Scaling 1
ggplot(ac.cid1, aes(x = RDA1, y = RDA2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid1, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*min(spec.cid1$RDA2),1.2*max(spec.cid1$RDA2)), x=c(1.2*min(spec.cid1$RDA1),1.3*max(spec.cid1$RDA1))) + 
                    #NB note that this is a convenience wrapper and may cut data out of your plot
                    #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar1[-c(1:4),] , aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar1[-c(1:4),] , aes(x= (mult + mult/20) * RDA1, y = (mult + mult/20) * RDA2,  label = env.ar1[-c(1:4),]$txt), size = 5, hjust = 0.5) +
                               # we add 10% to the text to push it slightly out from arrows
                               # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores1, aes(x = RDA1, y = RDA2, shape=env.fin$age_class)) +
  geom_point(data=site_constraints1, aes(x = RDA1, y = RDA2, shape=env.fin$age_class)) +
  #geom_path(data=df_ell1.1, aes(x=RDA1, y=RDA2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell1.2, aes(x=RDA1, y=RDA2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = RDA1, y = RDA2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()

# use env.ar1[-c(1:4),]  to ommit factor levels



#Scaling 2
ggplot(ac.cid2, aes(x = RDA1, y = RDA2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid2, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*mult*min(env.ar2$RDA2),1.2*max(site_scores2$RDA2)), x=c(1.2*mult*min(env.ar2$RDA1),1.3*max(site_scores2$RDA1))) + 
  #NB note that this is a convenience wrapper and may cut data out of your plot
  #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar2[-c(1:4),], aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar2[-c(1:4),], aes(x= (mult + mult/20) * RDA1, y = (mult + mult/20) * RDA2,  label = env.ar2[-c(1:4),]$txt), size = 5, hjust = 0.5) +
  # we add 10% to the text to push it slightly out from arrows
  # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores2, aes(x = RDA1, y = RDA2, shape=env.fin$age_class)) +
  geom_point(data=site_constraints2, aes(x = RDA1, y = RDA2, shape=env.fin$age_class)) +
  #geom_path(data=df_ell2.1, aes(x=RDA1, y=RDA2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell2.2, aes(x=RDA1, y=RDA2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = RDA1, y = RDA2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()




# Scaling 1
par(mfrow=c(1,1))
plot(fam.repmes, scaling=1, display=c("sp", "lc","cn"))

# Scaling 2
plot(fam.repmes, display=c("sp", "lc","cn"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4. partial db-RDA ####
# to reduce the influence of repeated measurements and use the Bray Curtis distance 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env.fin.db <- env.fin[,-c(7,10,11)] 
fam.fin.db <- fam.fin 

fam.db.repmes <- capscale(fam.fin.db  ~ . -field.ID + Condition(field.ID), distance="bray", data=env.fin.db, add=T)
summary(fam.db.repmes)
stepping <- ordiR2step(capscale(fam.fin.db ~ 1 + Condition(field.ID),distance="bray", data=env.fin.db, add=T), scope=formula(fam.db.repmes),direction="forward",pstep=1000,trace=F)
anova(stepping)

fam.db.repmes <- capscale(fam.fin ~ age_class + n + clay + Condition(field.ID), distance="bray", data=env.fin.db, add=TRUE)
summary(fam.db.repmes, display=NULL)

anova(fam.db.repmes, step=1000, perm.max=1000)
anova(fam.db.repmes, by="axis", step=1000, perm.max=1000)
# 2 axes are significant
# However field.ID was not used as a factor!!



#Triplots:
x11()
par(mfrow=c(1,2))
plot(fam.db.repmes, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(fam.db.repmes, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

x11()
par(mfrow=c(1,2))
plot(fam.db.repmes, scaling=1, main="Triplot RDA - scaling 1 - wa scores")
plot(fam.db.repmes, scaling=2, main="Triplot RDA - scaling 2 - wa scores")



fam.rm.sc1 <- scores(fam.db.repmes, display=c("sp", "lc","cn", "bp"), scaling=1)
fam.rm.sc1
fam.rm.sc2 <- scores(fam.db.repmes, display=c("sp", "lc","cn", "bp"), scaling=2)

# age_class centroids
ac.cid1 <- data.frame(fam.rm.sc1$centroids)
ac.cid1$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")
ac.cid2 <- data.frame(fam.rm.sc2$centroids)
ac.cid2$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")


# species
spec.cid1 <- data.frame(fam.rm.sc1$species)
spec.cid1$txt <- rownames(spec.cid1)
spec.cid2 <- data.frame(fam.rm.sc2$species)
spec.cid2$txt <- rownames(spec.cid2) 
# you can also put label = rownames() etc in the ggplot code. 
# But sometimes you may want it in a column for later

# arrows for continous environmental
env.ar1 <- data.frame(fam.rm.sc1$biplot)
env.ar1$txt <- rownames(env.ar1)
env.ar2 <- data.frame(fam.rm.sc2$biplot)
env.ar2$txt <- rownames(env.ar2)

mult <- attributes(scores(fam.db.repmes))$const

#sites
site_scores1 <- data.frame(scores(fam.db.repmes, scaling=1)$sites)
site_scores1$age_class <- env.fin$age_class

site_scores2 <- data.frame(scores(fam.db.repmes, scaling=2)$sites)
site_scores2$age_class <- env.fin$age_class

# site constraints 
site_constraints1 <- data.frame(fam.rm.sc1$constraints)
site_constraints1$age_class <- env.fin$age_class
site_constraints2 <- data.frame(fam.rm.sc2$constraints)
site_constraints2$age_class <- env.fin$age_class


#sc.mean=aggregate(site_scores[,1:2],list(group=site_scores$age_class),mean)

df_ell1.1 <- data.frame()
for(g in levels(site_scores1$age_class)){
  df_ell1.1 <- rbind(df_ell1.1, cbind(as.data.frame(with(site_scores1[site_scores1$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}

df_ell1.2 <- data.frame()
for(g in levels(site_constraints1$age_class)){
  df_ell1.2 <- rbind(df_ell1.2, cbind(as.data.frame(with(site_constraints1[site_constraints1$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}


df_ell2.1 <- data.frame()
for(g in levels(site_scores2$age_class)){
  df_ell2.1 <- rbind(df_ell2.1, cbind(as.data.frame(with(site_scores2[site_scores2$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}

df_ell2.2 <- data.frame()
for(g in levels(site_constraints2$age_class)){
  df_ell2.2 <- rbind(df_ell2.2, cbind(as.data.frame(with(site_constraints2[site_constraints2$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}


# Scaling 1
ggplot(ac.cid1, aes(x = CAP1, y = CAP2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid1, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*min(spec.cid1$CAP2),1.2*max(spec.cid1$CAP2)), x=c(1.2*min(spec.cid1$CAP1),1.3*max(spec.cid1$CAP1))) + 
  #NB note that this is a convenience wrapper and may cut data out of your plot
  #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar1[-c(1:4),] , aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar1[-c(1:4),] , aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar1[-c(1:4),]$txt), size = 5, hjust = 0.5) +
  # we add 10% to the text to push it slightly out from arrows
  # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores1, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  geom_point(data=site_constraints1, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  #geom_path(data=df_ell1.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell1.2, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()

# use env.ar1[-c(1:4),]  to ommit factor levels



#Scaling 2
ggplot(ac.cid2, aes(x = CAP1, y = CAP2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid2, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*mult*min(env.ar2$CAP2),1.2*max(site_scores2$CAP2)), x=c(1.2*mult*min(env.ar2$CAP1),1.3*max(site_scores2$CAP1))) + 
  #NB note that this is a convenience wrapper and may cut data out of your plot
  #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar2[-c(1:4),], aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar2[-c(1:4),], aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar2[-c(1:4),]$txt), size = 5, hjust = 0.5) +
  # we add 10% to the text to push it slightly out from arrows
  # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores2, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  geom_point(data=site_constraints2, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  #geom_path(data=df_ell2.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell2.2, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



