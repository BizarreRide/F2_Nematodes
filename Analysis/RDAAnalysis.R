#%%%%%%%%%%%%%%%%%%%%%%%%%
# F2 Nematodes
# PCA and RDA Analysis
# Quentin Schorpp
# 06.08.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/RMAkeLikeFile.R")
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


## 1. Data Preparation ####
# after Ralf Schäfer 
# https://www.youtube.com/watch?v=gY_iktfpSpQ
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
env.fin <- subset(env1, select=c("field.ID","age_class","crop","samcam","age","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation"))
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
sort(faraway::vif(env.fin[,-c(1:4)]))
# Borcard et al 2011: 175 argue that in the case of RDA, VIFs > 10 should be avoided
# Theres a video of Ralf schäfer, that tells about how to deal with multicollinearity in multiple regression


# checking gradient length with detrended correspondence analysis
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

## 2. PCA ####


#******************************************************Repeated Measurements Data*******************************************************************####

## 3. RDA ####
### for each sampling campaign ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### SAMPLING CAMPAIGN 2 - AUTMUN 2012 ---------------------

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

### SAMPLING CAMPAIGN 4 - AUTMUN 2013  ---------------------

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

### Triplots ####

# For Sampling campaign 2 only, s.o.

par(mfrow=c(1,2))
plot(fam.rda1, scaling=1, main="Triplot RDA scaling 1 - wa scores")
fam.sc <- scores(fam.rda1, choices=1:2, scaling=1, display="sp")
arrows(0,0, fam.sc[,1], fam.sc[,2], length =0, lty=1, col="red")
# lines for species can be omitted

plot(fam.rda1, scaling=2, main="Triplot RDA scaling 2 - wa scores")
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

field.ID <- env$field.ID
env1 <- env.fin[,1]
env2 <- env.fin[,2]

x <- varpart(fam, ~ .,env2,data=env1, transfo = "hel")
showvarparts(x)
plot(x)
fam.hel1 <- fam.hel
fam.hel1 <- fam.hel1[!env.fin$crop=="Maize",]
env.fin1 <- env.fin[-c(3,8,11,12,15)]
env.fin1 <- env.fin1[!env.fin$crop=="Maize",]
field.ID <- env.fin$field.ID
field.ID <- env.fin[!env.fin$crop=="Maize",]$field.ID

fam.repmes1 <- rda(fam.hel ~ . - c + Condition(field.ID),env.fin) #+ Condition(field.ID) !env.fin$crop=="Maize",
fam.repmes0 <- rda(fam.hel1 ~ 1 + Condition(field.ID),env.fin1)

vif(fam.repmes1)

step.res <- ordiR2step(fam.repmes0, scope=formula(fam.repmes1), perm.max=200, trace=F)
step.res$anova
anova(step.res)

x11()
par(mfrow=c(1,2))
plot(fam.repmes1, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(fam.repmes1, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

stepping <- ordiR2step(rda(fam.hel ~ 1 + Condition(field.ID),data=env.fin), scope=formula(fam.repmes),direction="both",pstep=1000,trace=F)
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

# env.fin$age <- as.numeric(env.fin$age)
#env.fin$field.ID <- as.numeric(env.fin$field.ID)
#env.fin$field.ID <- as.factor(env.fin$field.ID)
#env.fin$field.ID <- as.factor(chartr("ABCDEFGHIA","1234567890",  env.fin$field.ID))
str(env.fin)

decorana(fam.fin)
fam.rel <- fam.fin/indices$N
decorana(fam.rel)
fam.hel <- decostand(fam.fin, "hel")
decorana(fam.hel)

fam.db.repmes1 <- capscale(fam.rel  ~ .   + Condition(samcam), distance="bray", data=env.fin[,-1], add=T);ordiplot(fam.db.repmes1,type="p", scaling=1)
fam.db.repmes0 <- capscale(fam.rel  ~ 1   + Condition(samcam), distance="bray", data=env.fin, add=T)

stepping <- ordiR2step(fam.db.repmes0, fam.db.repmes1,direction="forward", perm.max=200, pstep=1000,trace=F)
#stepping <- ordistep(fam.db.repmes0, fam.db.repmes1,direction="both", perm.max=200, pstep=1000,trace=F)
summary(stepping)
anova(stepping)

fam.db.repmes <- capscale(fam.rel ~ age_class + clay + Condition(samcam), distance="bray", data=env.fin[,-1], add=TRUE)
summary(fam.db.repmes)

anova(fam.db.repmes, step=1000, perm.max=1000)
anova(fam.db.repmes, by="axis", step=1000, perm.max=1000)
# 2 axes are significant

# Notes:
# If field.ID was taken as conditioning factor, then only the parameters: samcam, pH, mc, ata2, hum2, ata1, prec1, fertilisation and intensity, work.
# This rda()s with the other parameters (each by each) always create the error message:
#Error in cbind(x$CCA$v, x$CA$v) : number of rows of matrices must match (see arg 2)
# If field.ID is the Condition() argument, and no SIlphie fields are included, I get the same error message.

R2.all <- RsquareAdj(fam.db.repmes)$adj.r.squared
RsquareAdj(fam.db.repmes1)
RsquareAdj(fam.db.repmes0)
RsquareAdj(fam.db.repmes)

vif.cca(fam.db.repmes)

coef(fam.db.repmes)


## Ordinary Plots #####
colvec1 <- palette()[1:5]
#colvec1 <- brewer.pal(5, "Set1")
colvec <- rep(rep(colvec1, each=3),2)

par(mfrow=c(1,1))
plot1 <- ordiplot(fam.db.repmes,type="p", scaling=1, display=c("sp","lc","cn")) #, display=c("sp","lc","cn")
ordiequilibriumcircle(fam.db.repmes,plot1)
identify(plot1,"sp", labels=names(fam.fin), cex=1.0)
points(plot1, "constraints", pch=25, bg=colvec, cex=1) 
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pt.bg=colvec1, pch=25, pt.cex=0.7)



plot1 <- ordiplot(fam.db.repmes,type="p", scaling=2)
ordiequilibriumcircle(fam.db.repmes,plot1)
identify(plot1,"sp", labels=names(fam.fin), cex=1.0)
points(plot1, "sites", pch=25, bg=colvec, cex=1) 
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pt.bg=colvec1, pch=25, pt.cex=0.7)

cleanplot.pca(fam.db.repmes)

require(vegan3d)
plot3d <- ordiplot3d(fam.db.repmes, scaling=2, angle=50, type="n")
points(plot3d, "points", pch=25, bg=colvec, cex=0.7) 
text(plot3d, "arrows", col="blue", pos=3)
ordirgl(fam.db.repmes, size=2)

ordirgl(fam.db.repmes, display = "sites", type = "p",col=colvec, size=10)
orgltext(fam.db.repmes, display = "species", choices = 1:3,pch=25, col="orange", size=10)
orgltext(fam.db.repmes, text, display = "species", choices = 1:3, justify = "center", adj = 0.5)
rgl.quit()

# variation partitioning
field.ID <- groups$field.ID
age_class <- groups$age_class
numerics <- soil[,c(4,5,7,11)]
vp <- varpart(fam.hel,numerics,~age_class, ~field.ID)

plot(vp)

####
# Triplots:
x11()
par(mfrow=c(1,2))
plot(fam.db.repmes, scaling=1,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 1 - lc scores")
plot(fam.db.repmes, scaling=2,display=c("sp","lc", "cn"), main="Triplot RDA - scaling 2 - lc scores")

x11()
par(mfrow=c(1,2))
plot(fam.db.repmes, scaling=1, main="Triplot RDA - scaling 1 - wa scores")
plot(fam.db.repmes, scaling=2, main="Triplot RDA - scaling 2 - wa scores")

## ggplot2 Triplots ####

fam.rm.sc1 <- scores(fam.db.repmes, display=c("sp", "lc","cn", "bp"), scaling=1)
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
Scaling1 <- ggplot(ac.cid1, aes(x = CAP1, y = CAP2)) +
  coord_cartesian(y=c(1.2*min(spec.cid1$CAP2),1.3*max(spec.cid1$CAP2)), x=c(1.5*min(spec.cid1$CAP1),1.5*max(spec.cid1$CAP1))) + 
  geom_hline(aes(yintercept=0),lty=3) + 
  geom_vline(aes(xintercept=0),lty=3) + 
  
  ## Site Scores: ##
  #geom_path(data=df_ell1.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) + 
  #stat_ellipse(data=df_ell1.1, geom = "polygon", alpha = 1/2, aes(x=CAP1, y=CAP2, fill=age_class, color=age_class)) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+
  #geom_point(data=site_scores1, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  
  ## Site constraints: ##
  geom_path(data=df_ell1.2, aes(x=CAP1, y=CAP2, bg=age_class, color=age_class), size=0.5, linetype=2) +
  stat_ellipse(data=df_ell1.2, geom = "polygon", alpha = 1/2, aes(x=CAP1, y=CAP2, fill=age_class, color=age_class)) +
  geom_point(pch=13, size=3, col="red") + # points for age_class means
  geom_text(aes(label = txt), family="Arial Narrow",fontface="italic", size=5) + # text for age_class means
  geom_point(data=site_constraints1, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  
  # species scores:
  geom_segment(data = spec.cid1[-c(1:4),] , aes(x = 0, xend = 0.9*CAP1, y = 0, yend = 0.96*CAP2),
               arrow = arrow(angle=20, type="closed",length = unit(0.25, "cm")), lty=3,lwd=0.5, colour = "dark grey") +
  geom_text(data = spec.cid1, aes(label = txt), cex=3, colour = "black") +
  
  # Environmental constraints
  geom_segment(data = env.ar1[-c(1:4),] , aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm"), type="closed", angle=20), colour = "black") + #grid is required for arrow to work.
  geom_text(data = env.ar1[-c(1:4),] , aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar1[-c(1:4),]$txt), size = 3, hjust = 0.5) +
  
  mytheme + theme(legend.position="none")


ggsave(file="partial db-RDASacling-1.svg",Scaling1, width=11, height=9, units="cm", dpi=300)

# use env.ar1[-c(1:4),]  to ommit factor levels

Scaling2 <- ggplot(ac.cid2, aes(x = CAP1, y = CAP2)) +
  coord_cartesian(y=c(1.2*mult*min(env.ar2$CAP2),1.2*max(site_scores2$CAP2)), x=c(1.2*mult*min(env.ar2$CAP1),1.3*max(site_scores2$CAP1))) + 
  geom_hline(aes(yintercept=0),lty=3) + 
  geom_vline(aes(xintercept=0),lty=3) + 
  
  ## Site Scores: ##
  #geom_path(data=df_ell2.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) + 
  #stat_ellipse(data=df_ell2.1, geom = "polygon", alpha = 1/2, aes(x=CAP1, y=CAP2, fill=age_class, color=age_class)) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+
  #geom_point(data=site_scores2, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  
  ## Site constraints: ##
  stat_ellipse(data=df_ell2.2, geom = "polygon", alpha = 1/2, aes(x=CAP1, y=CAP2, fill=age_class, color=age_class)) +
  geom_path(data=df_ell2.2, aes(x=CAP1, y=CAP2, bg=age_class, color=age_class), size=0.5, linetype=2) +
  geom_point(pch=13, size=3, col="red") + # points for age_class means
  geom_text(aes(label = txt), family="Arial Narrow",fontface="italic", size=5) + # text for age_class means
  geom_point(data=site_constraints2, aes(x = CAP1, y = CAP2, shape=env.fin$age_class)) +
  
  # species scores:
  geom_segment(data = spec.cid2[-c(1:4),] , aes(x = 0, xend = 0.9*CAP1, y = 0, yend = 0.96*CAP2),
               arrow = arrow(angle=20, type="closed",length = unit(0.25, "cm")), lty=3,lwd=0.5, colour = "dark grey") +
  geom_text(data = spec.cid2, aes(label = txt), cex=3, colour = "black") +
  
  # Environmental constraints
  geom_segment(data = env.ar2[-c(1:4),] , aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm"), type="closed", angle=20), colour = "black") + #grid is required for arrow to work.
  geom_text(data = env.ar2[-c(1:4),] , aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar2[-c(1:4),]$txt), size = 3, hjust = 0.5) +
  
  mytheme + theme(legend.position="none")

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
  mytheme + legend.position("none")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Some testing rda() #####
env.fin.x <- env.fin[!env.fin$crop=="Silphie",] #1
env.fin.x <- env.fin[!env.fin$crop=="Maize",] #2
env.fin.x <- env.fin#3
fam.fin.x <- fam.fin[!env.fin$crop=="Silphie",] #1
fam.fin.x <- fam.fin[!env.fin$crop=="Maize",] #2
fam.fin.x <- fam.fin #3
condition <- as.factor(chartr("1234567890", "ABCDEFGHIA",  env.fin$field.ID))

### Condition: Age Class

fam.db.repmes <- capscale(fam.fin.x ~ . + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ crop + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ age_class + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ c + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ n + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ clay + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
# NOT: age_class, crop,   
# Not, but works alone: c,n,clay
# for maize fields only: subscripts out of bounds!

# All Others:
fam.db.repmes <- capscale(fam.fin.x ~ samcam + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ field.ID + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ pH + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ mc + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ ata1 + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ ata2 + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ hum2 + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ prec1 + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ intensity + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ fertilisation + Condition(age_class), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
# All work!
# for maize fields only: subscripts out of bounds!

### Condition field.ID

fam.db.repmes <- capscale(fam.fin.x ~ . + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ field.ID + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ age_class + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ crop + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ n + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ c + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ clay + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
# Not: field.ID,age_class,crop,n,c,clay,
# Not but work alone: NONE
# for maize fields only: number of rows of matrices must match!

# All Others:
fam.db.repmes <- capscale(fam.fin.x ~ samcam + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ pH + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ mc + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ ata1 + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ ata2 + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ hum2 + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ prec1 + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ intensity + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
fam.db.repmes <- capscale(fam.fin.x ~ fertilisation + Condition(field.ID), distance="bray", data=env.fin.x, add=TRUE);summary(fam.db.repmes)
# All work!
# for maize fields only: number of rows of matrices must match!


#******************************************************Averaged Data *******************************************************************************####


# db - RDA with averaged data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## All families ####

fam.db1 <- capscale(fam.av  ~ age_class + pH + n + mc + prec1 + intensity + ata1 , distance="bray", add=TRUE, data=env.av)
fam.db0 <- capscale(fam.av  ~ 1 , distance="bray", data=env.av, add=T)
summary(fam.db1)

stepping <- ordiR2step(fam.db0, fam.db1,direction="forward",pstep=1000,trace=F)
warnings()
anova(stepping)

fam.db <- capscale(fam.av ~ age_class + n + intensity, distance="bray", data=env.av, add=TRUE)
summary(fam.db, display=NULL)

anova(fam.db, step=1000, perm.max=1000)
anova(fam.db, by="axis", step=1000, perm.max=1000)
# Three significant Axes

# Triplots
colvec1 <- palette()[1:5]
#colvec1 <- brewer.pal(5, "Set1")
colvec <- c(rep(colvec1, each=3),rep(colvec1[5],3))
plot1 <- ordiplot(fam.db,type="p", scaling=1)
ordiequilibriumcircle(fam.db,plot1)
identify(plot1,"sp", labels=names(fam.av), cex=1.0)
points(plot1, "sites", pch=25, bg=colvec, cex=0.7) 
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pt.bg=colvec1, pch=25, pt.cex=0.7)




cleanplot.pca(fam.db)

dbRDA <- fam.db
dbenv <- env.av
source("Analysis/dbRDA_ggplots.R")
Scaling1
Scaling2

## selected families ####

# presence-absence transformation to calculate species number per site
fam.pa <- decostand(fam.av, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)

# remove species that occur at less than 5 sites
fam.av.fin <- fam.av[, !fam.sum<5]

fam.db1 <- capscale(fam.av.fin  ~ age_class + pH + n + mc + prec1 + intensity + ata1 , distance="bray", data=env.av, add=T)
fam.db0 <- capscale(fam.av.fin  ~ 1 , distance="bray", data=env.av, add=T)

summary(fam.db1)
stepping <- ordiR2step(fam.db0, fam.db1,direction="forward",pstep=1000,trace=F)
anova(stepping)

fam.db <- capscale(fam.av.fin ~ age_class + n + intensity, distance="bray", data=env.av, add=TRUE)
summary(fam.db, display=NULL)

anova(fam.db, step=1000, perm.max=1000)
anova(fam.db, by="axis", step=1000, perm.max=1000)

# Colored
require(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
dev.off()

colvec1 <- palette()[1:5]
#colvec1 <- brewer.pal(5, "Set1")
colvec <- c(rep(colvec1, each=3),rep(colvec1[5],3))
plot1 <- ordiplot(fam.db, scaling=2)
ordiequilibriumcircle(fam.db,plot1)
identify(plot1,"sp", labels=names(fam.av.fin), cex=1.0)
points(plot1, "sites", pch=25, bg=colvec, cex=0.7) 
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pt.bg=colvec1, pch=25, pt.cex=0.7)

# Grey
colvec1 <- gray.colors(5, start = 0, end = 1, gamma = 2.2, alpha = NULL)
colvec <- c(rep(colvec1, each=3),rep(colvec1[5],3))
plot1 <- ordiplot(fam.db, scaling=1)
ordiequilibriumcircle(fam.db,plot1)
identify(plot1,"sp", labels=names(fam.av.fin), cex=1.0)
points(plot1, "sites", pch=25, bg=colvec, cex=0.7) 
legend("bottomleft", legend=c("Sp_O", "SP_I2", "SP_I1", "Sp_Y", "Cm"), text.col="black", pt.bg=colvec1, pch=25, pt.cex=0.7) 

require(vegan3d)
plot3d <- ordiplot3d(fam.db, scaling=2, angle=50, type="n")
points(plot3d, "points", pch=25, bg=colvec, cex=0.7) 
text(plot3d, "arrows", col="blue", pos=3)
ordirgl(fam.db, size=2)
ordirgl(fam.db, display = "species", type = "t")
rgl.quit()

# I'd say:
# Maize fields are very similar, 
# Old SIlphie fields are very similar,
# Intermediate 1 stages are very similar,
# Young silphie fileds can be very different, so they can't be categorized, or said to be similar to one of the other groups
# Intermediate 2 are rahter similar to intermediate 1, however one seems to develop towards older stages. And they seem to be associated with 
# Tylenchidae

dbRDA <- fam.db
dbenv <- env.av
source("Analysis/dbRDA_ggplots.R")
Scaling1
Scaling2




