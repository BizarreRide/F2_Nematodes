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



