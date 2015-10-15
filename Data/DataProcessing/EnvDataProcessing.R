#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Processing of environmental data
# Quentin Schorpp
# 13.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$


# Environmental variables ####

# Slice and categorize env data, according to 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spa <- env1[,c(3,5,6)]
soil <- env1[,c(8:10)]
soil <- cbind(soil,env1[,which(names(env1)=="pH"):which(names(env1)=="clay")])
groups <- env1[,c(2,11,12,15)]
climate30 <- env1[,which(names(env1)=="ata1"):which(names(env1)=="rad1")]
climate1 <- env1[,which(names(env1)=="ata2"):which(names(env1)=="rad2")]
mngmnt <- env1[,which(names(env1)=="soil_coverage"):which(names(env1)=="compaction")]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Select orthogonal environmental variables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## PCA of environmental subsets ####
# Management
mngmnt.pca <- rda(mngmnt[,-1], scale=T)
summary(mngmnt.pca, display=NULL)
ev <- mngmnt.pca$CA$eig
evplot(ev)
cleanplot.pca(mngmnt.pca)
mngmnt$intensity <- scores(mngmnt.pca, display="sites")[,1] # fertilisation would be a uncorrelated addition to intensity
env1$intensity <- scores(mngmnt.pca, display="sites")[,1] # fertilisation would be a uncorrelated addition to intensity
# intensity summarizes all management parameters, except fertilisation
rm(mngmnt.pca)

#Soil
soil.pca <- rda(soil[,-c(1:3)], scale=T)
cleanplot.pca(soil.pca)
summary(soil.pca, display=NULL)
ev <- soil.pca$CA$eig
evplot(ev) # ommit non-othogonal: cn, silt, sand
rm(soil.pca)


# Climate
climate.pca <- rda(climate1, scale=T)
cleanplot.pca(climate.pca) # select:ata2 + hum2
rm(climate.pca)
climate30.pca <- rda(climate30, scale=T)
cleanplot.pca(climate30.pca) # select:ata1 + prec1
rm(climate30.pca)


## Subset of (mostly) orthogonal variables ####
env.fin <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","age","clay","ata2","hum2","ata1","prec1", "fertilisation", "intensity"))

env.fin.pca <- rda(env.fin[,-c(1:4)], scale=T)
cleanplot.pca(env.fin.pca)
rm(env.fin.pca)

par(mfrow=c(1,1))
dev.off()
rm(ev)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
