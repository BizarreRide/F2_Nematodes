###########################
# F2 Nematodes
# Data Processing
# Quentin Schorpp
# 06.08.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Community Data (Family level)
spe <- read.delim("Data/spe.txt")
str(spe)
spe$samcam <- as.factor(spe$samcam)
spe$field.ID <- as.factor(spe$field.ID)

indices <- spe[,which(names(spe)=="N"):47] # Nematode Indices - without gneral diversity indices (besides total abundance N)
fam <- spe[, which(names(spe)=="Aphelenchidae"):which(names(spe)=="Pratylenchidae")]
fety <- spe[, which(names(spe)=="carnivore"):which(names(spe)=="omnivore")]

env <- read.delim("Data/env.txt")
str(env)
env$samcam <- as.factor(env$samcam)
env$field.ID <- as.factor(env$field.ID)

env1 <- droplevels(env[16:45,])

spa <- env1[,c(3,5,6)]
soil <- env1[,c(8:10)]
soil <- cbind(soil,env1[,which(names(env1)=="pH"):which(names(env1)=="clay")])
groups <- env1[,c(2,11,12,15)]
climate30 <- env1[,which(names(env1)=="ata1"):which(names(env1)=="rad1")]
climate1 <- env1[,which(names(env1)=="ata2"):which(names(env1)=="rad2")]
mngmnt <- env1[,which(names(env1)=="soil_coverage"):which(names(env1)=="compaction")]


env2 <- env[rep(seq_len(nrow(env)), each=3),-1]

counts <- read.delim("Data/counts.txt")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculate biodiversity Indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# families Richness
indices$SR <- rowSums(fam >0)

# rarefaction
indices$rarefy <- vegan::rarefy(fam,90,se=F, MARGIN=1)

# Shannon entropy
indices$H <- vegan::diversity(fam, index="shannon")

# simpson dominance
indices$D <- vegan::diversity(fam, index="simpson")

# Pielou Evenness
indices$J <- indices$H/log(indices$SR)

# Hill's N1
indices$H1 <- exp(indices$H)

# camargo's diversity

# McIntosh dominance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Select environmental variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Management
mngmnt.pca <- rda(mngmnt[,-1], scale=T)
summary(mngmnt.pca, display=NULL)
ev <- mngmnt.pca$CA$eig
evplot(ev)
cleanplot.pca(mngmnt.pca)
mngmnt$intensity <- scores(mngmnt.pca, display="sites")[,1] # fertilisation would be a uncorrelated addition to intensity

#Soil
soil.pca <- rda(soil[,-c(1:3)], scale=T)
cleanplot.pca(soil.pca)
summary(soil.pca, display=NULL)
ev <- soil.pca$CA$eig
evplot(ev) # ommit: cn, silt, sand

# Climate
climate.pca <- rda(climate1, scale=T)
cleanplot.pca(climate.pca) # ata2 + hum2

climate30.pca <- rda(climate30, scale=T)
cleanplot.pca(climate30.pca) # ata1 + prec1

# Subset of (mostly) orthogonal variables
env.fin <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation"))
env.fin$intensity <- mngmnt$intensity

env.fin.pca <- rda(env.fin[,-c(1:4)], scale=T)
cleanplot.pca(env.fin.pca)

par(mfrow=c(1,1))
dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str(env1)

fam.av <- aggregate(fam[sapply(fam,is.numeric)], list(env1$field.ID),mean)
fam.av <- fam.av[,-1]
env.fac <- env.fin[!sapply(env.fin,is.numeric)]
env.fac <- env.fac[-c(16:27),]
env.av <- aggregate(env.fin[sapply(env.fin,is.numeric)], list(env1$field.ID),mean) 
env.av <- cbind(env.fac,env.av)




