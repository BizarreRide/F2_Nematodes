###########################
# F2 Nematodes
# Data Processing
# Quentin Schorpp
# 06.08.2015
###########################


# Load and slice Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Community Data (Family level)
nema <- read.delim("Data/spe.txt") # Not really species data
nema$samcam <- as.factor(nema$samcam)
nema$field.ID <- as.factor(nema$field.ID)

indices <- nema[,which(names(nema)=="N"):47] # Nematode Indices - without gneral diversity indices (besides total abundance N)
fam <- nema[, which(names(nema)=="Aphelenchidae"):which(names(nema)=="Pratylenchidae")]
fety <- nema[, which(names(nema)=="carnivore"):which(names(nema)=="omnivore")]

env <- read.delim("Data/env.txt")
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

counts <- read.delim("Data/counts.txt")
counts.env <- env[rep(seq_len(nrow(env)), each=3),-1]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Indices ####
# Calculate and add biodiversity Indices
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


# Environmental variables ####
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
env.fin <- subset(env1, select=c("field.ID","age_class","crop","samcam","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation", "intensity"))

env.fin.pca <- rda(env.fin[,-c(1:4)], scale=T)
cleanplot.pca(env.fin.pca)
rm(env.fin.pca)

par(mfrow=c(1,1))
dev.off()
rm(ev)
## Averaged environmental variables ####
env.fac <- env.fin[!sapply(env.fin,is.numeric)]
env.fac <- env.fac[-c(16:27),]
env.av <- aggregate(env.fin[sapply(env.fin,is.numeric)], list(env1$field.ID),mean) 
env.av <- cbind(env.fac,env.av)
rm(env.fac)

env.av.fin <- subset(env.av, select=c("field.ID","age_class","crop","pH","mc","c","n","clay","ata2","hum2","ata1","prec1", "fertilisation", "intensity"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Nematode Families ####
# Select orthogonal environmental variables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Averaged Data ####
fam.av <- aggregate(fam[sapply(fam,is.numeric)], list(env1$field.ID),mean)
fam.av <- fam.av[,-1]

## Differences between  Sampling Campaigns ####
fam.slope <- fam[!env1$crop=="Maize",]
fam.slope <- fam.slope[13:24,] - fam.slope[1:12,]

## Selected families ####

fam.pa <- decostand(fam, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)
# remove species that occur at less than 5 sites
fam.av.fin <- fam.av[, !fam.sum<5]

### Averaged Data
# presence-absence transformation to calculate species number per site
fam.pa <- decostand(fam.av, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)
# remove species that occur at less than 5 sites
fam.av.fin <- fam.av[, !fam.sum<5]

rm(fam.pa, fam.sum)
