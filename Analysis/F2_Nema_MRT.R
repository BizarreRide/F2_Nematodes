###########################
# F2 Nematodes
# MRT + ISA
# Quentin Schorpp
# 15.09.2015
###########################

# Note: The package randomforest could provide more interesting features!

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
#env.fin$c <-  env1$c

source("Data/DataProcessing/AverageData.R")

data <- fam.av
source("Data/DataProcessing/FamDatProcessing.R") 

fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dev.off()

# family Proportion ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env.x <- env.av[,c("age_class","crop","mc","n","pH","ata1","prec1")] 
env.x <- env.av[,c("age","crop","mc","n","pH","ata1","prec1")] 

fam.bc <- vegdist(fam.rel, "bray")
fam.hel <- decostand(fam.rel, "hel")

require(mvpart)
fam.mvpart <- mvpart(data.matrix(fam.bc) ~ . , env.x, margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
fam.mvpart <- mvpart(data.matrix(fam.hel) ~ ., env.x, margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
summary(fam.mvpart)
printcp(fam.mvpart)

par(mfrow=c(1,2))
hist(residuals(fam.mvpart), col="grey")
plot(predict(fam.mvpart), residuals(fam.mvpart), main="Residuals vs Predicted")
abline(h=0, lty=3, col="grey")

# group composition
fam.mvpart$where

# group identity
(groups.mrt  <- levels(as.factor(fam.mvpart$where)))

# Nematode composition of first leaf
fam.hel[which(fam.mvpart$where==groups.mrt[1]),]

# Environmental variables of first leaf
env.fin[which(fam.mvpart$where==groups.mrt[1]),]

# Table and pie charts of fish composition of leaves
leaf.sum  <-  matrix(0, length(groups.mrt), ncol(fam.hel))
colnames(leaf.sum)  <- colnames(fam.hel)

for(i in 1:length(groups.mrt)) {
  leaf.sum[i,] <- 
    apply(fam.hel[which(fam.mvpart$where==groups.mrt[i]),],2,sum)
}
leaf.sum

par(mfrow=c(1,2))
for(i in 1:length(groups.mrt)){
  pie(which(leaf.sum[i,]>0), raius=1, main=c("leaf #", groups.mrt[i]))
}

# Extracting MRT results from an mvpart object
# Pckages MVPARTwrap and rdaTest must have been loaded

# Packages not aivailable


# Indicator species search on the MRT result
require(labdsv)
fam.mvpart.indval <- indval(fam.hel, fam.mvpart$where)
fam.mvpart.indval$pval

# for each significant species find the leaf with the highest IndVal
fam.mvpart.indval$maxcls[which(fam.mvpart.indval$pval<=0.05)]

# Indval value in the best leaf for each significant species
fam.mvpart.indval$indcls[which(fam.mvpart.indval$pval<=0.05)]

# Table of the significant indicator families
gr <- fam.mvpart.indval$maxcls[fam.mvpart.indval$pval<=0.05]
iv <- fam.mvpart.indval$indcls[fam.mvpart.indval$pval<=0.05]
pv <- fam.mvpart.indval$pval[fam.mvpart.indval$pval<=0.05]
fr <- apply(fam.usc>0,2,sum)[fam.mvpart.indval$pval<=0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




# family abundance ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
env.x <- env.av[,c("age_class","crop","mc","n","pH","ata1","prec1")] 
env.x <- env.av[,c("age","crop","mc","n","pH","ata1","prec1")] 


fam.bc <- vegdist(fam.usc, "bray")
fam.hel <- decostand(fam.usc, "hel")

require(mvpart)
fam.mvpart <- mvpart(data.matrix(fam.bc) ~ ., env.x, margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
fam.mvpart <- mvpart(data.matrix(fam.hel) ~ ., env.x[,-1], margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
summary(fam.mvpart)
printcp(fam.mvpart)

par(mfrow=c(1,2))
hist(residuals(fam.mvpart), col="grey")
plot(predict(fam.mvpart), residuals(fam.mvpart), main="Residuals vs Predicted")
abline(h=0, lty=3, col="grey")

# group composition
fam.mvpart$where

# group identity
(groups.mrt  <- levels(as.factor(fam.mvpart$where)))

# Nematode composition of first leaf
fam.hel[which(fam.mvpart$where==groups.mrt[1]),]

# Environmental variables of first leaf
env.fin[which(fam.mvpart$where==groups.mrt[1]),]

# Table and pie charts of fish composition of leaves
leaf.sum  <-  matrix(0, length(groups.mrt), ncol(fam.hel))
colnames(leaf.sum)  <- colnames(fam.hel)

for(i in 1:length(groups.mrt)) {
  leaf.sum[i,] <- 
    apply(fam.hel[which(fam.mvpart$where==groups.mrt[i]),],2,sum)
}
leaf.sum

par(mfrow=c(2,2))
for(i in 1:length(groups.mrt)){
  pie(which(leaf.sum[i,]>0), raius=1, main=c("leaf #", groups.mrt[i]))
}

# Extracting MRT results from an mvpart object
# Pckages MVPARTwrap and rdaTest must have been loaded

# Packages not aivailable


# Indicator species search on the MRT result
require(labdsv)
fam.mvpart.indval <- indval(fam.hel, fam.mvpart$where)
fam.mvpart.indval$pval

# for each significant species find the leaf with the highest IndVal
fam.mvpart.indval$maxcls[which(fam.mvpart.indval$pval<=0.05)]

# Indval value in the best leaf for each significant species
fam.mvpart.indval$indcls[which(fam.mvpart.indval$pval<=0.05)]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





