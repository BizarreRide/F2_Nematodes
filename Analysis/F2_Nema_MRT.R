###########################
# F2 Nematodes
# MRT + ISA
# Quentin Schorpp
# 15.09.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# family composition ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env.x <- env.fin[,-c(3,7,8,9,15)]

fam.bc <- vegdist(fam, "bray")
fam.hel <- decostand(fam, "hel")

require(mvpart)
fam.mvpart <- mvpart(data.matrix(fam.bc) ~ ., env.x[,-1], margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




# family abundance ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
env.x <- env.fin[,-c(3,7,8,9,15)]
fam.count <- fam*counts$counts

fam.bc <- vegdist(fam.count, "bray")
fam.hel <- decostand(fam.count, "hel")

require(mvpart)
fam.mvpart <- mvpart(data.matrix(fam.bc) ~ ., env.x[,-1], margin=0.08, cp=0, xv="pick", xval=10, xvmult=100, which=4)
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





