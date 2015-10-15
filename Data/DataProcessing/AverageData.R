#§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Calculate Average Data
# Quentin Schorpp
# 13.10.2015
#$$$$$$$$$$$$$$$$$$$$$$$$$

## Averaged environmental variables ####
env.fac <- env1[!sapply(env1,is.numeric)]
env.fac <- env.fac[-c(16:27),]
env.av <- aggregate(env1[sapply(env1,is.numeric)], list(env1$field.ID),mean) 
env.av <- cbind(env.fac,env.av)
rm(env.fac)

#env.av.fin <- subset(env.av, select=c("field.ID","age_class","crop","pH","mc","c","age","clay","ata2","hum2","ata1","prec1", "fertilisation", "intensity"))


## Averaged families Data ####
fam.av <- aggregate(fam.org[sapply(fam.org,is.numeric)], list(env1$field.ID),mean)
fam.av <- fam.av[,-1]

## Differences in SIlphie fields between  Sampling Campaigns ####
fam.slope <- fam.org[!env1$crop=="Maize",]
fam.slope <- fam.slope[13:24,] - fam.slope[1:12,]

## Averaged families Data ####
counts.av <- aggregate(counts[sapply(counts,is.numeric)], list(env.org$field.ID),mean)
counts.av <- counts.av[,-1]
