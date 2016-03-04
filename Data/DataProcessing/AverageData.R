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
fam.av.org <- aggregate(. ~ env1$field.ID,fam.org,mean)
fam.av.org <- fam.av.org[,-1]

fam.av.usc <- aggregate(. ~ env1$field.ID,fam.usc,mean)
fam.av.usc <- fam.av.usc[,-1]

## Differences in SIlphie fields between  Sampling Campaigns ####
fam.slope.org <- fam.org[!env1$crop=="Maize",]
fam.slope.org <- fam.slope.org[13:24,] - fam.slope.org[1:12,]

fam.slope.usc <- fam.usc[!env1$crop=="Maize",]
fam.slope.usc <- fam.slope.usc[13:24,] - fam.slope.usc[1:12,]

## Averaged count Data ####
counts.av <- aggregate(. ~ field.ID, counts,mean)
counts.av <- counts.av[,-1]
