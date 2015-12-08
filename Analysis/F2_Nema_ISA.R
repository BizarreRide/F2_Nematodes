###########################
# F2 Nematodes
# ISA
# Quentin Schorpp
# 19.08.2015
###########################

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

# Faunal Profile on subsample abundances
fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)



# presence-absence transformation to calculate species number per site
fam.pa <- decostand(fam.usc, "pa")
# calculate sum per species
fam.sum <- apply(fam.pa,2,sum)
sort(fam.sum)

# remove species that occur at less than 5 sites
fam.fin <- fam.usc[, !fam.sum<5]

env.x <- env.av[,c("age","crop","mc","n","pH","ata1","prec1")] # try instead!
env.x <- env.av[,c("age_class","crop","mc","n","pH","ata1","prec1")] 
summary(env.x)

#env.x$field.ID <- as.factor(env.x$field.ID)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Indicator Species Analysis for Age Class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam.hel <- decostand(fam.fin, "hel")

library(labdsv)
iva <- indval(fam.fin, as.numeric(env.x$age_class)) # almost the same result
iva <- indval(fam.hel, as.numeric(env.x$age_class))

p = 0.1

# Table of the significant indicator families
gr <- iva$maxcls[iva$pval<=p]
iv <- iva$indcls[iva$pval<=p]
pv <- iva$pval[iva$pval<=p]
fr <- apply(fam.hel>0,2,sum)[iva$pval<=p]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg


# Indicator Species Analysis for Clusters based on environmental data (incl. Factors)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Gower dissimilarity matrix

env.gow <- cluster::daisy(env.x, "gower")


env.kmeans <- kmeans(env.gow, centers = 5, nstart = 100)
env.kmeans$cluster
iva <- indval(fam.hel,env.kmeans$cluster)

p = 0.1
# Table of the significant indicator families
gr <- iva$maxcls[iva$pval<=p]
iv <- iva$indcls[iva$pval<=p]
pv <- iva$pval[iva$pval<=p]
fr <- apply(fam.hel>0,2,sum)[iva$pval<=p]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg


# Indicator Species Analysis for Ward Clusters based on species data 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fam.bc <- vegdist(fam.fin, method="bray")
fam.ward <- hclust(fam.bc, method = "ward.D")
groups<- as.factor(cutree(fam.ward,k=5))
plot(fam.ward, labels=interaction(groups,env.av$field.ID))

library(labdsv)
p = 0.1
iva <- indval(fam.fin, groups)
# Table of the significant indicator families
gr <- iva$maxcls[iva$pval<=p]
iv <- iva$indcls[iva$pval<=p]
pv <- iva$pval[iva$pval<=p]
fr <- apply(fam.fin>0,2,sum)[iva$pval<=p]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Indicator Species analysis using package indispecies
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

require(indicspecies)

indval <- multipatt(fam.hel, as.numeric(env.x$age_class), control = how(nperm=999))
summary(indval)



data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

fam.hel <- decostand(fam.usc, "hel") 

indval <- multipatt(fam.rel, as.numeric(env1$age_class), #duleg=TRUE,
                    control = how(within = Within(type = "series", mirror = FALSE), 
                              plots = Plots(strata = env1$field.ID, type = "free"), 
                              #blocks = env1$samcam, 
                              nperm=999))
summary(indval, alpha=0.05, indvalcomp=TRUE)
indval$sign


fam.pa <- decostand(fam.usc, "pa")

phi <- multipatt(fam.pa, as.numeric(env1$age_class), func="r.g", 
                    control = how(within = Within(type = "series", mirror = FALSE),
                                  plots = Plots(strata = env1$field.ID, type = "free"), 
                                  #blocks = env1$samcam, 
                                  nperm=999))
summary(phi, alpha=1, indvalcomp=TRUE)
round(head(phi$str),3)


indval <- multipatt(fam.usc, as.numeric(env1$age_class), control = how(nperm=999))
summary(indval, indvalcomp=TRUE)

phi <- multipatt(fam.pa, as.numeric(env1$age_class), func="r.g", control = how(nperm=999))
summary(phi, alpha=1, indvalcomp=TRUE)
round(head(phi$str),3)


fam.comb = combinespecies(fam.usc, max.order = 2)$XC
dim(fam.comb)

indval <- multipatt(fam.comb, as.numeric(env1$age_class), #duleg=TRUE,
                    control = how(within = Within(type = "series", mirror = FALSE), 
                                  plots = Plots(strata = env1$field.ID, type = "free"), 
                                  #blocks = env1$samcam, 
                                  nperm=999))
summary(indval, alpha=0.05, indvalcomp=TRUE)

# A and B gives us additional information about why species can be used as in-
# dicators.  For example, Pratzlenchiae is a good indicator of Group 1 because it
# occurs almost exclusivelz in sites belonging to this group (i.e., A = 0.9575), although not
# all  sites  belonging  to  Group  1  include  the  species  (i.e.,  B  =  0.8333).
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

