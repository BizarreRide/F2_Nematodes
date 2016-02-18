#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Nematodes
# GLMM with Dredge
# Quentin Schorpp
# 12.02.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

asinTransform <- function(p) { asin(sqrt(p)) }

dispersion_glmer<- function(modelglmer)
{   
  # computing  estimated scale  ( binomial model)
  #following  D. Bates :
  #That quantity is the square root of the penalized residual sum of
  #squares divided by n, the number of observations, evaluated as:
  
  n <- length(modelglmer@resid)
  
  return(  sqrt( sum(c(modelglmer@resid, modelglmer@u) ^2) / n ) )
}

library(influence.ME)
library(lmerTest)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Data procsessing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant

env1$agsam <- with(env1, interaction(age_class,samcam))
#env.fin$c <-  env1$c # include carbon

#fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)

indices <- data.frame(ID = 1:nrow(env1),
                      age_class=env1$age_class,
                      samcam = env1$samcam,
                      nsamcam = as.numeric(factor(env1$samcam)),
                      agsam = env1$agsam,
                      location=env1$location,
                      field.ID = env1$field.ID,
                      cnratio = env1$cn,
                      mc = env1$mc,
                      pH = env1$pH,
                      ats1 = env1$ats1,
                      N=rowSums(fam.org))

indices$nsamcam <- as.numeric(factor(indices$samcam))

indices.backup <- indices
indices <- droplevels(indices.backup[!indices.backup$age_class %in% "A_Cm",])

explanatory <- c("age_class", "samcam", "age_class:samcam", "pH", "mc", "cnratio", "ats1") # include "intercept" when using Anova type III
q <- length(explanatory)#+1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis Nematode Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.org
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 
biodiv <- biodiv.fun(round(fam.usc,0))

df.response2 <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv)
df.response1 <- df.response2[!indices.backup$age_class %in% "A_Cm",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


df.response3 <- fety[!indices.backup$age_class %in% "A_Cm",-c(1:4)]
df.response4 <- df.response3/indices$N     # Percentage data
df.response5 <- round(df.response4*counts[16:45,][!indices.backup$age_class %in% "A_Cm", "counts"],0)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of selected Taxa ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.response6 <- round(fam.usc[,c("Tylenchidae", "Aphelenchidae", "Hoplolaimidae", "Cephalobidae", "Plectidae", "Telotylenchidae", "Rhabditidae", "Aporcelaimidae", "Aphelenchoididae", "Panagrolaimidae")],0)
df.response7 <- df.response6[!indices.backup$age_class%in% "A_Cm",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.responseX <- cbind(df.response1, df.response5, df.response7)

p <- ncol(df.responseX)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row.names(df.responseX) <- 1:nrow(df.responseX)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- df.responseX[,i]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$samcam)
  car::Boxplot(indices$y ~ indices$agsam)
  title(names(df.responseX)[i],outer=TRUE)
}

# outliers:
# fungivores: 2
# omnivores: 1


# age_class:
# bacterivores decrease
# herbivores increase
# fungivores have an polynomial relationship

# samcam:
# carnivores increase in the second year, only in Silphie!
# no big differences for the others
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# drop Outliers ####
outlier <- list(nema.BI <- -c(24,9),
                nema.SI <- 1:24,
                nema.EI <- 1:24,
                nema.CI <- -24,
                nema.MI <- -17,
                nema.PPI <- -c(6,19),
                nema.MI25 <- 1:24,
                nema.sigmaMI <- -c(3,17),
                nema.sigmaMI25 <- 1:24,
                nema.PPI1 <- 1:24,
                nema.SR <- 1:24,
                nema.rarefy <- 1:24,
                nema.H <- -23,
                nema.D <- -23,
                nema.J <- -23,
                nema.H1 <- 1:24,
                nema.N <- 1:24,
                fety.fungi <- -9,
                fety.bacti <- -20,
                fety.carni <- -20,
                fety.omni <- -c(22,17),
                fety.Tyli <- -13,
                fety.herbi <- -19,
                fety.herbi2 <- -21,
                fety.fungi2 <- 1:24,
                spec.Tyli <- -5,
                spec.Aph <- -2,
                spec.Hop <- -c(23,13),
                spec.Cph <- 1:24,
                spec.Plec <- -11,
                spec.Telo <- -c(24,9),
                spec.Rha <- -12,
                spec.Apc <- -c(22,20),
                spec.Aphdd <- -10,
                spec.Pan <- -8)



# change factor properties
indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2

indices$samcam2 <- indices$samcam
indices$samcam <- indices$nsamcam
indices$samcam <- indices$samcam2

str(indices)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Global Models
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# for (i in 1:ncol(df.responseX)) {
#   v1[i] <- c(paste("glbM",i,abbreviate(colnames(df.responseX),3)[i], " <- ", sep="."))
# }

ls.glbModels <- list()

fml.glb <- as.formula(y ~ age_class + samcam + age_class:samcam + (1|field.ID)) # + pH + cnratio + mc + ats1 
fml.glb.log <- as.formula(log1p(y) ~ age_class + samcam + age_class:samcam + (1|field.ID)) #+ pH + cnratio + mc + ats1 

con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = TRUE, check.conv.grad="ignore")


df.families <- read.delim("Data/F2_Nema_GLMMFamilies.txt", header=FALSE)  

for ( i in 1:ncol(df.responseX)) {
  indices$y <- df.responseX[,i]
  indices2 <- indices[outlier[[i]],]
  if(df.families[i,2] == "normal")
    G1 <- lmer(fml.glb, indices2)
  if(df.families[i,2] == "lognormal")
    G1 <- lmer(fml.glb.log, indices2)
  if(df.families[i,2] == "poisson")
    G1 <- glmer(fml.glb, family = "poisson", indices2, control = con)
  name <- c(paste("glbM",i,abbreviate(colnames(df.responseX),3)[i], sep="."))
  ls.glbModels[[i]] <- assign(name, G1)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Dredge Models # demo(dredge.subset) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate Cluster ####
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
#clusterExport(clust, varlist=c("dt.exp", "con"))

# Subsets of models excluded from dredge: ####
opo <- indices[,c("mc","pH","cnratio","ats1")]
opo <- as.data.frame(opo)


options(na.action = na.fail)

library(MuMIn)
options(na.action = na.fail)

# Dredge Abundance ####
p <- ncol(df.responseX)
ls.dredge <- list()
for(i in 1:p) {
  indices$y <- df.responseX[,i]
  indices2 <- indices[outlier[[i]],]
#   dt.exp$y <- dt.rsp.abn[,i, with=F]
  clusterExport(clust, varlist=c("indices2", "con"))
  GM.dredge <- pdredge(ls.glbModels[[i]], cluster=clust)
  #fixed=c("age_class", "samcam"),
  name <- c(paste("Dredge",i,abbreviate(colnames(df.responseX),3)[i], sep="."))
  assign(name, GM.dredge)
  ls.dredge[[i]] <- assign(name, GM.dredge)
}


# 1.a get the best models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls.bestmodels <- list()

for ( i in 1:p) {
  indices$y <- df.responseX[,i]
  indices2 <- indices[outlier[[i]],]
  #dt.exp$y <- dt.rsp[,i, with=F]
  M.best <- get.models(ls.dredge[[i]], 1)[[1]]
  name <- c(paste("BM",i,abbreviate(colnames(df.responseX),3)[i], sep="."))
  assign(name, M.best)
  ls.bestmodels[[i]] <- assign(name, M.best)
  names(ls.bestmodels)[[i]] <- name
}


for ( i in 1:p) {
  print(formula(ls.bestmodels[[i]]))
}


for ( i in 1:p) {
  indices$y <- df.responseX[,i]
  indices2 <- indices[outlier[[i]],]
  print(drop1(ls.glbModels[[i]]))
}

step(ls.glbModels[[1]])

df.compM2 <- vector()

for(i in 1:p){
  delta4 <- subset(ls.dredge[[i]], delta < 4)
  df.compM <- data.frame(delta4)
  df.compM[,ncol(df.compM)+1] <- rep(colnames(df.responseX)[i], length(rownames(df.compM)))
  df.compM2 <- rbind(df.compM2, df.compM)
}

