#§§§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Many alternative and wrong tests to figure out difference in p values
# Without Maize
# Quentin Schorpp
# 11.11.2015
#§§§§§§§§§§§§§§§§§§§§§§§§§§



# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

asinTransform <- function(p) { asin(sqrt(p)) }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# as needed
# source("Data/DataProcessing/MaturityIndices.R") 
# source("Data/DataProcessing/FaunalProfileIndices.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant

env.fin$agsam <- with(env.fin, interaction(age_class,samcam))
#env.fin$c <-  env1$c # include carbon

#fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)

indices <- droplevels(cbind(groups,env.fin, location=env1$location, N=rowSums(fam.org)))
indices$ID <- 1:nrow(indices)
indices.backup <- indices
indices <- indices.backup
indices$nsamcam <- as.numeric(factor(indices$samcam))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


fety.backup <- fety/indices$N
fety.backup <- cbind(fety.backup, group=indices$age_class, time=indices$samcam, subject=indices$field.ID, crop = indices$crop)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#fety3 <- fety.backup[,-c(1:4)]
fety3 <- fety.backup[!fety.backup$group %in% "A_Cm",-c(1:4)]


# PERMANOVA no strata ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fety.adis0 <- list()
for (i in 1:8) {
  fety3$response <- fety3[,i]
  adi1 <- adonis(response ~ group*time, method="euclidean", data=fety3, perm=999)
  name <- paste("adifety",i,names(fety3)[i], sep = ".")
  assign(name, adi1)
  fety.adis0[[i]] <- assign(name, adi1)
}

p.fety.adis0 <- matrix(NA,3,10)
colnames(p.fety.adis0) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.adis0 <- p.fety.adis0

for(i in 1:8) {
  f.fety.adis0[,i+2] <- round(fety.adis0[[i]]$aov.tab$"F.Model"[1:3],2)
  p.fety.adis0[,i+2] <- round(fety.adis0[[i]]$aov.tab$`Pr(>F)`[1:3],3)
}
f.fety.adis0[,1] <- p.fety.adis0[,1] <- c("group", "time", "group:time")
f.fety.adis0[,2] <- p.fety.adis0[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PERMANOVA strata subject ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fety.adis <- list()
for (i in 1:8) {
  fety3$response <- fety3[,i]
  adi1 <- adonis(response ~ group*time, method="euclidean", data=fety3, strata=fety3$subject, perm=999)
  name <- paste("adifety",i,names(fety3)[i], sep = ".")
  assign(name, adi1)
  fety.adis[[i]] <- assign(name, adi1)
}

p.fety.adis <- matrix(NA,3,10)
colnames(p.fety.adis) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.adis <- p.fety.adis

for(i in 1:8) {
  f.fety.adis[,i+2] <- round(fety.adis[[i]]$aov.tab$"F.Model"[1:3],2)
  p.fety.adis[,i+2] <- round(fety.adis[[i]]$aov.tab$`Pr(>F)`[1:3],3)
}
f.fety.adis[,1] <- p.fety.adis[,1] <- c("group", "time", "group:time")
f.fety.adis[,2] <- p.fety.adis[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PERMANOVA strata time####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fety.adissc <- list()
for (i in 1:8) {
  fety3$response <- fety3[,i]
  adi1 <- adonis(response ~ group*time, method="euclidean", data=fety3, strata=fety3$time, perm=999)
  name <- paste("adifety",i,names(fety3)[i], sep = ".")
  assign(name, adi1)
  fety.adissc[[i]] <- assign(name, adi1)
}

p.fety.adissc <- matrix(NA,3,10)
colnames(p.fety.adissc) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.adissc <- p.fety.adissc

for(i in 1:8) {
  f.fety.adissc[,i+2] <- round(fety.adissc[[i]]$aov.tab$"F.Model"[1:3],2)
  p.fety.adissc[,i+2] <- round(fety.adissc[[i]]$aov.tab$`Pr(>F)`[1:3],3)
}
f.fety.adissc[,1] <- p.fety.adissc[,1] <- c("group", "time", "group:time")
f.fety.adissc[,2] <- p.fety.adissc[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Simple ANOVA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.aov <- matrix(NA,3,10)
colnames(p.fety.aov) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.aov <- p.fety.aov

fety.aov <- list()

for (i in 1:8) {
  fety3$response <- fety3[,i]
  A1 <- aov(response ~ group*time, fety3)
  name <- paste("aovfety",i,names(fety3)[i], sep = ".")
  assign(name, A1)
  fety.aov[[i]] <- assign(name, A1)
}

for(i in 1:8) {
  f.fety.aov[,i+2] <- round(summary(fety.aov[[i]])[[1]][["F value"]][1:3],2)
  p.fety.aov[,i+2] <- round(summary(fety.aov[[i]])[[1]][["Pr(>F)"]][1:3],3)
}
f.fety.aov[,1] <- p.fety.aov[,1] <- c("group", "time", "group:time")
f.fety.aov[,2] <- p.fety.aov[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Simple ANOVA with arcsine transformation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.aovas <- matrix(NA,3,10)
colnames(p.fety.aovas) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.aovas <- p.fety.aovas

fety.aovas <- list()

for (i in 1:8) {
  fety3$response <- asinTransform(fety3[,i])
  A1 <- aov(response ~ group*time, fety3)
  name <- paste("aovfety",i,names(fety3)[i], sep = ".")
  assign(name, A1)
  fety.aovas[[i]] <- assign(name, A1)
}

for(i in 1:8) {
  f.fety.aovas[,i+2] <- round(summary(fety.aovas[[i]])[[1]][["F value"]][1:3],2)
  p.fety.aovas[,i+2] <- round(summary(fety.aovas[[i]])[[1]][["Pr(>F)"]][1:3],3)
}
f.fety.aovas[,1] <- p.fety.aovas[,1] <- c("group", "time", "group:time")
f.fety.aovas[,2] <- p.fety.aovas[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Repeated measures ANOVA Long format ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.rmaov <- matrix(NA,3,10)
colnames(p.fety.rmaov) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.rmaov <- p.fety.rmaov

fety.rmaov <- list()

for (i in 1:8) {
  fety3$response <- fety3[,i]
  rmA1 <- aov(response ~ group*time + Error(subject), fety3)
  name <- paste("aovfety",i,names(fety3)[i], sep = ".")
  assign(name, rmA1)
  fety.rmaov[[i]] <- assign(name, rmA1)
}

for(i in 1:8) {
  f.fety.rmaov[,i+2] <- round(c(summary(fety.rmaov[[i]])[[1]][[1]][1,"F value"], summary(fety.rmaov[[7]])[[2]][[1]][1:2,"F value"]),2)
  p.fety.rmaov[,i+2] <- round(c(summary(fety.rmaov[[i]])[[1]][[1]][1,"Pr(>F)"], summary(fety.rmaov[[7]])[[2]][[1]][1:2,"Pr(>F)"]),2)
}
f.fety.rmaov[,1] <- p.fety.rmaov[,1] <- c("group", "time", "group:time")
f.fety.rmaov[,2] <- p.fety.rmaov[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Repeated measures ANOVA with arcsine Transformation Long format ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.rmaovas <- matrix(NA,3,10)
colnames(p.fety.rmaovas) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.rmaovas <- p.fety.rmaovas

fety.rmaovas <- list()

for (i in 1:8) {
  fety3$response <- fety3[,i]
  rmA1 <- aov(response ~ group*time + Error(subject), fety3)
  name <- paste("aovfety",i,names(fety3)[i], sep = ".")
  assign(name, rmA1)
  fety.rmaovas[[i]] <- assign(name, rmA1)
}

for(i in 1:8) {
  f.fety.rmaovas[,i+2] <- round(c(summary(fety.rmaovas[[i]])[[1]][[1]][1,"F value"], summary(fety.rmaovas[[7]])[[2]][[1]][1:2,"F value"]),2)
  p.fety.rmaovas[,i+2] <- round(c(summary(fety.rmaovas[[i]])[[1]][[1]][1,"Pr(>F)"], summary(fety.rmaovas[[7]])[[2]][[1]][1:2,"Pr(>F)"]),2)
}
f.fety.rmaovas[,1] <- p.fety.rmaovas[,1] <- c("group", "time", "group:time")
f.fety.rmaovas[,2] <- p.fety.rmaovas[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Non parametric Test with repeated measurements ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


library(nparLD)
fety3 <- fety.backup[!fety.backup$group %in% "A_Cm",-c(1:4)]
p <- ncol(fety3)

#fety3$group <- as.ordered(fety3$group)



p.fety.nparLD <- matrix(NA,3,10)
colnames(p.fety.nparLD) <- c("Env", "DF", colnames(fety3)[1:8])
f.fety.nparLD <- p.fety.nparLD

fety.npmLDs <- list()

for (i in 1:8) {
  fety3$response <- fety3[,i]
  npM1 <- nparLD(response ~ group*time,  data = fety3,subject = "subject", description = FALSE)
  par(xpd=TRUE)
  plot(npM1)
  summary(npM1)
  name <- paste("npfety",i,names(fety3)[i], sep = ".")
  assign(name, npM1)
  fety.npmLDs[[i]] <- assign(name, npM1)
  f.fety.nparLD[,i+2] <- round(fety.npmLDs[[i]]$Wald.test[1:3,1],2)
  p.fety.nparLD[,i+2] <- round(fety.npmLDs[[i]]$Wald.test[1:3,3],2)
}
f.fety.nparLD[,1] <- p.fety.nparLD[,1] <- c("group", "time", "group:time")
f.fety.nparLD[,2] <- p.fety.nparLD[,2] <- c(3,1,3)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For bacterivores and Tylenchidae the mean RankMeans in both samcam's are the same, and factor time is NA


save(list=c("p.fety.aov" ,"p.fety.aovas" ,"p.fety.rmaov" ,"p.fety.rmaovas", "p.fety.nparLD"), file="pValuesNorm.rda")

