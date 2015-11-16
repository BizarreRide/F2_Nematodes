###########################
# F2 Nematodes
# Distribution check of Nematode Indices
# Quentin Schorpp
# 10.11.2015
###########################


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 
source("Data/DataProcessing/AverageData.R") 

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

# as needed
# source("Data/DataProcessing/MaturityIndices.R") 
# source("Data/DataProcessing/FaunalProfileIndices.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant

env.av$agsam <- with(env.av, interaction(age_class,samcam))
#env.av$c <-  env1$c # include carbon

#fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)

indices <- data.frame(age_class=env.av$age_class, 
                      field.ID=env.av$field.ID, 
                      crop=env.av$crop, 
                      trials=round(rowSums(fam.av),0), 
                      n=env.av$n, 
                      mc = env.av$mc, 
                      pH = env.av$pH, 
                      ata1 = env.av$ata1)
indices$ID <- 1:nrow(indices)

indices.backup <- indices
indices <- indices.backup
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis Nematode Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.av
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 
biodiv <- biodiv.fun(round(fam.av*counts.av$counts,0))


nema.backup <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv[,-7], N=counts.av$counts)
nema <- round(nema.backup,2)


p <- ncol(nema)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
  fk <-  fligner.test(nema[,i] ~ indices$crop)
  pvector[i] <- fk$p.value
  boxplot(nema[,i]~indices$crop)
}
any(pvector<0.05) # Here is sth. wrong!!
pvector

# There is no violation of heterosecedasticity for any of the indices
# Null Hypothesis: All variances are equal can't be abgelehnt at the 5% level
# The probability to judge real diffeences in varainces although the NULL is true is too high
# in more than 5 of a hundred cases we would judge wrong if we say there are differences
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data summary ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- nema[,i]
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(nema)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row.names(nema) <- 1:nrow(nema)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- nema[,i]
  car::Boxplot(indices$y ~ indices$crop)
  title(names(nema)[i],outer=TRUE)
}

# outliers:
# fungivores: 2
# omnivores: 1


# age_class:
# bacterivores decrease
# herbivores iindiease
# fungivores have an polynomial relationship

# samcam:
# carnivores iindiease in the second year, only in Silphie!
# no big differences for the others
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# drop Outliers ####
outlier <- list(nema.BI <- 1:18,
                nema.SI <- 1:18,
                nema.EI <- 1:18,
                nema.CI <- 1:18,
                nema.MI <- 1:18,
                nema.PPI <- -6,
                nema.MI25 <- 1:18,
                nema.sigmaMI <- 1:18,
                nema.sigmaMI25 <- 1:18,
                nema.PPI1 <- -6,
                nema.SR <- 1:18,
                nema.rarefy <- 1:18,
                nema.H <- 1:18,
                nema.D <- 1:18,
                nema.J <- 1:18,
                nema.H1 <- 1:18,
                nema.N <- -15)


# change factor properties
indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2

str(indices)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Distribution plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(MASS)
library(car)


for (i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  y <- nema[outlier[[i]],i] + 1
  qqp(y, "norm")
  qqp(y, "lnorm")
  
  #nbinom <- fitdistr(y, "negative binomial")
  #qqp(y, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  
  poisson <- fitdistr(y, "Poisson")
  qqp(y, "pois", poisson$estimate)
  
  #gamma <- fitdistr(y, "gamma")
  #qqp(y, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
  
  
  title(names(nema)[i],outer=TRUE)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# normal LMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p <- ncol(nema)


fp.indi.lmer.crop <- matrix(NA,2,2+(2*p))
colnames(fp.indi.lmer.crop) <- c("Env", "DF", rep(colnames(nema)[1:p], each=2))
fp.indi.lmer.crop[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))

indi.lmer.crop <- list()



for(i in 1:p) {
  indices2 <- indices[outlier[[17]],]
  indices2$y <- nema[outlier[[17]],17]
  model <- lm(y ~ age_class , indices2)
  model <- lmer(y ~ age_class + (1|crop) , indices2)
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("indi",i,names(nema)[i], sep = ".")
  assign(name, model)
  indi.lmer.crop[[i]] <- assign(name, model)
  fp.indi.lmer.crop[2,2+(i*2)] <- round(car::Anova(model, type="II")$"Chisq",2)
  fp.indi.lmer.crop[2,2+((i*2)-1)] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
fp.indi.lmer.crop[2,1]  <- row.names(Anova(model))
fp.indi.lmer.crop[2,2]  <- Anova(model)$"Df"


mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(nema)[i], sep = "."))}
names(indi.lmer.crop)[1:p] <- mod.names


indi.rsquared <- matrix(NA,2,2+2*p)
indi.rsquared[1:2,1] <- c("R2m", "R2c")

for(i in 1:p) {
  indi.rsquared[,2+2*i] <- round(MuMIn::r.squaredGLMM(indi.biglmer.crop[[i]]),2)
}
colnames(indi.rsquared) <- c("X", "X", rep(colnames(nema)[1:p],each=2))

fpr2.indi.lmer.crop <- rbind(fp.indi.lmer.crop, indi.rsquared, c("X", "X", rep("binomial", 2*p)))

# save(list=c("f.indi.lmer.crop","p.indi.lmer.crop", "r2.indi.lmer.crop"), file="Results/ANOVATables/CHi2+p_indi_LM_crop.rda")
# write.csv(fpr2.indi.lmer.crop, file="Results/ANOVATables/fpr2_indi_LM_crop.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detach("package:lmerTest", unload=TRUE)
library(lsmeans)

for( i in 1:p) {
lsm <- lsmeans(indi.lmer.crop[[i]], "age_class")
print(lsm)
lsm2 <- contrast(lsm, "trt.vs.ctrl", ref=c(2:5))
print(lsm2)
lsm2 <- contrast(lsm, "trt.vs.ctrl", ref=c(1))
print(lsm2)
Boxplot(y ~ age_class, indices2)
}

### normal LMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(indi.lmer.crop[[k]]),Anova(indi.lmer.crop[[k]], type="III")))
  #corvif(indi.lmer.crop[[k]])
  
  E1 <- resid(indi.lmer.crop[[k]], type="pearson")
  E2 <- resid(indi.lmer.crop[[k]], type="response")
  F1 <- fitted(indi.lmer.crop[[k]], type="response")
  P1 <- predict(indi.lmer.crop[[k]], type="response")
  
  par(mfrow=c(2,2),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  # Plot fitted vs. residuals
  scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(P1, E1, cex.lab = 1.5, xlab="Predicted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  #scatter.smooth(F1, P1, cex.lab = 1.5, xlab="Fitted values", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  # Histogram of Residuals
  hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
  lines(density(E1), col="light blue", lwd=3)
  lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
  
  title(names(indi.lmer.crop)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$samcam[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$agsam[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(indi.lmer.crop)[k], outer=TRUE)
  
  indices$y <- nema[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(indi.lmer.crop)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(indi.lmer.crop[[k]])
  
  E1 <- resid(indi.lmer.crop[[k]], type="pearson")
  E2 <- resid(indi.lmer.crop[[k]], type="response")
  F1 <- fitted(indi.lmer.crop[[k]], type="response")
  P1 <- predict(indi.lmer.crop[[k]], type="response")
  
  par(mfrow=c(2,2),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  
  # Plot fitted vs. residuals
  scatter.smooth(indices$n[outlier[[k]]], E1, cex.lab = 1.5, xlab="Nitrogen", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(indices$pH[outlier[[k]]], E1, cex.lab = 1.5, xlab="pH", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  scatter.smooth(indices$mc[outlier[[k]]], E1, cex.lab = 1.5, xlab="mc", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  # plot fitted vs. predicted
  scatter.smooth(indices$ata1[outlier[[k]]], E1, cex.lab = 1.5, xlab="Temperature", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  title(names(indi.lmer.crop)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(indi.lmer.crop[[k]], obs=TRUE)
}

for(k in 1:p) {
  plot1 <- plot(infl[[k]], which = "cook", main=mod.names[k])
  print(plot1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#####

# END





