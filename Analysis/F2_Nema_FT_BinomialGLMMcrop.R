###########################
# F2 Nematodes
# Binomial GLMM
# Quentin Schorpp
# 11.08.2015
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
                      N=round(rowSums(fam.av),0), 
                      n=env.av$n, 
                      mc = env.av$mc, 
                      pH = env.av$pH, 
                      ata1 = env.av$ata1)
indices$ID <- 1:nrow(indices)

indices.backup <- indices
indices <- droplevels(indices.backup)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.av
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


fety.backup <- fety/indices$N
fety2 <- round(fety[,-c(1:4)],0)

p <- ncol(fety2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
  fk <-  fligner.test(fety[,i] ~ env.fin$agsam)
  pvector[i] <- fk$p.value
  boxplot(fety[,i]~env.fin$agsam)
}
any(pvector<0.05)
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
  indices$y <- asinTransform(fety2[,i])
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(fety2)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- fety.backup[,i+4]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$crop)
  title(names(fety2)[i],outer=TRUE)
}


for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- fety2[,i]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$crop)
  title(names(fety2)[i],outer=TRUE)
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
outlier <- list(fety.fungi <- -15,
                fety.bacti <- -8,
                fety.carni <- 1:18,
                fety.omni <- 1:18,
                fety.Tyli <- 1:18,
                fety.herbi <- 1:18,
                fety.herbi2 <- 1:18,
                fety.fungi2 <- 1:18)



indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2


str(indices)

# Bi☺nomial GLMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fp.fety.biglmer.crop <- matrix(NA,2,2+(2*p))
colnames(fp.fety.biglmer.crop) <- c("Env", "DF", rep(colnames(fety2)[1:p], each=2))
fp.fety.biglmer.crop[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))

fety.biglmer.crop <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$scs <- fety2[outlier[[i]],i]
  indices2$fail <- indices2[,"N"] - fety2[outlier[[i]],i]
  model <- glmer(cbind(scs,fail) ~ crop + (1|ID), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  #model <- glmer(cbind(scs,fail) ~ crop - 1   + (1|ID) + (1|age_class), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("fety",i,names(fety2)[i], sep = ".")
  assign(name, model)
  fety.biglmer.crop[[i]] <- assign(name, model)
  fp.fety.biglmer.crop[2,2+(i*2)] <- round(car::Anova(model, type="II")$"Chisq",2)
  fp.fety.biglmer.crop[2,2+((i*2)-1)] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
fp.fety.biglmer.crop[2,1]  <- row.names(Anova(model))
fp.fety.biglmer.crop[2,2]  <- Anova(model)$"Df"

mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(fety2)[i], sep = "."))}
names(fety.biglmer.crop)[1:p] <- mod.names

fety.rsquared <- matrix(NA,2,2+2*p)
fety.rsquared[1:2,1] <- c("R2m", "R2c")

for(i in 1:p) {
  fety.rsquared[,2+2*i] <- round(MuMIn::r.squaredGLMM(fety.biglmer.crop[[i]]),2)
}
colnames(fety.rsquared) <- c("X", "X", rep(colnames(fety2)[1:p],each=2))

FpR2.frame <- rbind(fp.fety.biglmer.crop, fety.rsquared, c("X", "X", rep("binomial", 2*p)))

overdisp_fun(model)
# 
# save(list=c("f.fety.biglmer.crop","p.fety.biglmer.crop", "r2.fety.biglmer.crop"), file="Results/ANOVATables/CHi2+p_Fety_bnGLMM_crop.rda")
# write.csv(f.fety.biglmer.crop, file="Results/ANOVATables/Chi2_Fety_bnGLMM_crop.csv")
# write.csv(p.fety.biglmer.crop, file="Results/ANOVATables/p_Fety_bnGLMM_crop.csv")
# write.csv(r2.fety.biglmer.crop, file="Results/ANOVATables/r2_Fety_bnGLMM_crop.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Overdispersion ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


### binomial GLMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(fety.biglmer.crop[[k]]),Anova(fety.biglmer.crop[[k]], type="III")))
  #corvif(fety.biglmer.crop[[k]])
  
  E1 <- resid(fety.biglmer.crop[[k]], type="pearson")
  E2 <- resid(fety.biglmer.crop[[k]], type="response")
  F1 <- fitted(fety.biglmer.crop[[k]], type="response")
  P1 <- predict(fety.biglmer.crop[[k]], type="response")
  
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
  
  title(names(fety.biglmer.crop)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$crop[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(fety.biglmer.crop)[k], outer=TRUE)
  
  indices$y <- fety2[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(fety.biglmer.crop)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(fety.biglmer.crop[[k]])
  
  E1 <- resid(fety.biglmer.crop[[k]], type="pearson")
  E2 <- resid(fety.biglmer.crop[[k]], type="response")
  F1 <- fitted(fety.biglmer.crop[[k]], type="response")
  P1 <- predict(fety.biglmer.crop[[k]], type="response")
  
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
  
  title(names(fety.biglmer.crop)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(fety.biglmer.crop[[k]], obs=TRUE)
}

for(k in 1:p) {
  plot1 <- plot(infl[[k]], which = "cook", main=mod.names[k])
  print(plot1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#####

# END

