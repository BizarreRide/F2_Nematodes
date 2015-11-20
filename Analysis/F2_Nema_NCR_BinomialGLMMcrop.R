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
ncr <- as.data.frame(fety[,-1])

df.response <- round(ncr[,c(5,6)],0)

df.backup <- df.response/rowSums(df.response)

p <- ncol(df.response)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(1:p)
for (i in 1:p) {
  fk <-  fligner.test(df.response[,i] ~ indices$crop)
  pvector[i] <- fk$p.value
  boxplot(df.response[,i]~indices$crop)
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
  indices$y <- df.response[,i]
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(df.response)[i],outer=TRUE)
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
  indices$y <- df.response[,i]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$crop)
  title(names(df.response)[i],outer=TRUE)
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
outlier <- list(ncr.fungi <- -c(15,9),
                ncr.bacti <- -c(15,9))


indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2


str(indices)

# Bi☺nomial GLMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# F and p-value squared *****************************************************************


df.Fpvalue <- matrix(NA,2,2+(2*p))
colnames(df.Fpvalue) <- c("Env", "DF", rep(colnames(df.response)[1:p], each=2))
df.Fpvalue[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))

ls.models <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$scs <- df.response[outlier[[i]],i]
  indices2$fail <- rowSums(df.response[outlier[[i]],]) - df.response[outlier[[i]],i]
  model <- glmer(cbind(scs,fail) ~ age_class + (1|ID), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  #model <- glmer(cbind(scs,fail) ~ crop - 1   + (1|ID) + (1|age_class), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("ncr",i,names(df.response)[i], sep = ".")
  assign(name, model)
  ls.models[[i]] <- assign(name, model)
  df.Fpvalue[2,2+(i*2)] <- round(car::Anova(model, type="II")$"Chisq",2)
  df.Fpvalue[2,2+((i*2)-1)] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
df.Fpvalue[2,1]  <- row.names(Anova(model))
df.Fpvalue[2,2]  <- Anova(model)$"Df"

mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("ncr",i,names(df.response)[i], sep = "."))}
names(ls.models)[1:p] <- mod.names


# R² - R squared ***********************************************************************

df.rsquared <- matrix(NA,2,2+2*p)
df.rsquared[1:2,1] <- c("R2m", "R2c")

for(i in 1:p) {
  df.rsquared[,2+2*i] <- round(MuMIn::r.squaredGLMM(ls.models[[i]]),2)
}
colnames(df.rsquared) <- c("X", "X", rep(colnames(nema)[1:p],each=2))


# Summary *******************************************************************************

df.FpvalueR2 <- rbind(df.Fpvalue, df.rsquared, c("X", "X", rep("binomial", 2*p)))

save("df.FpvalueR2", file="Results/ANOVATables/FpR2_ncr_bnGLM_crop.rda")
write.csv(df.FpvalueR2, file="Results/ANOVATables/FpR2_ncr_bnGLM_crop.csv")

# p-values with afex ********************************************************************
df.FpvalueR2.1 <- df.FpvalueR2 
df.FpvalueR2.1[2,] <- "N"

for(i in 1:p){
  indices2 <- indices[outlier[[i]],]
  indices2$scs <- df.response[outlier[[i]],i]
  indices2$fail <- rowSums(df.response[outlier[[i]],]) - df.response[outlier[[i]],i]
  obj.afex <- afex::mixed(cbind(scs,fail) ~ age_class + (1|ID), family=binomial(link="logit"), indices2, method="LRT") 
  df.FpvalueR2.1[2,2+(i*2)-1] <- round(obj.afex[[1]]$"Chisq",2)
  df.FpvalueR2.1[2,2+((i*2))] <- round(obj.afex[[1]]$"Pr(>Chisq)",3)
}

write.csv(df.FpvalueR2.1, file="Results/ANOVATables/FpR2afex_ncr_bnGLM_crop.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Post Hoc ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detach("package:piecewiseSEM", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:afex", unload=TRUE)
detach("package:lsmeans", unload=TRUE)
library(lsmeans)

df.posthoc <- matrix(NA,1,1+2*p)
df.posthoc2 <- matrix(NA,4,1+2*p)

for(i in 1:p) {
  lsm <- lsmeans::lsmeans(ls.models[[i]], "age_class")
  lsm2 <- contrast(lsm, "trt.vs.ctrl", ref=c(2:5))
  print(lsm2)
  df.posthoc[,(2*i)] <- round(summary(lsm2)$estimate,3)
  df.posthoc[,(2*i)+1] <- round(summary(lsm2)$p.value,3)
  lsm2 <- contrast(regrid(lsm), "trt.vs.ctrl", ref=c(2:5)) # Why do the p-values differ if i use regrid() ???
  df.posthoc[,(2*i)] <- round(summary(lsm2)$estimate,3)
  lsm3 <- contrast(lsm, "trt.vs.ctrl", ref=c(1))
  print(lsm3)
  df.posthoc2[,(2*i)] <- round(summary(lsm3)$estimate,3)
  df.posthoc2[,(2*i)+1] <- round(summary(lsm3)$p.value,3)
  lsm3 <- contrast(regrid(lsm), "trt.vs.ctrl", ref=c(1))
  df.posthoc2[,(2*i)] <- round(summary(lsm3)$estimate,3)
  
  indices2 <- indices[outlier[[i]],]
  indices2$scs <- df.response[outlier[[i]],i]
  indices2$fail <- rowSums(df.response[outlier[[i]],]) - df.response[outlier[[i]],i]
  Boxplot(scs ~ age_class, indices2)
}

df.posthoc[,1] <- paste(summary(lsm2)$"contrast")
df.posthoc2[,1] <- paste(summary(lsm3)$"contrast")

colnames(df.posthoc) <- c("Contrast", rep(c("estimate", "p-value"), p))
colnames(df.posthoc2) <- c("Contrast", rep(c("estimate", "p-value"), p))


write.csv(df.posthoc, file="Results/ANOVATables/PostHocC_ncr_bnGLM_crop.csv")
write.csv(df.posthoc2, file="Results/ANOVATables/PostHocAC_ncr_bnGLM_crop.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




### binomial GLMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(ls.models[[k]]),Anova(ls.models[[k]], type="III")))
  #corvif(ls.models[[k]])
  
  E1 <- resid(ls.models[[k]], type="pearson")
  E2 <- resid(ls.models[[k]], type="response")
  F1 <- fitted(ls.models[[k]], type="response")
  P1 <- predict(ls.models[[k]], type="response")
  
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
  
  title(names(ls.models)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$crop[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(ls.models)[k], outer=TRUE)
  
  indices$y <- df.response[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(ls.models)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(ls.models[[k]])
  
  E1 <- resid(ls.models[[k]], type="pearson")
  E2 <- resid(ls.models[[k]], type="response")
  F1 <- fitted(ls.models[[k]], type="response")
  P1 <- predict(ls.models[[k]], type="response")
  
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
  
  title(names(ls.models)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(ls.models[[k]], obs=TRUE)
}

for(k in 1:p) {
  plot1 <- plot(infl[[k]], which = "cook", main=mod.names[k])
  print(plot1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#####

# END

