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


df.backup <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv[,-7], N=counts.av$counts)
df.response <- round(df.backup,2)


p <- ncol(df.response)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
  fk <-  fligner.test(df.response[,i] ~ indices$crop)
  pvector[i] <- fk$p.value
  boxplot(df.response[,i]~indices$crop)
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
row.names(df.response) <- 1:nrow(df.response)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- df.response[,i]
  car::Boxplot(indices$y ~ indices$crop)
  title(names(df.response)[i],outer=TRUE)
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
  y <- df.response[outlier[[i]],i] + 1
  qqp(y, "norm")
  qqp(y, "lnorm")
  
  #nbinom <- fitdistr(y, "negative binomial")
  #qqp(y, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  
  poisson <- fitdistr(y, "Poisson")
  qqp(y, "pois", poisson$estimate)
  
  #gamma <- fitdistr(y, "gamma")
  #qqp(y, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
  
  
  title(names(df.response)[i],outer=TRUE)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# normal LM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# F and p-value squared *****************************************************************

p <- ncol(df.response)

df.Fpvalue <- matrix(NA,2,2+(2*p))
colnames(df.Fpvalue) <- c("Env", "DF", rep(colnames(df.response)[1:p], each=2))
df.Fpvalue[1,] <- c("X", "X", rep(c("F", "p-value"),p))

ls.models <- list()

vf1 <- nlme::varIdent(form= ~ 1|age_class)

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$y <- df.response[outlier[[i]],i]
  #model <- lm(y ~ age_class , indices2)
  model <- nlme::gls(y ~ age_class , indices2, weights=vf1)
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("indi",i,names(df.response)[i], sep = ".")
  assign(name, model)
  ls.models[[i]] <- assign(name, model)
  #df.Fpvalue[2,2+((i*2)-1)] <- round(car::Anova(model, type="II")$"F value"[1],2)
  #df.Fpvalue[2,2+((i*2)-0)] <- round(car::Anova(model)$"Pr(>F)"[1],3)
  df.Fpvalue[2,2+((i*2)-1)] <- round(car::Anova(model, type="II")$"Chisq"[1],2)
  df.Fpvalue[2,2+((i*2)-0)] <- round(car::Anova(model)$"Pr(>Chisq)"[1],3)
}
df.Fpvalue[2,1]  <- row.names(Anova(model))[1]
df.Fpvalue[2,2]  <- Anova(model)$"Df"[1]


mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(df.response)[i], sep = "."))}
names(ls.models)[1:p] <- mod.names


# R² - R squared ***********************************************************************

df.rsquared <- matrix(NA,2,2+2*p)
df.rsquared[1:2,1] <- c("R2m", "R2c")

for(i in 1:p) {
  df.rsquared[,2+2*i] <- round(MuMIn::r.squaredGLMM(ls.models[[i]]),2)
}
colnames(df.rsquared) <- c("X", "X", rep(colnames(df.response)[1:p],each=2))


# Summary *******************************************************************************

df.FpvalueR2 <- rbind(df.Fpvalue, df.rsquared, c("X", "X", rep("normal", 2*p)))

# save("df.FpvalueR2", file="Results/ANOVATables/FpR2_indi_LM_crop.rda")
# write.csv(df.FpvalueR2, file="Results/ANOVATables/FpR2_indi_LM_crop.csv")

# p-values with afex ********************************************************************
df.FpvalueR2.1 <- df.FpvalueR2 
df.FpvalueR2.1[2,] <- "NA"

for(i in 1:p){
  indices2 <- indices[outlier[[i]],]
  indices2$y <- df.response[outlier[[i]],i]
  obj.afex <- afex::mixed(y ~ age_class , indices2,  method="LRT") 
  df.FpvalueR2.1[2,2+(i*2)-1] <- round(obj.afex[[1]]$"Chisq",2)
  df.FpvalueR2.1[2,2+((i*2))] <- round(obj.afex[[1]]$"Pr(>Chisq)",3)
}

# write.csv(df.FpvalueR2.1, file="Results/ANOVATables/FpR2afex_indi_LM_crop.csv")

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
  indices2$y <- df.response[outlier[[i]],i]
  Boxplot(y ~ age_class, indices2)
}

df.posthoc[,1] <- paste(summary(lsm2)$"contrast")
df.posthoc2[,1] <- paste(summary(lsm3)$"contrast")

colnames(df.posthoc) <- c("Contrast", rep(c("estimate", "p-value"), p))
colnames(df.posthoc2) <- c("Contrast", rep(c("estimate", "p-value"), p))

# write.csv(df.posthoc, file="Results/ANOVATables/PostHocC_indi_LM_crop.csv")
# write.csv(df.posthoc2, file="Results/ANOVATables/PostHocAC_indi_LM_crop.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detach("package:lmerTest", unload=TRUE)
library(lsmeans)

for( i in 1:p) {
lsm <- lsmeans(ls.models[[i]], "age_class")
print(lsm)
lsm2 <- contrast(lsm, "trt.vs.ctrl", ref=c(2:5))
print(lsm2)
lsm2 <- contrast(lsm, "trt.vs.ctrl", ref=c(1))
print(lsm2)
Boxplot(y ~ age_class, indices2)
}

### normal LM - Model Validation ####
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
  scatter.smooth(F1, E1, span=4/5, cex.lab = 1.5, xlab="Fitted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(P1, E1, span=4/5, cex.lab = 1.5, xlab="Predicted values", ylab=" Residuals")
  abline(h = 0, v=0, lty=2)
    
  # Histogram of Residuals
  hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
  lines(density(E1), col="light blue", lwd=3)
  lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
  
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  title(names(ls.models)[k], outer=TRUE)
  

  
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$crop[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  title(names(ls.models)[k], outer=TRUE)
  
  
  indices$y <- df.response[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,df.response[outlier[[k]],k], span=4/5, cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
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





