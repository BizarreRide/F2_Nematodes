###########################
# F2 Nematodes
# Normal LMM for Nematode Indices
# Quentin Schorpp
# 12.11.2015
###########################


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

# as needed
# source("Data/DataProcessing/MaturityIndices.R") 
# source("Data/DataProcessing/FaunalProfileIndices.R") 
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
                      n = env1$n,
                      mc = env1$mc,
                      pH = env1$pH,
                      ata1 = env1$ata1,
                      N=rowSums(fam.org))

indices$nsamcam <- as.numeric(factor(indices$samcam))

indices.backup <- indices
indices <- droplevels(indices.backup[!indices.backup$age_class %in% "A_Cm",])

explanatory <- c("age_class", "samcam", "age_class:samcam") # include "intercept" when using Anova type III
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


p <- ncol(df.response1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
  fk <-  fligner.test(df.response1[,i] ~ indices$agsam)
  pvector[i] <- fk$p.value
  boxplot(df.response1[,i]~indices$agsam)
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
  indices$y <- df.response1[,i]
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(df.response1)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row.names(df.response1) <- 1:nrow(df.response1)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- df.response1[,i]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$samcam)
  car::Boxplot(indices$y ~ indices$agsam)
  title(names(df.response1)[i],outer=TRUE)
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
                nema.N <- 1:24)
                


# change factor properties
indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2

indices$samcam2 <- indices$samcam
indices$samcam <- indices$nsamcam
indices$samcam <- indices$samcam2

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
  y <- df.response1[outlier[[i]],i] + 1
  qqp(y, "norm")
  qqp(y, "lnorm")
  
  #nbinom <- fitdistr(y, "negative binomial")
  #qqp(y, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  
  poisson <- fitdistr(y, "Poisson")
  qqp(y, "pois", poisson$estimate)
  
  gamma <- fitdistr(y, "gamma")
  qqp(y, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
  
  
  title(names(df.response1)[i],outer=TRUE)
}

# nema.norm <- subset(nema, select=c("SI","EI","sigmaMI","sigmaMI25", "SR","rarefy","H","D","J","H1"))
df.response1 <- subset(df.response1, select=c("BI","CI", "PPI", "PPI.1", "MI", "MI25", "N"))
# nema.gamma <- subset(nema, select=c("MI", "MI25", "N"))
# 
# outlier.norm  <- list(nema.SI <- 1:24,
#                       nema.EI <- 1:24,
#                       nema.sigmaMI <- -c(3,17),
#                       nema.sigmaMI25 <- 1:24,
#                       nema.SR <- 1:24,
#                       nema.rarefy <- 1:24,
#                       nema.H <- 1:24,
#                       nema.D <- 1:24,
#                       nema.J <- 1:24,
#                       nema.H1 <- 1:24)
outlier <- list(nema.BI <- -c(24,9),
                      nema.CI <- -c(11,20,24),
                      nema.PPI <- -c(6,19),
                      nema.PPI1 <- 1:24,
                      nema.MI <- -17,
                      nema.MI25 <- 1:24,
                      nema.N <- 1:24)
p <- ncol(df.response1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#df.response1 <- df.response1[,-4]
#outlier <- outlier[-4]


# log-normal LMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df.Fpvalue <- matrix(NA,1+q,2+(2*p))
colnames(df.Fpvalue) <- c("Env", "DF", rep(colnames(df.response1)[1:p], each=2))
df.Fpvalue[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))

ls.models <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$y <- df.response1[outlier[[i]],i]
  model <- glmer(y ~ age_class*samcam +(1|field.ID), family=gaussian(link="log"), indices2)
  #model <- glmer(y ~ age_class*samcam  +  (1|field.ID), family=gaussian(link="log"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("indi",i,names(df.response1)[i], sep = ".")
  assign(name, model)
  ls.models[[i]] <- assign(name, model)
  df.Fpvalue[2:(1+q),2+((i*2)-1)] <- round(car::Anova(model, type="II")$"Chisq",2)
  df.Fpvalue[2:(1+q),2+(i*2)] <- round(car::Anova(model, type="II")$"Pr(>Chisq)",3)
}
df.Fpvalue[2:(1+q),1]  <- row.names(Anova(model, type="II"))
df.Fpvalue[2:(1+q),2]  <- Anova(model, type="II")$"Df"


mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(df.response1)[i], sep = "."))}
names(ls.models)[1:p] <- mod.names

df.rsquared <- matrix(NA,2,2+2*p)
df.rsquared[1:2,1] <- c("R2m", "R2c")

for(i in 1:p) {
  df.rsquared[,2+2*i] <- round(piecewiseSEM::sem.model.fits(ls.models[[i]]),2)
}
colnames(df.rsquared) <- c("X", "X", rep(colnames(df.response1)[1:p],each=2))

df.FpvalueR2 <- rbind(df.Fpvalue, df.rsquared, c("X", "X", rep("log-normal", 2*p)))

#write.csv(df.FpvalueR2, file="Results/ANOVATables/FpR2_Indi_logLMMT2.csv")
write.csv(df.FpvalueR2, file="Results/ANOVATables/FpR2_Indi_logLMMT3.csv")
write.table(df.FpvalueR2, file="Results/ANOVATables/FpR2_Indi_logLMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)


# p-values with afex ********************************************************************
df.FpvalueR2.1 <- df.FpvalueR2 
df.FpvalueR2[2:(1+q),] <- "NA"

for(i in 1:p){
  indices2 <- indices[outlier[[i]],]
  indices2$y <- df.response1[outlier[[i]],i]
  obj.afex <- afex::mixed(y ~ age_class*samcam  + (1|field.ID), family=gaussian(link="log"), indices2,  method="LRT") 
  df.FpvalueR2[2:(1+q),2+(i*2)-1] <- round(obj.afex[[1]]$"Chisq",2)
  df.FpvalueR2[2:(1+q),2+((i*2))] <- round(obj.afex[[1]]$"Pr(>Chisq)",3)
}

# write.csv(df.FpvalueR2.1, file="Results/ANOVATables/FpR2afex_Indi_logLMM.csv")
# write.csv(df.FpvalueR2.1, file="Results/ANOVATables/FpR2afex_Indi_logLMM2.csv")
# write.table(df.FpvalueR2.1, file="Results/ANOVATables/FpR2_Indi_logLMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Post Hoc data inspection with lsmeans package ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detach("package:piecewiseSEM", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:afex", unload=TRUE)
detach("package:lsmeans", unload=TRUE)
library(lsmeans)
library(multcompView)

ls.lsm<- list()

df.posthoc <- matrix(NA,8,2+(2*p))

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans::lsmeans(ls.models[[i]],  ~ age_class*samcam, contr= "cld")
  x <- cld(lsm, type = "response", sort=FALSE)
  df.posthoc[,2+((2*i)-1)] <- x$"response"
  df.posthoc[,2+((2*i)-0)] <- x$".group"
  print(x)
  name <- paste("lsm",i,names(df.response1)[i], sep = ".")
  ls.lsm[[i]] <- assign(name, lsm)
  # to see the results graphically
  p1 <- plot(lsm, by = "samcam", intervals = TRUE, type = "response")
  print(p1)
  #title(names(ncr.biglmer)[i], outer=TRUE)
}


for(i in 1:p){
  lsmip(ls.lsm[[i]], age_class ~ samcam, type = "response")
  summary(pairs(ls.lsm[[i]]), type = "response")
  summary(pairs(regrid(ls.lsm[[i]])), type = "response")
}

df.posthoc[,1] <- paste(x$"age_class")
df.posthoc[,2] <- paste(x$"samcam")

colnames(df.posthoc) <- c("Factor1", "Factor1", rep(colnames(df.response1),each=2))
df.posthoc <- rbind(c("age_class", "samcam", rep(c("response", "group"), p)), df.posthoc)


#************************************************************************


ls.lsmAC <- list()
df.posthocAC <- matrix(NA,12,2+(2*p))
colnames(df.posthocAC) <- c("contrast", "samcam", rep(colnames(df.response1),each=2))

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(ls.models[[i]],  ~ age_class|samcam, at=list(samcam=c("2","4")))
  x <- contrast(lsm, "pairwise" , type = "response")
  print(x)
  xx <- summary(x)
  df.posthocAC[,2+((2*i)-1)] <- round(xx$"estimate",2)
  df.posthocAC[,2+((2*i)-0)] <- round(xx$"p.value",3)
  name <- paste("lsm",i,names(df.response1)[i], sep = ".")
  ls.lsmAC[[i]] <- assign(name, lsm)
}

df.posthocAC[,1] <- paste(xx$"contrast")
df.posthocAC[,2] <- paste(xx$"samcam")


#************************************************************************


ls.lsmSC <- list()
df.posthocSC <- matrix(NA,4,2+(2*p))
colnames(df.posthocSC) <- c("contrast", "age_class", rep(colnames(df.response1),each=2))

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(ls.models[[i]], pairwise ~ samcam|age_class)
  x <- contrast(lsm, "pairwise" , type = "response")
  print(x)
  xx <- summary(x)
  df.posthocSC[,2+((2*i)-1)] <- round(xx$"estimate",2)
  df.posthocSC[,2+((2*i)-0)] <- round(xx$"p.value",3)
  name <- paste("lsm",i,names(df.response1)[i], sep = ".")
  ls.lsmSC[[i]] <- assign(name, lsm)
}

df.posthocSC[,1] <- paste(xx$"contrast")
df.posthocSC[,2] <- paste(xx$"age_class")

# save(list=c("ls.models","ls.lsm", "ls.lsmAC", "ls.lsmSC",  "df.FpvalueR2", "df.posthoc", "df.posthocAC", "df.posthocSC"), file="Results/ANOVATables/Indi_logLMM.rda")
# write.csv(df.posthoc, file="Results/ANOVATables/PostHoc_Indi_logLMM.csv")
# write.csv(df.posthocAC, file="Results/ANOVATables/PostHocAC_Indi_logLMM.csv")
# write.csv(df.posthocSC, file="Results/ANOVATables/PostHocSC_Indi_logLMM.csv")
# write.table(df.posthoc, file="Results/ANOVATables/FpR2_Indi_logLMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
# write.table(df.posthocAC, file="Results/ANOVATables/FpR2_Indi_logLMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
# write.table(df.posthocSC, file="Results/ANOVATables/FpR2_Indi_logLMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



### log normal LMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(ls.models[[k]]),Anova(ls.models[[k]], type="II")))
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
  boxplot(E1 ~ indices$samcam[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$agsam[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(ls.models)[k], outer=TRUE)
  
  indices$y <- df.response1[,k]
  
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


# Post Hoc Tukey Tests  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source("Analysis/F2_Nema_ContrastMatrix.R")

indices <- droplevels(indices)
PH.list <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$y <- df.response1[outlier[[i]],i]
  model <- glmer(y ~ agsam + (1|field.ID), family=gaussian(link="log"), indices2) # No Difference
  posthoc <- glht(model, linfct=mcp(agsam=cm1))
  #   posthoc.ci <- confint(posthoc)
  #   posthoc.sig <- which(posthoc.ci$confint[,2]>0)
  #   data.frame(names(potshoc.ci))
  posthoc <- summary(posthoc, test = adjusted(type = "bonferroni"))
  print(posthoc)
  nam <- paste("glht",i,names(df.response1[i]), sep = ".")
  PH.list[[i]] <- assign(nam, posthoc)
}

#####

# Plot Errorbars
phfig1 <- ggplot(endad.pw.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prediction plots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


testdata = expand.grid(age_class=unique(indices$age_class),
                       samcam = unique(indices$samcam))

test.list <- list()

for(i in 1:p) {
  X <- model.matrix(~ age_class*samcam, data = testdata)
  testdata$fit <- X %*% fixef(ls.models[[i]])
  testdata$SE <- sqrt(  diag(X %*%vcov(ls.models[[i]]) %*% t(X))  )
  testdata$upr=testdata$fit+1.96*testdata$SE
  testdata$lwr=testdata$fit-1.96*testdata$SE
  nam <- paste("tdata",i,names(df.response1[i]), sep = ".")
  test.list[[i]] <- assign(nam, testdata)
}


for (i in 1:p) {
  print(ggplot(test.list[[i]], aes(x = age_class, y = fit, ymin = lwr, ymax = upr)) + 
          geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
          geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
          facet_grid(.~samcam) +
          geom_hline(xintercept = 1, size=0.15) +
          ylab("Nematodes?") +
          xlab("Age Class") +
          scale_x_discrete(labels=c("Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
          mytheme +
          theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)))
}

#nema <- fety.backup[!indices.backup$age_class %in% "A_Cm",-c(1:4)]
df.response1$age_class <- indices$age_class

for (i in 1:p) {
  df.response1$response <- df.response1[,i]
  print(ggplot(test.list[[i]], aes(x = age_class, y = fit)) + 
          #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
          geom_point(aes(x=as.numeric(age_class)+0.3),pch=23, bg="aquamarine2") + 
          geom_errorbar(aes(x=as.numeric(age_class)+0.3, ymin = lwr, ymax = upr),position = position_dodge(1),col="black",width=0.15, size=0.15) + 
          geom_boxplot(aes(y=response), data=df.response1[outlier[[i]],]) +
          facet_grid(.~samcam) +
          geom_hline(xintercept = 1, size=0.15) +
          ylab("Nematodes?") +
          xlab("Age Class") +
          scale_x_discrete(labels=c("Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
          mytheme +
          theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)))
}

df.response1 <- df.response2[!indices.backup$age_class %in% "A_Cm",]

#####

# END