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

p.fety.biglmer.crop <- matrix(NA,3,2+p)
colnames(p.fety.biglmer.crop) <- c("Env", "DF", colnames(fety2)[1:p])

f.fety.biglmer.crop <- p.fety.biglmer.crop

fety.biglmer.crop <- list()
i=1
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
  f.fety.biglmer.crop[,i+2] <- round(car::Anova(model)$"Chisq",2)
  p.fety.biglmer.crop[,i+2] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
f.fety.biglmer.crop[,1] <- p.fety.biglmer.crop[,1]  <- row.names(Anova(model))
f.fety.biglmer.crop[,2] <- p.fety.biglmer.crop[,2]  <- Anova(model)$"Df"

mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(fety2)[i], sep = "."))}
names(fety.biglmer.crop)[1:p] <- mod.names

dispersion_glmer(model)

#save(list=c("f.fety.biglmer.crop","p.fety.biglmer.crop"), file="Results/CHi2+p_Fety_bnGLMM_crop.rda")
#write.csv(f.fety.biglmer.crop, file="Results/Chi2_Fety_bnGLMM_crop.csv")
#write.csv(p.fety.biglmer.crop, file="Results/p_Fety_bnGLMM_crop.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Post Hoc data inspection with lsmeans package ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detach("package:lmerTest", unload=TRUE)
require(lsmeans)

fety.lsm <- list()

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(fety.biglmer.crop[[i]], ~ crop)
  summary(lsm, type = "response")
  name <- paste("lsm",i,names(fety2)[i], sep = ".")
  fety.lsm[[i]] <- assign(name, lsm)
  # to see the results graphically
  p1 <- plot(lsm,  intervals = TRUE, type = "response")
  print(p1)
  #title(names(fety.biglmer.crop)[i], outer=TRUE)
}


for(i in 1:p){
  lsmip(fety.lsm[[i]], age_class ~ samcam, type = "response")
  summary(pairs(fety.lsm[[i]]), type = "response")
  summary(pairs(regrid(fety.lsm[[i]])), type = "response")
}
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


# Post Hoc Tukey Tests  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indices <- droplevels(indices)

LettersC <- matrix(NA,2,p)
PH.list <- list()

for(i in 1:p) {
  posthoc <- glht(fety.biglmer.crop[[i]], linfct=mcp(crop="Tukey"))
  posthoc <- summary(posthoc, test = adjusted(type = "bonferroni"))
  nam <- paste("glht",i,names(fety2[i]), sep = ".")
  PH.list[[i]] <- assign(nam, posthoc)
  x <- cld(posthoc)$mcletters
  xx <- x$Letters
  LettersC[,i] <- xx
  plot(cld(posthoc), main=nam)
}
colnames(LettersC) <- colnames(fety2)
row.names(LettersC) <- levels(indices$crop)

#####
# Pairwise comparisons (with interaction term)
endad.pairwise <- glht(endad.tukey, linfct=mcp(ia.acl.smc = cm1))
endad.pw.ci <- confint(endad.pairwise)
summary(endad.pairwise, test=adjusted(type="fdr"))

# Confidence intervals including 0
endad.pw.sig <- which(endad.pw.ci$confint[,2]>0)
data.frame(names(endad.pw.sig))

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

indices <- droplevels(indices)

testdata = expand.grid(crop=unique(indices$crop))

test.list <- list()

for(i in 1:p) {
  X <- model.matrix(~ crop, data = testdata)
  testdata$fit <- X %*% fixef(fety.biglmer.crop[[i]])
  testdata$SE <- sqrt(diag(X %*%vcov(fety.biglmer.crop[[i]]) %*% t(X)))
  testdata$upr=testdata$fit+1.96*testdata$SE
  testdata$lwr=testdata$fit-1.96*testdata$SE
  nam <- paste("tdata",i,names(fety2[i]), sep = ".")
  test.list[[i]] <- assign(nam, testdata)
}

fety2 <- fety.backup[,-c(1:4)]
fety2$crop <- indices$crop

for (i in 1:p) {
  fety2$response <- fety2[,i]
  print(ggplot(test.list[[i]], aes(x = crop, y = exp(fit))) + 
          #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
          geom_point(aes(x=as.numeric(crop)+0.5),pch=23, size=5, bg="aquamarine2") + 
          geom_errorbar(aes(x=as.numeric(crop)+0.5, ymin = exp(lwr), ymax = exp(upr)),position = position_dodge(1),col="black",width=0.15, size=0.15) + 
          geom_boxplot(aes(y=response, fill=crop), data=fety2[outlier[[i]],]) +
          geom_hline(xintercept = 1, size=0.15) +
          ylab("Nematodes?") +
          xlab("Age Class") +
          scale_x_discrete(labels=c("Maize", "Silphie")) +
          theme_bw() +
          theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)))
}

fety2 <- round(fety[,-c(1:4)],0)

#####

# END

