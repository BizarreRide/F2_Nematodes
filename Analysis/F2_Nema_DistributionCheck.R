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
indices <- indices.backup[!indices.backup$age_class %in% "A_Cm",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis Nematode Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.org
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 

nema.backup <- cbind(FaPro[,-c(1:5)], MaturityIndices, N=counts$counts[16:45])
nema <- nema.backup[!indices.backup$age_class %in% "A_Cm",]


p <- ncol(nema)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(2:10)
for (i in 2:10) {
  fk <-  fligner.test(nema[,i] ~ indices$agsam)
  pvector[i] <- fk$p.value
  boxplot(nema[,i]~indices$agsam)
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
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$samcam)
  car::Boxplot(indices$y ~ indices$agsam)
  title(names(nema)[i],outer=TRUE)
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
                nema.N <- -20)


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
  y <- nema[outlier[[i]],i] + 1
  qqp(y, "norm")
  qqp(y, "lnorm")
  
  #nbinom <- fitdistr(y, "negative binomial")
  #qqp(y, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
  
  poisson <- fitdistr(y, "Poisson")
  qqp(y, "pois", poisson$estimate)
  
  gamma <- fitdistr(y, "gamma")
  qqp(y, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
  
  
  title(names(nema)[i],outer=TRUE)
}


# BI <- normal or gamma
# SI <- normal, poisson or gamma
# EI <- normal, poisson or gamma
# CI <- normal, all others
# MI <- normal or gamma
# PPI <- lnorm
# MI25 <- lnorm
# sigmaMI <- normal or gamma
# sigmaMI25 <- normal or gamma
# PPI.1 <- lnorm

# nema.norm <- subset(nema, select=c("BI"))
# nema.gamma <- subset(nema, select=c("SI","EI","MI","sigmaMI","sigmaMI25"))
# nema.lnorm <- subset(nema, select=c("CI", "PPI","MI25","PPI.1"))
# 
# outlier.norm <- list(nema.BI <- -c(24,9))
# outlier.gamma <- list(nema.SI <- 1:24,
#                       nema.EI <- 1:24,
#                       nema.MI <- -17,
#                       nema.sigmaMI <- -c(3,17),
#                       nema.sigmaMI25 <- 1:24)
# outlier.lnorm <- list(nema.CI <- -24,
#                       nema.PPI <- -c(6,19),
#                       nema.MI25 <- 1:24,
#                       nema.PPI1 <- 1:24 )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#nema <- nema[,-4]
#outlier <- outlier[-4]


# normal LMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p <- ncol(nema)

p.indi.lmer <- matrix(NA,3,2+p)
colnames(p.indi.lmer) <- c("Env", "DF", colnames(nema)[1:p])

f.indi.lmer <- p.indi.lmer

indi.lmer <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$y <- nema[outlier[[i]],i]
  model <- lmer(y ~ age_class*samcam  +  (1|field.ID), indices2)
  #model <- glmer(y ~ age_class*samcam  +  (1|field.ID), family=gaussian(link="log"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("fety",i,names(nema)[i], sep = ".")
  assign(name, model)
  indi.lmer[[i]] <- assign(name, model)
  f.indi.lmer[,i+2] <- round(car::Anova(model)$"Chisq",2)
  p.indi.lmer[,i+2] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
f.indi.lmer[,1] <- p.indi.lmer[,1]  <- row.names(Anova(model))
f.indi.lmer[,2] <- p.indi.lmer[,2]  <- Anova(model)$"Df"


mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(nema)[i], sep = "."))}
names(indi.lmer)[1:p] <- mod.names
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### normal LMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
  # print(list(summary(indi.lmer[[k]]),Anova(indi.lmer[[k]], type="III")))
  #corvif(indi.lmer[[k]])
  
  E1 <- resid(indi.lmer[[k]], type="pearson")
  E2 <- resid(indi.lmer[[k]], type="response")
  F1 <- fitted(indi.lmer[[k]], type="response")
  P1 <- predict(indi.lmer[[k]], type="response")
  
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
  
  title(names(indi.lmer)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$samcam[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$agsam[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(indi.lmer)[k], outer=TRUE)
  
  indices$y <- nema[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(indi.lmer)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(indi.lmer[[k]])
  
  E1 <- resid(indi.lmer[[k]], type="pearson")
  E2 <- resid(indi.lmer[[k]], type="response")
  F1 <- fitted(indi.lmer[[k]], type="response")
  P1 <- predict(indi.lmer[[k]], type="response")
  
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
  
  title(names(indi.lmer)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(indi.lmer[[k]], obs=TRUE)
}

for(k in 1:p) {
  plot1 <- plot(infl[[k]], which = "cook", main=mod.names[k])
  print(plot1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Post Hoc Tukey Tests  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indices <- droplevels(indices)

LettersAS <- matrix(NA,8,p)
PH.list <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$y <- nema[outlier[[i]],i]
  #model <- lm(x ~ agsam + pH + c + mc + ata1, data=indices2)
  #model <- lme(x ~ agsam + pH + c + mc + clay + intensity + fertilisation + ata1 + prec1, random = ~ 1 | field.ID, data=indices2)
  model <- lmer(y ~ agsam + (1|field.ID), indices2) # No Difference
  posthoc <- glht(model, linfct=mcp(agsam="Tukey"))
  posthoc <- summary(posthoc, test = adjusted(type = "bonferroni"))
  nam <- paste("glht",i,names(nema[i]), sep = ".")
  PH.list[[i]] <- assign(nam, posthoc)
  x <- cld(posthoc)$mcletters
  xx <- x$Letters
  LettersAS[,i] <- xx
  plot(cld(posthoc), main=nam)
}
colnames(LettersAS) <- colnames(nema)
row.names(LettersAS) <- levels(indices2$agsam)

#**************************************************************************

LettersAC <- matrix(NA,4,p)
PH.list <- list()

for(i in 1:p) {
  posthoc <- glht(indi.lmer[[i]], linfct=mcp(age_class="Tukey"))
  posthoc <- summary(posthoc, test = adjusted(type = "bonferroni"))
  nam <- paste("glht",i,names(nema[i]), sep = ".")
  PH.list[i] <- assign(nam, posthoc)
  x <- cld(posthoc)$mcletters
  xx <- x$Letters
  LettersAC[,i] <- xx
  plot(cld(posthoc), main=nam)
}
colnames(LettersAC) <- colnames(nema)
row.names(LettersAC) <- levels(indices$age_class)

#**************************************************************************

LettersSC <- matrix(NA,2,p)
PH.list <- list()

for(i in 1:p) {
  posthoc <- glht(indi.lmer[[i]], linfct=mcp(samcam="Tukey"))
  posthoc <- summary(posthoc, test = adjusted(type = "bonferroni"))
  nam <- paste("glht",i,names(nema[i]), sep = ".")
  PH.list[i] <- assign(nam, posthoc)
  x <- cld(posthoc)$mcletters
  xx <- x$Letters
  LettersSC[,i] <- xx
  plot(cld(posthoc), main=nam)
}
colnames(LettersSC) <- colnames(nema)
row.names(LettersSC) <- levels(indices$samcam)

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


testdata = expand.grid(age_class=unique(indices$age_class),
                       samcam = unique(indices$samcam))

test.list <- list()

for(i in 1:p) {
  X <- model.matrix(~ age_class*samcam, data = testdata)
  testdata$fit <- X %*% fixef(indi.lmer[[i]])
  testdata$SE <- sqrt(  diag(X %*%vcov(indi.lmer[[i]]) %*% t(X))  )
  testdata$upr=testdata$fit+1.96*testdata$SE
  testdata$lwr=testdata$fit-1.96*testdata$SE
  nam <- paste("tdata",i,names(nema[i]), sep = ".")
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
nema$age_class <- indices$age_class

for (i in 1:p) {
  nema$response <- nema[,i]
  print(ggplot(test.list[[i]], aes(x = age_class, y = fit)) + 
          #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
          geom_point(aes(x=as.numeric(age_class)+0.3),pch=23, bg="aquamarine2") + 
          geom_errorbar(aes(x=as.numeric(age_class)+0.3, ymin = lwr, ymax = upr),position = position_dodge(1),col="black",width=0.15, size=0.15) + 
          geom_boxplot(aes(y=response), data=nema[outlier[[i]],]) +
          facet_grid(.~samcam) +
          geom_hline(xintercept = 1, size=0.15) +
          ylab("Nematodes?") +
          xlab("Age Class") +
          scale_x_discrete(labels=c("Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
          mytheme +
          theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)))
}

nema <- nema.backup[!indices.backup$age_class %in% "A_Cm",]

#####

# END





