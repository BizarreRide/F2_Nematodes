###########################
# F2 Nematodes
# Binomial GLMM for Nematode Feeding Types
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

#library(influence.ME)
#library(lmerTest)
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


# 1. Analysis Zof FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


fety.backup <- fety/indices$N
fety2 <- fety[!indices.backup$age_class %in% "A_Cm",-c(1:4)]

p <- ncol(fety2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Fligner Killeen Test for Heteroscedasticity: ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvector <- c(1:p)
for (i in 1:p) {
  fk <-  fligner.test(fety2[,i] ~ env.fin$agsam)
  pvector[i] <- fk$p.value
  boxplot(fety2[,i]~env.fin$agsam)
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
  indices$y <- fety2[,i]
  plot(indices$y)
  boxplot(indices$y)
  hist(indices$y, main="")
  plot(indices$y^2)
  title(names(fety2)[i],outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Detecting Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row.names(fety2) <- 1:nrow(fety2)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  indices$y <- fety[!indices.backup$age_class %in% "A_Cm",i+4]
  car::Boxplot(indices$y ~ indices$age_class)
  car::Boxplot(indices$y ~ indices$samcam)
  car::Boxplot(indices$y ~ indices$agsam)
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
outlier <- list(fety.fungi <- -9,
                fety.bacti <- -20,
                fety.carni <- -20,
                fety.omni <- -c(22,17),
                fety.Tyli <- -13,
                fety.herbi <- -19,
                fety.herbi2 <- -21,
                fety.fungi2 <- 1:24)


# change factor properties
indices$age_class2 <- indices$age_class
indices$age_class <- as.ordered(indices$age_class)
indices$age_class <- indices$age_class2

indices$samcam2 <- indices$samcam
indices$samcam <- indices$nsamcam
indices$samcam <- indices$samcam2

str(indices)

# Bi☺nomial GLMM ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.fety.biglmer <- matrix(NA,3,2+p)
colnames(p.fety.biglmer) <- c("Env", "DF", colnames(fety2)[1:p])

f.fety.biglmer <- p.fety.biglmer

fety.biglmer <- list()

for(i in 1:p) {
  indices2 <- indices[outlier[[i]],]
  indices2$scs <- fety2[outlier[[i]],i]
  indices2$fail <- indices2[,"N"] - fety2[outlier[[i]],i]
  model <- glmer(cbind(scs,fail) ~ age_class*samcam  + (1|ID) + (1|field.ID), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  # I set REML to FALSE since m random factors are nested and i have only one random factor, and the data are balanced
  # if it is "disregarded in glmer() it is OK
  print(summary(model))
  name <- paste("fety",i,names(fety2)[i], sep = ".")
  assign(name, model)
  fety.biglmer[[i]] <- assign(name, model)
  f.fety.biglmer[,i+2] <- round(car::Anova(model)$"Chisq",2)
  p.fety.biglmer[,i+2] <- round(car::Anova(model)$"Pr(>Chisq)",3)
}
f.fety.biglmer[,1] <- p.fety.biglmer[,1]  <- row.names(Anova(model))
f.fety.biglmer[,2] <- p.fety.biglmer[,2]  <- Anova(model)$"Df"

mod.names <- c(1:p)
for(i in 1:p) { mod.names[i] <- c(paste("fety",i,names(fety2)[i], sep = "."))}
names(fety.biglmer)[1:p] <- mod.names


fety.rsquared <- matrix(NA,2,p)
row.names(fety.rsquared) <- c("R2m", "R2c")
colnames(fety.rsquared) <- colnames(fety2)

for(i in 1:p) {
  fety.rsquared[,i] <- MuMIn::r.squaredGLMM(fety.biglmer[[i]])
}

#dispersion_glmer(model)

# save(list=c("f.fety.biglmer","p.fety.biglmer","fety.rsquared"), file="Results/CHi2+p_Fety_bnGLMM.rda")
# write.csv(f.fety.biglmer, file="Results/Chi2_Fety_bnGLMM.csv")
# write.csv(p.fety.biglmer, file="Results/p_Fety_bnGLMM.csv")
# write.csv(fety.rsquared, file="Results/r2_Fety_bnGLMM.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Post Hoc data inspection with lsmeans package ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detach("package:lmerTest", unload=TRUE)
require(lsmeans)
require(multcompView)

fety.lsm <- list()
fety.LettersAS <- matrix(NA,8,p)
fety.estimates <- matrix(NA,8,p)
fety.pairs <- matrix(NA,8,(4*p))
colnames(fety.LettersAS) <-  colnames(fety2)[1:p]
colnames(fety.estimates) <-  colnames(fety2)[1:p]
colnames(fety.pairs) <- rep(colnames(fety2)[1:p], each=4)

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(fety.biglmer[[i]],  ~ age_class*samcam, contr= "cld")
  x <- cld(lsm, type = "response")
  fety.LettersAS[,i] <- x$".group"
  fety.estimates[,i] <- x$"prob"
  fety.pairs[,((4*i)-3)] <- x$"age_class"
  fety.pairs[,((4*i)-2)] <- x$"samcam"
  fety.pairs[,((4*i)-1)] <- x$"prob"
  fety.pairs[,((4*i)-0)] <- x$".group"
  print(cld(lsm, type = "response"))
  name <- paste("lsm",i,names(fety2)[i], sep = ".")
  fety.lsm[[i]] <- assign(name, lsm)
  # to see the results graphically
  p1 <- plot(lsm, by = "samcam", intervals = TRUE, type = "response")
  print(p1)
  #title(names(ncr.biglmer)[i], outer=TRUE)
}

#  save(list=c("f.fety.biglmer","p.fety.biglmer","fety.rsquared","fety.LettersAS", "fety.estimates", "fety.pairs"), file="Results/CHi2+p_Fety_bnGLMM.rda")
#  write.csv(fety.LettersAS, file="Results/letters_Fety_bnGLMM.csv")
#  write.csv(fety.estimates, file="Results/estim_Fety_bnGLMM.csv")
#  write.csv(fety.pairs, file="Results/pairs_Fety_bnGLMM.csv")

for(i in 1:p){
  lsmip(fety.lsm[[i]], age_class ~ samcam, type = "response")
  summary(pairs(fety.lsm[[i]]), type = "response")
  summary(pairs(regrid(fety.lsm[[i]])), type = "response")
}

#************************************************************************


fety.lsmAC <- list()

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(fety.biglmer[[i]], pairwise ~ age_class|samcam)
  print(summary(lsm, type = "response"))
  name <- paste("lsm",i,names(fety2)[i], sep = ".")
  fety.lsmAC[[i]] <- assign(name, lsm)
}

#************************************************************************


fety.lsmSC <- list()

for (i in 1:p) {
  # get the results on a back transformed scale:
  lsm <- lsmeans(fety.biglmer[[i]], pairwise ~ samcam|age_class)
  summary(lsm, type = "response")
  name <- paste("lsm",i,names(fety2)[i], sep = ".")
  fety.lsmSC[[i]] <- assign(name, lsm)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Overdispersion ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


### binomial GLMM - Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(k in 1:p){ 
 # print(list(summary(fety.biglmer[[k]]),Anova(fety.biglmer[[k]], type="III")))
  #corvif(fety.biglmer[[k]])
  
  E1 <- resid(fety.biglmer[[k]], type="pearson")
  E2 <- resid(fety.biglmer[[k]], type="response")
  F1 <- fitted(fety.biglmer[[k]], type="response")
  P1 <- predict(fety.biglmer[[k]], type="response")
  
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
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
  # Normal QQ Plots
  qqnorm(E2)
  qqline(E2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ indices$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$samcam[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ indices$agsam[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
  indices$y <- fety2[,k]
  
  par(mfrow=c(1,1))
  scatter.smooth(F1,indices$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,indices$y,label=interaction(indices$field.ID,indices$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  title(names(fety.biglmer)[k], outer=TRUE)
  
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Residuals against variables not in the model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(k in 1:p){ 
  #corvif(fety.biglmer[[k]])
  
  E1 <- resid(fety.biglmer[[k]], type="pearson")
  E2 <- resid(fety.biglmer[[k]], type="response")
  F1 <- fitted(fety.biglmer[[k]], type="response")
  P1 <- predict(fety.biglmer[[k]], type="response")
  
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
  
  title(names(fety.biglmer)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Influence measures ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infl <- list()
for (k in 1:p) {
  infl[[k]] <- influence(fety.biglmer[[k]], obs=TRUE)
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
  indices2$scs <- fety2[outlier[[i]],i]
  indices2$fail <- indices2[,"N"] - fety2[outlier[[i]],i]
  model <- glmer(cbind(scs, fail) ~ agsam + (1|ID) + (1|field.ID), family=binomial(link="logit"), indices2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) # No Difference
  posthoc <- glht(model, linfct=mcp(agsam=cm1))
#   posthoc.ci <- confint(posthoc)
#   posthoc.sig <- which(posthoc.ci$confint[,2]>0)
#   data.frame(names(potshoc.ci))
  posthoc <- print(summary(posthoc, test = adjusted(type = "bonferroni")))
  nam <- paste("glht",i,names(fety2[i]), sep = ".")
  PH.list[[i]] <- assign(nam, posthoc)
}

#####

# Plot Errorbars
phfig1 <- ggplot(posthoc.ci, aes(y = lhs, x = exp(estimate), xmin = exp(lwr), xmax = exp(upr))) + 
  geom_errorbarh() + 
  geom_point() + 
  geom_vline(xintercept = 1) +
  mytheme
phfig1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prediction plots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indices <- droplevels(indices)

testdata = expand.grid(age_class=unique(indices$age_class),
                       samcam = unique(indices$samcam))

test.list <- list()

for(i in 1:p) {
X <- model.matrix(~ age_class*samcam, data = testdata)
testdata$fit <- X %*% fixef(fety.biglmer[[i]])
testdata$SE <- sqrt(  diag(X %*%vcov(fety.biglmer[[i]]) %*% t(X))  )
testdata$upr=testdata$fit+1.96*testdata$SE
testdata$lwr=testdata$fit-1.96*testdata$SE
nam <- paste("tdata",i,names(fety2[i]), sep = ".")
test.list[[i]] <- assign(nam, testdata)
}

fety2 <- fety.backup[!indices.backup$age_class %in% "A_Cm",-c(1:4)]
fety2$age_class <- indices$age_class

for (i in 1:p) {
  fety2$response <- fety2[,i]
  print(ggplot(test.list[[i]], aes(x = age_class, y = exp(fit))) + 
          #geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
          geom_point(aes(x=as.numeric(age_class)+0.5),pch=23, size=5, bg="aquamarine2") + 
          geom_errorbar(aes(x=as.numeric(age_class)+0.5, ymin = exp(lwr), ymax = exp(upr)),position = position_dodge(1),col="black",width=0.15, size=0.15) + 
          geom_boxplot(aes(y=response, fill=age_class), data=fety2[outlier[[i]],]) +
          facet_grid(.~samcam) +
          geom_hline(xintercept = 1, size=0.15) +
          ylab("Nematodes?") +
          xlab("Age Class") +
          scale_x_discrete(labels=c("Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
          theme_bw() +
          theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1)))
}

fety2 <- fety[!indices.backup$age_class %in% "A_Cm",-c(1:4)]

#####

# END